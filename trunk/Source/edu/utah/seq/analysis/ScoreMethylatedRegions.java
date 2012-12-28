package edu.utah.seq.analysis;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.data.Point;
import edu.utah.seq.data.PointData;
import edu.utah.seq.data.SmoothingWindow;
import edu.utah.seq.useq.data.Region;

import trans.anno.*;
import trans.misc.*;
import trans.roc.ParseSgrsForParticularRegions;
import trans.roc.Positive;
import trans.roc.PositiveComparator;
import util.bio.parsers.*;
import util.bio.seq.Seq;
import util.gen.*;
import util.bio.annotation.*;

/**
 * Returns values that overlap particular regions, can generate p-values associated with random back ground model 
 */
public class ScoreMethylatedRegions {

	//fields
	private File[] convertedPointDirs;
	private File[] nonConvertedPointDirs;
	private int minimumNumberObservations = 10;
	private int minimumReadCoverage = 8;
	private int minimumNumberBaseFractions = 1;
	private boolean printAll = true;
	private Positive[] regions;
	private String currentChrom = "";
	private int numberRandom = 1000;
	private boolean makeRandom = false;
	private HashMap<String, File> chromGCFileMap;
	private File gcGenomeDir;
	private boolean[] currentGCContent;
	private double fractionGCTolerance = 0.1;
	private double fractionObsTolerance = 0.3;
	private HashMap<String, Region[]> interrogatedRegions;
	private Region[] currentInterrogatedRegions;
	private int totalBPInterrogatedRegions;
	private int[] startsBPInterrogatedRegions;
	double totalNumberRegions = 1;
	private RegionResult[] regionResults;
	private PrintWriter out;

	//from BisStat
	private HashMap<String,PointData[]> convertedPlusPointData;
	private HashMap<String,PointData[]> convertedMinusPointData;
	private HashMap<String,PointData[]> nonConvertedPlusPointData;
	private HashMap<String,PointData[]> nonConvertedMinusPointData;

	private PointData convertedMergedChromPlus = null;
	private PointData convertedMergedChromMinus = null;	
	private PointData nonConvertedMergedChromPlus = null;
	private PointData nonConvertedMergedChromMinus = null;
	private MethylatedBaseObservationOneSample[] methylatedBases;
	private int[] methylatedBasePositions;


	public void printDefaultValues(){
		System.out.println();
		System.out.println(numberRandom+"\tNumber Random Sets");
		System.out.println(fractionGCTolerance+"\tFraction GC Tolerance");
		System.out.println(fractionObsTolerance+"\tFraction Obs Tolerance");
	}
	
	public ScoreMethylatedRegions(String[] args){

		processArgs(args);

		System.out.println();
		System.out.println(minimumNumberObservations + "\tMinimum Number Observations");
		System.out.println(minimumReadCoverage + "\tMinimum Read Coverage");
		System.out.println(minimumNumberBaseFractions + "\tMinimum Number Cs passing Read Coverage in region");
		if(makeRandom) printDefaultValues();
		System.out.println();

		//print header ave score, fold enrich, pvalUp, pvalDown
		if (makeRandom) out.println("#Chrom\tStart\tStop\tNumBaseFractions\tMeanBaseFractions\tMedianBaseFractions\tRegionSumFractionMeth\tNumNonCon\tNumTotal\tNumNonConPlus\tNumConPlus\tNumNonConMinus\tNumConMinus\tAveRandomFractionMeth\tFoldDiff(Real/Random)\tBonCorrPValEnriched\tBonCorrPValReduced");
		else out.println("#Chrom\tStart\tStop\tNumBaseFractions\tMeanBaseFractions\tMedianBaseFractions\tRegionSumFractionMeth\tNumNonCon\tNumTotal\tNumNonConPlus\tNumConPlus\tNumNonConMinus\tNumConMinus");

		//for each region
		regionResults = new RegionResult[regions.length];
		
		ArrayList<ArrayList<Float>> randomResultsAgg = new ArrayList<ArrayList<Float>>();
		ArrayList<Float> regionResultsAgg = new ArrayList<Float>();
		
		for (int x=0; x< regions.length; x++){
			
			//new chromosome? load data!
			if (regions[x].getChromosome().equals(currentChrom) == false) {
				currentChrom = regions[x].getChromosome();
				System.out.print("\n"+currentChrom+" ");

				//load PointData
				fetchDataAndRemove();
				//check they all exist
				if (convertedMergedChromPlus == null || convertedMergedChromMinus == null || nonConvertedMergedChromPlus == null || nonConvertedMergedChromMinus == null) {
					Misc.printErrAndExit("\nError: missing all four PointDatasets for -> "+currentChrom +"\n");
				}
				//making random?
				if (makeRandom) loadRandomData();
				
				//merge strands
				ArrayList<PointData> al = new ArrayList<PointData>();
				al.add(convertedMergedChromMinus);
				al.add(convertedMergedChromPlus);
				
				PointData mergedCon = PointData.mergePointData(al, false, false);
				al.clear();
				al.add(nonConvertedMergedChromMinus);
				al.add(nonConvertedMergedChromPlus);
				PointData mergedNonCon = PointData.mergePointData(al, false, false);
				methylatedBases = MethylatedBaseObservationOneSample.fetchCommonBasesWithMinimumObservations(mergedNonCon, mergedCon, minimumReadCoverage);
				methylatedBasePositions = MethylatedBaseObservationOneSample.fetchPositions(methylatedBases);
			}

			//load select region with PointData
			loadRegion(regions[x]);
			regionResults[x] = new RegionResult(regions[x].getIndex());

			//fetch scores numberBaseFractions, mean, median, fractionNonConverted, numberNonCon, numberNonCon+numberCon
			float[] scores = (float[])(regions[x].getScores().get(0));
			int numberObservations = (int)scores[5];

			//enough observations?
			if (numberObservations < minimumNumberObservations || scores[0] < minimumNumberBaseFractions){
				regionResults[x].results = regions[x].toStringSimpleFloat()+"\tToo few observations";
				continue;
			}
			

			//make random matched regions?
			if (makeRandom) {
				ArrayList<Float> f = makeRandom(regions[x],numberObservations);
				regionResults[x] = makeRandomInd(regions[x],numberObservations,f);
				
				if (f.size() == numberRandom) {
					randomResultsAgg.add(f);
					regionResultsAgg.add(scores[3]);
				}
			} else {
				regionResults[x].results = regions[x].toStringSimpleFloat();
				regionResultsAgg.add(scores[3]);
			}

		}
		
		System.out.println();
		
		
		
		//sort and print results
		Arrays.sort(regionResults);
		for (RegionResult rr : regionResults) {
			if (printAll || rr.results.contains("Too few") == false) out.println(rr.results);
		}
		out.close();
		
		this.printAggregateData(randomResultsAgg, regionResultsAgg, regions.length);

		System.out.println("\nDone!\n");
	}

	
	private class RegionResult implements Comparable<RegionResult>{
		int index;
		String results;
		
		public RegionResult (int index){
			this.index = index;
		}

		public int compareTo(RegionResult other) {
			if (other.index< this.index) return 1;
			if (other.index> this.index) return -1;
			return 0;
		}
	}

	public void loadRandomData(){
		//load gc
		File file = chromGCFileMap.get(currentChrom);
		if (file == null ) Misc.printExit("\nError: cannot find a gc binary file for "+currentChrom);
		currentGCContent = (boolean[])IO.fetchObject(file);
		//load interrogatedRegions
		currentInterrogatedRegions = interrogatedRegions.get(currentChrom);
		if (currentInterrogatedRegions == null) Misc.printExit("\nError: cannot find interrogated regions for "+currentChrom);
		totalBPInterrogatedRegions = Region.totalBP(currentInterrogatedRegions);
		startsBPInterrogatedRegions = Region.startsInBases(currentInterrogatedRegions);
	}
	
	public double getMedian(ArrayList<Float> values) {
		int num = values.size();
		if (num == 0) return Float.NaN;
		if (num == 1) return values.get(0);
		if (values.size() == 2) {
			float total = values.get(0) + values.get(1);
			return total/ 2.0f;
		}
		//calc median
		float[] f = Num.arrayListOfFloatToArray(values);
		Arrays.sort(f);
		return Num.median(f);
	}
	
	private void printAggregateData(ArrayList<ArrayList<Float>> randomData, ArrayList<Float> regionData,int length) {
		System.out.println(String.format("Aggregating data for %d of %d regions\n",regionData.size(),length));
		
		double realMed = getMedian(regionData);
		ArrayList<ArrayList<Float>> realData = new ArrayList<ArrayList<Float>>();
		int maxSize = Integer.MIN_VALUE;
		
		if (randomData != null) {
			int numGreaterThan = 0;
			int numLessThan = 0;
			double average = 0;
			double med;
			
			//Find the largest row in the data
			for (ArrayList<Float> row: randomData) {
				if (row.size() > maxSize) {
					maxSize = row.size();
				}
			}
			
			//Transpose random data
			for (int i=0; i<maxSize;i++) {
				ArrayList<Float> col = new ArrayList<Float>();
				for (ArrayList<Float> rowdata: randomData) {
					try {
						col.add(rowdata.get(i));
					} catch (IndexOutOfBoundsException iobex) {
						//ArrayLists might not all have the same size, so this is here to make sure
						//errors aren't thrown.
					}
				}
				realData.add(col);
			}
			
			for (ArrayList<Float> rd: realData) {
				med = getMedian(rd);
				average += med;
				
				if (med > realMed) {
					numGreaterThan += 1;
				} else if (med < realMed) {
					numLessThan += 1;
				} else {
					numGreaterThan += 1;
					numLessThan += 1;
				}
			}
			

			average = average / ((double)realData.size());
			double numGreaterThanP = numGreaterThan/ ((double)realData.size());
			if (numGreaterThanP > 1) numGreaterThanP = 1;
			double numLessThanP = numLessThan/ ((double)realData.size());
			if (numLessThanP > 1) numLessThanP = 1;
			double foldEnrich = Math.log(realMed/average)/Math.log(2);
			System.out.println(String.format("Median methylation in the targed regions is: %.3f",realMed));
			System.out.println(String.format("Median methylation in randomly generated regions: %.3f.\n" +
											 "P-value: Targeted region more enriched than random regions: %.3f (%d/%d).\n" + 
											 "P-value: Targeted region less enriched than random regions: %.3f (%d/%d).\n" +
											 "Fold-Enrichment target regions / random regions: %.3f.",
					average,numGreaterThanP,numGreaterThan,realData.size(),numLessThanP,numLessThan,realData.size(),foldEnrich));
		}
		
	}
	
	private RegionResult makeRandomInd(Positive region, double numberObservations, ArrayList<Float> fractionsNonConverted) {
		//ArrayList<Float> fractionsNonConverted = makeRandom(region,numberObservations);
		
		RegionResult rr = new RegionResult(region.getIndex());
		
		if (fractionsNonConverted.size() != numberRandom) {
			rr.results = region.toStringSimpleFloat()+"\tSkipping, too few random regions found";
		}
		else {
			ArrayList regionScoresAL = region.getScores();
			float[] regionScores = (float[])regionScoresAL.get(0);
			double real = regionScores[3];

			//calculate p-value up and down. 
			//count number > and <
			double numGreaterThan = 0;
			double numLessThan = 0;
			double average =0;
			for (int i=0; i<fractionsNonConverted.size(); i++){
				average += fractionsNonConverted.get(i);
				if (fractionsNonConverted.get(i) > real) numGreaterThan++;
				else if (fractionsNonConverted.get(i) < real ) numLessThan++;
				//same score
				else{
					numGreaterThan++;
					numLessThan++;
				}
			}

			//calc ave, and pvalues, bonferroni corrected
			average = average/ ((double)numberRandom);
			numGreaterThan =  numGreaterThan/ ((double)numberRandom);
			if (numGreaterThan > 1) numGreaterThan = 1;
			numLessThan = numLessThan/ ((double)numberRandom);
			if (numLessThan > 1) numLessThan = 1;
			double foldEnrich = Math.log(real/average)/Math.log(2);
			//save results 
			//coord, ave score, fold enrich, pvalUp, pvalDown
			rr.results = region.toStringSimpleFloat()+"\t"+average+"\t"+foldEnrich+"\t"+numGreaterThan+"\t"+numLessThan;
		}
		
		return rr;
		
	}
	
	
	/**Takes a Positive and finds 1000 chrom, length, gc, and number scores matched random regions.*/
	private ArrayList<Float> makeRandom(Positive region, double numberObservations){
		//calculate gc content of real region and set min max
		double realGC = calculateFractionGCContent(region.getStart(), region.getStop());		
		int sizeRealRegionMinOne = region.getStop() - region.getStart();
		double minGC = realGC - fractionGCTolerance;
		double maxGC = realGC + fractionGCTolerance;

		//calculate min max number of scores
		double diff = numberObservations * fractionObsTolerance;
		int minObs = (int)Math.round(numberObservations - diff);
		if (minObs < minimumNumberObservations) minObs = minimumNumberObservations; 
		
		//try 100,000 times to find a gc and min num matched random region
		Random randomGenerator = new Random();
		ArrayList<Float> fractionsNonConverted = new ArrayList<Float>();
		
		for (int x=0; x<numberRandom; x++){
			for (int y=0; y< 10000; y++){
				//randomly pick an interrogatedRegion unbiased by size
				//pick a base from total
				int randomBase = randomGenerator.nextInt(totalBPInterrogatedRegions);
				//find it's index
				int index = Arrays.binarySearch(startsBPInterrogatedRegions, randomBase);
				//pull the region
				int indexInterrogatedRegion;
				if (index >=0) indexInterrogatedRegion = index;
				else indexInterrogatedRegion = -index -1;
				Region interrogatedRegion = currentInterrogatedRegions[indexInterrogatedRegion];
				//get it's size
				int lengthIntReg = interrogatedRegion.getLength();
				//check size to be sure it's big enough
				if (lengthIntReg < sizeRealRegionMinOne) continue;
				//generate a random start stop matching the length of the region
				int start = randomGenerator.nextInt(lengthIntReg) + interrogatedRegion.getStart();
				//check to see if start is < 0
				if (start < 0) continue;
				int stop = start + sizeRealRegionMinOne;
				//check to see if stop goes past the stop of the interrogated region
				if (stop >= interrogatedRegion.getStop()) continue;
				//check gc
				double testGC = calculateFractionGCContent(start, stop);
				if (testGC < minGC || testGC > maxGC) continue;
				//check num scores
				Positive testRegion = new Positive(currentChrom, start, stop);
				//load region with scores
				loadRegion(testRegion);
				//correct number?
				ArrayList testAL = testRegion.getScores();
				float[] scores = (float[])testAL.get(0);
				int numTestScores = (int)scores[5];
				if (numTestScores < minObs ) continue;
				//otherwise add fractionNonConverted

				fractionsNonConverted.add(scores[3]);
				break;
			}
		}
		//System.out.print(".");
		return fractionsNonConverted;
		
	}


	/**Calculates the gc content.*/
	public double calculateFractionGCContent(int start, int stop){
		double ave = 0;
		int realStop = stop +1;
		if (realStop >= currentGCContent.length){
			realStop = currentGCContent.length -1;
		}
		for (int i = start; i< realStop; i++){
			if (currentGCContent[i]) ave++;
		}
		double length = realStop - start;
		return ave/ length;
	}



	/**Fetchs the data for a particular chromosome.*/
	public void fetchDataAndRemove(){
		ArrayList<PointData> al = null;
		PointData[] pd;
		//merge converted
		convertedMergedChromPlus = null;
		if (convertedPlusPointData.containsKey(currentChrom)) {
			pd = convertedPlusPointData.remove(currentChrom);
			al = PointData.convertArray2ArrayList(pd);
			convertedMergedChromPlus = PointData.mergePointData(al, false, false);
		}
		convertedMergedChromMinus = null;
		if (convertedMinusPointData.containsKey(currentChrom)) {
			pd = convertedMinusPointData.remove(currentChrom);
			al = PointData.convertArray2ArrayList(pd);
			convertedMergedChromMinus = PointData.mergePointData(al, false, false);
		}
		//merge nonConverted
		nonConvertedMergedChromPlus = null;
		if (nonConvertedPlusPointData.containsKey(currentChrom)) {
			pd = nonConvertedPlusPointData.remove(currentChrom);
			al = PointData.convertArray2ArrayList(pd);
			nonConvertedMergedChromPlus = PointData.mergePointData(al, false, false);
		}
		nonConvertedMergedChromMinus = null;
		if (nonConvertedMinusPointData.containsKey(currentChrom)) {
			pd = nonConvertedMinusPointData.remove(currentChrom);
			al = PointData.convertArray2ArrayList(pd);
			nonConvertedMergedChromMinus = PointData.mergePointData(al, false, false);
		}
		pd = null;
		al = null;
	}

	/**Loads a Positive with fractionNonConverted.*/
	public void loadRegion(Positive region){
		int start = region.getStart();
		int stop = region.getStop();
		//nonConPlus, conPlus, nonConMinus, conMinus
		double[] counts = scoreRegion(start, stop);
		//calc fraction methylated, numNonCon/total
		double numberNonConPositions = counts[0]+ counts[2];
		double total = Num.sumArray(counts);
		double fractionNonConverted;
		if (total == 0) fractionNonConverted = 0;
		else fractionNonConverted = numberNonConPositions/ total;
		
		//calculate mean base fraction methylation
		float[] med = scoreBaseFractionMethylation(start, stop);
		
		float[] scores = new float[]{med[0], med[1], med[2], (float)fractionNonConverted,(float)numberNonConPositions, (float)total, (float)counts[0], (float)counts[1], (float)counts[2], (float)counts[3]};
		
		//set scores
		ArrayList al = region.getScores();
		al.add(scores);
		region.setScores(al);
	}

	/**Returns number of fractions in region and the median*/
	private float[] scoreBaseFractionMethylation(int start, int stop) {
		int[] indexes = Num.findIndexes(start, stop, methylatedBasePositions);

		//for each index, calculate fraction methylated
		float[] fractions = new float[indexes[1]- indexes[0]];
		int f = 0;
		for (int x=indexes[0]; x< indexes[1]; x++) {
			fractions[f++] = methylatedBases[x].getFractionMethylation();
		}
		float median = -1;
		float mean = -1;
		if (fractions.length == 1) {
			median = fractions[0];
			mean = median;
		}
		else if (fractions.length == 2) {
			median = Num.mean(fractions);
			mean = median;
		}
		else if (fractions.length > 2) {
			Arrays.sort(fractions);
			median = (float) Num.median(fractions);
			mean = Num.mean(fractions);
		}
		return new float[]{fractions.length, mean, median};
	}

	/**Returns nonConPlus, conPlus, nonConMinus, conMinus*/
	private double[] scoreRegion (int start, int stop){
		//treatment
		double conPlus = convertedMergedChromPlus.sumScoreBP(start, stop);
		double nonConPlus = nonConvertedMergedChromPlus.sumScoreBP(start, stop);
		double conMinus = convertedMergedChromMinus.sumScoreBP(start, stop); 
		double nonConMinus = nonConvertedMergedChromMinus.sumScoreBP(start, stop);
		return new double[]{nonConPlus, conPlus, nonConMinus, conMinus};
	}

	public static void main(String[] args) {
		if (args.length==0){
			printDocs();
			System.exit(0);
		}	
		new ScoreMethylatedRegions(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File interrogatedRegionsFile = null;
		File regionFile = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': convertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'n': nonConvertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'r': regionFile = new File(args[i+1]); i++; break;
					case 'g': gcGenomeDir = new File(args[i+1]); i++; makeRandom = true; break;
					case 'u': numberRandom = Integer.parseInt(args[i+1]); i++; break;
					case 'm': minimumNumberObservations = Integer.parseInt(args[i+1]); i++; break;
					case 'o': minimumReadCoverage = Integer.parseInt(args[++i]); break;
					case 'b': minimumNumberBaseFractions = Integer.parseInt(args[++i]); break;
					case 'p': printAll = false; break;
					case 'i': interrogatedRegionsFile = new File(args[++i]); makeRandom = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (convertedPointDirs == null || nonConvertedPointDirs == null || regionFile == null || regionFile.exists() == false) Misc.printExit("\nError: Please set the -c, -n, and or -r parameters.\n");

		//only one directory look deeper
		if (convertedPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(convertedPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) convertedPointDirs = otherDirs;
		}
		if (nonConvertedPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(nonConvertedPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) nonConvertedPointDirs = otherDirs;
		}
		
		
		//load PointData maps
		//fetch converted PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (convertedPointDirs);
		convertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		convertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
		combo = PointData.fetchStrandedPointDataNoMerge (nonConvertedPointDirs);
		nonConvertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		nonConvertedMinusPointData = PointData.convertArrayList2Array(combo[1]);

		//load and sort Regions File
		regions = ParseSgrsForParticularRegions.parseRegionFile(regionFile);
		Arrays.sort(regions, new PositiveComparator());

		//make hash of chromosome text and gc boolean file?
		if (makeRandom){
			File[] gcFiles = IO.extractFiles(gcGenomeDir, ".gc");
			chromGCFileMap = Seq.makeChromosomeNameFileHash( gcFiles );
			if (chromGCFileMap == null) Misc.printExit("\nError parsing xxx.gc files.\n");
			if (interrogatedRegionsFile == null) Misc.printExit("\nError, cannot find your interrogated regions file!\n");

			//find smallest region
			int shortestLength = Integer.MAX_VALUE;
			for (int i=0; i< regions.length; i++){
				int size = regions[i].getLength();
				if (size < shortestLength) shortestLength = size;
			}
			totalNumberRegions = regions.length;
			interrogatedRegions = Region.parseStartStops(interrogatedRegionsFile, 0, 0, shortestLength);
		}
		
		//make print writer for results
		try {
			File results = new File (Misc.removeExtension(regionFile.getCanonicalPath())+"_ScrMethReg.txt");
			System.out.println("\nSaving results to -> "+results);
			out = new PrintWriter( new FileWriter (results));
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}	


	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Score Methylated Regions: March 2012                      **\n" +
				"**************************************************************************************\n" +
				"For each region finds the underlying methylation data. A p-value (Bon Corr) for each\n" +
				"region's fraction methylated (# nonConObs/ # totalObs) as well as a fold enrichment\n" +
				"can be calculated using regions randomly drawn matched by chromosome, region length,\n" +
				"# obs, and GC content.\n\n" +

				"Options:\n"+
				"-c Converted PointData directories, full path, comma delimited. These should\n" +
				"       contain stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories. See the NovoalignBisulfiteParser app.\n" +
				"-n Non-converted PointData directories, ditto. \n" +
				"-r Full path file text for your region of interest file (tab delim: chr start stop).\n" +
				"-g To calculate p-values for methylation enrichment/ reduction,  provide a full path\n"+
				"         directory containing for chromosome specific gc content boolean arrays. See\n"+
				"         the ConvertFasta2GCBoolean app. Complete option -i\n"+
				"-i Likewise, to calculate p-values, also provide a full path file text containing the\n" +
				"         interrogated regions (tab delim: chr start stop ...) to use in drawing\n" +
				"         random regions.\n"+
				"-u Number of random region sets, defaults to 1000.\n" +
				"-m Minimum number of observations in a region to score, defaults to 10.\n"+
				"-o Minimum read coverage to count mC fraction, defaults to 8\n"+
				"-b Minimum number of Cs passing read coverage in region to score, defaults to 1\n"+
				"-p Print only regions that pass thresholds, defaults to all\n"+

				"\nExample: java -jar pathTo/Apps/ScoreMethylatedRegions -c /Data/Sperm/Converted -n \n" +
				"      /Data/Sperm/NonConverted -r /Res/miRNARegions.bed -i /Res/interrRegions.bed\n" +
				"       -g /Genomes/Hg18/GCBooleanArrays/\n\n" +

		"**************************************************************************************\n");		
	}

}
