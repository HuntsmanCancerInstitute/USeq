package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import org.apache.commons.math3.stat.inference.ChiSquareTest;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import edu.utah.seq.useq.apps.Bar2USeq;
import edu.utah.seq.useq.data.Region;
import trans.tpmap.*;

/**Takes PointData from the RNAEditingPileupParser and generates a variety of statistics for user defined regions.
 * @author Nix
 * */
public class DefinedRegionRNAEditing {

	//user defined fields
	private File[] convertedPointDirs;
	private File[] nonConvertedPointDirs;
	private double errorRateMinOne;
	private int minimumBaseCoverage = 5;
	private boolean runStrandedAnalysis = false;

	//internal fields
	private boolean plusStrand = true;
	private boolean minusStrand =  true;
	private int numberRandomPermutations = 10;
	private int targetNumberRandomRegions = 1000000;
	private ChiSquareTest chiSquare = new ChiSquareTest();
	private HashMap<String,PointData[]> convertedPlusPointData;
	private HashMap<String,PointData[]> convertedMinusPointData;
	private HashMap<String,PointData[]> nonConvertedPlusPointData;
	private HashMap<String,PointData[]> nonConvertedMinusPointData;
	private LinkedHashMap<String,SmoothingWindow[]> chromosomeRegion;
	private File regionsFile;
	private RandomScoreArray[] randomScores;
	private SmoothingWindow[] regions;
	private String genomeVersion = null;
	private int numberRegions;
	private float[][] regionScoresPlus = null;

	//by chromosome
	private String chromosome;
	private PointData convertedMergedChromPlus = null;
	private PointData convertedMergedChromMinus = null;	
	private PointData nonConvertedMergedChromPlus = null;
	private PointData nonConvertedMergedChromMinus = null;	
	private PointData convertedChrom = null;
	private PointData nonConvertedChrom = null;

	MethylatedBaseObservationOneSample[] editedBases;
	int[] editedBasePositions;
	private int maxSize;
	private int numberRegionsWithData;

	//constructor
	public DefinedRegionRNAEditing(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);
		
		//load data pointers
		loadPointDataArrays();
		
		if (runStrandedAnalysis){
			//scan plus
			plusStrand = true;
			minusStrand = false;
			doWork();
			collectScoresPlus();
			//scan minus
			plusStrand = false;
			minusStrand = true;
			doWork();
		}
		else doWork();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	public void doWork(){
		numberRegionsWithData = 0;

		//make containers for randomScores
		randomScores = new RandomScoreArray[maxSize+10];

		//for each chromosome
		System.out.print("Scanning Chromosomes:\n\t");

		for (String chr : chromosomeRegion.keySet()){
			chromosome = chr;
			System.out.print(chromosome + " ");
			regionScanChromosome();
		}
		System.out.println();

		//calculate p-values and FDRs
		calculateRegionPValues();

		//print spreadsheet
		printSpreadSheet();

	}
	
	public void collectScoresPlus(){
		regionScoresPlus = new float[numberRegions][];
		int index = 0;
		for (String chr: chromosomeRegion.keySet()){
			for (SmoothingWindow sw : chromosomeRegion.get(chr)){
				regionScoresPlus[index++] = sw.getScores();
				sw.setScores(null);
			}
		}

	}

	/**Writes out an excel compatible tab delimited spreadsheet with hyperlinks for IGB.*/
	public void printSpreadSheet(){
		try{
			String strand = "";
			if (regionScoresPlus != null) strand = "Stranded";
			File file = new File(regionsFile.getParentFile(), Misc.removeExtension(regionsFile.getName())+"_DRRE"+ strand+ ".xls");
			PrintWriter out = new PrintWriter (new FileWriter (file));
			//print header line
			if (regionScoresPlus == null ) out.println("#"+genomeVersion+"_IGBHyperLinks\tChr\tStart\tStop\tNumberObservations\tPseudoMedianBaseFractionEditing\t-10Log10(FDR)");
			else out.println("#"+genomeVersion+"_IGBHyperLinks\tChr\tStart\tStop\tNumberObservationsPlus\tPseudoMedianBaseFractionEditingPlus\t-10Log10(FDR)Plus\tNumberObservationsMinus\tPseudoMedianBaseFractionEditingMinus\t-10Log10(FDR)Minus");
			String url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";
			String tab = "\t";
			//for each region
			int index = 0;
			int regionIndex = 0;
			for (String chr: chromosomeRegion.keySet()){
				for (SmoothingWindow sw : chromosomeRegion.get(chr)){
					//url
					int winStart = sw.getStart() - 10000;
					if (winStart < 0) winStart = 0;
					int winEnd = sw.getStop() + 10000;
					out.print(url+chr+"&start="+winStart+"&end="+winEnd+"\",\""+(++index)+"\")\t");
					//chrom
					out.print(chr); out.print(tab);
					//start
					out.print(sw.getStart()); out.print(tab);
					//stop
					out.print(sw.getStop()); out.print(tab);
					//print plus strand scores?
					if (regionScoresPlus !=null){
						float[] scores = regionScoresPlus[regionIndex++];
						if (scores == null) {
							out.print("0\t0\t0\t"); 
						}
						else {
							//	num
							out.print((int)scores[0]); out.print(tab);
							//	pse
							out.print(scores[1]); out.print(tab);
							//	FDR
							out.print(scores[2]); out.print(tab);
						}
					}
					//scores = numberObs, pseMedian, pval
					float[] scores = sw.getScores();
					if (scores == null) {
						out.println("0\t0\t0"); 
					}
					else {
						//	num
						out.print((int)scores[0]); out.print(tab);
						//	pse
						out.print(scores[1]); out.print(tab);
						//	FDR
						out.println(scores[2]);
					}
				}
			}
			out.close();
		}catch (Exception e){
			System.out.println("\nError: problem printing spreadsheet report");
			e.printStackTrace();
		}
	}



	private class RandomScoreArray{
		int index = 0;
		float[] scores;

		public RandomScoreArray(int maxNumberScores){
			scores = new float[maxNumberScores];
		}
		/**Returns true if counted, false if full.*/
		public boolean countScore(double score){
			if (index == scores.length) return false;
			scores[index++] = (float)score;
			return true;
		}
		public void sortScores() {
			//clip at index, it's always one more than what is in the array
			float[] clipped = new float[index];
			System.arraycopy(scores, 0, clipped, 0, index);
			scores = clipped;
			Arrays.sort(scores);
		}
		public float calculateTransformedPValue (float realScore){
			int i = Num.findClosestIndexToValue(scores, realScore);
			double numScores = index - i;
			if (numScores == 0) return 1+ Num.minus10log10Float(1.0/(double)index);
			else return Num.minus10log10Float(numScores/(double)index);
		}

		public boolean isMaxed() {
			if (index == scores.length) return true;
			return false;
		}
	}

	/**Fetchs the data for a particular chromosome. Returns true if all four datasets were found.*/
	public boolean fetchData(){
		ArrayList<PointData> al = null;
		PointData[] pd;
		//merge converted
		convertedMergedChromPlus = null;
		if (plusStrand == true && convertedPlusPointData.containsKey(chromosome)) {
			pd = convertedPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			convertedMergedChromPlus = PointData.mergePointData(al, false, true);
		}
		convertedMergedChromMinus = null;
		if (minusStrand == true && convertedMinusPointData.containsKey(chromosome)) {
			pd = convertedMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			convertedMergedChromMinus = PointData.mergePointData(al, false, true);
		}
		//merge nonConverted
		nonConvertedMergedChromPlus = null;
		if (plusStrand == true && nonConvertedPlusPointData.containsKey(chromosome)) {
			pd = nonConvertedPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			nonConvertedMergedChromPlus = PointData.mergePointData(al, false, true);
		}
		nonConvertedMergedChromMinus = null;
		if (minusStrand == true && nonConvertedMinusPointData.containsKey(chromosome)) {
			pd = nonConvertedMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			nonConvertedMergedChromMinus = PointData.mergePointData(al, false, true);
		}
		pd = null;
		al = null;
		if (plusStrand){
			if (convertedMergedChromPlus == null || nonConvertedMergedChromPlus == null) return false;
			if (genomeVersion == null) genomeVersion = convertedMergedChromPlus.getInfo().getVersionedGenome();
		}
		if (minusStrand){
			if (convertedMergedChromMinus == null || nonConvertedMergedChromMinus == null) return false;
			if (genomeVersion == null) genomeVersion = convertedMergedChromMinus.getInfo().getVersionedGenome();
		}
		return true;
	}



	/**Region scans a chromosome collecting read count data and calculating binomial p-values.*/
	public void regionScanChromosome(){
		//fetch data
		if (fetchData() == false) {
			System.out.println("\n\tSkipping "+chromosome+". Failed to find all required PointData sets.");
			return;
		}
			
		//stranded analysis?
		if (runStrandedAnalysis){
			if (plusStrand){
				convertedChrom = convertedMergedChromPlus;
				nonConvertedChrom = nonConvertedMergedChromPlus;
			}
			else {
				convertedChrom = convertedMergedChromMinus;
				nonConvertedChrom = nonConvertedMergedChromMinus;
			}
		}
		else {
			convertedChrom = PointData.mergePairedPointDataNoSumming(convertedMergedChromPlus, convertedMergedChromMinus);
			nonConvertedChrom = PointData.mergePairedPointDataNoSumming(nonConvertedMergedChromPlus, nonConvertedMergedChromMinus);
		}

		editedBases = MethylatedBaseObservationOneSample.fetchCommonBasesWithMinimumObservations(nonConvertedChrom, convertedChrom, minimumBaseCoverage);

		//remove those 100% edited (likely snvs)
		ArrayList<MethylatedBaseObservationOneSample> good = new ArrayList<MethylatedBaseObservationOneSample>();
		for (int i=0; i< editedBases.length; i++){
			if (editedBases[i].getCon() != 0) good.add(editedBases[i]);
		}
		editedBases = new MethylatedBaseObservationOneSample[good.size()];
		good.toArray(editedBases);
		editedBasePositions = MethylatedBaseObservationOneSample.fetchPositions(editedBases);

		//calculate expect
		double edited = 0;
		double reference = 0;
		for (MethylatedBaseObservationOneSample ob: editedBases){
			edited+= ob.getNonCon();
			reference+= ob.getCon();
		}
		errorRateMinOne = 1.0 - (edited/ (edited+reference));

		//fetch regions
		regions = chromosomeRegion.get(chromosome);

		//scan for real scores
		scanRegionsReal();

		//scan permutated 
		scanPermutedRegions();

	}


	private void calculateRegionPValues(){
		//sort RandomScoreArray
		for (int x=0; x< randomScores.length; x++){
			if (randomScores[x]!= null) {
				randomScores[x].sortScores();
			}
		}
		//container for b&h correction
		Point[] pvals = new Point[numberRegionsWithData];
		int index = 0;
		//for each chromosome
		for (String chr : chromosomeRegion.keySet()){
			regions = chromosomeRegion.get(chr);
			//for each 
			for (int i=0; i< regions.length; i++){
				//scores = chi, mean, startIndex, stopIndex
				float[] scores = regions[i].getScores();
				if (scores == null) continue;
				int regionLength = (int)(scores[3]- scores[2]);

				//pval
				float mean = scores[1];
				float transPVal = randomScores[regionLength].calculateTransformedPValue(scores[0]);
				pvals[index] = new Point(index, transPVal);
				index++;

				//final scores = numberObs, pseMedian, pval
				scores = new float[]{regionLength, mean, transPVal};
				regions[i].setScores(scores);	
			}
		}

		//sort by score
		Arrays.sort(pvals, new ComparatorPointAscendingScore());
		//correct
		Point.benjaminiHochbergCorrect(pvals, 0);
		//sort back to original position
		Arrays.sort(pvals, new ComparatorPointPosition());
		//assign FDRs to pVals
		index = 0;
		//for each chromosome
		for (String chr : chromosomeRegion.keySet()){
			regions = chromosomeRegion.get(chr);
			//for each 
			for (int i=0; i< regions.length; i++){
				//scores = chi, mean, startIndex, stopIndex
				float[] scores = regions[i].getScores();
				if (scores == null) continue;
				//scores = numberObs, pseMedian, FDR
				scores[2] = pvals[index++].getScore();
			}
		}
	}

	private void scanPermutedRegions(){
		for (int i=0; i<numberRandomPermutations; i++){
			//randomize chromPointData values
			Misc.randomize(editedBases, System.currentTimeMillis());
			//scan it
			scanRegionsRandom();
		}
	}

	private void scanRegionsRandom(){
		//for each region 
		HashSet<Integer> scannedSizes = new HashSet<Integer>();
		for (int i=0; i< regions.length; i++){
			//get num obs
			//scores = chi, mean, startIndex, stopIndex
			float[] winScores = regions[i].getScores();
			//if null then no data in region
			if (winScores == null) continue;
			int regionStart = (int)winScores[2];
			int regionStop = (int) winScores[3];
			int numObs = (int)(regionStop - regionStart);
			if (scannedSizes.contains(numObs)) continue;
			scannedSizes.add(numObs);

			//check that randomScoreArray exists and isn't maxed
			if (randomScores[numObs] == null) randomScores[numObs] = new RandomScoreArray(targetNumberRandomRegions);


			//walk the dataset, this is ok since we're assuming the regions aren't overlapping
			int start =0;
			while (true){
				if (randomScores[numObs].isMaxed()) break;
				int stop = start + numObs;
				if (stop >= editedBases.length) break;
				//skip real region?
				if (start >= regionStop || stop < regionStart){ 
					long[] observed = new long[numObs *2];
					double[] expect = new double[observed.length];
					double totalObs;
					int index = 0;
					for (int j= start; j< stop; j++){
						MethylatedBaseObservationOneSample base = editedBases[j];
						observed[index] = (long)base.getNonCon();
						totalObs = observed[index];
						index++;
						observed[index] = (long)base.getCon();
						totalObs+= observed[index];
						//set errorCounts
						expect[index] = totalObs*errorRateMinOne;
						expect[index-1] = totalObs - expect[index]; 
						index++;	
					}
					//calc chiSquare
					double chi = chiSquare.chiSquare(expect, observed);
					if (randomScores[numObs].countScore(chi) == false) break;
				}
				start = stop;
			}
		}
	}

	private void scanRegionsReal(){
		//for each region 
		for (int i=0; i< regions.length; i++){
			//any data? if not return leaving scores = null;
			int[] indexes = Num.findIndexes(regions[i].getStart(), regions[i].getStop(), editedBasePositions);
			int numBases = indexes[1] - indexes[0];
			if (numBases == 0) continue;
			numberRegionsWithData++;

			//make array for chiSquare goodness of fit test
			long[] observed = new long[numBases *2];
			double[] expect = new double[observed.length];
			float[] fractions = new float[numBases];
			int index = 0;
			int counter = 0;
			double totalObs;
			for (int j= indexes[0]; j< indexes[1]; j++){
				MethylatedBaseObservationOneSample base = editedBases[j];
				fractions[counter++]= base.getFractionMethylationNoAddOne();
				observed[index] = (long)base.getNonCon();
				totalObs = observed[index];
				index++;
				observed[index] = (long)base.getCon();
				totalObs+= observed[index];
				//set errorCounts
				expect[index] = totalObs*errorRateMinOne;
				expect[index-1] = totalObs - expect[index]; 
				index++;
			}

			//calc chiSquare
			double chi = chiSquare.chiSquare(expect, observed);

			float pseMean = (float)Num.pseudoMedian(fractions);
			//scores = chi, mean, startIndex, stopIndex
			float[] scores = new float[]{(float)chi, pseMean ,  (float)indexes[0], (float)indexes[1]};
			regions[i].setScores(scores);
		}
	}


	private void loadPointDataArrays(){
		//fetch converted PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (convertedPointDirs);
		convertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		convertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
		combo = PointData.fetchStrandedPointDataNoMerge (nonConvertedPointDirs);
		nonConvertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		nonConvertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new DefinedRegionRNAEditing(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'r': convertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'e': nonConvertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'b': regionsFile = new File(args[++i]); break;
					case 't': runStrandedAnalysis = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for point directories
		if (convertedPointDirs == null || convertedPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your converted PointData directories(s)!\n");
		//only one directory look deeper
		if (convertedPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(convertedPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) convertedPointDirs = otherDirs;
		}

		//nonConverted data
		if (nonConvertedPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your nonConverted PointData directories(s)!\n");
		//only one directory look deeper
		if (nonConvertedPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(nonConvertedPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) nonConvertedPointDirs = otherDirs;
		}

		//load regions
		if (regionsFile == null) Misc.printExit("\nError: please provide a text bed file (tab delimited: chr start stop) of regions to score for RNA editing.\n");
		HashMap<String, Region[]>cr = Region.parseStartStops(regionsFile, 0, 0, 0);
		chromosomeRegion = new LinkedHashMap<String, SmoothingWindow[]>();
		for (String chr: cr.keySet()){
			Region[] regions = cr.get(chr);
			SmoothingWindow[] sw = new SmoothingWindow[regions.length];
			numberRegions+= sw.length;
			for (int i=0; i< regions.length; i++){
				sw[i] = new SmoothingWindow (regions[i].getStart(), regions[i].getStop(), null);
				int size = regions[i].getLength();
				if (size > maxSize) maxSize = size;
			}
			chromosomeRegion.put(chr, sw);
		}


	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Defined Region RNA Editing: April 2013                 **\n" +
				"**************************************************************************************\n" +
				"DRRE scores regions for the pseudomedian of the base fraction edits as well as the\n" +
				"probability that the observations occured by chance using a permutation test based on\n" +
				"the chiSquare goodness of fit statistic. \n\n" +

				"Options:\n"+
				"-b A bed file of regions to score (tab delimited: chr start stop ...)\n"+
				"-e Edited PointData directory from the RNAEditingPileUpParser.\n" +
				"       These should contain stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories. These will be merged when scanning.\n" +
				"-r Reference PointData directory from the RNAEditingPileUpParser. Ditto.\n" +
				"-t Run a stranded analysis, defaults to non-stranded.\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/DefinedRegionRNAEditing -b hg19UTRs.bed\n" +
				"-e /PointData/Edited -r /PointData/Reference \n\n" +

		"**************************************************************************************\n");

	}
}
