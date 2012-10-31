package trans.roc;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.useq.data.Region;

import trans.anno.*;
import trans.misc.*;
import util.bio.parsers.*;
import util.bio.seq.Seq;
import util.gen.*;
import util.bio.annotation.*;

/**
 * Returns values that overlap particular regions, can generate p-values associated with random back ground model
 */
public class ScoreParsedBars {
	private Positive[] regions;
	private File barDirectory;
	private HashMap chromBarFileMap;
	private String currentChrom = "";
	private int[] positions;
	private float[] scores;
	private int numberRandom = 1000;
	private boolean makeRandom = false;
	private HashMap chromGCFileMap;
	private File gcGenomeDir;
	private boolean[] currentGCContent;
	private double fractionGCTolerance = 0.1;
	private double fractionScoreTolerance = 0.2;
	private StringBuffer randomResults = new StringBuffer();
	private int bpPositionOffSetBar = 0;
	private int bpPositionOffSetRegion = 0;
	private HashMap<String, Region[]> interrogatedRegions;
	private Region[] currentInterrogatedRegions;
	private int totalBPInterrogatedRegions;
	private int[] startsBPInterrogatedRegions;
	private boolean unlog = false;
	private boolean printBaseScoresToScreen = true;
	
	public ScoreParsedBars(File barDirectory, int bpPositionOffSetBar, int bpPositionOffSetRegion){
		this.barDirectory = barDirectory;
		makeBarFileHashMap();
		this.bpPositionOffSetBar = bpPositionOffSetBar;
		this.bpPositionOffSetRegion = bpPositionOffSetRegion;
	}
	
	public void removeRegion(int index){
		ArrayList<Positive> subAL = new ArrayList<Positive>();
		for (int i=0; i< regions.length; i++){
			if (i != index) subAL.add(regions[i]);
		}
		regions = new Positive[regions.length-1];
		subAL.toArray(regions);		
	}


	public ScoreParsedBars(String[] args){
		processArgs(args);

		//for each region
		if (printBaseScoresToScreen) System.out.println("\nIndividual hits to each region with their scores:");
		ArrayList allScores = new ArrayList(regions.length*4);
		for (int x=0; x< regions.length; x++){
			//load new bar file?
			if (regions[x].getChromosome().equals( currentChrom) == false) {
				currentChrom = regions[x].getChromosome();
				loadChromosomeBarFile();
				//making random?
				if (makeRandom){
					//load gc
					System.out.println("\tLoading serialized GC boolean file...");
					Object obj = chromGCFileMap.get(currentChrom);
					if (obj == null ) Misc.printExit("\nError: cannot find a gc binary file for "+currentChrom);
					currentGCContent = (boolean[])IO.fetchObject((File)obj);
					//load interrogatedRegions
					currentInterrogatedRegions = interrogatedRegions.get(currentChrom);
					if (currentInterrogatedRegions == null) Misc.printExit("\nError: cannot find interrogated regions for "+currentChrom);
					totalBPInterrogatedRegions = Region.totalBP(currentInterrogatedRegions);
					startsBPInterrogatedRegions = Region.startsInBases(currentInterrogatedRegions);
				}
			}
			//load region with scores
			loadRegion(regions[x]);
			
			//any scores
			if (regions[x].getScores().size() == 0){
				System.out.println("\tSkipping region "+regions[x].toStringSimple()+" no associated scores.\n");
				randomResults.append(regions[x].toStringSimple()+"\tSkipped, no scores\n");
				removeRegion(x);
				x--;
				continue;
			}
			//print region
			System.out.print(regions[x].toStringSimple()+"\t");
			//print scores, calculate mean
			double mean = 0;
			ArrayList grsAL = regions[x].getScores();
			double numberScores = grsAL.size();
			if (numberScores !=0){
				//double[] forMedian = new double[(int)numberScores];
				for (int i=0; i< numberScores; i++){
					Gr gr = (Gr)grsAL.get(i);
					if (printBaseScoresToScreen) {
						//System.out.println("\t"+gr.toString());
						System.out.print(gr.getScore()+" ");
					}
					double antiLog;
					if (unlog) antiLog = Num.antiLog(gr.getScore(), 2);
					else antiLog = gr.getScore();
					Double antiLogScore = new Double (antiLog);
					allScores.add(antiLogScore);
					mean += antiLog;
					//forMedian[i] = antiLog;
				}
				mean = mean/numberScores;
				//Arrays.sort(forMedian);
				//mean = Num.median(forMedian);
			}
			System.out.println();
			//set mean as first score in scores
			grsAL.clear();
			grsAL.add(new Double(mean));

			//make random matched regions?
			if (makeRandom){
				makeRandom(regions[x], numberScores);
				//calculate p-value up and down. 
				//first value is the real all subsequent are random
				double[] scores = Num.arrayListOfDoubleToArray(regions[x].getScores());
				
				double real = scores[0];
				//count number > and <
				double numGreaterThan = 0;
				double numLessThan = 0;
				double average =0;
				for (int i=1; i<scores.length; i++){
					average += scores[i];
					if (scores[i] > real) numGreaterThan++;
					else if (scores[i] < real ) numLessThan++;
					//same score
					else{
						numGreaterThan++;
						numLessThan++;
					}
				}
				
				//calc ave, and pvalues
				average = average/ ((double)numberRandom);
				numGreaterThan = numGreaterThan/ ((double)numberRandom);
				numLessThan = numLessThan/ ((double)numberRandom);
				double foldEnrich = real/average;
				//save results in StringBuffer
				//coord, real score, ave score, fold enrich, pvalUp, pvalDown
				randomResults.append(regions[x].toStringSimple()+"\t"+real+"\t"+average+"\t"+foldEnrich+"\t"+numGreaterThan+"\t"+numLessThan+"\n");
			}
		}

		//any regions?
		if (regions.length == 0) Misc.printExit("\nNo remaining regions? Aborting!");
		
		//print stats on all scores
		double[] scores = Num.arrayListOfDoubleToArray(allScores);
		float[] scoresF = new float[scores.length];
		for (int k=0; k< scores.length; k++) scoresF[k] = new Double (scores[k]).floatValue();
		if (unlog) System.out.println("Stats on all antiLog2ed scores:");
		else System.out.println("Stats on all scores:");		
		Num.statFloatArray(scoresF, false);
		
		if (makeRandom){
			System.out.println("\nCoord\tReal_Score\tMean_Rnd_Score\tFold_enrichment\tPvalUp\tPvalDown\n"+randomResults);
			//calc summary stats on list
			double[] scoreTotals = new double[numberRandom +1];
			for (int r=0; r< regions.length; r++){
				//first value is the real all subsequent are random
				double[] s = Num.arrayListOfDoubleToArray(regions[r].getScores());			
				for (int z=0; z< s.length; z++) scoreTotals [z] += s[z];
			}			
			double pval = 0;
			double totalRandomScores = 0;
			//calc averages
			double numRegions = regions.length;
			double aveRealScore = scoreTotals[0]/numRegions;
			for (int r=1; r< scoreTotals.length; r++) {
				totalRandomScores += scoreTotals[r];
				scoreTotals[r] = scoreTotals[r]/numRegions; 
				if (scoreTotals[r] >= aveRealScore) pval++;
			}
			double aveRandomScore = totalRandomScores/(numRegions * (double)numberRandom);
			System.out.println("\nMean real scores: "+aveRealScore);
			System.out.println("Mean random scores: "+aveRandomScore);
			System.out.println("Fold over random: "+aveRealScore/aveRandomScore);
			System.out.println("PVal: "+pval+"/"+numberRandom);
		}
		else {
			//print means
			System.out.println("Individual region summaries:\n\nCoordinates\tMean");
			//System.out.println("Individual region summaries:\n\nCoordinates\tlog2(Median)");
			int numGreaterThan1 = 0;
			int numLessThan1=0;
			for (int x=0; x< regions.length; x++){
				ArrayList grsAL = regions[x].getScores();
				double mean = ((Double)grsAL.get(0)).doubleValue();
				double log2;
				if (unlog) log2 = Num.log2(mean);
				else log2 = mean;
				StringBuilder scoresSB = new StringBuilder();
				for (int a=1; a< grsAL.size(); a++){
					scoresSB.append(((Double)grsAL.get(a)).toString());
					scoresSB.append(" ");
				}
				System.out.println(regions[x].toStringSimple()+"\t"+log2+"\t"+scoresSB);
				if (log2 >=1 )numGreaterThan1 ++;
				else if (log2 <= -1)numLessThan1 ++;
			}
			System.out.println("Total\t>1\t<-1");
			System.out.println(regions.length+"\t"+numGreaterThan1+"\t"+numLessThan1);
		}
		
		System.out.println("\nDone!\n");
	}
	
	/**Takes a Positive and finds 1000 chrom, length, gc, and number scores matched random regions.
	 * Calculates each mean and sets in the Positive.*/
	public int makeRandom(Positive region, double numberScores){
		
		//calculate gc content of real region and set min max
		double realGC = calculateFractionGCContent(region.getStart(), region.getStop());		
		int sizeRealRegionMinOne = region.getStop() - region.getStart();
		double minGC = realGC - fractionGCTolerance;
		double maxGC = realGC + fractionGCTolerance;
		//System.out.println("\t\tGC real "+ realGC+" minGC "+minGC+" maxGC "+maxGC);
		
		//calculate min max number of scores
		double diff = numberScores * fractionScoreTolerance;
		int minScores = (int)Math.round(numberScores - diff);
		//int maxScores = (int)Math.round(numberScores + diff);
		//System.out.println("\t\tNum scores "+numberScores+" min Scores "+minScores+" max Scores"+maxScores);
		
		//try 10000 times to find a gc and num oligo matched random region
		Random randomGenerator = new Random();
		int numberRandomFound = 0;
		ArrayList regionScores = region.getScores();
		System.out.print("\tPicking random regions ");		
		for (int x=0; x<numberRandom; x++){
			System.out.print(".");
			for (int y=0; y< 100000; y++){
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
				if (lengthIntReg < sizeRealRegionMinOne) {
					//System.out.println("Fail too small "+sizeRealRegionMinOne+" min "+lengthIntReg);
					continue;
				}
				//generate a random start stop matching the length of the region
				int start = randomGenerator.nextInt(lengthIntReg) + interrogatedRegion.getStart();
				//check to see if start is < 0
				if (start < 0) continue;
				int stop = start + sizeRealRegionMinOne;
				//check to see if stop goes past the stop of the interrogated region
				if (stop >= interrogatedRegion.getStop()){
					//System.out.println("Fail size stop greater than stop "+stop+" IR stop "+interrogatedRegion.getEnd());
					continue;
				}

				//check gc
				double testGC = calculateFractionGCContent(start, stop);
				if (testGC < minGC || testGC > maxGC) {
					//System.out.println("Fail GC "+testGC+" min "+minGC+" max "+maxGC);
					continue;
				}
				
				//check num scores			
				Positive testRegion = new Positive(currentChrom, start, stop);
				//load region with scores
				loadRegion(testRegion);
				//correct number?
				ArrayList testAL = testRegion.getScores();
				int numTestScores = testAL.size();				
				if (numTestScores < minScores ) { //|| numTestScores > maxScores) {	
					//System.out.println("Fail scores "+numTestScores+" min "+minScores);
					continue;
				}
				//otherwise calc mean and add!
				numberRandomFound++;
				double mean = 0;
				for (int i=0; i< numTestScores; i++){
					Gr gr = (Gr)testAL.get(i);
					if (unlog) mean += Num.antiLog(gr.getScore(), 2);
					else mean += gr.getScore();
				}
				mean = mean/(double)numTestScores;
				regionScores.add(new Double(mean));
				break;
			}
		}
		System.out.println("\n");
		if (numberRandomFound != numberRandom) {
			System.out.println("\nError: could not generate the requested number of random regions, found "+numberRandomFound+ " for "+region.toStringSimple());
		}
		return numberRandomFound;
	}
	

	/**Calculates the gc content.*/
	public double calculateFractionGCContent(int start, int stop){
		double ave = 0;
		int realStop = stop +1;
		for (int i = start; i< realStop; i++){
			if (currentGCContent[i]) ave++;
		}
		double length = realStop - start;
		return ave/ length;
	}
	
	public void loadChromosomeBarFile(){
		Object obj = chromBarFileMap.get(currentChrom);
		if (obj == null) Misc.printExit("\nError: cannot find a chromosome bar file for "+currentChrom + "\n\t"+chromBarFileMap);
		File chromBarFile = (File)obj;
		GrGraph currentGrGraph = Util.readSimpleGrBarFile(chromBarFile);
		//use Gr[] to sort
		Gr[] grs = Gr.parseGrGraph(currentGrGraph);
		Arrays.sort(grs, new GrComparator());
		//break out positions and scores
		positions = new int[grs.length];
		scores = new float[grs.length];
		for (int i=0; i< grs.length; i++){
			positions[i] = grs[i].getPosition() + bpPositionOffSetBar;
			scores[i] = grs[i].getScore();
		}
		grs = null;
		currentGrGraph = null;
	}
	
	/**Checks to see if an sgr object is in one of the positives, if not returns -1, 
	 * otherwise returns index number.*/
	public static int inRegion(Sgr sgr, Positive[] pos){
		int numPos = pos.length;
		for (int i=0; i< numPos; i++){
			if (pos[i].matches(sgr)) return i;
		}
		return -1;
	}
	
	/**Loads a Positive with Gr objects.*/
	public void loadRegion(Positive region){
		//use fast lookup to find closest
		int closestIndex = Num.findClosestIndexToValue(positions, region.getStart()) - 1;
		if (closestIndex < 0) closestIndex = 0;
		//run through positions
		boolean inside = false;
		ArrayList grs = new ArrayList();
		for (int i=closestIndex; i< positions.length; i++){
			if (region.contains(positions[i])) {
				grs.add(new Gr(positions[i], scores[i]));
				inside = true;
			}
			else if (inside) break;
		}
		region.setScores(grs);
	}
	
	/**Fetches the Grs under the region,
	 * Uses the bpPositionOffSets so set appropriately!*/
	public ArrayList fetchRegionGrs(String chr, int start, int stop){
		Positive region = new Positive(chr, start, stop + bpPositionOffSetRegion);
		//can it hold an oligo?
		if (region.getStop()- region.getStart() < 1) return null;
		//load new bar file?
		if (region.getChromosome().equals( currentChrom) == false) {
			currentChrom = region.getChromosome();
			loadChromosomeBarFile();
		}
		//load region with scores
		loadRegion(region);
		//save scores
		return region.getScores();
	}
	
	/**Calculates the median intensity for the oligos that lie under a gene.*/
	public void calculateMedians(UCSCGeneLine[] lines){
		ArrayList allGrs = new ArrayList();
		for (int x=0; x< lines.length; x++){
			//get exons
			ExonIntron[] exons = lines[x].getExons();
			allGrs.clear();
			//for each exon
			for (int y=0; y< exons.length; y++){
				Positive region = new Positive(lines[x].getChrom(), exons[y].getStart(), exons[y].getEnd() + bpPositionOffSetRegion);
				//can it hold an oligo?
				if (region.getStop()- region.getStart() < 1) continue;
				//load new bar file?
				if (region.getChromosome().equals( currentChrom) == false) {
					currentChrom = region.getChromosome();
					loadChromosomeBarFile();
				}
				//load region with scores
				loadRegion(region);
				//save scores
				allGrs.addAll(region.getScores());
			}

			//calculate median
			float median = 0;
			int numberScores = allGrs.size();
			if (numberScores !=0){
				double[] values = new double[numberScores];
				for (int i=0; i< numberScores; i++){
					Gr gr = (Gr)allGrs.get(i);
					if (unlog) values[i] = (float)Num.antiLog(gr.getScore(), 2);
					else values[i] = gr.getScore();
				}
				Arrays.sort(values);
				median = (float)Num.median(values);
			}
			if (unlog) lines[x].setScores(new float[]{Num.log2(median)});
			else lines[x].setScores(new float[]{median});
		}
	}
	
	public static void main(String[] args) {
		if (args.length==0){
			printDocs();
			System.exit(0);
		}	
		new ScoreParsedBars(args);
	}		
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File interrogatedRegionsFile = null;
		barDirectory = null;
		File regionFile = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': barDirectory = new File(args[i+1]); i++; break;
					case 'r': regionFile = new File(args[i+1]); i++; break;
					case 'g': gcGenomeDir = new File(args[i+1]); i++; makeRandom = true; break;
					case 'n': numberRandom = Integer.parseInt(args[i+1]); i++; break;
					case 'o': bpPositionOffSetBar = Integer.parseInt(args[i+1]); i++; break;
					case 's': bpPositionOffSetRegion = Integer.parseInt(args[i+1]); i++; break;
					case 'u': unlog = true; break;
					case 'd': printBaseScoresToScreen = false; break;
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
		if (barDirectory == null || barDirectory.isDirectory()== false || regionFile == null || regionFile.exists() == false) Misc.printExit("\nError: Please set the -b and or -r parameters.\n");
		
		//load and sort Regions File
		System.out.println("\nLoading region file "+regionFile.getName());
		regions = ParseSgrsForParticularRegions.parseRegionFile(regionFile);
		if (bpPositionOffSetRegion !=0 ) regions = Positive.filter(bpPositionOffSetRegion, regions);
		Arrays.sort(regions, new PositiveComparator());
		//System.out.println("Pos "+regions[0].toStringSimple());
		
		if (bpPositionOffSetBar !=0 )System.out.println("BP Offset Bar Position= "+bpPositionOffSetBar);
		if (bpPositionOffSetRegion !=0 )System.out.println("BP Offset GenomicRegion= "+bpPositionOffSetRegion);
		
		//make hash of chromosome text and barFile
		makeBarFileHashMap();
		
		//make hash of chromosome text and gc boolean file?
		if (makeRandom){
			File[] gcFiles = IO.extractFiles(gcGenomeDir, ".gc");
			chromGCFileMap = Seq.makeChromosomeNameFileHash( gcFiles );
			if (chromGCFileMap == null) Misc.printExit("\nError parsing xxx.gc files.\n");
			if (interrogatedRegionsFile == null) Misc.printExit("\nError, cannot find your interrogated regions file!\n");
			interrogatedRegions = Region.parseStartStops(interrogatedRegionsFile, 0, 0, 0);
		}
	}	
	
	public void makeBarFileHashMap(){
		File[] barFiles = IO.extractFiles(barDirectory, ".bar");
		chromBarFileMap = new HashMap();
		for (int i=0; i< barFiles.length; i++){
			String truncName = barFiles[i].getName().replaceFirst(".bar","");
			chromBarFileMap.put(truncName, barFiles[i]);
		}
		if (chromBarFileMap == null) Misc.printExit("\nError parsing xxx.bar files.\n");
	}
	
	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           ScoreParsedBars: Sept 2008                             **\n" +
				"**************************************************************************************\n" +
				"For each region finds the underlying scores from the chromosome specific bar files.\n" +
				"Prints the scores as well as their mean . A p-value for each region's score can be\n" +
				"calculated using chromosome, interrogated region, length, # scores, and gc matched\n" +
				"random regions. Be sure to set the -u flag if your scores are log2 values.\n\n" +
				
				"-r Full path file text for your region file (tab delimited: chr start stop(inclusive)).\n" +
				"-b Full path directory text for the chromosome specific data xxx.bar files.\n" +
				"-o Bp offset to add to the position coordinates, defaults to 0.\n"+
				"-s Bp offset to add to the stop of each region, defaults to 0.\n"+
				"-u Unlog the bar values, set this flag if your scores are log2 transformed.\n"+
				"-g Estimate a p-value for the score associated with each region. Provide a full path\n"+
				"         directory text for chromosome specific gc content boolean arrays. See\n"+
				"         ConvertFasta2GCBoolean app. Complete option -i\n"+
				"-i If estimating p-values, provide a full path file text containing the interrogated\n"+
				"         regions (tab delimited: chr start stop ...) to use in drawing random regions.\n"+
				"-n Number of random region sets, defaults to 1000.\n" +
				"-d Don't print individual scores to screen.\n"+
				
				"\nExample: java -jar pathTo/Apps/ScoreParsedBars -b /BarFiles/Oligos/\n" +
				"       -r /Res/miRNARegions.bed -o -30 -s -60 -i /Res/interrRegions.bed\n" +
				"       -g /Genomes/Hg18/GCBooleans/\n\n" +
				
		"**************************************************************************************\n");		
	}
	
}
