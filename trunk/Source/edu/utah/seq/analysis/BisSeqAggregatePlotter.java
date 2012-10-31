package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.bio.annotation.*;
import util.bio.parsers.*;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import trans.tpmap.*;

/**Generates aggregate plots/ class averages given a list of regions and point data.
 * @author nix
 * */
public class BisSeqAggregatePlotter {

	//user defined fields
	private File[] convertedPointDirs;
	private File[] nonConvertedPointDirs;
	private File bedFile;
	private int scaleSize = 0;
	private int minimumObservations = 8;
	private boolean useMean = false;

	//strand usage
	boolean invertNegativeRegions = true;

	//internal fields
	private HashMap<String,PointData[]> convertedPlusPointData;
	private HashMap<String,PointData[]> convertedMinusPointData;
	private HashMap<String,PointData[]> nonConvertedPlusPointData;
	private HashMap<String,PointData[]> nonConvertedMinusPointData;

	//by chromosome
	private String chromosome;
	private String strand;
	private PointData convertedMergedChromPlus = null;
	private PointData convertedMergedChromMinus = null;	
	private PointData nonConvertedMergedChromPlus = null;
	private PointData nonConvertedMergedChromMinus = null;
	private BisSeqRegion mergedBisSeqPoints = null;
	private HashMap<String,ArrayListStartStop[]> chromRegions;
	private HashSet<String> chromosomes = new HashSet<String>();
	private int maxSize = 0;
	private float[][] regionScores = null;
	private float[] summedRegionScores = null;

	//constructors
	public BisSeqAggregatePlotter(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//load regions
		chromRegions = ArrayListStartStop.parseArrayListStartStops(bedFile, 0, 0, true);
		if (chromRegions == null | chromRegions.size() == 0) Misc.printExit("\nProblem loading your bed file.\n");
		
		//get chromosomes without strand info
		loadChromosomes();

		//load hashes
		loadDataHashes();

		//for each region load data by chromosome
		System.out.print("Loading regions:\n\t");
		Iterator<String> it = chromosomes.iterator();

		ArrayList<String> toRemove = new ArrayList<String>();
		while (it.hasNext()){
			chromosome = it.next();
			registerData();
			
			//check datasets exist
			if (convertedMergedChromPlus == null || convertedMergedChromMinus == null || nonConvertedMergedChromPlus == null || nonConvertedMergedChromMinus == null) {
				toRemove.add(chromosome);
			}
			else {
				System.out.print(chromosome+" ");
				//load regions
				strand = ".";
				loadRegions();
				strand = "+";
				loadRegions();
				strand = "-";
				loadRegions();
			}

		}
		System.out.println();

		//remove regions that were skipped due to lack of all 4 PointData sets
		if (toRemove.size() != 0){
			System.err.println("\nWARNING: regions on the following chromosomes were ignored due to missing one or more of the four required PointData: "+toRemove);
			for (String chrStr: toRemove) {
				chromRegions.remove(chrStr+".");
				chromRegions.remove(chrStr+"-");
				chromRegions.remove(chrStr+"+");
			}
		}

		//scale?
		if (scaleSize != 0){
			System.out.println("Scaling regions to "+scaleSize+"bp");
			scaleToSize(scaleSize);
		}
		else if (sameSizes() != true) {
			System.out.println("Regions of different size found. Scaling to max size "+maxSize);
			scaleToSize(maxSize);
			scaleSize = maxSize;
		}
		else scaleSize = maxSize;

		//invert - stranded regions? only works if their regions are stranded
		if (invertNegativeRegions) invertRegions();

		//merge the regions into the mergedBisSeqPoints object
		mergeRegionPoints();

		//fetch the counts nonConPlus[0][] conPlus[1][] nonConMinus[2][] conMinus[3][]
		float[][] counts = mergedBisSeqPoints.fetchCounts(scaleSize);
		
		//print actual counts
		System.out.println("\nMerged count data:");
		System.out.println("Index\tNonConPlus\tConPlus\tNonConMinus\tConMinus");
		System.out.println(BisSeqRegion.printCounts(counts));

		//print window smoothed arrays as well as the original
		System.out.println("\nPrinting region size window fraction data:\n");
		
		float[] allNonCon = Num.sum(counts[0], counts[2]);
		float[] allCon = Num.sum(counts[1], counts[3]);

		if (useMean){
			System.out.println("Unstranded Median Window Fraction");
			printMeanSmoothedArrays(allNonCon, allCon);

			System.out.println("\nSense Median Window Fraction");
			printMeanSmoothedArrays(counts[0], counts[1]);

			System.out.println("\nAntisense Median Window Fraction");
			printMeanSmoothedArrays(counts[2], counts[3]);
		}
		else{
			System.out.println("Unstranded Summed Window Fraction");
			printSmoothedArrays(allNonCon, allCon);

			System.out.println("\nSense  Summed Window Fraction");
			printSmoothedArrays(counts[0], counts[1]);

			System.out.println("\nAntisense Summed Window Fraction");
			printSmoothedArrays(counts[2], counts[3]);
		}
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	public void loadChromosomes(){
		for (String cs: chromRegions.keySet()){
			String chr = cs.substring(0,cs.length()-1);
			chromosomes.add(chr);
		}
	}


	/**Loads data hashes without actually loading data.*/
	private void loadDataHashes(){
		//fetch converted PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (convertedPointDirs);
		convertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		convertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
		combo = PointData.fetchStrandedPointDataNoMerge (nonConvertedPointDirs);
		nonConvertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		nonConvertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
	}

	/**Fetchs the data for a particular chromosome.*/
	public void registerData(){
		
		ArrayList<PointData> al = null;
		PointData[] pd = null;
		//merge converted
		convertedMergedChromPlus = null;
		if (convertedPlusPointData.containsKey(chromosome)) {
			pd = convertedPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			convertedMergedChromPlus = PointData.mergePointData(al, false, true);
		}

		convertedMergedChromMinus = null;
		if (convertedMinusPointData.containsKey(chromosome)) {
			pd = convertedMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			convertedMergedChromMinus = PointData.mergePointData(al, false, true);
		}
		//merge nonConverted
		nonConvertedMergedChromPlus = null;
		if (nonConvertedPlusPointData.containsKey(chromosome)) {
			pd = nonConvertedPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			nonConvertedMergedChromPlus = PointData.mergePointData(al, false, true);
		}
		nonConvertedMergedChromMinus = null;
		if (nonConvertedMinusPointData.containsKey(chromosome)) {
			pd = nonConvertedMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			nonConvertedMergedChromMinus = PointData.mergePointData(al, false, true);
		}
		pd = null;
		al = null;
	}

	/**Prints smoothed arrays of 5, 10, and 20% of region size.*/
	private void printSmoothedArrays(float[] nonCon, float[] con){
		float[] fractionsbp = calculateFractions(nonCon, con);
		int fivePercent = Math.round((float)scaleSize * 0.05f);
		float[] smoothedFive = windowFractionNoShrink(nonCon, con, fivePercent);
		int tenPercent = Math.round((float)scaleSize * 0.10f);
		float[] smoothedTen = windowFractionNoShrink(nonCon, con, tenPercent);
		int twentyFivePercent = Math.round((float)scaleSize * 0.25f);
		float[] smoothedTwentyFive = windowFractionNoShrink(nonCon, con, twentyFivePercent);

		System.out.println("Index\t0%\t5%("+fivePercent+")\t10%("+tenPercent+")\t25%("+twentyFivePercent+")");
		for (int i=0; i< fractionsbp.length; i++){
			System.out.println((i+1)+"\t"+fractionsbp[i]+"\t"+smoothedFive[i]+"\t"+smoothedTen[i]+"\t"+smoothedTwentyFive[i]);
		}
	}

	/**Prints smoothed arrays of 5, 10, and 20% of region size.*/
	private void printMeanSmoothedArrays(float[] nonCon, float[] con){
		float[] fractionsbp = calculateFractions(nonCon, con);
		int fivePercent = Math.round((float)scaleSize * 0.05f);
		float[] smoothedFive = meanWindowFractionNoShrink(nonCon, con, fivePercent);
		int tenPercent = Math.round((float)scaleSize * 0.10f);
		float[] smoothedTen = meanWindowFractionNoShrink(nonCon, con, tenPercent);
		int twentyFivePercent = Math.round((float)scaleSize * 0.25f);
		float[] smoothedTwentyFive = meanWindowFractionNoShrink(nonCon, con, twentyFivePercent);

		System.out.println("Index\t0%\t5%("+fivePercent+")\t10%("+tenPercent+")\t25%("+twentyFivePercent+")");
		for (int i=0; i< fractionsbp.length; i++){
			System.out.println((i+1)+"\t"+fractionsbp[i]+"\t"+smoothedFive[i]+"\t"+smoothedTen[i]+"\t"+smoothedTwentyFive[i]);
		}
	}

	/** Slides a window along an array, one index at a time, averaging the contents. 
	 * Scans all windows including start and ends. Note, if windowSize is odd then the scan size is actually windowSize-1*/
	public float[] meanWindowFractionNoShrink(float[] nonCon, float[] con, int windowSize) {
		double[] bpFractions =  calculateFractionsWithMinObs(nonCon, con);
		//double[] smoothed = Num.windowMedianScoresNoShrinkIgnoreZeros(bpFractions, windowSize);
		double[] smoothed = Num.windowAverageScoresNoShrinkIgnoreZeros(bpFractions, windowSize);
		return Num.doubleArrayToFloatArray(smoothed);

	}

	/** Slides a window along an array, one index at a time, averaging the contents. 
	 * Scans all windows including start and ends. Note, if windowSize is odd then the scan size is actually windowSize-1*/
	public static float[] windowFractionNoShrink(float[] nonCon, float[] con, int windowSize) {
		if (windowSize == 1) return calculateFractions(nonCon, con);
		double window = windowSize;
		int halfWindow = (int)(window/2);
		float[] means = new float[nonCon.length];
		for (int i=0; i<nonCon.length; i++){
			//set start and stop
			int start = i-halfWindow;
			if (start < 0) start = 0;
			int stop = i+halfWindow;
			if (stop > nonCon.length) stop = nonCon.length;
			//sum window
			double totalNonCon = 0;
			double totalCon = 0;
			for (int j=start; j< stop; j++){
				totalNonCon += nonCon[j];
				totalCon += con[j];
			}
			double num = stop - start;
			double fraction = totalNonCon/(totalNonCon+totalCon);
			means[i] = (float)fraction;
			if (Float.isNaN(means[i])) means[i] = 0;
		}
		return means;
	}

	private static float[] calculateFractions(float[] nonCon, float[] con){
		float[] fractions = new float[nonCon.length];
		for (int i=0; i< nonCon.length; i++) {
			fractions[i] = nonCon[i]/(nonCon[i]+con[i]);
			if (Float.isNaN(fractions[i])) fractions[i] = 0;
		}
		return fractions;
	}

	private double[] calculateFractionsWithMinObs(float[] nonCon, float[] con){
		double[] fractions = new double[nonCon.length];
		for (int i=0; i< nonCon.length; i++) {
			int total = (int)(nonCon[i]+con[i]);
			if (total >= minimumObservations) fractions[i] = nonCon[i]/(nonCon[i]+con[i]);
			else fractions[i] = 0;
		}
		return fractions;
	}

	/**Sums the columns/ bp for all the regions*/
	private void sumRegionScoreColumns(){
		summedRegionScores = new float[regionScores[0].length];
		//for each bp
		for (int i=0; i< regionScores[0].length; i++){
			//for each region
			for (int j=0; j< regionScores.length; j++){
				summedRegionScores[i] += regionScores[j][i];
			}
		}
	}

	/**Loads the scores and merges.*/
	private void mergeRegionPoints(){
		Iterator<String> it = chromRegions.keySet().iterator();
		//for each chromosome
		while (it.hasNext()){

			//fetch chromosome and strand
			String chromStrand = it.next();

			ArrayListStartStop[] alss = chromRegions.get(chromStrand);
			
			//fetch first region with data?
			int i=0;
			if (mergedBisSeqPoints == null){
				for (;i<alss.length; i++){
					mergedBisSeqPoints = (BisSeqRegion)alss[i].getArrayList().get(0);
					if (mergedBisSeqPoints.dataPresent) {
						i++;
						break;
					}
				}
			}

			//for each subsequent region
			for (; i< alss.length; i++){
				BisSeqRegion pts = (BisSeqRegion)alss[i].getArrayList().get(0);
				mergedBisSeqPoints.merge(pts);
			}
		}
	}

	/**Scales the Point positions to the defined size.*/
	private void scaleToSize(int size){

		Iterator<String> it = chromRegions.keySet().iterator();
		double sizeMinOne = size -1;

		while (it.hasNext()){
			//fetch chromosome and strand
			String chromStrand = it.next();
			ArrayListStartStop[] alss = chromRegions.get(chromStrand);

			//for each region
			for (int i=0; i< alss.length; i++){
				//calculate length
				double length = alss[i].getStop() - alss[i].getStart();

				//OK length skip
				if (length == sizeMinOne) continue;

				//calculate scalar
				double scalar = sizeMinOne/length;

				//fetch Point container
				BisSeqRegion pts = (BisSeqRegion)alss[i].getArrayList().get(0);

				//scale and sum identical positions
				pts.scale(scalar);
				pts.sumIdentialPositions();
			}
		}		
	}

	/**Inverts the negative strand region data based on the scaleSize.*/
	private void invertRegions(){
		Iterator<String> it = chromRegions.keySet().iterator();
		while (it.hasNext()){
			//fetch chromosome and strand
			String chromStrand = it.next();
			if (chromStrand.endsWith("-") == false) continue;
			ArrayListStartStop[] alss = chromRegions.get(chromStrand);

			//for each region with data
			for (int i=0; i< alss.length; i++){
				//fetch Point container
				BisSeqRegion pts = (BisSeqRegion)alss[i].getArrayList().get(0);
				pts.invert(scaleSize);
			}
		}		
	}


	private boolean sameSizes(){
		Iterator<String> it = chromRegions.keySet().iterator();
		int size = -1;
		while (it.hasNext()){
			//fetch chromosome and strand
			String chromStrand = it.next();
			ArrayListStartStop[] alss = chromRegions.get(chromStrand);
			for (int i=0; i< alss.length; i++){
				int length = alss[i].getStop() - alss[i].getStart();
				if (size == -1) size = length;
				else if (size != length) return false;
			}
		}
		return true;
	}


	/**Load Point[] into each region paying attention to strand.*/
	private void loadRegions(){

		//get regions for this particular chromosome
		ArrayListStartStop[] regions = chromRegions.get(chromosome+strand);
		
		//any regions
		if (regions == null) return;

		//for each region
		for (int i=0; i< regions.length; i++){

			//fetch Point[], may be null
			int start = regions[i].getStart();
			int stop = regions[i].getStop();
			int length = stop-start;
			if (length > maxSize) maxSize = length;

			regions[i].getArrayList().add(new BisSeqRegion(start, stop, nonConvertedMergedChromPlus, convertedMergedChromPlus, nonConvertedMergedChromMinus, convertedMergedChromMinus));
		}
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BisSeqAggregatePlotter(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': convertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'n': nonConvertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'i': invertNegativeRegions = false; break;
					case 'm': useMean = true; break;
					case 'o': minimumObservations = Integer.parseInt(args[++i]); useMean = true; break;
					case 's': scaleSize = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for bed file
		if (bedFile == null) Misc.printExit("\nPlease enter a regions file to use in scoring regions.\n");

		//print params
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " "));

		System.out.println();

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Bis Seq Aggregate Plotter: October 2012                    **\n" +
				"**************************************************************************************\n" +
				"BSAP merges bisulfite data over equally sized regions to generate data for class\n" +
				"average agreggate plots of fraction methylation.  A smoothing window is also applied.\n" +
				"Data for unstranded, sense, and antisense are produced.\n\n" +

				"Options:\n"+
				"-c Converted PointData directories, full path, comma delimited. These should\n" +
				"       contain stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories. See the NovoalignBisulfiteParser app.\n" +
				"-n Non-converted PointData directories, ditto. \n" +
				"-b Bed file (tab delim: chr start stop name score strand(+/-/.)), full path.\n"+
				"-i Don't invert - stranded regions, defaults to inverting.\n"+
				"-s Scale all regions to a particular size. Defaults to scaling to max region size.\n"+
				"-m Calculate individual base fractions and then take a mean, ignoring zeros, over\n" +
				"       the window, instead of summing the obs in the window and taking the fraction.\n"+
				"-o Minimum number of observations before scoring base fraction methylation, defaults\n"+
				"       to 8.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/BisSeqAgregatePlotter -c\n" +
				"      /NBP/Con -n /NBP/NonCon -b /Anno/tssSites.bed -m\n\n" +

		"**************************************************************************************\n");

	}
}
