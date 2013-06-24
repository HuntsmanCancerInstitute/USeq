package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import edu.utah.seq.useq.apps.Bar2USeq;
import trans.tpmap.*;

/**Takes PointData from the RNAEditingPileupParser and generates sliding window scan statistics for identifying clustered edited regions. Estimates PVals using random permutation.
 * @author Nix
 * */
public class RNAEditingScanSeqs {

	//user defined fields
	private File[] convertedPointDirs;
	private File[] nonConvertedPointDirs;
	private File saveDirectory;
	private int windowSize = 50;
	private int minimumNumberObservationsInWindow = 3;
	private double errorRateMinOne;
	private int minimumBaseCoverage = 5;
	private float minimumPseMedian = 0.005f;
	private float minimumBaseFractionEdited = 0.01f;
	private boolean runStrandedAnalysis = false;

	//internal fields
	private int numberRandomPermutations = 100;
	private int targetNumberRandomWindows = 1000000;
	private ChiSquareTest chiSquare = new ChiSquareTest();
	private File meanDir;
	private File pvalDir;
	private boolean plusStrand = true;
	private boolean minusStrand =  true;
	private HashMap<String,PointData[]> convertedPlusPointData;
	private HashMap<String,PointData[]> convertedMinusPointData;
	private HashMap<String,PointData[]> nonConvertedPlusPointData;
	private HashMap<String,PointData[]> nonConvertedMinusPointData;
	private ArrayList<SmoothingWindowInfo> smiAL = new ArrayList<SmoothingWindowInfo>();
	private RandomScoreArray[] randomScores;

	private WindowMaker windowMaker; 
	private int[][] windows;
	private SmoothingWindow[] smoothingWindow;
	private SmoothingWindowInfo[] smoothingWindowInfo;
	private File swiFile;
	private String[] scoreNames;
	private String[] scoreDescriptions;
	private String[] scoreUnits;

	//by chromosome
	private String chromosome;
	private PointData convertedMergedChromPlus = null;
	private PointData convertedMergedChromMinus = null;	
	private PointData nonConvertedMergedChromPlus = null;
	private PointData nonConvertedMergedChromMinus = null;	
	private PointData convertedChrom = null;
	private PointData nonConvertedChrom = null;

	MethylatedBaseObservationOneSample[] methylatedBases;
	int[] methylatedBasePositions;

	//constructor
	public RNAEditingScanSeqs(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);
		
		//set initial score items
		setScoreStrings(false);
		
		//make containers for randomScores
		randomScores = new RandomScoreArray[windowSize+10];
		
		//load data pointers
		loadPointDataArrays();

		//stranded analysis?
		if (runStrandedAnalysis){
			meanDir = new File(saveDirectory, "MeanPlus");
			meanDir.mkdir();
			pvalDir = new File(saveDirectory, "PValPlus");
			pvalDir.mkdir();
			plusStrand = true;
			minusStrand =  false;
			doWork();
			//reset
			setScoreStrings(false);
			randomScores = new RandomScoreArray[windowSize+10];
			smiAL.clear();
			meanDir = new File(saveDirectory, "MeanMinus");
			meanDir.mkdir();
			pvalDir = new File(saveDirectory, "PValMinus");
			pvalDir.mkdir();
			plusStrand = false;
			minusStrand =  true;
			doWork();
		}
		//combine analysis
		else {
			meanDir = new File(saveDirectory, "Mean");
			meanDir.mkdir();
			pvalDir = new File(saveDirectory, "PVal");
			pvalDir.mkdir();
			doWork();
		}


		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	//methods
	
	public void doWork(){

		//for each chromosome
		System.out.print("Scanning Chromosomes");
		
		String[] chromosomes = fetchAllChromosomes();
			for (int i=0; i< chromosomes.length; i++){
				chromosome = chromosomes[i];
				windowScanChromosome();
			}
		
		//any data, make arrays
		if (smiAL.size() == 0) Misc.printExit("\nNo windows found. Thresholds too strict?\n");
		smoothingWindowInfo = new SmoothingWindowInfo[smiAL.size()];
		smiAL.toArray(smoothingWindowInfo);

		//calculate p-values
		calculateWindowPValues();
		
		//write out graph data
		System.out.println("\nSaving graph data...");
		writeBarFileGraphs();
		
		//convert graph data to useq
		new Bar2USeq(meanDir, true);
		new Bar2USeq(pvalDir, true);
		
		//multiple test correct the pvalues
		System.out.println("Converting pvalues to FDRs...");
		convertPValuesToFDRs();
		setScoreStrings(true);

		//save window data 
		System.out.println("Writing window objects for the EnrichedRegionMaker...");
		String strand = "";
		if (runStrandedAnalysis){
			if (plusStrand) strand = "Plus";
			else strand = "Minus";
		}
		swiFile = new File (saveDirectory, "windowData"+windowSize+"bp"+strand+".swi");
		IO.saveObject(swiFile, smoothingWindowInfo);
	}

	private void convertPValuesToFDRs() {
		//remove low scoring windows
		//for each chromosome
		int numGoodWin = 0;
		for (int i=0; i< smoothingWindowInfo.length; i++){
			smoothingWindow = smoothingWindowInfo[i].getSm();
			//for each window
			ArrayList<SmoothingWindow> goodWin = new ArrayList<SmoothingWindow>((int)(((double)smoothingWindow.length) * 0.1));
			for (int j=0; j< smoothingWindow.length; j++){
				//scores = numberObs, pseMedian, pval
				float[] scores = smoothingWindow[j].getScores();
				if (scores[1]>= minimumPseMedian) goodWin.add(smoothingWindow[j]);
			}
			smoothingWindow = new SmoothingWindow[goodWin.size()];
			goodWin.toArray(smoothingWindow);
			smoothingWindowInfo[i].setSm(smoothingWindow);
			numGoodWin += smoothingWindow.length;
		}
		//collect pvals for B&H corr
		Point[] pvals = new Point[numGoodWin];
		int index = 0;
		for (int i=0; i< smoothingWindowInfo.length; i++){
			smoothingWindow = smoothingWindowInfo[i].getSm();
			//for each window
			for (int j=0; j< smoothingWindow.length; j++){
				//scores = numberObs, pseMedian, pval
				float[] scores = smoothingWindow[j].getScores();
				pvals[index] = new Point(index, scores[2]);
				index++;
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
		for (int i=0; i< smoothingWindowInfo.length; i++){
			smoothingWindow = smoothingWindowInfo[i].getSm();
			//for each window
			for (int j=0; j< smoothingWindow.length; j++){
				//scores = numberObs, pseMedian, pval
				float[] scores = smoothingWindow[j].getScores();
				scores[2] = pvals[index++].getScore();
			}
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

	
	/**Fetches the names of all the chromosomes in the data*/
	public String[] fetchAllChromosomes(){
		HashSet<String> c = new HashSet<String>();
		Iterator<String> it;
		if (plusStrand){
			it = convertedPlusPointData.keySet().iterator();
			while (it.hasNext()) c.add(it.next());
			it = nonConvertedPlusPointData.keySet().iterator();
			while (it.hasNext()) c.add(it.next());
		}
		if (minusStrand){
			it = convertedMinusPointData.keySet().iterator();
			while (it.hasNext()) c.add(it.next());
			it = nonConvertedMinusPointData.keySet().iterator();
			while (it.hasNext()) c.add(it.next());
		}
		return Misc.hashSetToStringArray(c);
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
		}
		if (minusStrand){
			if (convertedMergedChromMinus == null || nonConvertedMergedChromMinus == null) return false;
		}
		return true;
	}



	/**Window scans a chromosome collecting read count data and calculating binomial p-values.*/
	public void windowScanChromosome(){
		//fetch data
		if (fetchData() == false) {
			System.out.println("\n\tSkipping "+chromosome+". Failed to find required PointData sets.");
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
		methylatedBases = MethylatedBaseObservationOneSample.fetchCommonBasesWithMinimumObservations(nonConvertedChrom, convertedChrom, minimumBaseCoverage);
		
		//remove those 100% edited (likely snvs) and fail the minimumBaseFractionEdited
		ArrayList<MethylatedBaseObservationOneSample> good = new ArrayList<MethylatedBaseObservationOneSample>();
		for (int i=0; i< methylatedBases.length; i++){
			float bfe = methylatedBases[i].getFractionMethylationNoAddOne();
			if (bfe != 1.0f && bfe >= minimumBaseFractionEdited) good.add(methylatedBases[i]);
		}
		methylatedBases = new MethylatedBaseObservationOneSample[good.size()];
		good.toArray(methylatedBases);
		methylatedBasePositions = MethylatedBaseObservationOneSample.fetchPositions(methylatedBases);
		
		//calculate expect
		double edited = 0;
		double reference = 0;
		for (MethylatedBaseObservationOneSample ob: methylatedBases){
			edited+= ob.getNonCon();
			reference+= ob.getCon();
		}
		errorRateMinOne = 1.0 - (edited/ (edited+reference));
		
		//make windows using all of the reads
		makeWindows(methylatedBasePositions);
		
		//any windows?
		if (windows.length == 0){
			System.out.println("\n\tSkipping "+chromosome+". No windows found with minimum observations of "+minimumNumberObservationsInWindow+" within a window size of "+windowSize);
			return;
		}
		System.out.print(".");
		//make SmoothingWindow[] container
		smoothingWindow = new SmoothingWindow[windows.length];
		
		//scan for real scores
		scanWindowsReal();
		
		//scan permutated 
		scanPermutedWindows();
		
		//save results
		addWindowDataToSMIArrayList();
	}

	/**Saves window data to SMI Array.*/
	public void addWindowDataToSMIArrayList(){
		Info info = convertedChrom.getInfo();
		HashMap<String,String> notes = new HashMap<String,String>();
		notes.put(BarParser.WINDOW_SIZE, windowSize+"");
		notes.put(BarParser.DESCRIPTION_TAG, Misc.stringArrayToString(scoreNames, ","));
		notes.put(BarParser.UNIT_TAG, Misc.stringArrayToString(scoreUnits, ","));
		info.setNotes(notes);	
		if (plusStrand && minusStrand == false) info.setStrand("+");
		else if (minusStrand && plusStrand == false) info.setStrand("-");
		else info.setStrand(".");
		smiAL.add(new SmoothingWindowInfo(smoothingWindow, info));
	}

	private void calculateWindowPValues(){
		//sort RandomScoreArray
		for (int x=0; x< randomScores.length; x++){
			if (randomScores[x]!= null) {
				randomScores[x].sortScores();
			}
		}

		//for each chromosome
		for (int i=0; i< smoothingWindowInfo.length; i++){
			smoothingWindow = smoothingWindowInfo[i].getSm();
			//for each window
			for (int j=0; j< smoothingWindow.length; j++){
				//scores = chi, mean, startIndex, stopIndex
				float[] scores = smoothingWindow[j].getScores();
				int windowLength = (int)(scores[3]- scores[2]);
				
				//pval
				float mean = scores[1];
				float transPVal = randomScores[windowLength].calculateTransformedPValue(scores[0]);
				
				//final scores = numberObs, pseMedian, pval
				scores = new float[]{windowLength, mean, transPVal};
				smoothingWindow[j].setScores(scores);	
			}
		}
	}

	private void scanPermutedWindows(){
		for (int i=0; i<numberRandomPermutations; i++){
			//randomize chromPointData values
			Misc.randomize(methylatedBases, System.currentTimeMillis());
			//scan it
			scanWindowsRandom();
		}
	}
	
	private void scanWindowsRandom(){
		//calculate random pseudomedians?
		//float[] fractions = null;
		//if (randomPseMedianIndex < randomPseMedians.length) fractions = new float[0];
		//for each window 
		for (int i=0; i< windows.length; i++){
			//get start stop
			//scores = chi, mean, startIndex, stopIndex
			float[] winScores = smoothingWindow[i].getScores();
			int startIndex = (int) winScores[2];
			int stopIndex = (int) winScores[3];
			
			//num bases
			int numBases = stopIndex - startIndex;
			//check that randomScoreArray exists and isn't maxed
			if (randomScores[numBases] == null) randomScores[numBases] = new RandomScoreArray(targetNumberRandomWindows);
			if (randomScores[numBases].isMaxed()) continue;
			
			long[] observed = new long[numBases *2];
			double[] expect = new double[observed.length];
			//boolean addFraction = (fractions!= null && (randomPseMedianIndex < randomPseMedians.length));
			//if (addFraction) fractions = new float[numBases];
			//int counter = 0;
			double totalObs;
			int index = 0;
			for (int j= startIndex; j< stopIndex; j++){
				MethylatedBaseObservationOneSample base = methylatedBases[j];
				//if (addFraction) fractions[counter++]= base.getFractionMethylationNoAddOne();
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
			//calc pseMedian
			//if (addFraction) randomPseMedians[randomPseMedianIndex++] = (float)Num.pseudoMedian(fractions);
			
			//calc chiSquare stat
			double chi = chiSquare.chiSquare(expect, observed);
			randomScores[numBases].countScore(chi);
		}
	}

	private void scanWindowsReal(){
		
		//for each window 
		for (int i=0; i< windows.length; i++){
			int[] indexes = Num.findIndexes(windows[i][0], windows[i][1], methylatedBasePositions);
			//make array for chiSquare goodness of fit test
			int numBases = indexes[1] - indexes[0];
			long[] observed = new long[numBases *2];
			double[] expect = new double[observed.length];
			float[] fractions = new float[numBases];
			int index = 0;
			int counter = 0;
			double totalObs;
			for (int j= indexes[0]; j< indexes[1]; j++){
				MethylatedBaseObservationOneSample base = methylatedBases[j];
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
			
			//calc chiSquare stat
			double chi = chiSquare.chiSquare(expect, observed);
			
			float pseMean = (float)Num.pseudoMedian(fractions);
			//scores = chi, mean, startIndex, stopIndex
			float[] scores = new float[]{(float)chi, pseMean ,  (float)indexes[0], (float)indexes[1]};
			
			//make window
			smoothingWindow[i] = new SmoothingWindow (windows[i][0], windows[i][1], scores);
		}
	}




	
	private void loadPointDataArrays(){
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (convertedPointDirs);
		convertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		convertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
		combo = PointData.fetchStrandedPointDataNoMerge (nonConvertedPointDirs);
		nonConvertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		nonConvertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
	}

	/**Makes a common set of windows using merged positions.*/
	public void makeWindows(int[] positions){
		windows = windowMaker.makeWindows(positions);
		//assign bp positions
		for (int i=0; i< windows.length; i++){
			windows[i][0] = positions[windows[i][0]];
			windows[i][1] = positions[windows[i][1]]+1;	//last base isn't included
		}
	}

	/**Writes stair step window bar graph files*/
	public void writeBarFileGraphs(){

			//for each chromosome
			for (int i=0; i< smoothingWindowInfo.length; i++){
				Info info = smoothingWindowInfo[i].getInfo().copy();
				SmoothingWindow[] sm = smoothingWindowInfo[i].getSm();
				if (plusStrand && minusStrand == false) info.setStrand("+");
				else if (minusStrand && plusStrand == false) info.setStrand("-");
				else info.setStrand(".");
				
				//final scores = numberObs, mean, PVal			
				saveSmoothedHeatMapData (1, sm, info, meanDir, "#FF00FF"); //magenta
				saveSmoothedHeatMapData (2, sm, info, pvalDir, "#00FF00"); //green
			}	
	}

	/**Saves bar heatmap/ stairstep graph files*/
	public void saveSmoothedHeatMapData (int scoreIndex, SmoothingWindow[] sm, Info info, File dir, String color){
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = new HashMap<String,String>();		
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		//color red
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, color);
		//what's the source
		String fileNames = Misc.stringArrayToString(IO.fetchFileNames(nonConvertedPointDirs),",");
		map.put(BarParser.SOURCE_TAG, fileNames);
		//what's window size
		map.put(BarParser.WINDOW_SIZE, windowSize+"");
		//what's the unit on the scores
		map.put(BarParser.UNIT_TAG, scoreUnits[scoreIndex]);
		//description
		map.put(BarParser.DESCRIPTION_TAG, scoreDescriptions[scoreIndex]);
		//save in info
		info.setNotes(map);
		//get heatmap positions and values
		HeatMapMaker hm = new HeatMapMaker(scoreIndex,0);
		PointData pd = hm.makeHeatMapPositionValues(sm, false);
		pd.setInfo(info);
		pd.writePointData(dir);
		//clean up
		pd.nullPositionScoreArrays();
	}

	/**Sets score names/ descriptions/ units base.*/
	public void setScoreStrings(boolean fdrPresent){
		//final scores = numberObs, pseMedian, FDR
			scoreNames = new String[]{
					"NumberObs",
					"PseMedianBaseFracEdits",
					"-10Log10(FDR)",
			};
			scoreDescriptions = new String[]{
					"# base fraction edits in window",
					"Pseudo median of the base fraction edits in the window",
					"B&H corrected permutation chiSquare -10Log10(PVal)",
			};
			scoreUnits = new String[]{
					"Count",
					"Fraction",
					"-10Log10(FDR)",
			};
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new RNAEditingScanSeqs(args);
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
					case 's': saveDirectory = new File(args[++i]); break;
					case 'w': windowSize = Integer.parseInt(args[++i]); break;
					case 'm': minimumNumberObservationsInWindow = Integer.parseInt(args[++i]); break;
					case 'p': minimumPseMedian = Float.parseFloat(args[++i]); break;
					case 'b': minimumBaseFractionEdited = Float.parseFloat(args[++i]); break;
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


		if (windowSize == 0 ) Misc.printExit("\nPlease enter a positive length for the window size if you are going to set the peak shift to 0\n");

		//make window maker 
		windowMaker = new WindowMaker(windowSize,minimumNumberObservationsInWindow);

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		if (saveDirectory.exists() == false) saveDirectory.mkdir();

	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           RNA Editing Scan Seqs: June 2013                       **\n" +
				"**************************************************************************************\n" +
				"RESS attempts to identify clustered editing sites across a genome using a sliding\n" +
				"window approach.  Each window is scored for the pseudomedian of the base fraction\n" +
				"edits as well as the probability that the observations occured by chance using a\n" +
				"permutation test based on the chiSquare goodness of fit statistic. \n\n" +

				"Options:\n"+
				"-s Save directory, full path.\n"+
				"-e Edited PointData directory from the RNAEditingPileUpParser.\n" +
				"       These should contain stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories. These will be merged when scanning.\n" +
				"-r Reference PointData directory from the RNAEditingPileUpParser. Ditto.\n" +
				
				"\nAdvanced Options:\n"+
				"-p Minimum pseudomedian, defaults to 0.005.\n"+
				"-b Minimum base fraction edited to use in analysis, defaults to 0.01\n"+
				"-w Window size, defaults to 50.\n"+
				"-m Minimum number observations in window, defaults to 3. \n" +
				"-t Run a stranded analysis, defaults to non-stranded.\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/RNAEditingScanSeqs -s /Results/RESS -p 0.01\n" +
				"-e /PointData/Edited -r /PointData/Reference \n\n" +

				"**************************************************************************************\n");

	}

	public SmoothingWindowInfo[] getSmoothingWindowInfo() {
		return smoothingWindowInfo;
	}

	public File getSwiFile() {
		return swiFile;
	}
}
