package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import trans.tpmap.*;

/**Takes PointData and generates sliding window scan statistics. Estimates FDRs using two different methods.
 * @author Nix
 * */
public class ScanSeqs {

	//user defined fields
	private File[] treatmentPointDirs;
	private File[] controlPointDirs;
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private int windowSize = -1;
	private int peakShift = -1;
	private float numberStandardDeviations = 4;
	private int minimumNumberReadsInWindow = 2;
	private boolean printGraphs = true;
	private boolean printPointGraphs = false;
	private boolean stripHotWindows = false;
	private boolean useWeightedReads = false;
	private boolean findReducedRegions = false;
	private boolean filterWindowsViaQValue = true;
	private boolean plusStrand = true;
	private boolean minusStrand =  true;

	//internal fields
	private int halfPeakShift = 0;
	private double numberTreatmentObservations = 0;
	private double numberControlObservations = 0;
	private double expectedFractionUp;
	private double expectedFractionDown;
	private double scalarTC;
	private double scalarCT;
	private HashMap<String,PointData[]> treatmentPlusPointData;
	private HashMap<String,PointData[]> treatmentMinusPointData;
	private HashMap<String,PointData[]> controlPlusPointData;
	private HashMap<String,PointData[]> controlMinusPointData;
	private ArrayList<SmoothingWindowInfo> smiAL = new ArrayList<SmoothingWindowInfo>();
	private int numberCachedBinPVals = 500;
	private float[][] pValsDataUp;
	private float[][] pValsDataDown;
	private float[][] pValsSkew;
	private WindowMaker windowMaker; 
	private int[][] windows;
	private SmoothingWindow[] smoothingWindow;
	private SmoothingWindowInfo[] smoothingWindowInfo;
	private File swiFile;
	private ArrayList<SmoothingWindow> allSmoothingWindowAL = new ArrayList<SmoothingWindow>();
	private float[] controlThresholds;
	private int[] controlNumberERs;
	private float qValueThresholdForSaving = 4;
	private float pValueThresholdForEmpFDR = 7;
	private String[] scoreNames;
	private String[] scoreDescriptions;
	private String[] scoreUnits;
	private String adapterName = "chrAdapter";
	private boolean verbose = true;

	//by chromosome
	private String chromosome;
	private PointData treatmentChromPlus = null;
	private PointData treatmentChromMinus = null;	
	private PointData controlChromPlus = null;
	private PointData controlChromMinus = null;	

	//for Ken's weighted reads
	private double totalTreatmentScores;
	private double totalControlScores;
	private float defaultProportionT;
	private float defaultProportionC;

	//constructors

	/**For integration with the RNASeq app.*/
	public ScanSeqs(File[] treatmentPointDirs, File[] controlPointDirs, File saveDirectory, File fullPathToR, int windowSize, int peakShift, int minimumNumberReadsInWindow, boolean findReducedRegions, boolean verbose){
		//set params
		this.treatmentPointDirs = treatmentPointDirs;
		this.controlPointDirs = controlPointDirs;
		this.saveDirectory = saveDirectory;
		this.fullPathToR = fullPathToR;
		this.windowSize = windowSize;
		this.peakShift = peakShift;
		halfPeakShift = (int)Math.round( ((double)peakShift)/2 );
		this.minimumNumberReadsInWindow = minimumNumberReadsInWindow;
		this.findReducedRegions = findReducedRegions;
		this.verbose = verbose;

		//make window maker 
		windowMaker = new WindowMaker(windowSize,minimumNumberReadsInWindow);

		//set score items
		setScoreStrings();

		scan();
	}

	/**Stand alone.*/
	public ScanSeqs(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		scan();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	//methods

	public void scan(){
		//info
		if (verbose){
			if (useWeightedReads) System.out.println("Using read score probabilities...");
			if (findReducedRegions) System.out.println("Scanning for enriched and reduced regions, no EmpFDR estimation...");
			else System.out.println("Scanning for enriched regions...");
			System.out.println("\t"+minimumNumberReadsInWindow+ "\tMinium number reads in window");
			System.out.println("\t"+peakShift+ "\tPeak shift");
			System.out.println("\t"+windowSize+ "\tWindow size");
		}
		
		//fetch counts
		if (verbose) System.out.println("Calculating read count stats and binomial p-value look up tables...");
		calculateReadCountStatistics();

		//for each chromosome
		if (verbose) System.out.print("Scanning Chromosomes");
		String[] chromosomes = fetchAllChromosomes();
		for (int i=0; i< chromosomes.length; i++){
			chromosome = chromosomes[i];
			windowScanChromosome();
		}
		if (verbose) System.out.println();

		//any data, make arrays
		if (smiAL.size() == 0) Misc.printExit("\nNo windows found. Thresholds too strict?\n");
		smoothingWindowInfo = new SmoothingWindowInfo[smiAL.size()];
		smiAL.toArray(smoothingWindowInfo);

		//control data? run FDR estimations
		if (controlPointDirs != null){
			//strip windows with high control counts?
			if (stripHotWindows){
				if (verbose) System.out.println("Removing windows with high number of control reads...");
				removeWindowsWithHighControlReadCounts();
			}		

			//collect remaining SmoothingWindows
			ArrayList<SmoothingWindow> smAL = new ArrayList<SmoothingWindow>();
			for (int i=0; i< smoothingWindowInfo.length; i++){
				SmoothingWindow[] win = smoothingWindowInfo[i].getSm();
				for (int j = 0; j< win.length; j++){
					smAL.add(win[j]);
				}
			}
			smoothingWindow = new SmoothingWindow[smAL.size()];
			smAL.toArray(smoothingWindow);

			//clean up!
			System.gc(); System.gc(); System.gc();
			
			//convert binomial pvalues to FDRs
			if (verbose) System.out.println("Converting fuzzy binomial p-values to q-value FDRs in R...");

			//calculate reduce regions and EmpFDRs?
			if (findReducedRegions == false){
				convertPValuesToQValuesFuzzedUp();

				if (verbose) System.out.println("Estimating empirical FDRs...");
				//process controls for empFDRs
				processControlsForEmpFDR();

				//first for regions with enrichment, 
				//scores = pVal,upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus, cSumMinus,realScore,mockScore,upEmpFDR
				int scoreIndex = 8;
				int fdrIndex = 10;			
				estimateEmpiricalFDRs(smoothingWindow, scoreIndex, fdrIndex);

			}
			else convertPValuesToQValuesFuzzedUpDown();

			//collapse scores
			collapseScores();

			//filter windows ?
			if (filterWindowsViaQValue){
				if (verbose) System.out.println("Removing windows with a qvalue FDR threshold of > 40% ");
				filterWindows();
			}
		}

		//save window data 
		if (verbose) System.out.println("Saving serialized window data...");
		String ext = "";
		if (plusStrand == true && minusStrand == true) ext = "";
		else if (plusStrand) ext = "Plus";
		else if (minusStrand) ext = "Minus";

		swiFile = new File (saveDirectory, "windowData"+windowSize+"bp"+ext+".swi");
		IO.saveObject(swiFile, smoothingWindowInfo);

		//write out bar file graphs, call after saving since Info is modified
		if (printGraphs){
			if (verbose) System.out.println("Writing bar file graphs...");
			writeBarFileGraphs();
		}
	}

	/**Collapse scores from:
	 * scores = pVal, upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
	 * to:
	 * scores = pVal, qValFDR, eFDR, pValSkew, Log2((sumT+1)/(sumC+1)), tSumPlus, tSumMinus, cSumPlus, cSumMinus*/
	public void collapseScores(){
		for (int i=0; i< smoothingWindow.length; i++){
			float[] scores = smoothingWindow[i].getScores();
			float[] col = new float[9];
			//pval
			col[0] = scores[0];
			//qValueFDR
			if (scores[1]>=scores[2]) col[1] = scores[1];
			else col[1] = scores[2] * -1;
			//eFDR
			if (findReducedRegions) col[2] = 0;
			else col[2] = scores[10];	
			//skew
			col[3] = scores[3];
			//log2(ratio)
			col[4] = calculateLog2Ratio(scores[4]+scores[5], scores[6]+scores[7]);
			//sums
			col[5] = scores[4];
			col[6] = scores[5];
			col[7] = scores[6];
			col[8] = scores[7];
			smoothingWindow[i].setScores(col);
		}
	}

	/**Fetches the names of all the chromosomes in the data.*/
	public String[] fetchAllChromosomes(){
		HashSet<String> c = new HashSet<String>();
		Iterator<String> it = treatmentPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = treatmentMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		if(controlPointDirs != null) {
			it = controlPlusPointData.keySet().iterator();
			while (it.hasNext()) c.add(it.next());
			it = controlMinusPointData.keySet().iterator();
			while (it.hasNext()) c.add(it.next());
		}
		return Misc.hashSetToStringArray(c);
	}

	/**Fetchs the data for a particular chromosome.*/
	public boolean fetchData(){

		//fetch treatment
		treatmentChromPlus = null;
		if (treatmentPlusPointData.containsKey(chromosome) && plusStrand) treatmentChromPlus = PointData.combinePointData(treatmentPlusPointData.get(chromosome), true);
		treatmentChromMinus = null;
		if (treatmentMinusPointData.containsKey(chromosome) && minusStrand) treatmentChromMinus = PointData.combinePointData(treatmentMinusPointData.get(chromosome), true);
		//fetch control
		if (controlPointDirs != null){
			controlChromPlus = null;
			if (controlPlusPointData.containsKey(chromosome) && plusStrand) controlChromPlus = PointData.combinePointData(controlPlusPointData.get(chromosome), true);
			controlChromMinus = null;
			if (controlMinusPointData.containsKey(chromosome) && minusStrand) controlChromMinus = PointData.combinePointData(controlMinusPointData.get(chromosome), true);

			//check that datasets are present when processing two sample comparison
			if (plusStrand){
				if (treatmentChromPlus==null || controlChromPlus==null) return false;
			}
			if (minusStrand){
				if (treatmentChromMinus==null || controlChromMinus==null) return false;
			}


		}
		return true;
	}

	/**Window scans a chromosome collecting read count data and calculating binomial p-values.*/
	public void windowScanChromosome(){
		//fetch data
		if (fetchData() == false) {
			if (verbose) System.out.println("\n\tSkipping "+chromosome+". Failed to find paired PointData datasets.");
			return;
		}
		//fetch, shift, and merge all positions from the treatment and control
		int[] positions = fetchShiftStripPointData();
		//make windows using all of the reads
		makeWindows(positions);
		//any windows?
		if (windows.length == 0){
			if (verbose) System.out.println("\n\tSkipping "+chromosome+". No windows found with minimum reads of "+minimumNumberReadsInWindow+" within a window size of "+windowSize);
			return;
		}
		if (verbose) System.out.print(" "+chromosome);
		//make SmoothingWindow[] container
		smoothingWindow = new SmoothingWindow[windows.length];
		//scan
		if (controlPointDirs !=null) {
			calculateBinomialPValues();
			if (findReducedRegions == false) calculateBinomialPValsForEmpiricalFDRs();
		}
		else sumTreatmentReads();
		//save results
		saveWindowDataToSMIArrayList();
	}

	/**Saves window data to SMI Array.*/
	public void saveWindowDataToSMIArrayList(){
		Info info;
		if (treatmentChromPlus != null) info= treatmentChromPlus.getInfo();
		else info = treatmentChromMinus.getInfo();
		HashMap<String,String> notes = new HashMap<String,String>();
		notes.put("totalTreatmentObservations",""+(int)numberTreatmentObservations);
		if (controlPointDirs !=null) notes.put("totalControlObservations",""+(int)numberControlObservations);
		notes.put(BarParser.WINDOW_SIZE, windowSize+"");
		notes.put(BarParser.BP_3_PRIME_SHIFT, halfPeakShift+"");
		notes.put(BarParser.DESCRIPTION_TAG, Misc.stringArrayToString(scoreNames, ","));
		notes.put(BarParser.UNIT_TAG, Misc.stringArrayToString(scoreUnits, ","));
		info.setNotes(notes);
		info.setStrand(".");		
		smiAL.add(new SmoothingWindowInfo(smoothingWindow, info));
	}

	/**Calculates a log2( (tSum+1)/(cSum+1) ) on linearly scaled tSum and cSum based on the total observations.*/
	public float calculateLog2Ratio( double tSum, double cSum){
		double t;
		double c;
		if (tSum !=0 ) {
			t = tSum * scalarCT;
			c = cSum;
		}
		else {
			c = cSum * scalarTC;
			t = tSum;
		}
		double ratio = (t+1)/(c+1);
		return (float)Num.log2(ratio);
	}

	/**Keeps windows with a pVal score that is less than or greater than the pValueThreshold.*/
	private int filterWindows(){
		int numberWindows = 0;
		float negPValueThreshold = -1 * qValueThresholdForSaving;
		ArrayList<SmoothingWindowInfo> al = new ArrayList<SmoothingWindowInfo>();
		for (int i=0; i< smoothingWindowInfo.length; i++){
			ArrayList<SmoothingWindow> posWinAL = new ArrayList<SmoothingWindow>();
			SmoothingWindow[] allWin = smoothingWindowInfo[i].getSm();
			for (int j=0; j< allWin.length; j++){
				float pVal = allWin[j].getScores()[0];
				if (pVal > qValueThresholdForSaving || pVal < negPValueThreshold) {
					posWinAL.add(allWin[j]);
				}
			}
			if (posWinAL.size() !=0){
				numberWindows += posWinAL.size();
				SmoothingWindow[] posWin = new SmoothingWindow[posWinAL.size()];
				posWinAL.toArray(posWin);
				al.add(new SmoothingWindowInfo(posWin, smoothingWindowInfo[i].getInfo()));
			}
		}
		//convert
		smoothingWindowInfo = new SmoothingWindowInfo[al.size()];
		al.toArray(smoothingWindowInfo);		
		return numberWindows;
	}

	/**Calculates the number of EnrichedRegions are made at all thresholds using the control vs control data.*/
	private void processControlsForEmpFDR (){

		//scores = pVal, upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
		//loop through SmoothingWindows saving those with a control pval > pValueThresholdForEmpFDR 
		HashSet<Float> controlScoresHM = new HashSet<Float>();
		SmoothingWindowInfo[] controlWindowInfo = new SmoothingWindowInfo[smoothingWindowInfo.length];
		for (int i=0; i< smoothingWindowInfo.length; i++){
			ArrayList<SmoothingWindow> posWinAL = new ArrayList<SmoothingWindow>();
			SmoothingWindow[] allWin = smoothingWindowInfo[i].getSm();
			for (int j=0; j< allWin.length; j++){
				float mockPVal = allWin[j].getScores()[9];
				if (mockPVal > pValueThresholdForEmpFDR) {
					posWinAL.add(allWin[j]);
					controlScoresHM.add(new Float(mockPVal));
				}
			}
			SmoothingWindow[] posWin = new SmoothingWindow[posWinAL.size()];
			posWinAL.toArray(posWin);
			controlWindowInfo[i] = new SmoothingWindowInfo(posWin, smoothingWindowInfo[i].getInfo());
		}
		//Make and sort scores
		controlThresholds = Num.hashSetToFloat(controlScoresHM);
		Arrays.sort(controlThresholds);
		controlNumberERs = new int[controlScoresHM.size()];
		controlScoresHM = null;

		//Make ERM
		EnrichedRegionMaker controlERM = new EnrichedRegionMaker(windowSize, controlWindowInfo);
		//for each score calculate number of ERs and save
		for (int i=0; i< controlThresholds.length; i++){
			controlNumberERs[i] = controlERM.countEnrichedRegions(true, 9, controlThresholds[i]);
			if (controlNumberERs[i] == 0) break;
		}
		//clean up 
		controlWindowInfo = null;
		controlScoresHM = null;
	}

	/**Removes windows that exceed 4 standard deviations.*/
	private void removeWindowsWithHighControlReadCounts(){
		//calculate the standard deviation and median
		ArrayList<Float> countsAL = new ArrayList<Float>();
		for (int i=0; i< smoothingWindowInfo.length; i++){
			SmoothingWindow[] win = smoothingWindowInfo[i].getSm();
			for (int j = 0; j< win.length; j++){
				//scores = pVal, upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
				float[] scores = win[j].getScores();
				countsAL.add(new Float (scores[6]+ scores[7]));
			}
		}
		//calc median
		float[] controlCounts = Num.arrayListOfFloatToArray(countsAL);
		Arrays.sort(controlCounts);
		double controlCountMedian = Num.median(controlCounts);
		double stndDev = Num.standardDeviation(controlCounts);
		float controlWindowThreshold = (float)(controlCountMedian + stndDev * numberStandardDeviations);
		if (verbose) System.out.println("\tMedian "+controlCountMedian+"\tStndDev "+Num.formatNumberOneFraction(stndDev)+"\tCntrlWinThres "+(int)controlWindowThreshold);

		//strip windows
		double numPassingWindows = 0;
		double totalWindows = 0;
		for (int i=0; i< smoothingWindowInfo.length; i++){
			SmoothingWindow[] win = smoothingWindowInfo[i].getSm();
			totalWindows+= win.length;
			ArrayList<SmoothingWindow> goodWin = new ArrayList<SmoothingWindow>();
			for (int j = 0; j< win.length; j++){
				//scores = pVal, upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
				float[] scores = win[j].getScores();
				if ((scores[6]+scores[7]) < controlWindowThreshold) {
					goodWin.add(win[j]);
					numPassingWindows++;
				}
			}
			win = new SmoothingWindow[goodWin.size()];
			goodWin.toArray(win);
			smoothingWindowInfo[i].setSm(win);
		}
		if (verbose) System.out.println("\t"+Num.formatNumber(numPassingWindows/totalWindows, 3)+" windows passed control sum threshold");
	}

	/**Converts binomial pvalues in SmoothingWindowInfo to qvalues using Storey's method.
	 * Adds random noise to windows with counts <21
	 **/
	private void convertPValuesToQValuesFuzzedUp(){
		//make up and down fuzzy pvalue generators
		FuzzyBinomialPValueGenerator upFuzz = new FuzzyBinomialPValueGenerator( saveDirectory, fullPathToR, expectedFractionUp);
		//collect -10Log10(pvalues)
		float[] up = new float[smoothingWindow.length+1];
		float maxUpPVal = 0;
		for (int i=0; i< smoothingWindow.length; i++){
			//scores = pVal, upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
			float[] scores = smoothingWindow[i].getScores();
			//low number of obs? then fuzz pval
			int totalObs = Math.round(scores[4]+scores[5]+scores[6]+scores[7]);
			if (totalObs < 21){
				up[i] = (float)upFuzz.fetchFuzzyPVal((int)(scores[4]+scores[5]), (int)(scores[6]+scores[7]), true);
				if (scores[1]> maxUpPVal) maxUpPVal = scores[1];
			}
			else{
				up[i] = scores[1];
			}
		}

		//set max
		up[smoothingWindow.length] = maxUpPVal;
		//calc qvals
		float[] upQVals = Num.qValueFDR(saveDirectory, up, fullPathToR, true, true);

		//transform total
		float transTotal = (float)Num.minus10log10(smoothingWindow.length);

		//assign to windows and multiple test correct the skew pvalue
		ArrayList<SmoothingWindow> toSetQVal = new ArrayList<SmoothingWindow>();
		for (int i=0; i< smoothingWindow.length; i++){
			//scores = pVal, upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
			float[] scores = smoothingWindow[i].getScores();
			int totalObs = Math.round(scores[4]+scores[5]+scores[6]+scores[7]);
			if (totalObs < 21) {
				toSetQVal.add(smoothingWindow[i]);
				upFuzz.addAdjPValQVal(up[i],upQVals[i]);
			}
			else {
				scores[1] = upQVals[i];
				scores[2] = 0;
			} 
			//multiple test correct skew binomial pvalue
			scores[3] = scores[3] + transTotal;
			if (scores[3] < 0) scores[3] = 0;
		}		
		//set max in FuzzyBinomialPValGenerators
		upFuzz.addAdjPValQVal(up[smoothingWindow.length],upQVals[smoothingWindow.length]);

		//make qvalue calculators
		upFuzz.generateSortedPAndQValArrays();
		//convert toSeqQval
		SmoothingWindow[] toFixQVal = new SmoothingWindow[toSetQVal.size()];
		toSetQVal.toArray(toFixQVal);
		//load qVals
		upFuzz.loadQValues(1, toFixQVal);
	}

	/**Converts binomial pvalues in SmoothingWindowInfo to qvalues using Storey's method.
	 * Adds random noise to windows with counts <21
	 **/
	private void convertPValuesToQValuesFuzzedUpDown(){
		//make up and down fuzzy pvalue generators
		FuzzyBinomialPValueGenerator upFuzz = new FuzzyBinomialPValueGenerator( saveDirectory, fullPathToR, expectedFractionUp);
		FuzzyBinomialPValueGenerator downFuzz = new FuzzyBinomialPValueGenerator( saveDirectory, fullPathToR, expectedFractionDown);

		//collect -10Log10(pvalues)
		float[] up = new float[smoothingWindow.length+1];
		float[] down = new float[smoothingWindow.length+1];
		float maxPValUp = 0;
		float maxPValDown = 0;
		for (int i=0; i< smoothingWindow.length; i++){
			//scores = pVal, upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
			float[] scores = smoothingWindow[i].getScores();
			//low number of obs? then fuzz pval
			int totalObs = Math.round(scores[4]+scores[5]+scores[6]+scores[7]);
			if (totalObs < 21){
				up[i] = (float)upFuzz.fetchFuzzyPVal((int)(scores[4]+scores[5]), (int)(scores[6]+scores[7]), true);
				down[i] = (float)downFuzz.fetchFuzzyPVal((int)(scores[6]+scores[7]), (int)(scores[4]+scores[5]), true);
				//set max
				if (scores[1] > maxPValUp) maxPValUp = scores[1];
				if (scores[2] > maxPValDown) maxPValDown = scores[2];
			}
			else{
				up[i] = scores[1];
				down[i] = scores[2];
			}
		}

		//set max as last
		up[smoothingWindow.length] = maxPValUp;
		down[smoothingWindow.length] = maxPValDown;

		//calc qvals
		float[] upQVals = Num.qValueFDR(saveDirectory, up, fullPathToR, true, true);
		float[] downQVals = Num.qValueFDR(saveDirectory, down, fullPathToR, true, true);

		//transform total
		float transTotal = (float)Num.minus10log10(smoothingWindow.length);

		//assign to windows and multiple test correct the skew pvalue
		ArrayList<SmoothingWindow> toSetQVal = new ArrayList<SmoothingWindow>();
		for (int i=0; i< smoothingWindow.length; i++){
			//scores = pVal, upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
			float[] scores = smoothingWindow[i].getScores();
			int totalObs = Math.round(scores[4]+scores[5]+scores[6]+scores[7]);
			if (totalObs < 21) {
				toSetQVal.add(smoothingWindow[i]);
				upFuzz.addAdjPValQVal(up[i],upQVals[i]);
				downFuzz.addAdjPValQVal(down[i],downQVals[i]);
			}
			else {
				scores[1] = upQVals[i];
				scores[2] = downQVals[i];
			} 
			//multiple test correct skew binomial pvalue
			scores[3] = scores[3] + transTotal;
			if (scores[3] < 0) scores[3] = 0;
		}
		//set max in FuzzyBinomialPValGenerators
		upFuzz.addAdjPValQVal(up[smoothingWindow.length],upQVals[smoothingWindow.length]);
		downFuzz.addAdjPValQVal(down[smoothingWindow.length],downQVals[smoothingWindow.length]);
		//make qvalue calculators
		upFuzz.generateSortedPAndQValArrays();
		downFuzz.generateSortedPAndQValArrays();
		//convert toSeqQval
		SmoothingWindow[] toFixQVal = new SmoothingWindow[toSetQVal.size()];
		toSetQVal.toArray(toFixQVal);
		//load qVals
		upFuzz.loadQValues(1, toFixQVal);
		downFuzz.loadQValues(2, toFixQVal);
	}


	/**Calculates binomial p-values for emp FDR estimation.
	 * Only call this after you are done with the PointData for a particular chromosome!
	 * This method subsamples and replaces the data!*/
	private void calculateBinomialPValsForEmpiricalFDRs(){
		//subsample PointData so T: Ca: Cb, on a chromosome basis not genome basis, theoretical problem here?! 
		PointData[] pdTCC = subSamplePointDataForEmpFDR();
		ArrayList<SmoothingWindow> forBinom = new ArrayList<SmoothingWindow>();
		//for each window
		for (int i=0; i< smoothingWindow.length; i++){
			//fetch bases
			int bpStart = smoothingWindow[i].getStart();
			int bpStop = smoothingWindow[i].getStop();
			//fetch window sums
			float tSum = pdTCC[0].sumScoreBP(bpStart, bpStop);
			float c1Sum = pdTCC[1].sumScoreBP(bpStart, bpStop);
			float c2Sum = pdTCC[2].sumScoreBP(bpStart, bpStop);
			//calc pval from cache?
			//yes
			if (tSum < numberCachedBinPVals && c1Sum < numberCachedBinPVals && c2Sum< numberCachedBinPVals){
				//calc pvals
				int t = Math.round(tSum);
				int c1 = Math.round(c1Sum);
				int c2 = Math.round(c2Sum);
				//calc realPVal
				float up = pValsSkew[t][c2];
				float down = pValsSkew[c2][t];
				float realPVal;
				if (up > down) realPVal = up;
				else realPVal = -1 * down;
				//calc mockPVal
				up = pValsSkew[c1][c2];
				down = pValsSkew[c2][c1];
				float mockPVal;
				if (up > down) mockPVal = up;
				else mockPVal = -1 * down;
				//trim pVals sig fig to one, xxx.x
				realPVal = pseudoRound(realPVal);
				mockPVal = pseudoRound(mockPVal);
				//scores = pVal, upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
				float[] scores = smoothingWindow[i].getScores();
				scores[8] = realPVal;
				scores[9] = mockPVal;
			}
			//no so save for calculating with a call to R
			else {
				float[] scores = smoothingWindow[i].getScores();
				scores[8] = tSum;
				scores[9] = c1Sum;
				scores[10] = c2Sum;
				forBinom.add(smoothingWindow[i]);
			}
		}
		if (forBinom.size() != 0) calculatePValsFromRForEmpFDR(forBinom);
	}
	/**Calcs p-vals from R for EmpFDR that aren't in the look up tables.*/
	private void calculatePValsFromRForEmpFDR (ArrayList<SmoothingWindow> al){
		//for each window calculate 4 pvals, up/down real, up/down mock
		int num = al.size();
		int[][] obs = new int[num*4][2];
		int index = 0;
		for (int i=0; i< num; i++){
			float[] scores = al.get(i).getScores();
			//scores = pVal,upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus, tSum , c1Sum, c2Sum 
			int t = Math.round(scores[8]);
			int c1 = Math.round(scores[9]);
			int c2 = Math.round(scores[10]);
			//t vs c2
			obs[index++] = new int[]{t,c2};
			//c2 vs t
			obs[index++] = new int[]{c2,t};
			//c1 vs c2
			obs[index++] = new int[]{c1,c2};
			//c2 vs c1
			obs[index++] = new int[]{c2,c1};
		}
		//fetch pvalues
		double[] pVals = Num.binomialPValues(0.5, saveDirectory, obs, fullPathToR, true, false, false);
		if (pVals == null || pVals.length != (4*num)) Misc.printErrAndExit("\nProblem fetching binomial pvalues from R for EmpFDR. Check that R really is at "+fullPathToR);
		//assign
		index = 0;
		for (int i=0; i< num; i++){
			float realPVal;
			float upReal = (float) pVals[index++];
			float downReal = (float) pVals[index++];
			if (upReal > downReal) realPVal = upReal;
			else realPVal = -1* downReal;
			float mockPVal;
			float upMock = (float) pVals[index++];
			float downMock = (float) pVals[index++];
			if (upMock > downMock) mockPVal = upMock;
			else mockPVal = -1* downMock;
			//trim pVals sig fig to one, xxx.x
			realPVal = pseudoRound(realPVal);
			mockPVal = pseudoRound(mockPVal);
			//scores = pVal, upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
			float[] scores = al.get(i).getScores();
			scores[8] = realPVal;
			scores[9] = mockPVal;
			scores[10] = 0;
		}
	}

	/**Subsamples point data for empirical FDR calculations. t=c1=c2*/
	private PointData[] subSamplePointDataForEmpFDR(){
		//count observations
		int halfCtrlObs = (int)Math.round(numberControlObservations/2);
		//merge PointData
		PointData treatmentPD = PointData.combinePointData(new PointData[]{treatmentChromPlus,treatmentChromMinus}, true);
		PointData controlPD = PointData.combinePointData(new PointData[]{controlChromPlus,controlChromMinus}, true);
		//make c 2x t or t 1/2 c

		//match number of treatments to half the number of controls
		//subsample treatment?
		if (halfCtrlObs < numberTreatmentObservations) {
			double totalNumberToMatch = halfCtrlObs;
			double totalNumberInChrom = treatmentChromPlus.getInfo().getNumberObservations() + treatmentChromMinus.getInfo().getNumberObservations();
			int numToFetchForChrom = (int)Math.round(totalNumberToMatch * totalNumberInChrom/numberTreatmentObservations);
			if (numToFetchForChrom < treatmentPD.getInfo().getNumberObservations())treatmentPD = PointData.fetchRandomObservations(treatmentPD, numToFetchForChrom);
		}
		//subsample controls?
		else if (halfCtrlObs > numberTreatmentObservations) {
			double totalNumberToMatch = numberTreatmentObservations *2;
			double totalNumberInChrom = controlChromPlus.getInfo().getNumberObservations() + controlChromMinus.getInfo().getNumberObservations();
			int numToFetchForChrom = (int)Math.round(totalNumberToMatch * totalNumberInChrom/numberControlObservations);
			if (numToFetchForChrom < controlPD.getInfo().getNumberObservations()) controlPD = PointData.fetchRandomObservations(controlPD, numToFetchForChrom);
		}
		//split controls
		Point[] controlPts = Point.makePoints(controlPD.getPositions(), controlPD.getScores());
		Point[][] splitPts = Point.split(controlPts);
		//set
		PointData c1 = Point.extractPositionScores(splitPts[0]);
		PointData c2 = Point.extractPositionScores(splitPts[1]);
		return new PointData[]{treatmentPD,c1,c2};
	}

	/**For cases when just interested in treatment data.*/
	private void sumTreatmentReads(){
		//for each window 
		for (int i=0; i< windows.length; i++){
			//fetch scores
			float tSumPlus = 0; 
			float tSumMinus = 0;
			if (treatmentChromPlus != null) tSumPlus = treatmentChromPlus.sumScoreBP(windows[i][0], windows[i][1]); 
			if (treatmentChromMinus != null)tSumMinus = treatmentChromMinus.sumScoreBP(windows[i][0], windows[i][1]);
			float tSum = tSumPlus+ tSumMinus;
			//scores = tSum, tSumPlus,tSumMinus
			float[] scores = new float[]{tSum,tSumPlus,tSumMinus};
			//make window
			smoothingWindow[i] = new SmoothingWindow (windows[i][0], windows[i][1], scores);
			allSmoothingWindowAL.add(smoothingWindow[i]);
		}
	}

	/**Main window scanner using binomial p-value as the score based on read counts.
	 * Be sure to replace the PointData scores with 1 if you want to use the real bin p-val.*/
	private void calculateBinomialPValues(){
		//make arrays to calculate binomial p-vals with R
		ArrayList<SmoothingWindow> forBinom = new ArrayList<SmoothingWindow>();

		//for each window 
		for (int i=0; i< windows.length; i++){
			//fetch scores
			float tSumPlus = 0; 
			float tSumMinus = 0;
			float cSumPlus = 0;
			float cSumMinus = 0;
			if (treatmentChromPlus != null) tSumPlus = treatmentChromPlus.sumScoreBP(windows[i][0], windows[i][1]); 
			if (treatmentChromMinus != null)tSumMinus = treatmentChromMinus.sumScoreBP(windows[i][0], windows[i][1]);
			if (controlChromPlus != null)cSumPlus = controlChromPlus.sumScoreBP(windows[i][0], windows[i][1]);
			if (controlChromMinus != null)cSumMinus = controlChromMinus.sumScoreBP(windows[i][0], windows[i][1]);

			float tSum = tSumPlus+ tSumMinus;
			float cSum = cSumPlus+ cSumMinus;
			//scores =                 pVal,upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus, realND,mockND,upEmpFDR
			float[] scores = new float[]{0,   0,     0,       0, tSumPlus,tSumMinus,cSumPlus, cSumMinus,  0,     0,      0};

			//make window
			smoothingWindow[i] = new SmoothingWindow (windows[i][0], windows[i][1], scores);
			allSmoothingWindowAL.add(smoothingWindow[i]);
			//calc diff binomial from cache?
			float max = tSum;
			if (max < cSum) max = cSum;
			//yes
			if (max < numberCachedBinPVals){
				//calc pvals
				int t = Math.round(tSum);
				int c = Math.round(cSum);
				float up = pValsDataUp[t][c];
				float down = pValsDataDown[c][t];
				//only test for skew in one direction, that minus is way less than plus
				float skew = pValsSkew[Math.round(tSumPlus)][Math.round(tSumMinus)];
				//set scores
				//scores = pVal,upPVal,downPVal,skew,tSumPlus,tSumMinus,cSum
				if (up > down) scores[0] = up;
				else scores[0] = -1* down;
				scores[1] = up;
				scores[2] = down;
				scores[3] = skew;
				smoothingWindow[i].setScores(scores);
			}
			//no
			else forBinom.add(smoothingWindow[i]);
		}
		//calc binom from R
		if (forBinom.size()!=0) calculateBinomialsFromR(forBinom);
	}

	/**Takes SmoothingWindows that had observations that exceeded cached pvalues and calculates
	 * directly from R.*/
	private void calculateBinomialsFromR(ArrayList<SmoothingWindow> sm){
		int num = sm.size();
		//for each win, want to calculate up, down, skew
		int[][] skewObs = new int[num][2];
		int[][] upObs = new int[num][2];
		int[][] downObs = new int[num][2];
		for (int i=0; i< num; i++){
			//scores = pVal,upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus, realND,mockND,upEmpFDR
			float[] scores = sm.get(i).getScores();
			int tSumPlus = Math.round(scores[4]);
			int tSumMinus = Math.round(scores[5]);
			int tSum = tSumPlus + tSumMinus;
			int cSum = Math.round(scores[6]+scores[7]);
			//set obs
			skewObs[i] = new int[]{tSumPlus, tSumMinus};
			upObs[i] = new int[]{tSum,cSum};
			downObs[i] = new int[]{cSum,tSum};
		}
		//fetch pvalues
		double[] skewPvals = Num.binomialPValues(0.5, saveDirectory, skewObs, fullPathToR, true, false, false);
		double[] upPvals = Num.binomialPValues(expectedFractionUp, saveDirectory, upObs, fullPathToR, true, false, false);
		double[] downPvals = Num.binomialPValues(expectedFractionDown, saveDirectory, downObs, fullPathToR, true, false, false);
		if (skewPvals == null || upPvals == null || downPvals == null || skewPvals.length != num || upPvals.length != num || downPvals.length != num) {
			Misc.printErrAndExit("\nProblem fetching binomial pvalues from R. Check that R really is at "+fullPathToR);
		}
		//assign scores
		for (int i=0; i< num; i++){
			SmoothingWindow s = sm.get(i);
			//scores = pVal,upPVal,downPVal,skew,tSumPlus,tSumMinus,cSum,realND,mockND,upEmpFDR
			float[] scores = s.getScores();
			scores[1] = (float)upPvals[i];
			scores[2] = (float)downPvals[i];
			scores[3] = (float)skewPvals[i];
			//assign pVal score
			if (scores[1]> scores[2]) scores[0] = scores[1];
			else scores[0] = -1 * scores[2];
		}
	}

	/**Uses weighted reads in window scanning to generate binomial p-values.*/
	/*
	private void calculateWeightedReadBinomialPValues(){
		//internal fields
		BinomialDistributionImpl b;
		float numTObs = (float)numberTreatmentObservations;
		float numCObs = (float)numberControlObservations;

		//make arrays to calculate binomial p-vals with R
		ArrayList<SmoothingWindow> forBinom = new ArrayList<SmoothingWindow>();
		//for each window 
		for (int i=0; i< windows.length; i++){
			//fetch float[]{numReads, sumScores}
			float[] tPlus = null; 
			float[] tMinus = null;
			float[] cPlus = null;
			float[] cMinus = null;
			//summary stats
			float tNumReads = 0;
			float cNumReads = 0;
			float tSumScores = 0;
			float cSumScores = 0;
			//might be null due to strandedness
			if (treatPlusPD != null) {
				tPlus = treatPlusPD.sumScoresPositionsBP(windows[i][0], windows[i][1]);
				tNumReads += tPlus[0];
				tSumScores += tPlus[1];
			}
			if (treatMinusPD != null) {
				tMinus = treatMinusPD.sumScoresPositionsBP(windows[i][0], windows[i][1]);
				tNumReads += tMinus[0];
				tSumScores += tMinus[1];
			}
			if (ctrlPlusPD != null) {
				cPlus = ctrlPlusPD.sumScoresPositionsBP(windows[i][0], windows[i][1]);
				cNumReads += cPlus[0];
				cSumScores += cPlus[1];
			}
			if (ctrlMinusPD != null) {
				cMinus = ctrlMinusPD.sumScoresPositionsBP(windows[i][0], windows[i][1]);
				cNumReads += cMinus[0];
				cSumScores += cMinus[1];
			}
			//calc binomial pvalue 
			float tProportion = defaultProportionT;
			float cProportion = defaultProportionC;
			if (tNumReads !=0) tProportion = tSumScores/tNumReads;
			if (cNumReads !=0) cProportion = cSumScores/cNumReads;

			float cProbNumT = cProportion * numTObs;
			float probOfSuccess = cProbNumT / (cProbNumT  + (tProportion * numCObs));
			int numTrials = (int) (tNumReads + cNumReads);
			b = new BinomialDistributionImpl(numTrials, probOfSuccess);
			double pvalUp = 0;

			try {
				pvalUp = 1-b.cumulativeProbability(tNumReads-1);
			} catch (MathException e){
				e.printStackTrace();
			}
			if (pvalUp==0) pvalUp = 300;
			else pvalUp = Num.minus10log10(pvalUp);
			float[] scores = new float[]{(float) pvalUp};

			 //if (verbose) System.out.println(tNumReads +"\t"+ tSumScores +"\t"+ tProportion +"\t"+ cNumReads +"\t"+ cSumScores +"\t"+ cProportion +"\t"+ numTrials +"\t"+ probOfSuccess+"\t"+pvalUp);

			//scores =                 pVal,upPVal,downPVal,skew,tSumPlus,tSumMinus,cSumPlus,cSumMinus, realND,mockND,upEmpFDR
			//float[] scores = new float[]{0,   0,     0,       0, tSumPlus,tSumMinus,cSumPlus, cSumMinus,  0,     0,      0,       0};

			//make window
			smoothingWindow[i] = new SmoothingWindow (windows[i][0], windows[i][1], scores);
			allSmoothingWindowAL.add(smoothingWindow[i]);
			//calc diff binomial from cache?
			/*
			float max = tSum;
			if (max < cSum) max = cSum;
			//yes
			if (max < numberCachedBinPVals){
				//calc pvals
				int t = Math.round(tSum);
				int c = Math.round(cSum);
				float up = pValsDataUp[t][c];
				float down = pValsDataDown[c][t];
				//only test for skew in one direction, that minus is way less than plus
				float skew = pValsSkew[Math.round(tSumPlus)][Math.round(tSumMinus)];
				//set scores
				//scores = pVal,upPVal,downPVal,skew,tSumPlus,tSumMinus,cSum
				if (up > down) scores[0] = up;
				else scores[0] = -1* down;
				scores[1] = up;
				scores[2] = down;
				scores[3] = skew;
				smoothingWindow[i].setScores(scores);
			}
			//no
			else forBinom.add(smoothingWindow[i]);
	 */
	/*
		}
		//calc binom from R
		if (forBinom.size()!=0) calculateBinomialsFromR(forBinom);


	}*/

	/**Collects and calculates a bunch of stats re the PointData.*/
	private void calculateReadCountStatistics(){
		//fetch treatment PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (treatmentPointDirs);
		//remove chrAdapter
		combo[0].remove(adapterName);
		combo[1].remove(adapterName);
		treatmentPlusPointData = PointData.convertArrayList2Array(combo[0]);
		treatmentMinusPointData = PointData.convertArrayList2Array(combo[1]);
		if (numberTreatmentObservations == 0){
			if (plusStrand) numberTreatmentObservations = PointData.totalObservationsMultiPointData(treatmentPlusPointData);
			if (minusStrand) numberTreatmentObservations += PointData.totalObservationsMultiPointData(treatmentMinusPointData);
		}
		if (verbose) System.out.println("\t"+(int)numberTreatmentObservations+"\tTreatment Observations");

		//likewise control data
		if (controlPointDirs != null){
			combo = PointData.fetchStrandedPointDataNoMerge (controlPointDirs);
			//remove chrAdapter
			combo[0].remove(adapterName);
			combo[1].remove(adapterName);
			controlPlusPointData = PointData.convertArrayList2Array(combo[0]);
			controlMinusPointData = PointData.convertArrayList2Array(combo[1]);
			if (numberControlObservations == 0){
				if (plusStrand) numberControlObservations = PointData.totalObservationsMultiPointData(controlPlusPointData);
				if (minusStrand) numberControlObservations += PointData.totalObservationsMultiPointData(controlMinusPointData);
			}
			if (verbose) System.out.println("\t"+(int)numberControlObservations+"\tControl Observations");
			//some stats for binomial pvalue calc and scaling the ratios
			expectedFractionUp = numberTreatmentObservations/(numberTreatmentObservations+numberControlObservations);
			expectedFractionDown = 1- expectedFractionUp;		
			scalarTC = numberTreatmentObservations/ numberControlObservations;
			scalarCT = numberControlObservations/numberTreatmentObservations;

			//make look up tables
			pValsDataUp = Num.convertToFloat(Num.binomialPValMatrix(numberCachedBinPVals+1, expectedFractionUp, saveDirectory, fullPathToR, true));
			pValsDataDown = Num.convertToFloat(Num.binomialPValMatrix(numberCachedBinPVals+1, expectedFractionDown, saveDirectory, fullPathToR, true));
			pValsSkew = Num.convertToFloat(Num.binomialPValMatrix(numberCachedBinPVals+1, 0.5, saveDirectory, fullPathToR, true));		

			if (useWeightedReads) {
				totalTreatmentScores = PointData.totalScoreMultiPointData(treatmentPlusPointData);
				//check to see if scores present
				if (totalTreatmentScores == -1) Misc.printExit("\nCannot use weighted reads. Cannot calculate total scores.  Reparse your data using an updated USeq package.\n");
				totalTreatmentScores += PointData.totalScoreMultiPointData(treatmentMinusPointData);
				totalControlScores = PointData.totalScoreMultiPointData(controlPlusPointData);
				totalControlScores += PointData.totalScoreMultiPointData(controlMinusPointData);
				defaultProportionT = (float)(totalTreatmentScores/numberTreatmentObservations);
				defaultProportionC = (float)(totalControlScores/numberControlObservations);
				if (verbose) System.out.println("Default proportions "+defaultProportionT+"\t"+ defaultProportionC);
			}
		}
	}

	/**Shifts the positions halfPeakShift (+ for sense, - for antisense) sets the positions into the data
	 * returns all of the positions after sorting. May replace all scores with 1 if stripScores == true.*/
	private int[] fetchShiftStripPointData(){
		ArrayList<int[]> posAL = new ArrayList<int[]>();
		//fetch data from treatments
		if (treatmentChromPlus !=null) {
			int[] p = treatmentChromPlus.getPositions();
			addShift(p,halfPeakShift);
			posAL.add(p);
			if (useWeightedReads == false) treatmentChromPlus.stripScores();
		}
		if (treatmentChromMinus!=null) {
			int[] p = treatmentChromMinus.getPositions();
			addShift(p, -1*halfPeakShift);
			posAL.add(p);
			if (useWeightedReads == false) treatmentChromMinus.stripScores();
		}
		//fetch data from controls
		if (controlPointDirs != null) {
			if (controlChromPlus !=null) {
				int[] p = controlChromPlus.getPositions();
				addShift(p, halfPeakShift);
				posAL.add(p);
				if (useWeightedReads == false) controlChromPlus.stripScores();
			}
			if (controlChromMinus!=null){
				int[] p = controlChromMinus.getPositions();
				addShift(p, -1*halfPeakShift);
				posAL.add(p);
				if (useWeightedReads == false) controlChromMinus.stripScores();
			}
		}
		//merge
		int[][] toMerge = new int[posAL.size()][];
		for (int i=0; i< posAL.size(); i++) toMerge[i] = posAL.get(i);
		int[] merged = Num.collapseIntArray(toMerge);
		Arrays.sort(merged);
		//return
		return merged;
	}

	/**Adds the toAdd to each int.*/
	public static void addShift(int[] positions, int toAdd){
		for (int i=0; i< positions.length; i++){
			positions[i] += toAdd;
			if (positions[i]<0) positions[i] = 0;
		}
	}

	/**Given a score threshold, returns the number of EnrichedRegions that would be generated using the
	 * precomputed values.  This is a way of minimizing the memory requirements.*/
	private int calculateNumberControlERs(float threshold){
		int index = Arrays.binarySearch(controlThresholds, threshold);
		if (index < 0){
			int conv = -index -1;
			if (conv >= controlThresholds.length) return 0;
			else return controlNumberERs[conv];
		}
		else return controlNumberERs[index];
	}


	/**Calculates FDRs based on the number of EnrichedRegions.*/
	private void estimateEmpiricalFDRs(SmoothingWindow[] allSM, int scoreIndex, int fdrIndex){

		//sort windows by score
		Arrays.sort(allSM, new ComparatorSmoothingWindowScore(scoreIndex));

		//find first window with normScore greater than pValueThresholdForEmpFDR, set prior eFDRs to 0
		int i=0;
		for (; i< allSM.length; i++){
			float[] scores = allSM[i].getScores();
			if (scores[scoreIndex] > pValueThresholdForEmpFDR )  break;
			else scores[fdrIndex] = 0;
		}

		//make ER makers
		EnrichedRegionMaker realERM = new EnrichedRegionMaker(windowSize, smoothingWindowInfo);

		//for each score calc fdr = # control ERs/ #real ERs
		float testScore = 0;
		double numControlERs = 0;
		float maxFDR = 0;
		//if (verbose) System.out.println("NormDiffIndex "+scoreIndex+" Placing FDR at "+fdrIndex);
		for (; i< allSM.length; i++){
			//calc num control windows
			testScore = allSM[i].getScores()[scoreIndex];
			numControlERs = calculateNumberControlERs(testScore);
			//any control windows?
			//yes, controls, calc FDR and assign to maxFDR
			if (numControlERs !=0) {				
				//calc num real
				double numReal = realERM.countEnrichedRegions(false, scoreIndex, testScore);
				//calc fdr
				float fdr = 0;
				if (numReal !=0) fdr = new Double(-10 * Num.log10(numControlERs/numReal)).floatValue();
				if (fdr > maxFDR) maxFDR = fdr;
				//if (verbose) System.out.println(i+ "index: testScore-> "+testScore+" numControl-> "+numControlERs+" numReal-> "+numReal+" fdr-> "+fdr+" maxFDR-> "+maxFDR);				
				//advance until score changes and assign max
				for (; i< allSM.length; i++){
					//different so break
					if (allSM[i].getScores()[scoreIndex] != testScore) break;
					//not different so set maxFDR
					else {						
						float[] s = allSM[i].getScores();
						s[fdrIndex] = maxFDR;
						//if (verbose) System.out.println("\t"+i+" Setting MaxFDR-> "+maxFDR+" score-> "+s[scoreIndex]);						
					}
				}
				//back off one so calc numcontrol on diff score
				if (i < allSM.length) i--;
			}
			//no  more controls so assign maxFDR to remainder
			else {
				for (; i< allSM.length; i++){
					float[] s = allSM[i].getScores();
					s[fdrIndex] = maxFDR;
					//if (verbose) System.out.println("\t"+i+" Setting Final "+maxFDR+" score-> "+s[scoreIndex]);					
				}
			}
		}		
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
		if (controlPointDirs != null){
			//printin point data graphs?
			String type = "";
			File pValPt = null;
			File qValFDRPt = null;
			File empFDRPt = null;
			File log2RatioPt = null;
			if (printPointGraphs) {
				type = "_StairStep";
				pValPt = new File(saveDirectory, "BinPVal_Point");
				qValFDRPt = new File(saveDirectory, "QValFDR_Point");
				empFDRPt = new File(saveDirectory, "EmpFDR_Point");
				log2RatioPt = new File(saveDirectory, "Log2Ratio_Point");
				pValPt.mkdir();
				qValFDRPt.mkdir();
				log2RatioPt.mkdir();
				if (findReducedRegions == false) empFDRPt.mkdir();
			}
			//make directories
			File pVal = new File(saveDirectory, "BinPVal"+type);
			pVal.mkdir();
			File qValFDR = new File(saveDirectory, "QValFDR"+type);
			qValFDR.mkdir();
			File log2Ratio = new File(saveDirectory, "Log2Ratio"+type);
			log2Ratio.mkdir();
			File empFDR = new File(saveDirectory, "EmpFDR"+type);
			boolean upDown = true;
			if (findReducedRegions == false) {
				empFDR.mkdir();
				upDown = false;
			}
			//scores = pVal, qValFDR, eFDR, pValSkew, log2(tSum/cSum), tSumPlus,tSumMinus,cSumPlus,cSumMinus
			//for each chromosome
			for (int i=0; i< smoothingWindowInfo.length; i++){
				Info info = smoothingWindowInfo[i].getInfo();
				SmoothingWindow[] sm = smoothingWindowInfo[i].getSm();
				saveSmoothedHeatMapData (0, sm, info, pVal, "#FFFFFF", upDown); //white
				saveSmoothedHeatMapData (1, sm, info, qValFDR, "#FF0000", upDown); //red
				saveSmoothedHeatMapData (4, sm, info, log2Ratio, "#FF00FF", upDown); //purple
				if (findReducedRegions == false) saveSmoothedHeatMapData (2, sm, info, empFDR, "#00FF00", false); //green
				//save point graphs
				if (printPointGraphs){
					//(int scoreIndex, SmoothingWindow[] sm, Info info, File dir, String color)
					saveSmoothedPointData (0,sm,info,pValPt,"#E8E8E8"); //off white
					saveSmoothedPointData (1,sm,info,qValFDRPt,"#C00000"); //off red
					saveSmoothedPointData (4,sm,info,log2RatioPt,"#CC9900"); //mustard
					if (findReducedRegions == false) saveSmoothedPointData(2, sm, info, empFDRPt,"#00CC00"); //off green
				}
			}
		}
		else {
			//make directories
			File sum = new File(saveDirectory, "Sum");
			sum.mkdir();
			File sumPlus = new File(saveDirectory, "Sum+");
			sumPlus.mkdir();
			File sumMinus = new File(saveDirectory, "Sum-");
			sumMinus.mkdir();
			//scores = tSum, tSumPlus,tSumMinus
			//for each chromosome
			for (int i=0; i< smoothingWindowInfo.length; i++){
				Info info = smoothingWindowInfo[i].getInfo();
				SmoothingWindow[] sm = smoothingWindowInfo[i].getSm();
				saveSmoothedHeatMapData (0, sm, info, sum, "#FF00FF", false); //magenta
				info.setStrand("+");
				saveSmoothedHeatMapData (1, sm, info, sumPlus, "#00FFFF", false); //cyan
				info.setStrand("-");
				saveSmoothedHeatMapData (2, sm, info, sumMinus, "#FFFF00", false); //yellow
			}
		}
	}

	/**Saves bar heatmap/ stairstep graph files*/
	public void saveSmoothedHeatMapData (int scoreIndex, SmoothingWindow[] sm, Info info, File dir, String color, boolean posNeg){
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = info.getNotes();		
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		//color red
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, color);
		//what's the source
		String fileNames = Misc.stringArrayToString(IO.fetchFileNames(treatmentPointDirs),",");
		if (controlPointDirs !=null) fileNames = fileNames+" vs "+Misc.stringArrayToString(IO.fetchFileNames(controlPointDirs),",");
		map.put(BarParser.SOURCE_TAG, fileNames);
		//what's window size
		map.put(BarParser.WINDOW_SIZE, windowSize+"");
		//what's the peak shift
		map.put(BarParser.BP_3_PRIME_SHIFT, halfPeakShift+"");
		//what's the unit on the scores
		map.put(BarParser.UNIT_TAG, scoreUnits[scoreIndex]);
		//description
		map.put(BarParser.DESCRIPTION_TAG, scoreDescriptions[scoreIndex]);
		//save in info
		info.setNotes(map);
		//get heatmap positions and values
		PointData pd;		
		if (posNeg){
			HeatMapMakerPosNeg hm = new HeatMapMakerPosNeg(scoreIndex, 0, 0);
			pd = hm.makeHeatMapPositionValues(sm);
		}
		else {
			HeatMapMaker hm = new HeatMapMaker(scoreIndex,0);
			pd = hm.makeHeatMapPositionValues(sm, false);
		}

		pd.setInfo(info);
		pd.writePointData(dir);
		//clean up
		pd.nullPositionScoreArrays();
	}

	/**Saves bar point graph files*/
	public void saveSmoothedPointData (int scoreIndex, SmoothingWindow[] sm, Info info, File dir, String color){
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = info.getNotes();
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
		//color
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, color);
		//what's the source
		String fileNames = Misc.stringArrayToString(IO.fetchFileNames(treatmentPointDirs),",");
		if (controlPointDirs !=null) fileNames = fileNames+" vs "+Misc.stringArrayToString(IO.fetchFileNames(controlPointDirs),",");
		map.put(BarParser.SOURCE_TAG, fileNames);
		//what's window size
		map.put(BarParser.WINDOW_SIZE, windowSize+"");
		//what's the peak shift
		map.put(BarParser.BP_3_PRIME_SHIFT, halfPeakShift+"");
		//what's the unit on the scores
		map.put(BarParser.UNIT_TAG, scoreUnits[scoreIndex]);
		//description
		map.put(BarParser.DESCRIPTION_TAG, scoreDescriptions[scoreIndex]);
		//save in info
		info.setNotes(map);
		//make centered position and
		int[] positions = new int[sm.length];
		float[] scores = new float[sm.length];
		for (int i=0; i< sm.length; i++){
			positions[i] = Num.calculateMiddleIntergenicCoordinates(sm[i].getStart(), sm[i].getStop());
			scores[i] = sm[i].getScores()[scoreIndex];
		}
		//write bar
		PointData pd = new PointData();
		pd.setInfo(info);
		pd.setPositions(positions);
		pd.setScores(scores);
		pd.writePointData(dir);
		//clean up
		pd.nullPositionScoreArrays();
		positions = null;
		scores = null;
		sm = null;
	}

	/**Sets score names/ descriptions/ units base on whether control data is present.*/
	public void setScoreStrings(){
		if (controlPointDirs != null){
			scoreNames = new String[]{
					"BinPVal",
					"QValueFDR",
					"EmpFDR",
					"SkewPVal",
					"Log2((sumT+1)/(sumC+1))",
					"SumT+",
					"SumT-",
					"SumC+",
					"SumC-",
			};
			scoreDescriptions = new String[]{
					"Binomial p-value",
					"Q-value FDR estimation from binomial p-value",
					"Empirical FDR estimation",
					"Bonferroni corrected binomial p-value for strand skew",
					"Log2((sumT+1)/(sumC+1)) of linearly scaled sumT and sumC",
					"Sum Treatment sense strand",
					"Sum Treatment antisense strand",
					"Sum Control sense strand",
					"Sum Control minus strand"
			};
			scoreUnits = new String[]{
					"-10Log10(p-value)",	
					"-10Log10(q-valueFDR)",
					"-10Log10(empFDR)",
					"-10Log10(bonCorr p-value)",
					"Log2((sumT+1)/(sumC+1))",
					"Sum",
					"Sum",
					"Sum",
					"Sum"
			};
		}
		else{
			scoreNames = new String[]{
					"Sum",
					"Sum+",
					"SumT-"
			};
			scoreDescriptions = new String[]{
					"Sum Treatment",
					"Sum Treatment sense strand",
					"Sum Treatment antisense strand",
			};
			scoreUnits = new String[]{
					"Sum",
					"Sum",
					"Sum",
			};
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ScanSeqs(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		String strand = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': treatmentPointDirs = IO.extractFiles(args[++i]); break;
					case 'c': controlPointDirs = IO.extractFiles(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'p': peakShift = Integer.parseInt(args[++i]); break;
					case 'w': windowSize = Integer.parseInt(args[++i]); break;
					case 'm': minimumNumberReadsInWindow = Integer.parseInt(args[++i]); break;
					case 'i': printGraphs = false; break;
					case 'n': printPointGraphs = true; break;
					case 'q': filterWindowsViaQValue = false; break;
					case 'f': stripHotWindows = true; break;
					case 'e': findReducedRegions = true; break;
					case 'u': useWeightedReads = true; break;
					case 'j': strand = args[++i]; break;
					case 'a': numberTreatmentObservations = Double.parseDouble(args[++i]); break;
					case 'b': numberControlObservations = Double.parseDouble(args[++i]); break;
					case 'g': numberStandardDeviations = Float.parseFloat(args[++i]); stripHotWindows = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//dumb users?
		if (stripHotWindows && findReducedRegions) Misc.printErrAndExit("\nDon't elect to remove windows with high control read counts " +
		"(option -f and -g) when looking for reduced regions! You'll be tossing your best hits.\n");

		//look for point directories
		if (treatmentPointDirs == null || treatmentPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your treatment PointData directories(s)!\n");
		//only one directory look deeper
		if (treatmentPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(treatmentPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) treatmentPointDirs = otherDirs;
		}

		//control data
		if (controlPointDirs != null){
			if (controlPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your control PointData directories(s)!\n");
			//only one directory look deeper
			if (controlPointDirs.length == 1){
				File[] otherDirs = IO.extractOnlyDirectories(controlPointDirs[0]);
				if (otherDirs != null && otherDirs.length > 0) controlPointDirs = otherDirs;
			}
		}
		//set half peak shift and windowSize
		if (peakShift == -1) Misc.printExit("\nPlease enter a peak shift, the results from running the PeakShiftFinder application.\n");
		halfPeakShift = (int)Math.round( ((double)peakShift)/2 );
		if (windowSize == -1) windowSize = peakShift;
		if (windowSize == 0 ) Misc.printExit("\nPlease enter a positive length for the window size if you are going to set the peak shift to 0\n");

		//strand?
		if (strand != null && (strand.equals("+") == false && strand.equals("-") == false)){
			Misc.printExit("\nError: Enter either + or - for a stranded scan.\n");
		}
		plusStrand =  (strand == null || strand.equals("+"));
		minusStrand =  (strand == null || strand.equals("-"));

		//make window maker 
		windowMaker = new WindowMaker(windowSize,minimumNumberReadsInWindow);

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		if (saveDirectory.exists() == false) saveDirectory.mkdir();

		//check for R and required libraries
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}
		else {
			/*String errors = IO.runRCommandLookForError("library(DESeq)", fullPathToR, saveDirectory);
			if (errors == null || errors.length() !=0){
				Misc.printExit("\nError: Cannot find the required R library.  Did you install DESeq " +
						"(http://www-huber.embl.de/users/anders/DESeq/)?  See the author's websites for installation instructions. Once installed, " +
						"launch an R terminal and type 'library(DESeq)' to see if it is present. R error message:\n\t\t"+errors+"\n\n");
			}*/
			String errors = IO.runRCommandLookForError("library(qvalue)", fullPathToR, saveDirectory);
			if (errors == null || errors.length() !=0){
				Misc.printExit("\nError: Cannot find the required R library.  Did you install qvalue " +
						"(http://genomics.princeton.edu/storeylab/qvalue/)?  See the author's websites for installation instructions. Once installed, " +
						"launch an R terminal and type 'library(qvalue)' to see if it is present. R error message:\n\t\t"+errors+"\n\n");
			}
		} 

		//set score items
		setScoreStrings();



	}	

	public static float pseudoRound(float d){
		float comp = d*10.0f;
		comp = Math.round(comp);
		return comp;
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                  Scan Seqs: Feb 2012                             **\n" +
				"**************************************************************************************\n" +
				"Takes unshifted stranded chromosome specific PointData and uses a sliding window to\n" +
				"calculate several smoothed window statistics. These include a binomial p-value, a\n" +
				"q-value FDR, an empirical FDR, and a Bonferroni corrected binomial p-value for peak\n" +
				"shift strand skew. These are saved as heat map/ stairstep xxx.bar graph files for\n" +
				"direct viewing in the Integrated Genome Browser. The empFDR is only calculated when\n" +
				"scanning for enriched regions. Provide >2x the # of control reads relative to\n" +
				"treatment to prevent significant sub sampling when calculating the empFDR. If control\n" +
				"data is not provided, simple window sums are calculated.\n\n" +

				"Options:\n"+
				"-s Save directory, full path.\n"+
				"-t Treatment PointData directories, full path, comma delimited. These should\n" +
				"       contain unshifted stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories.\n" +
				"-c Control PointData directories, ditto. \n" +
				"-p Peak shift, see the PeakShiftFinder app. Average distance between + and - strand\n" +
				"       peaks. Will be used to shift the PointData and set the window size.\n"+
				"-r Full path to R loaded with Storey's q-value library, defaults to '/usr/bin/R'\n" +
				"       file, see http://genomics.princeton.edu/storeylab/qvalue/\n"+

				"\nAdvanced Options:\n"+
				"-w Window size, defaults to peak shift. A good alternative window size is the\n" +
				"       peak shift plus the standard deviation, see the PeakShiftFinder app.\n"+
				"-e Scan for both reduced and enriched regions, defaults to look for only enriched\n" +
				"       regions. This turns off the empFDR estimation.\n"+
				"-j Scan only one strand, defaults to both, enter either + or - \n"+
				"-q Don't filter windows using q-value FDR threshold, save all to bar graphs,\n" +
				"       defaults to saving those with a q-value < 40%.\n"+
				"-m Minimum number reads in window, defaults to 2. Increasing this threshold will\n" +
				"       speed up processing considerably but compromises the q-value estimation.\n"+
				"-f Filter windows with high read control read counts. Don't use if looking for\n" +
				"       reduced regions.\n"+
				"-g Control window read count threshold, # stnd devs off median, defaults to 4.\n"+
				"-n Print point graph window representation xxx.bar files.\n"+
				"-a Number treatment observations to use in defining expect and ratio scalars.\n"+
				"-b Number control observations to use in defining expect and ratio scalars.\n"+
				"-u Use read score probabilities (assumes scores are > 0 and <= 1), defaults to\n" +
				"       assigning 1 to each read score. Experimental.\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/ScanSeqs -t\n" +
				"      /Data/PolIIRep1/,/Data/PolIIRep2/ -c /Data/Input1/,Data/Input2/ -s\n" +
				"      /Data/PolIIResults -w 200 -p 100 -f -g 5 \n\n" +

		"**************************************************************************************\n");

	}

	public SmoothingWindowInfo[] getSmoothingWindowInfo() {
		return smoothingWindowInfo;
	}

	public File getSwiFile() {
		return swiFile;
	}
}
