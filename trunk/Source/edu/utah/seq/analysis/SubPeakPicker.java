package edu.utah.seq.analysis;

import java.io.*;
import java.util.*;

import edu.utah.seq.data.*;
import trans.tpmap.WindowMaker;
import util.gen.*;

/**Use to find the peaks within EnrichedRegions*/
public class SubPeakPicker {
	//fields
	private EnrichedRegion[] enrichedRegions;
	private File[] treatmentPointDirs;
	private File[] controlPointDirs;
	private File tempDirectory;
	private File fullPathToR;
	private double expectedFractionUp;
	private double expectedFractionDown;
	private double scalarTC;
	private double scalarCT;
	private HashMap<String, PointData[]> treatments;
	private HashMap<String, PointData[]> controls;
	private int windowSize = 50;
	private int halfPeakShift;
	private ArrayList<SmoothingWindow> allSmoothingWindows = new ArrayList<SmoothingWindow>();


	//constructors
	/**Picks the best peak based on rescanning each enriched region using the defined window size.*/
	public SubPeakPicker (EnrichedRegion[] enrichedRegions, File[] treatmentPointDirs, File[] controlPointDirs, int windowSize, int halfPeakShift, File tempDirectory, File fullPathToR){
		this.enrichedRegions = enrichedRegions;
		this.treatmentPointDirs = treatmentPointDirs;
		this.controlPointDirs = controlPointDirs;
		this.windowSize = windowSize;
		this.halfPeakShift = halfPeakShift;
		this.tempDirectory = tempDirectory;
		this.fullPathToR = fullPathToR;
		pickPeaks();
	}


	//methods
	private void fetchData(){
		treatments = PointData.fetchStrandedCombinePointData(treatmentPointDirs);
		double numberTreatmentObservations = PointData.totalObservationsMultiPointData(treatments);
		//System.out.println("\t\t"+(int)numberTreatmentObservations+" Treatment Observations");
		if (controlPointDirs != null){
			controls = PointData.fetchStrandedCombinePointData(controlPointDirs);
			double numberControlObservations = PointData.totalObservationsMultiPointData(controls);
			//System.out.println("\t\t"+(int)numberControlObservations+" Control Observations");
			//some stats for binomial pvalue calc and scaling the ratios
			expectedFractionUp = numberTreatmentObservations/(numberTreatmentObservations+numberControlObservations);		
			scalarTC = numberTreatmentObservations/ numberControlObservations;
			scalarCT = numberControlObservations/numberTreatmentObservations;
		}
	}

	/**Main method of class. This method will then re window scan each EnrichedRegion and set the subWindows and bestSubWindow.*/
	public void pickPeaks(){
		System.out.println("PickingPeaks");
		//fetch data 
		fetchData();

		if (controlPointDirs == null) scoreTreatments();

		else {
			//window score each EnrichedRegion
			scoreTreatmentsControls();

			//calculate binomial p-values for each subwindow
			calculateSubWindowBinomialPValues();

			//find and assign best peak for each EnrichedRegion
			findPeaks();

			//add stats to best peak for whole ER 
			calculateStatsForEnrichedRegions();
		}

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

	public void calculateStatsForEnrichedRegions(){
		//bsw 
		for (int x=0; x< enrichedRegions.length; x++){
			SmoothingWindow bestSubWindow = findHighestPeak(enrichedRegions[x].getSubWindows(), 0);
			enrichedRegions[x].setBestSubWindow(bestSubWindow);
		}
		//fetch observations, also calculate number uniques and log2Ratio
		int[][] treatContObs = new int[enrichedRegions.length][2];
		for (int i=0; i< enrichedRegions.length; i++){
			Point[] t = enrichedRegions[i].getTreatmentPoints();
			int numT = 0;
			if (t != null) numT = t.length;
			Point[] c = enrichedRegions[i].getControlPoints();
			int numC = 0;
			if (c != null) numC = c.length;
			treatContObs[i] = new int[]{numT, numC};
			enrichedRegions[i].setLog2Ratio(calculateLog2Ratio(numT, numC));
		}
		//calculate binomial pvalues
		double[] pvals = Num.binomialPValues(expectedFractionUp, tempDirectory, treatContObs, fullPathToR, true, false, false);
		if (pvals == null || pvals.length != enrichedRegions.length) Misc.printErrAndExit("\nProblem calculating binomial p-values for enriched regions.\n");
		//assign
		for (int i=0; i< enrichedRegions.length; i++) enrichedRegions[i].setBinomialPValue((float)pvals[i]);
	}

	public void calculateSubWindowBinomialPValues(){
		SmoothingWindow[] sm = new SmoothingWindow[allSmoothingWindows.size()];
		allSmoothingWindows.toArray(sm);
		//fetch observations
		int[][] treatContObs = new int[sm.length][2];
		for (int i=0; i< sm.length; i++){
			float[] scores = sm[i].getScores();
			treatContObs[i] = new int[]{(int)scores[0], (int)scores[1]};
		}
		//calculate binomial pvalues
		double[] pvals = Num.binomialPValues(expectedFractionUp, tempDirectory, treatContObs, fullPathToR, true, false, false);
		if (pvals == null || pvals.length != sm.length) Misc.printErrAndExit("\nProblem calculating a binomial p-value for sub windows.\n");
		//assign
		for (int i=0; i< sm.length; i++) sm[i].setScores(new float[]{(float)pvals[i]});
	}

	/**Finds the best peak within each loaded enriched region.  Call scoreTreatmentsControls() or scoreTreatments() first.*/
	public void findPeaks(){
		for (int x=0; x< enrichedRegions.length; x++){
			SmoothingWindow bestSubWindow = findHighestPeak(enrichedRegions[x].getSubWindows(), 0);
			enrichedRegions[x].setBestSubWindow(bestSubWindow);
		}
	}

	/**Find the window with the largest score.  If multiple windows have the same best score then the start of the first and stop
	 * of the last are used to create a new window and are returned.*/
	public static SmoothingWindow findHighestPeak(SmoothingWindow[] sm, int scoreIndex){
		if (sm == null || sm.length ==0) {
			return null;
		}
		//find highest score
		float bestScore = sm[0].getScores()[scoreIndex];
		for (int i=1; i< sm.length; i++){
			float testScore = sm[i].getScores()[scoreIndex];
			if (testScore > bestScore) bestScore = testScore;
		}
		//find smallest start and biggest stop
		int start = sm[sm.length-1].getStart();
		int end = 0;
		for (int i=0; i< sm.length; i++){
			float testScore = sm[i].getScores()[scoreIndex];
			if (testScore == bestScore) {
				//set start?
				if (sm[i].getStart() < start) start = sm[i].getStart();
				//set stop?
				if (sm[i].getStop() > end) end = sm[i].getStop();
			}
		}
		return new SmoothingWindow(start, end, new float[]{bestScore});
	}

	public void scoreTreatments(){

		//fetch PointData
		PointData[] pdT = null;

		//for each enriched region, these are sorted by chromosome
		String currentChromosome = "";
		for (int x=0; x< enrichedRegions.length; x++){
			//check chromosome and load PointData
			String testChromosome = enrichedRegions[x].getChromosome();
			if (currentChromosome.equals(testChromosome) == false){
				currentChromosome = testChromosome;
				//fetch PointData
				pdT = treatments.get(currentChromosome);
				//check for plus and minus
				if (pdT[0] == null || pdT[1] == null) {
					Misc.printErrAndExit("\nError: one of your PointData is not paired with + and - stranded data.  " +
							"Check your treatment PointData files for "+currentChromosome+"\n");
				}
				//shift	
				if (halfPeakShift != 0){
					pdT[0].shiftPositions(halfPeakShift);
					pdT[1].shiftPositions(halfPeakShift);
				}
			}

			enrichedRegions[x].setNumberPlusObservations(pdT[0].countPoints(enrichedRegions[x].getStart(), enrichedRegions[x].getStop()));
			enrichedRegions[x].setNumberMinusObservations(pdT[1].countPoints(enrichedRegions[x].getStart(), enrichedRegions[x].getStop()));
		}

	}

	public void scoreTreatmentsControls(){
		//make window maker
		WindowMaker wm = new WindowMaker(windowSize,1);

		//fetch PointData
		PointData[] pdT = null;
		PointData[] pdC = null;
		PointData treatPD = null;
		PointData ctrlPD = null;
		int[] treatPos = null;
		int[] ctrlPos = null;

		//for each enriched region, these are sorted by chromosome
		String currentChromosome = "";
		for (int x=0; x< enrichedRegions.length; x++){
			//check chromosome and load PointData
			String testChromosome = enrichedRegions[x].getChromosome();
			if (currentChromosome.equals(testChromosome) == false){
				currentChromosome = testChromosome;
				//fetch PointData
				pdT = treatments.get(currentChromosome);
				pdC = controls.get(currentChromosome);
				//check for plus and minus
				if (pdT[0] == null || pdT[1] == null || pdC[0] == null || pdC[1] == null) {
					Misc.printErrAndExit("\nError: one of your PointData is not paired with + and - stranded data.  " +
							"Check your treatment or control PointData files for "+currentChromosome+"\n");
				}
				//shift				
				pdT[0].shiftPositions(halfPeakShift);
				pdT[0].stripScores();
				pdT[1].shiftPositions(halfPeakShift);
				pdT[1].stripScores();
				pdC[0].shiftPositions(halfPeakShift);
				pdC[0].stripScores();
				pdC[1].shiftPositions(halfPeakShift);
				pdC[1].stripScores();
				//merge
				treatPD = PointData.combinePointData(pdT, false);
				ctrlPD = PointData.combinePointData(pdC, false);
				//get positions
				treatPos = treatPD.getPositions();
				ctrlPos = ctrlPD.getPositions();
			}
			//find indexes		
			int[] treatIndexes = treatPD.findIndexes(enrichedRegions[x].getStart(), enrichedRegions[x].getStop());
			int[] controlIndexes = ctrlPD.findIndexes(enrichedRegions[x].getStart(), enrichedRegions[x].getStop());

			//fetch positions
			ArrayList<Integer> positions = new ArrayList<Integer>();
			for (int i=treatIndexes[0]; i< treatIndexes[1]; i++) positions.add(new Integer(treatPos[i]));
			for (int i=controlIndexes[0]; i< controlIndexes[1]; i++) positions.add(new Integer(ctrlPos[i]));
			int[] pos = Misc.integerArrayListToIntArray(positions);
			if (pos == null){
				System.out.println(currentChromosome+"\t"+enrichedRegions[x].getStart()+"\t"+enrichedRegions[x].getStop());
				Misc.printExit("\nNo reads found within region? Are you sure you have entered the same PointData as what was used in ScanSeqs?\n");
			}
			Arrays.sort(pos);

			//create windows
			int[][] windows = wm.makeWindows(pos);		
			//assign bp positions
			for (int i=0; i< windows.length; i++){
				windows[i][0] = pos[windows[i][0]];
				windows[i][1] = pos[windows[i][1]]+1;	//last base isn't included
			}

			//for each window
			ArrayList<SmoothingWindow> smAL = new ArrayList<SmoothingWindow>();
			for (int i=0; i< windows.length; i++){
				//fetch scores
				float tSum = treatPD.sumPositionBP(windows[i][0], windows[i][1]);
				float cSum = ctrlPD.sumPositionBP(windows[i][0], windows[i][1]);
				//calculate sums, save window
				SmoothingWindow sm = new SmoothingWindow (windows[i][0], windows[i][1], new float[]{tSum,cSum});
				smAL.add(sm);
				//add to all for bin pval
				allSmoothingWindows.add(sm);
			}
			//convert and add to enriched region
			SmoothingWindow[] sm = new SmoothingWindow[smAL.size()];
			smAL.toArray(sm);			
			enrichedRegions[x].setSubWindows(sm);		

			//for rescore entire regions
			Point[] tPoints = treatPD.fetchPoints(enrichedRegions[x].getStart(), enrichedRegions[x].getStop());
			Point[] cPoints = ctrlPD.fetchPoints(enrichedRegions[x].getStart(), enrichedRegions[x].getStop());
			enrichedRegions[x].setTreatmentPoints(tPoints);
			enrichedRegions[x].setControlPoints(cPoints);

			//calculate uniques
			Point[] tPlus = null;
			int numTPlusUni = 0;
			Point[] tMinus = null;
			int numTMinusUni = 0;

			tPlus = pdT[0].fetchPoints(enrichedRegions[x].getStart(), enrichedRegions[x].getStop());
			if (tPlus!=null) numTPlusUni = Point.sumIdenticalPositionScores(tPlus).length;
			tMinus = pdT[1].fetchPoints(enrichedRegions[x].getStart(), enrichedRegions[x].getStop());
			if (tMinus!=null) numTMinusUni = Point.sumIdenticalPositionScores(tMinus).length;
			enrichedRegions[x].setNumberUniqueTreatmentObservations(numTPlusUni+ numTMinusUni);

			Point[] cPlus = null;
			int numCPlusUni = 0;
			Point[] cMinus = null;
			int numCMinusUni = 0;

			cPlus = pdC[0].fetchPoints(enrichedRegions[x].getStart(), enrichedRegions[x].getStop());
			if (cPlus != null) numCPlusUni = Point.sumIdenticalPositionScores(cPlus).length;
			cMinus = pdC[1].fetchPoints(enrichedRegions[x].getStart(), enrichedRegions[x].getStop());
			if (cMinus != null) numCMinusUni = Point.sumIdenticalPositionScores(cMinus).length;
			enrichedRegions[x].setNumberUniqueControlObservations(numCPlusUni+ numCMinusUni);
		}

	}

}
