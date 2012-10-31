package trans.main;
import java.io.*;
import java.util.*;
import java.util.regex.*;

import trans.misc.Util;
import trans.tpmap.MapSplitter;
import util.gen.*;

/**
 * Scores windows of intensity values using several tests. 
 */
public class ScanChip {
	
	//fields
	private String tpmapFile; 
	private String infoFile;
	private String resultsFile; 
	private String treatmentDir; 
	private String controlDir; 
	private Window[] windows;
	private ArrayList totalWindowsAL = new ArrayList(10000);
	private int numWindows;
	private int sizeOfOligoMinusOne = 24;
	private int halfOligoLength;
	private boolean convertScoresToQValues = false;
	private boolean pmOnly = true;
	private boolean printOligoRatios = false;
	private boolean printPointSgrs = false;
	private boolean printHeatMapSgrs = false;
	private boolean convertSummaries2Log2 = false;		//still need to fix Heat maps
	private boolean usePseudoMedian = false;
	private boolean averageGroups = false;
	private static String fullPathToPValue = "/home/BioApps/T2/OthersCode/RBSymPTest/symmetric_p_test";	//old version
	private static String fullPathToR = "/usr/local/R/bin/R"; //must have the qvalue package loaded!
	private int minOligos = 10;
	private int maxScores = 1000;
	private float[][] treatmentIntensities;
	private float[][] controlIntensities;
	private String uniqueIdT = "";
	private String uniqueIdC = "";
	private String extension = "";		
	private 	float singleScores[] = null;
	private float[][] layeredScores = null;
	private 	String chromosome;
	private 	int[] positions;
	private int[][] testWindows;
	private PrintWriter outIndividualRatios = null;
	private PrintWriter outSum = null;
	private File outSumFile = null;
	private PrintWriter outScore = null;
	private File outScoreFile = null;
	private PrintWriter tempRatios = null;
	private File tempRatioFile = null;
	private File oligoFile = null;
	//Randomization of labels
	private boolean randomizeLabels = false;
	private int numPermutations = 10;
	//histogram
	private Histogram histogram;  
	private int scoreIndexForHistogram = 1;
	private int pValIndexForRandomPermutation = 2;
	
	/**Array describing the scores in the Window.Scores.  It changes frequently!*/
	private static String[] SCORE_DESCRIPTION = {
		"Wilcoxon Rank Sum, -10Log10(p-val)", 
		"Trimmed Mean or Layered Pseudo Median Relative Difference",
		"-10Log10(uncorrected p-value) based on rnd perm of cel file labels (if set)",
		"-10Log10(FDR) based on rnd perm of cel file labels (if set)",
		"-10Log10(q-valFDR) based on symmetric null distribution, (if set)"};
	
	
	public ScanChip(String[] args){	
		//process args
		processArgs(args);
		
		//build extension set scoreThreshold
		System.out.println();
		buildExtensionSetThreshold();
		
		//get .tpmap file info
		ArrayList info = (ArrayList)IO.fetchObject(new File(infoFile));
		
		//fetch float[replicas][oligoIntensities] of processed intensity values for treatment and controls
		//break both apart by chromosome using the bpmapInfo ArrayList, save to disk
		treatmentIntensities = fetchFloatArrays(treatmentDir,false);
		System.out.println("\nProcessing "+treatmentIntensities.length+" treatment arrays...");
		uniqueIdT = MapSplitter.breakSaveIntensityValues(info, treatmentIntensities,treatmentDir);
		treatmentIntensities = null;
		
		if (pmOnly){
			controlIntensities = fetchFloatArrays(controlDir, false);
			System.out.println("Processing "+controlIntensities.length+" control arrays...");
			uniqueIdC = MapSplitter.breakSaveIntensityValues(info, controlIntensities,controlDir);
			controlIntensities = null;
		}
		
		//calculate number of additional score place holders needed
		int numPlaceHolders = 0;
		if (randomizeLabels) {
			numPlaceHolders+=2;
			//instantiate histogram for random lable permutation
			histogram = new Histogram(-2,2,5000);
		}
		if (convertScoresToQValues) numPlaceHolders++;
		if (pmOnly) numPlaceHolders+=2;
		else numPlaceHolders++;
		
		try{
			//write scores to .sgr file?
			if (printPointSgrs){
				if (pmOnly){
					outSumFile = new File (resultsFile+"Sum.sgr");
					outSum = new PrintWriter(new FileWriter(outSumFile));
				}
				outScoreFile = new File (resultsFile+extension+".sgr");
				outScore = new PrintWriter(new FileWriter(outScoreFile));
			}
			if (convertScoresToQValues) {
				tempRatioFile = new File(resultsFile+"tmpRatioFile");
				tempRatios = new PrintWriter(new FileWriter(tempRatioFile));
			}
			if (printOligoRatios){
				oligoFile = new File (resultsFile+"Oligo.sgr");
				outIndividualRatios = new PrintWriter(new FileWriter(oligoFile));
			}
			
			//for each chromosome	
			System.out.println();
			for (int i=1; i<info.size(); i+=4){
				ArrayList chromWindowAL = new ArrayList(10000);
				//initialize tester
				WilcoxonRankSumTest wt = new WilcoxonRankSumTest();
				wt.setTwoTailed(true);
				
				//fetch testWindows to scan, the int[window index][start, stop indexes] made in the WindowMaker
				chromosome = (String)info.get(i);
				File chromFile = new File(tpmapFile+chromosome+"Win");
				if (chromFile.exists() == false) continue;
				else System.out.println("Testing chromosome: "+chromosome);
				testWindows = (int[][])IO.fetchObject(chromFile);
				
				//fetch treatment and control normalized, transformed, and chromosome divided, float[replicas][oligoIntensities]
				treatmentIntensities = (float[][])IO.fetchObject(new File(treatmentDir+File.separator+chromosome+uniqueIdT));
				if (pmOnly) controlIntensities = (float[][])IO.fetchObject(new File(controlDir+File.separator+chromosome+uniqueIdC));
				
				//make chromosome wide relative differences, either layered or mean
				makeChromosomeWideScores();
				
				//fetch bp positions to use in converting window indexes to real base pairs
				positions = (int[])IO.fetchObject(new File(tpmapFile+chromosome));
				
				//print out individual oligo scores?
				if (printOligoRatios) printOligoRatios();
				
				//run the window tests
				numWindows = testWindows.length;
				for (int j=0; j<numWindows; j++){
					
					//get transformed values to test
					int[] startStop = testWindows[j];
					int windowSize = 1+startStop[1]-startStop[0];
					
					//score window 
					double score = scoreWindow(j);
					
					//wilcoxon rank sum test? positive or negative
					double wrs = 0;
					if (windowSize > 3 && pmOnly) wrs = wt.test(treatmentIntensities, controlIntensities, startStop[0], startStop[1]);
					
					//make Window?, need to make all if using Bourgon's Symmetric Null P-Test.
					if (windowSize >= minOligos  || convertScoresToQValues){
						//IF scores IS MODIFIED CORRECT SCORE_DESCRIPTION above and index of qvalue calculation below
						double[] scores = new double[numPlaceHolders];
						if (pmOnly == false) {
							scores[0] = score;
						}
						else{
							scores[0] = wrs;
							scores[1] = score;
							//set index of test window for later use
							if (randomizeLabels) scores[2] = j;
						}					
						//chromosome, bpStart, bpStop, numberOligos, scores
						Window win = new Window(chromosome, positions[startStop[0]],positions[startStop[1]], windowSize, scores);
						chromWindowAL.add(win);
						//save score for sym p test?
						if (convertScoresToQValues ) tempRatios.println(score);
					}
					
					//print sgr lines?
					if (printPointSgrs && windowSize >= minOligos){
						//print .sgr files for sum test and ave intensity value? centered in window
						double diff = sizeOfOligoMinusOne + positions[startStop[1]] - positions[startStop[0]];
						String pos = chromosome +"\t"+ ( (int)Math.round(diff/2.0) + positions[startStop[0]] )+"\t";
						if (printPointSgrs){
							if (pmOnly) outSum.println(pos+wrs);
							Double modScore;
							if (convertSummaries2Log2) modScore = new Double (Num.relativeDifferenceToLog2Ratio(score));
							else modScore = new Double (score);
							outScore.println(pos+ modScore.floatValue());
						}
					}
				}
			
			
				//run random label permutations for the chromosome, call last!!!
				if (randomizeLabels){
					//add permutated scores to histogram for this set of windows
					addPermutatedScoresToHistogram(chromWindowAL);
					System.out.println();
				}
				
				//add chromosome specific windows to total windows
				totalWindowsAL.addAll(chromWindowAL);
				
			}
			
			//close files 
			if (printPointSgrs){
				if (pmOnly) {
					outSum.close();
					IO.zipAndDelete(outSumFile);
				}
				outScore.close();
				IO.zipAndDelete(outScoreFile);
			}
			if (convertScoresToQValues ) tempRatios.close();
			if (printOligoRatios){
				outIndividualRatios.close();
				IO.zipAndDelete(oligoFile);
			}
			
			//delete Tmp files
			IO.deleteFiles(treatmentDir, uniqueIdT);
			if (pmOnly) IO.deleteFiles(controlDir, uniqueIdC);
			
			//convert ArrayList to Window[] 
			numWindows = totalWindowsAL.size();
			windows = new Window[numWindows];
			totalWindowsAL.toArray(windows);
			totalWindowsAL = null;
			System.out.println("\nNumber of windows screened: "+numWindows+"\n");
			
			//convert ratios to random permutation p-values
			if (randomizeLabels) {
				calculateRandomLabelPermutationConfidenceScores();
				//print point gr files?
				if (printPointSgrs){
					File sgr = new File( resultsFile+extension+ "RndPVal.sgr");
					makeSgrFile(windows, sgr, pValIndexForRandomPermutation);
					sgr = new File(resultsFile+extension+ "RndFDR.sgr");
					makeSgrFile(windows, sgr, pValIndexForRandomPermutation+1);
				}
			}
			
			//convert ratios to q values? and strip windows with too few oligos
			if (convertScoresToQValues ) {
				convertRatiosToPValues();
				//prune windows
				System.out.println("Pruning windows with with < "+minOligos+" oligo positions....");
				//instatiate new ArrayList with approximate number of expected windows
				totalWindowsAL = new ArrayList(numWindows/2);
				for (int z=0; z<numWindows; z++){
					if (windows[z].getNumberOligos()>=minOligos) totalWindowsAL.add(windows[z]);
				}
				numWindows = totalWindowsAL.size();
				windows = new Window[numWindows];
				totalWindowsAL.toArray(windows);
				totalWindowsAL = null;
			}
			
			//print window scores as heat map sgrs
			if (printHeatMapSgrs){
				System.out.println("Making window heat map sgrs.");
				WindowBlockMaker bm = new WindowBlockMaker(sizeOfOligoMinusOne +1);
				int index = 0;
				
				//Wilcoxon Rank Sum or if pmOnly == false then index 0 is actually the pmmm value
				bm.setScoreIndex(index++);
				if (pmOnly){
					bm.makeHeatMapSgrFile(windows, new File(resultsFile+"SumHM.sgr.zip"), true);
					//Ratio
					bm.setScoreIndex(index++);
				}
				bm.makeHeatMapSgrFile(windows, new File(resultsFile+extension+"HM.sgr.zip"), true);
				
				//Randomized PValues?
				if (randomizeLabels){
					//pvalues
					bm.setScoreIndex(index++);
					File sgr = new File(resultsFile+extension+ "RndPValHM.sgr.zip");
					bm.makeHeatMapSgrFile(windows, sgr, true);
					//fdr
					bm.setScoreIndex(index++);
					sgr = new File(resultsFile+extension+ "RndFDRHM.sgr.zip");
					bm.makeHeatMapSgrFile(windows, sgr, true);
				}
				
				//SymQValue is always last
				if (convertScoresToQValues ) {
					bm.setScoreIndex(windows[0].getScores().length-1);
					bm.makeHeatMapSgrFile(windows, new File(resultsFile+extension+"QValHM.sgr.zip"), true);
				}
			}

			//save Window[]
			if (numWindows != 0){
				File serWinIntArray = new File(resultsFile);
				System.out.print ("Saving "+numWindows+" windows "+ serWinIntArray+"\n");
				IO.saveObject(serWinIntArray, windows);
			}
			else System.out.println("No windows to save!");
			System.out.println("\nDone!");
			
		}catch (IOException e){
			e.printStackTrace();
			IO.deleteFiles(treatmentDir, uniqueIdT);
			if (pmOnly) IO.deleteFiles(controlDir, uniqueIdC);
		}
	}
	
	/**Uses random permutation of the cel file labels to estimate a null distribution, 
	 * store random scores in a histogram for p-value and FDR estimations.*/
	public void addPermutatedScoresToHistogram(ArrayList windowsAL){
		//for each permutation
		System.out.print("\t");
		//for each permutation compare the window score to the real and
		for (int i=0; i< numPermutations; i++){
			System.out.print((numPermutations-i)+" ");
			//shuffle intensities
			shuffleIntensities();
			//make chromosome wide ratios, layered, logged, etc.
			makeChromosomeWideScores();
			//get transformed values to test
			//for each window calculate the permuted window score 
			int numWindows = windowsAL.size();
			for (int j=0; j< numWindows; j++){
				Window win = (Window)windowsAL.get(j);
				//check to see if it is a good window
				if (win.getNumberOligos() >= minOligos){
					double[] scores = win.getScores();
					//find appropriate test window
					int testWindowIndex = (int)scores[pValIndexForRandomPermutation]; 
					//add score to histogram
					histogram.count(scoreWindow(testWindowIndex));
				}
			}
		}
		System.out.println();
	}
	
	/**Combines and shuffles the treatment and control intensities.*/
	public void shuffleIntensities(){
		//combine treat and control
		float[][] combine = new float[treatmentIntensities.length+ controlIntensities.length][];
		int index = 0;
		for (int i=0; i< treatmentIntensities.length; i++) combine[index++] = treatmentIntensities[i];
		for (int i=0; i< controlIntensities.length; i++) combine[index++] = controlIntensities[i];
		//shuffle
		Num.randomizeFirstArray(combine, System.currentTimeMillis());
		//split
		index = 0;
		for (int i=0; i< treatmentIntensities.length; i++) treatmentIntensities[i] = combine[index++];
		for (int i=0; i< controlIntensities.length; i++) controlIntensities[i] = combine[index++];
	}
	
	/**Takes the histogram of random permutation scores and calculates an uncorrected p-value and FDR for each real score.*/
	public void calculateRandomLabelPermutationConfidenceScores(){
		//check if hits found in max or min
		if (histogram.getNumLessThanMin() !=0 || histogram.getNumMoreThanMax() !=0){
			System.out.println("\nWarning! Histogram shows overrun, increase size.\n");
			histogram.printScaledHistogram();
			System.out.println("\nWarning! Histogram shows overrun, increase size.\n");
		}
		//calc p-val 
		double log10 = Math.log(10);
		for (int j=0; j<windows.length; j++){
			double[] scores = windows[j].getScores();
			if (windows[j].getNumberOligos() >= minOligos){
				double pval = histogram.pValue(scores[scoreIndexForHistogram]);
				//histo based p-value, make neg for reduced, pos for accum
				if (scores[scoreIndexForHistogram] >= 0) scores[pValIndexForRandomPermutation] = -10 * Math.log(pval)/log10;
				else scores[pValIndexForRandomPermutation] = 10 * Math.log(pval)/log10;
				
			}
			else scores[pValIndexForRandomPermutation] = 0;
		}
		//calc fdr, can get odd behaviour is some of the random regions have very high scores.
		//sort windows by score, smallest to largest
		for (int j=0; j<windows.length; j++){
			double[] scores = windows[j].getScores();
			windows[j].setSortBy(-1.0* scores[scoreIndexForHistogram]);
		}
		Arrays.sort(windows);
		//for each window calculate an fdr 
		double maxFdrRight = -100;
		double maxFdrLeft = 100;
		int indexMaxRight =0;
		int indexMaxLeft = 0;
		int fdrIndex = pValIndexForRandomPermutation+1; 
		ArrayList toSetMax = new ArrayList();
		
		//for each window
		for (int y=0; y< windows.length; y++){
			double[] windowScores = windows[y].getScores();
			//correct number of oligos?
			if (windows[y].getNumberOligos() >= minOligos){
				
				//find number of mock hits and real hits 
				double numMockWindows;
				//positive score, go right side
				if (windowScores[scoreIndexForHistogram] >= 0.0){ 
					numMockWindows = histogram.numberBinHitsToRightAndIncludingValue( windowScores[scoreIndexForHistogram] ) / numPermutations;
				}
				//negative score go left side
				else numMockWindows = histogram.numberBinHitsToLeftAndIncludingValue( windowScores[scoreIndexForHistogram] ) / numPermutations;
				
				//any mock windows?
				if (numMockWindows != 0 ) {
					
					//count number of real windows
					double numRealWindows;
					if (windowScores[scoreIndexForHistogram] >= 0.0) numRealWindows = windows.length - y; 
					else numRealWindows = 1+ y - 0;
					
					//calc fdr = num mock/ num real
					double fdr = numMockWindows/numRealWindows;
					//do number of mock exceed number of real?
					if (fdr > 1) windowScores[fdrIndex] = 0;
					else {
						//set transformed fdr in window scores
						if (windowScores[scoreIndexForHistogram] >= 0.0) {
							windowScores[fdrIndex] = (-10.0) * (Math.log(fdr)/log10);
							//new right side max
							if (windowScores[fdrIndex] > maxFdrRight) {
								maxFdrRight = windowScores[fdrIndex];
								indexMaxRight = y;
							}
						}
						else {
							windowScores[fdrIndex] = (10.0) * (Math.log(fdr)/log10);
							//new left side max
							if (windowScores[fdrIndex] < maxFdrLeft) {
								maxFdrLeft = windowScores[fdrIndex];
								indexMaxLeft = y;
							}
						}
					}
				}
				//save window index to set max FDR
				else {
					toSetMax.add(new Integer(y));
				}
			}
			//too few oligos set fdr to zero
			else{
				windowScores[fdrIndex] = 0;
			}
		}
		
		//set max FDR
		int length = toSetMax.size();
		for (int x=0; x < length; x++){
			int index= ((Integer)toSetMax.get(x)).intValue();
			double[] scores = windows[index].getScores();
			if (scores[scoreIndexForHistogram] >=0 ) scores[fdrIndex] = maxFdrRight;
			else scores[fdrIndex] = maxFdrLeft;
		}
		
		/*Altering the FDR here, in some cases where some random permutations 
		have very high scores the fdr becomes funky at higher real scores and will drop 
		 despite the increased confidence, thus to "fix" resetting all fdrs after max
		 to the max*/
		//set max left side
		for (int i=0; i< indexMaxLeft; i++){
			double[] scores = windows[i].getScores();
			scores[fdrIndex] = maxFdrLeft;
		}
		//set max right side
		for (int i=indexMaxRight+1; i< windows.length; i++){
			double[] scores = windows[i].getScores();
			scores[fdrIndex] = maxFdrRight;
		}
		
		//resort windows by chrom and position
		Arrays.sort(windows, new WindowComparator());
	} 

	
	/**Scores a window using a pseudo median or trimmed mean on oligo relative differences.
	 * If the number of relative differences exceed the maxScores then they are randomly 
	 * sampled to the number of maxScores prior to running the pseudo median.
	 * Returns 0 if using pseudoMedian and window size is < 2 or if using trm meand and window size < 3.
	 * @return trimmed mean average or pseudo median layered relative difference score -2 to 2.*/
	public double scoreWindow(int testWindowIndex){
		int[] startStop = testWindows[testWindowIndex];
		int windowSize = 1+startStop[1]-startStop[0];
		float[] windowRelDiffs = null;
		int counter = 0;
		double score = 0;
		
		//pseudo median on layered relative differences?
		if (usePseudoMedian && windowSize > 1) {
			//using averaged treatment and averaged control ratios?
			if (averageGroups){
				windowRelDiffs = new float[windowSize];
				System.arraycopy(singleScores, startStop[0], windowRelDiffs, 0, windowSize);
			}
			//no, using all layer
			else {
				windowRelDiffs = new float[windowSize * layeredScores[0].length];
				int size = layeredScores[0].length;
				for (int k=startStop[0]; k<= startStop[1]; k++){
					System.arraycopy(layeredScores[k], 0, windowRelDiffs, counter, size);
					counter += size; 
				}
				//watch out for too big arrays
				if (windowRelDiffs.length > maxScores) {
					System.out.print(" Rndm Smpling PseMed:"+windowRelDiffs.length);
					windowRelDiffs = Num.randomSample(windowRelDiffs, maxScores);
				}
			}
			score = Num.pseudoMedian(windowRelDiffs);
		}
		//trimmed mean relative difference (-2 to 2)?
		else if (windowSize > 2){
			//collect scores
			windowRelDiffs = new float[windowSize];
			System.arraycopy(singleScores, startStop[0], windowRelDiffs, 0, windowSize);
			Arrays.sort(windowRelDiffs);
			//trim 10% or min of 1 off ends
			int trimNumber = 1;
			if (windowSize > 10) {
				double num = 0.1 * ((double)windowSize);
				trimNumber = (int)num;
			}
			score = Num.trimmedMean(windowRelDiffs,trimNumber);
		}
		return score;
	}
	
	public void convertRatiosToPValues () {
		try {
			System.out.println("Launching Richard Bourgon's symmetric null p-value estimator...");
			//run Richards script
			String[] command = {
					fullPathToPValue,
					tempRatioFile.getCanonicalPath(),
					"50",	//number of bootstraps
					""+numWindows
			};
			
			System.out.println("\t"+Misc.stringArrayToString(command," "));
			
			String[] results = IO.executeCommandLine(command);
			//check that the correct number of result lines were found
			if (results.length != numWindows){
				Misc.printExit("\nProblem with sym p test! The number of result lines ("+results.length+
						") from the sym p test does not match the number of windows ("+numWindows+")?!\n");
			}
			//write output to a new tmpFile after converting to real pValues
			File tempSymPValueFile = new File(resultsFile+"tmpFileSymPValueFile");
			PrintWriter out = new PrintWriter(new FileWriter(tempSymPValueFile));
			for (int z=0; z<numWindows; z++){
				double score = Double.parseDouble(results[z]);
				if (score == 1) score = 0.999999;
				score = 1- score;
				out.println(score);
			}
			out.close();
			//run John Storeys qvalue from R, saving results in tempRatioFile, overwrites?
			System.out.println("Launching John Storey's q-value app in R...");
			//write temp .R file
			StringBuffer rScript =new StringBuffer("library(qvalue)\np <- scan (\"");
			rScript.append(tempSymPValueFile.getCanonicalPath());
			rScript.append("\")\nqobj <- qvalue (p)\nqwrite (qobj, filename=\"");
			rScript.append(tempRatioFile.getCanonicalPath());
			rScript.append("\")\n");
			File tempScriptFile = new File(resultsFile+"tmpFileRScriptFile");
			IO.writeString(rScript.toString(), tempScriptFile);
			//make command
			command = new String[] {
					fullPathToR,
					"CMD",
					"BATCH",
					tempScriptFile.getCanonicalPath()};			
			//execute
			IO.executeCommandLine(command);
			System.out.println("\t"+Misc.stringArrayToString(command," "));
			//assign results to windows, also write an sgr file for q-values
			System.out.println("Assigning q-value scores...");			
			importScores(windows, tempRatioFile);
			//kill temp files
			tempRatioFile.delete();
			tempScriptFile.delete();
			tempSymPValueFile.delete();
		} catch (IOException e){
			e.printStackTrace();
		}
	}
	
	/**Prints log2 oligo ratios.*/
	public void printOligoRatios(){
		double[] meanT = Num.averageFloatArraysFlipped(treatmentIntensities);
		int numberOligoMeasurements = meanT.length;
		double[] meanC = null;
		if (pmOnly) {
			meanC = Num.averageFloatArraysFlipped(controlIntensities);
			//for every measurement
			for (int j=0; j<numberOligoMeasurements; j++){
				float ratio = new Double (Math.log(meanT[j]/meanC[j])/ Num.log2).floatValue();
				//print
				outIndividualRatios.println(chromosome+"\t"+(positions[j]+halfOligoLength)+"\t"+ratio);
			}
		}
		else {
			//for every measurement
			for (int j=0; j<numberOligoMeasurements; j++){
				float ratio = new Double (Math.log(meanT[j])/ Num.log2).floatValue();
				//print
				outIndividualRatios.println(chromosome+"\t"+(positions[j]+halfOligoLength)+"\t"+ratio);
			}
		}
		meanC = null;
		meanT = null;
	}
	
	/**Makes chromosome wide scores.*/
	public void makeChromosomeWideScores(){
		//make chromosome wide layers float[oligoIndex][relative differences] for pseudoMedian?
		if (usePseudoMedian) {
			if (pmOnly) {
				if (averageGroups) singleScores= Num.relativeDifferences(treatmentIntensities, controlIntensities);
				else layeredScores= Num.layeredRelativeDifferencesSeperate(treatmentIntensities, controlIntensities);
			}
			else layeredScores = treatmentIntensities;
		}
		//make chromosome wide oligo relative differences
		else {
			if (pmOnly) {
				singleScores = Num.relativeDifferences(treatmentIntensities, controlIntensities);
			}
			else singleScores = Num.averageFloatArraysFlippedToFloats(treatmentIntensities);
		}
	}
	
	/**Reads in a file of scores from a qvalue output file
	 * Also writes and zips an .sgr file for the qValues.
	 * Saves score in last position of the window.scores[].*/
	public void importScores(Window[] windows, File importFile){
		int counter = 0;
		try {
			BufferedReader in = new BufferedReader( new FileReader (importFile));
			PrintWriter outQValues = null;
			File pointSgrQValuesFile = null;
			if (printPointSgrs){
				pointSgrQValuesFile  = new File(Misc.replaceEnd(IO.getFullPathName(importFile),"tmpRatioFile", extension+"QVal.sgr"));
				outQValues = new PrintWriter( new FileWriter (pointSgrQValuesFile));
			}
			String line;
			String[] tokens;
			double log10 = Math.log(10.0);
			//skip first three lines
			in.readLine();in.readLine();in.readLine();
			//start
			while ((line= in.readLine()) != null){
				tokens = line.split("\\s+");
				if (tokens.length!=2) continue;
				double score = Double.parseDouble(tokens[1]);
				score = (-10.0) * (Math.log(score)/log10);
				double[] scores= windows[counter].getScores();
				scores[scores.length-1] = score;
				windows[counter].setScores(scores);
				//make sgr line?
				if (printPointSgrs){
					double diff = sizeOfOligoMinusOne + windows[counter].getStartLastOligo() - windows[counter].getStart1stOligo();
					String pos = windows[counter].getChromosome() +"\t"+ ( (int)Math.round(diff/2.0) + windows[counter].getStart1stOligo() )+"\t";
					outQValues.println(pos+score);
				}
				counter++;
			}
			in.close();
			
			if (printPointSgrs){
				outQValues.close();
				IO.zipAndDelete(pointSgrQValuesFile);
			}
		} catch (Exception e){
			e.printStackTrace();
		}
		//check counter
		if (counter != windows.length){
			System.out.println("\nProblem with QValue R script: the number of lines in the import file ("+counter+
					") do not match the number of Windows ("+windows.length+")!\n");
			System.exit(1);
		}
	}
	
	/**Given a directory (or a file), fetches serialized float[] from all files with the .celp extension.
	 * Will randomize the float[]s if so indicated.*/
	public static float[][] fetchFloatArrays(String dir, boolean randomize){
		File[] files = IO.extractFiles(new File(dir),".celp");
		int numFiles = files.length;		
		float[][] inties = new float[numFiles][];
		//force a different seed
		long time = System.currentTimeMillis();
		for (int i=0; i<numFiles; i++){
			inties[i] = (float[])IO.fetchObject(files[i]);
			if (randomize) Num.randomize(inties[i], time+i);
		}
		return inties;
	}
	
	public static void main(String[] args) {
		if (args.length<4){
			printDocs();
			System.exit(0);
		}	
		new ScanChip(args);
	}		
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': treatmentDir = args[i+1]; i++; break;
					case 'c': controlDir = args[i+1]; i++; break;
					case 'r': resultsFile = args[i+1]; i++; break;
					case 's': fullPathToPValue = args[i+1]; i++; convertScoresToQValues = true;break;
					case 'q': fullPathToR = args[i+1]; i++; convertScoresToQValues = true; break;
					case 'p': convertScoresToQValues = true; break;
					case 'a': usePseudoMedian = true; break;
					case 'y': usePseudoMedian = true; averageGroups = true;
					case 'i': printOligoRatios = true; break;
					case 'd': printHeatMapSgrs = true; break;
					case 'e': printPointSgrs = true; break;
					case 'l': randomizeLabels = true;break;
					case 'f': convertSummaries2Log2 = true;break;
					case 'n': numPermutations = Integer.parseInt(args[i+1]); i++; randomizeLabels = true; break;
					case 'm': minOligos = Integer.parseInt(args[i+1]); i++; break;
					case 'z': sizeOfOligoMinusOne =Integer.parseInt(args[i+1])-1;i++; break;
					case 'x': maxScores = Integer.parseInt(args[i+1]); i++; break;
					case 'b': tpmapFile = new File(args[i+1],"tpmap.fa").toString(); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check if directories are directories
		File treat = new File(treatmentDir);
		if (treat==null || treat.isDirectory()==false){
			Misc.printExit("\nCheck your treatment directory, it can't be found or read!\n");
		}	
		//look for control directory
		if (controlDir == null) pmOnly = false;
		
		infoFile = tpmapFile+"Info";
		if (new File(infoFile).exists()==false){
			Misc.printExit("\nCheck your TPMapFiles directory, it appears to be incomplete!\n");
		}		
		
		halfOligoLength = (int)Math.round(((double)sizeOfOligoMinusOne)/2.0);
	}	
	
	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Scan Chip: March 2007                              **\n" +
				"**************************************************************************************\n" +
				"SC uses a sliding window to calculate window level statistics on processed cel\n" +
				"files (xxx.celp) using a Wilcoxon Rank Sum test and either a trimmed mean (10%) or\n" +
				"pseudo median on the windowed relative differences. To estimate confidence, SC wraps\n" +
				"Bourgon's Symmetric Null P-value estimator and John Storey's QValue R package. SC can\n" +
				"also use random label permutation to estimate a p-value and FDR. A variety of sgr\n" +
				"files can be saved for direct viewing in Affymetrix' IGB. xxxHM.sgr should be viewed\n" +
				"as heat maps or stair-step graphs. Negative scores/ p-values/ FDRs indicate reduction,\n" +
				"positive accumulation.\n\n" +
				
				"Note, the order of scores associated with each window are:" +
				"\n      [0]= "+SCORE_DESCRIPTION[0]+
				"\n      [1]= "+SCORE_DESCRIPTION[1]+ 
				"\n      [2]= "+SCORE_DESCRIPTION[2]+
				"\n      [3]= "+SCORE_DESCRIPTION[3]+
				"\n      [4]= "+SCORE_DESCRIPTION[4]+
				"\n\n"+
				
				"-b The full path 'TPMapFiles' directory text generated by the TPMapProcessor program\n" +
				"-r The full path file text to use in saving the results\n" +
				"-t The full path directory text containing the serialized float[] treatment file(s)\n" +
				"-c The full path directory text containing the serialized flaot[] control file(s)\n" +
				"-m Minimal number of oligo positions required in each window. Defaults to 10. \n" +
				"-z Size of oligo, defaults to 25.\n"+
				"-a Use a pseudo median of all pair relative differences instead of a trimmed mean.\n" +
				"       This is very slow but much more robust.\n"+
				"-x Maximum number of window scores before random sampling for pseudo median\n" +
				"       calculation, defaults to 1000. (#T x #C x #oligos in window).\n"+
				"-y Use average treatment/ average control relative differences instead of all layer in\n"+
				"       pseudo median calculation.\n"+
				"-p Convert window scores to q-value FDRs, -10log10(multiple test corr pval).\n" +
				"       Requires config R, and possibly compiling the SymPTest.\n"+
				"-s The full path to the symmetric_p_test, defaults to\n" +
				"      '"+fullPathToPValue+"'\n"+
				"-q The full path to qvalue loaded R, defaults to\n"+
				"      '"+fullPathToR+"'\n"+
				"-l Use random permutation of the chip labels to estimate a p-value and FDR for the\n" +
				"       window scores.\n"+
				"-n Number of random permutations, defaults to 10.\n"+
				"-i Print an individual oligo log2 ratio (aveT/aveC) sgr file.\n"+
				"-d Print heat map window summary sgr files.\n"+
				"-e Print point window summary sgr files.\n"+
				//"-f Print log2 ratios instead of relative differences for point and heat map window\n"+
				//"       summaries.\n"+
				
				"\nExample: java -Xmx1500M -jar pathTo/T2/Apps/ScanChip -b /affy/TPMapFiles -r\n" +
				"      /affy/res/zeste.res -t /affy/tCels -c /affy/cCels -m 5 -i -e\n\n" +
				
		"**************************************************************************************\n");		
	}
	public static String[] getSCORE_DESCRIPTION() {
		return SCORE_DESCRIPTION;
	}
	
	public void buildExtensionSetThreshold(){
		//which test?
		System.out.print("Using ");
		if (usePseudoMedian) {
			if (averageGroups) System.out.print("averaged ");
			else System.out.print("all layer ");
			System.out.print("replica pseudo median oligo relative differences.");
			extension = "Pse";
		}
		else {
			extension = "Trm";
			System.out.print("a trimmed mean on average replica oligo relative differences.");
		}
		
		//set params based on PMMM or just PM
		//PMMM?
		if (pmOnly == false) {
			extension = extension + "PmMm";
			convertScoresToQValues = false;
			System.out.print(" using max(PM-MM,1)) transformation");
		}
		System.out.println();
	}
	
	/**Makes a xxx.sgr.zip file from a sorted Window[].*/
	public void makeSgrFile(Window[] windows, File sgrFile, int scoreIndex){
		try{
			PrintWriter out = new PrintWriter(new FileWriter(sgrFile));
			for (int i=0; i< windows.length; i++){
					double diff = sizeOfOligoMinusOne + windows[i].getStartLastOligo() - windows[i].getStart1stOligo();
					String pos = windows[i].getChromosome() + "\t" +( (int)Math.round(diff/2.0) + windows[i].getStart1stOligo() )+"\t";
					float trunkScore = new Double(windows[i].getScores()[scoreIndex]).floatValue();
					out.println(pos+trunkScore);
			}
			out.close();
			IO.zipAndDelete(sgrFile);	
		} catch (IOException e){
			e.printStackTrace();
		}
	}
	
}
