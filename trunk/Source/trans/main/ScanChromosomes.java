package trans.main;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import edu.utah.seq.parsers.*;
import trans.misc.Util;
import trans.tpmap.*;
import util.gen.*;

/**
 * Scores windows of intensity values using several tests. 
 */
public class ScanChromosomes {

	//fields
	private File resultsDirectory;
	private File[] treatmentDirectories = null;
	private File[] controlDirectories = null;
	private File[] chromosomeOligoPositions = null;
	private Window[] windows;
	private ArrayList totalWindowsAL = new ArrayList(10000);
	private int sizeOfOligoMinusOne = 24;
	private int halfOligoLength;
	private boolean convertScoresToQValues = false;
	private boolean pmOnly = true;
	private String genomeVersion;
	private static String fullPathToPValue = "/home/BioApps/T2/OthersCode/RBSymPTest/symmetric_p_test";	//old version
	private static String fullPathToR = "/usr/local/R/bin/R"; //must have the qvalue package loaded!
	private boolean randomizeData = false;
	private boolean randomizeLabels = false;
	private boolean randomizeIntensities = false;
	private int numPermutations = 10;
	private boolean printOligoRatios = false;
	private boolean printPointBars = false;
	private boolean printHeatMapBars = false;
	private boolean convertRelDiff2Log2 = false;
	private boolean usePseudoMedian = false;
	private double minWindowFilter = 0;
	private double maxWindowFilter = 0;
	private int minOligos = 10;
	private int maxScores = 1000;
	private int windowSize = 675;
	private float[][] treatmentIntensities;
	private float[][] controlIntensities;
	private String extension = "";
	private boolean averageGroups = false;
	private float singleScores[] = null;
	private float[][] layeredScores = null;
	private String chromosome = null;
	private String[] particularChromosomes;
	private int[] positions;
	private int[][] testWindows;
	private PrintWriter tempRatios = null;
	private File tempRatioFile = null;
	private File sumDirectory;
	private File relDiffDirectory;
	private File oligoDirectory;
	private File windowDirectory;
	private String randomWord = null;
	private String split = null;
	private int splitUnit = 0;
	private int splitBy = 0;
	private String strand = "."; //default or + or -
	//histogram
	private Histogram histogram; 
	private int scoreIndexForHistogram = 1;
	private int pValIndexForRandomPermutation = 2;

	/**Array describing the scores in the Window.Scores.  It changes frequently!*/
	private static String[] SCORE_DESCRIPTION = {
		"Wilcoxon Rank Sum, -10Log10(p-val)", 
		"Trimmed Mean or Layered Pseudo Median Relative Difference",
		"-10Log10(uncorrected p-value) based on rnd perm of cel file labels (if set)",
		"-10Log10(FDR) based on random permutation  (if set)",
		"-10Log10(q-val FDR) based on symmetric null distribution, (if set)"};

	public ScanChromosomes(String[] args){
		//process args
		processArgs(args);

		//build extension 
		System.out.println();
		buildExtensionSetThreshold();

		//randomizing data?
		if (randomizeData) System.out.println("Warning, randomizing data, all significant regions are now false positives.");

		//calculate number of additional score place holders needed
		int numPlaceHolders = 0;
		if (randomizeLabels || randomizeIntensities) {
			numPlaceHolders+=2;
			//instantiate histogram for random permutations
			histogram = new Histogram(-2,2,5000);
		}
		if (convertScoresToQValues) numPlaceHolders++;
		if (pmOnly) numPlaceHolders+=2;
		else numPlaceHolders++;


		//make directories to hold point bar files and temp scores when so indicated
		makeDirectories();

		//make a WindowMaker to make windows for each chromosome from oligo positions
		WindowMaker windowMaker;
		if (convertScoresToQValues) windowMaker = new WindowMaker(windowSize, 1);
		else windowMaker = new WindowMaker(windowSize, minOligos);

		//for each chromosome	
		System.out.println();
		for (int i=0; i<chromosomeOligoPositions.length; i++){

			chromosome = chromosomeOligoPositions[i].getName();
			System.out.println("Testing chromosome: "+chromosome);			

			//get positions
			positions = (int[])IO.fetchObject(chromosomeOligoPositions[i]);
			//more than 1?
			if (positions.length < 2){
				System.out.println("Too few oligos found, skipping.");
				continue;
			}

			//make windows int[windowNumber][start index, stop index]
			testWindows = windowMaker.makeWindows(positions);
			int numTestWindows = testWindows.length;
			if (numTestWindows == 0) {
				System.out.println("No Windows found, skipping.");
				continue;
			}

			//make array list to hold chromosome specific windows
			ArrayList chromWindowsAL = new ArrayList(numTestWindows/2);

			//get treatment intensities float[replicas][values]
			treatmentIntensities = fetchIntensities (treatmentDirectories, chromosome, randomizeData);
			if (treatmentIntensities == null) Misc.printExit("\nProblem with fetching treatment intensities. Any missing or extra files?\n");
			if (checkSizes(treatmentIntensities) == false) Misc.printExit("\nYour treatment intensity arrays are of different lengths?!\n");

			//get control intensities float[replicas][values]
			if (pmOnly) {
				controlIntensities = fetchIntensities (controlDirectories, chromosome, randomizeData);
				if (controlIntensities == null) Misc.printExit("\nProblem with fetching control intensities. Any missing or extra files?\n");
				if (checkSizes(controlIntensities) == false) Misc.printExit("\nYour control intensity arrays are of different lengths?!\n");
				if (treatmentIntensities[0].length != controlIntensities[0].length) Misc.printExit("\nYour treatment and control intensity arrays are of different lengths?!\n");
			}

			//quick fix for messed up sex?
			/*
			if (chromosome.equals("chrX") || chromosome.equals("chrY")){
				System.out.println("\tMedian scaling to 100");
				treatmentIntensities = Num.medianNormalize(treatmentIntensities, 100);
				controlIntensities = Num.medianNormalize(controlIntensities, 100);
			}
			 */
			//initialize tester
			WilcoxonRankSumTest wt = new WilcoxonRankSumTest();
			wt.setTwoTailed(true);

			//make chromosome wide scores for window scanning.
			makeChromosomeWideScores();

			//print out individual oligo scores?
			if (printOligoRatios) printOligoRatios();

			//process split?
			int start = 0;
			int stop = numTestWindows;
			if (split!=null){
				double numToProc =  ((double)numTestWindows)/ ((double)splitBy);
				start = (int)Math.round((numToProc * ((double)splitUnit)) - numToProc);
				stop = (int)Math.round(numToProc * ((double)splitUnit));
			}

			//make bar summary arrays and Files
			File sumBarFile = null;
			File relDiffBarFile = null;
			int[] chrPos = null;
			float[] chrSumVal = null;
			float[] chrScoreVal = null;
			int barIndex =0;
			if (printPointBars){
				int num = stop - start;
				if (pmOnly) {
					chrSumVal = new float[num];
					sumBarFile = new File (sumDirectory, chromosome+".bar");
				}
				relDiffBarFile = new File (relDiffDirectory, chromosome+".bar");
				chrScoreVal = new float[num];
				chrPos = new int[num];
			}

			//run the window tests
			for (int j=start; j<stop; j++){

				//get transformed values to test
				int[] startStop = testWindows[j];
				int windowSize = 1+startStop[1]-startStop[0];
				//score window 
				double score = scoreWindow(j);

				//wilcoxon rank sum test?
				double wrs = 0;
				if (windowSize > 3 && pmOnly) {
					wrs = wt.test(treatmentIntensities, controlIntensities, startStop[0], startStop[1]);
					if (score < 0) wrs = wrs * -1;
				}

				//make Window?, need to make all if converting to PValues.
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
						if (randomizeLabels || randomizeIntensities) scores[2] = j;
					}					
					//chromosome, bpStart, bpStop, numberOligos, scores
					Window win = new Window(chromosome, positions[startStop[0]],positions[startStop[1]], windowSize, scores);
					chromWindowsAL.add(win);
				}

				//save scores?
				if (convertScoresToQValues || printPointBars){

					double diff = sizeOfOligoMinusOne + positions[startStop[1]] - positions[startStop[0]];
					if (convertScoresToQValues ) tempRatios.println(score);
					//print .bar files for sum test and ave intensity value? centered in window
					if (printPointBars){
						if (pmOnly) chrSumVal[barIndex] = new Double (wrs).floatValue();
						//convert score to log2?
						if (convertRelDiff2Log2) chrScoreVal[barIndex] = new Double (Num.relativeDifferenceToLog2Ratio(score)).floatValue();
						else chrScoreVal[barIndex] = new Double (score).floatValue();
						chrPos[barIndex] = (int)Math.round(diff/2.0) + positions[startStop[0]];
						barIndex++;
					}
				}
			}

			if (randomizeLabels || randomizeIntensities){
				//add permutated scores to histogram for this set of windows
				addPermutatedScoresToHistogram(chromWindowsAL);
			}

			//write bar files?
			if (printPointBars){
				BarParser bp = new BarParser();
				bp.setZipCompress(true);
				HashMap<String,String> tagVals = new HashMap<String,String>();
				tagVals.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
				if (pmOnly)	bp.writeBarFile(sumBarFile, chromosome, genomeVersion, strand.charAt(0), chrPos, chrSumVal, tagVals);
				bp.writeBarFile(relDiffBarFile, chromosome, genomeVersion, strand.charAt(0), chrPos, chrScoreVal, tagVals);
			}

			//add chromosome specific windows to total windows
			totalWindowsAL.addAll(chromWindowsAL);
		}

		//End chromosome specific scanning
		//Finish window processing

		//convert ArrayList to Window[] 
		int numWindows = totalWindowsAL.size();
		windows = new Window[numWindows];
		totalWindowsAL.toArray(windows);
		totalWindowsAL = null;

		//convert scores to random permutation p-values
		if (randomizeLabels || randomizeIntensities) {
			if (randomizeLabels) System.out.println("Calculating random label permutation p-values and FDRs...");
			else System.out.println("Calculating random intensity position permutation p-values and FDRs...");
			calculateRandomPermutationConfidenceScores();
			//print point gr files?
			if (printPointBars){
				File grDir = new File( resultsDirectory, "RndPVal");
				makeBarFiles(windows, grDir, 2, false);
				grDir = new File( resultsDirectory, "RndFDR");
				makeBarFiles(windows, grDir, 3, false);
			}
		}

		//convert scores to q values?
		if (convertScoresToQValues ) {
			tempRatios.close();
			System.out.println("Launching Richard Bourgon's symmetric null p-value estimator...");
			convertRatiosToPValues();
			//prune windows
			System.out.println("Pruning windows with with < "+minOligos+" oligo positions...");
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

		//save Window[]? This should preceed printHeatMapBars since the window rel diff score could be changed to log2
		if (numWindows !=0 ){
			File serWinIntArray;
			if (particularChromosomes == null) serWinIntArray = new File(windowDirectory,"all_Win");
			else {
				String concat = Misc.stringArrayToString(particularChromosomes,"_").replaceAll("chr","");
				serWinIntArray = new File(windowDirectory,concat+"_Win");
			}
			//convert relative difference to log2?
			if (convertRelDiff2Log2){
				System.out.println("Converting relative differences to log2 ratios.");
				for (int k=0; k< windows.length; k++){
					double [] scores = windows[k].getScores();
					scores[scoreIndexForHistogram] = Num.relativeDifferenceToLog2Ratio(scores[scoreIndexForHistogram]);
				}
			}
			System.out.println("Saving "+numWindows+" windows "+ serWinIntArray);
			IO.saveObject(serWinIntArray, windows);
		}
		else System.out.println("No windows to save!");

		//print window scores as heat map bar files?
		if (printHeatMapBars){
			System.out.println("Making window heat map bar files...");
			WindowBlockMakerTwoColor bm = new WindowBlockMakerTwoColor(sizeOfOligoMinusOne +1, genomeVersion, strand, minWindowFilter, maxWindowFilter);
			int index = 0;
			File barDirectory;

			//Wilcoxon Rank Sum or if pmOnly == false then index 0 is actually the pmmm value
			if (pmOnly){
				barDirectory = new File(resultsDirectory, "SumHM");
				if (barDirectory.exists() == false) barDirectory.mkdir();
				bm.makeHeatMapBarFiles(windows, barDirectory, index);
				index++;
			}
			barDirectory = new File(resultsDirectory, extension+"HM");
			if (barDirectory.exists() == false) barDirectory.mkdir();
			bm.makeHeatMapBarFiles(windows, barDirectory, index);
			index++;

			//Randomized scores?
			if (randomizeLabels || randomizeIntensities){
				//pvalues
				barDirectory = new File(resultsDirectory, "RndPValHM");
				if (barDirectory.exists() == false) barDirectory.mkdir();
				bm.makeHeatMapBarFiles(windows, barDirectory, index);
				index++;
				//fdr
				barDirectory = new File(resultsDirectory, "RndFDRHM");
				if (barDirectory.exists() == false) barDirectory.mkdir();
				bm.makeHeatMapBarFiles(windows, barDirectory, index);
				index++;
			}

			//SymQValue is always last
			if (convertScoresToQValues ) {
				index = windows[0].getScores().length-1;
				barDirectory = new File(resultsDirectory, "QValHM");
				if (barDirectory.exists() == false) barDirectory.mkdir();
				bm.makeHeatMapBarFiles(windows, barDirectory, index);
			}
		}
		System.out.println("\nDone!");
	}
	/**Takes the histogram of random permutation scores and calculates a p-value and FDR for each real score.*/
	public void calculateRandomPermutationConfidenceScores(){
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
		double maxFdrRight = 0;
		double maxFdrLeft = 0;
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
		System.out.println("\tMaxFDRRight "+maxFdrRight+"  MaxFDRLeft "+maxFdrLeft);
		//set max FDR
		int length = toSetMax.size();
		for (int x=0; x < length; x++){
			int index= ((Integer)toSetMax.get(x)).intValue();
			double[] scores = windows[index].getScores();
			if (scores[scoreIndexForHistogram] >=0 ) scores[fdrIndex] = maxFdrRight;
			else scores[fdrIndex] = maxFdrLeft;
		}

		/**Altering the FDR here, in some cases where some random permutations 
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

	/**Uses random permutation of the cel file labels to estimate a null distribution, 
	 * store random scores in a histogram for p-value and FDR estimations.*/
	public void addPermutatedScoresToHistogram(ArrayList windowsAL){
		//for each permutation
		System.out.print("\t");
		//for each permutation compare the window score to the real and
		for (int i=0; i< numPermutations; i++){
			System.out.print((numPermutations-i)+" ");
			//shuffle intensities by label or position
			if (randomizeLabels) shuffleIntensityLabels();
			else shuffleIntensityPositions();
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

	/**Shuffles the treatment and control intensities (seperately).*/
	public void shuffleIntensityPositions(){
		//float[replicas][values]
		Num.randomize(treatmentIntensities);
		Num.randomize(controlIntensities);
	}

	/**Combines and shuffles the treatment and control intensities.*/
	public void shuffleIntensityLabels(){
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

	/**Scores a window using a pseudo median or trimmed mean on oligo relative differences.
	 * If the number of relative differences exceed the maxScores then they are randomly 
	 * sampled to the number of maxScores prior to running the pseudo median.
	 * Returns 0 if using trm mean and window size < 3.
	 * @return trimmed mean average or pseudo median layered relative difference score -2 to 2.*/
	public double scoreWindow(int testWindowIndex){
		int[] startStop = testWindows[testWindowIndex];
		int windowSize = 1+startStop[1]-startStop[0];
		float[] windowRelDiffs = null;
		int counter = 0;
		double score = 0;

		//pseudo median on layered relative differences?
		if (usePseudoMedian) {
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
		else {
			//collect scores
			windowRelDiffs = new float[windowSize];
			System.arraycopy(singleScores, startStop[0], windowRelDiffs, 0, windowSize);
			if (windowSize > 2){
				Arrays.sort(windowRelDiffs);
				//trim 10% or min of 1 off ends
				int trimNumber = 1;
				if (windowSize > 10) {
					double num = 0.1 * ((double)windowSize);
					trimNumber = (int)num;
				}
				score = Num.trimmedMean(windowRelDiffs,trimNumber);
			}
			else score = Num.mean(windowRelDiffs);
		}
		return score;
	}

	/**Makes xxx.bar.zip files from a Window[], one per chromosome.*/
	public void makeBarFiles(Window[] windows, File directory, int scoreIndex, boolean stairStep){

		if (directory.exists() == false) directory.mkdir();
		String chrom = windows[0].getChromosome();
		File gr = new File(directory, chrom+".bar");

		ArrayList positionsAL = new ArrayList();
		ArrayList valuesAL = new ArrayList();
		for (int i=0; i< windows.length; i++){
			//same chromosome?
			if (chrom.equals(windows[i].getChromosome())){
				double diff = sizeOfOligoMinusOne + windows[i].getStartLastOligo() - windows[i].getStart1stOligo();
				int pos = (int)Math.round(diff/2.0) + windows[i].getStart1stOligo();
				float trunkScore = new Double(windows[i].getScores()[scoreIndex]).floatValue();
				positionsAL.add(new Integer (pos));
				valuesAL.add(new Float(trunkScore));
			}
			else{
				//close old
				int[] positions = Num.arrayListOfIntegerToInts(positionsAL);
				positionsAL.clear();
				float[] values = Num.arrayListOfFloatToArray(valuesAL);
				valuesAL.clear();
				Util.writeSimpleBarFile(chrom, genomeVersion, strand, positions, values, gr);

				//start new
				chrom = windows[i].getChromosome();
				gr = new File(directory, chrom+".bar");

				//print gr line
				double diff = sizeOfOligoMinusOne + windows[i].getStartLastOligo() - windows[i].getStart1stOligo();
				int pos = (int)Math.round(diff/2.0) + windows[i].getStart1stOligo();
				float trunkScore = new Double(windows[i].getScores()[scoreIndex]).floatValue();
				positionsAL.add(new Integer (pos));
				valuesAL.add(new Float(trunkScore));		
			}
		}
		//write last
		int[] positions = Num.arrayListOfIntegerToInts(positionsAL);
		positionsAL = null;
		float[] values = Num.arrayListOfFloatToArray(valuesAL);
		valuesAL = null;


		BarParser bp = new BarParser();
		bp.setZipCompress(true);
		HashMap<String,String> tagVals = new HashMap<String,String>();
		if (stairStep) tagVals.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		else tagVals.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
		bp.writeBarFile(gr, chromosome, genomeVersion, strand.charAt(0), positions, values, tagVals);

		positions = null;
		values = null;
	}

	public void convertRatiosToPValues () {
		try {
			String randomFileName= Passwords.createRandowWord(10);
			int numWindows = windows.length;
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
			File tempSymPValueFile = new File(resultsDirectory,randomFileName+"tmpFileSymPValueFile");
			PrintWriter out = new PrintWriter(new FileWriter(tempSymPValueFile));
			for (int z=0; z<numWindows; z++){
				double score = Double.parseDouble(results[z]);
				if (score == 1) score = 0.999999;
				score = 1- score;
				out.println(score);
			}
			out.close();
			//run John Storeys qvalue from R, saving results in new file?
			File qValueResultsFile = new File(resultsDirectory,randomFileName+"tmpFileQValResFile");
			System.out.println("Launching John Storey's q-value app in R...");
			//write temp .R file
			StringBuffer rScript =new StringBuffer("library(qvalue)\np <- scan (\"");
			rScript.append(tempSymPValueFile.getCanonicalPath());
			rScript.append("\")\nqobj <- qvalue (p)\nqwrite (qobj, filename=\"");
			rScript.append(qValueResultsFile.getCanonicalPath());
			rScript.append("\")\n");
			File tempScriptFile = new File(resultsDirectory,randomFileName+"tmpFileRScriptFile");
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
			importScores(windows, qValueResultsFile);
			//kill temp files
			qValueResultsFile.delete();
			tempRatioFile.delete();
			tempScriptFile.delete();
			tempSymPValueFile.delete();
		} catch (IOException e){
			e.printStackTrace();
		}
	}

	/**Prints log2 oligo ratios as a chromosome specific xxx.bar file.
	 * Convert to print a bar file directly!*/
	public void printOligoRatios(){
		File oligoBarFile = new File(oligoDirectory, chromosome+".bar");
		double[] meanT = Num.averageFloatArraysFlipped(treatmentIntensities);
		int numberOligoMeasurements = meanT.length;
		int[] posChr = new int[numberOligoMeasurements];
		float[] valChr = new float[numberOligoMeasurements];

		double[] meanC = null; 
		if (pmOnly) {
			meanC = Num.averageFloatArraysFlipped(controlIntensities);
			//for every measurement
			for (int j=0; j<numberOligoMeasurements; j++){
				double ratio = Math.log(meanT[j]/meanC[j]) / Num.log2;
				posChr[j] = positions[j]+halfOligoLength;
				valChr[j] = new Double(ratio).floatValue();	
			}
		}
		else {
			//for every measurement
			for (int j=0; j<numberOligoMeasurements; j++){
				double ratio = Math.log(meanT[j])/ Num.log2;
				posChr[j] = positions[j]+halfOligoLength;
				valChr[j] = new Double(ratio).floatValue();	
			}
		}
		meanC = null;
		meanT = null;
		BarParser bp = new BarParser();
		bp.setZipCompress(true);
		HashMap<String,String> tagVals = new HashMap<String,String>();
		tagVals.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
		bp.writeBarFile(oligoBarFile, chromosome, genomeVersion, strand.charAt(0), posChr, valChr, tagVals);
		posChr = null;
		valChr = null;
	}

	/**Makes directories and files to hold point gr files and tempRatios*/
	public void makeDirectories(){
		if (printPointBars){
			if (pmOnly) {
				sumDirectory = new File(resultsDirectory, "Sum");
				if (sumDirectory.exists() == false) sumDirectory.mkdir();

			}
			relDiffDirectory = new File(resultsDirectory, extension);
			if (relDiffDirectory.exists() == false) relDiffDirectory.mkdir();

		}
		if (convertScoresToQValues) {
			//make random word for unique temp files;
			randomWord = Passwords.createRandowWord(10);
			tempRatioFile = new File(resultsDirectory, randomWord+"tmpRatioFile");
			try{
				tempRatios = new PrintWriter(new FileWriter(tempRatioFile));
			} catch (IOException e){
				e.printStackTrace();
			}
		}
		if (printOligoRatios) {
			oligoDirectory = new File(resultsDirectory, "Oligos");
			if (oligoDirectory.exists() == false) oligoDirectory.mkdir();

		}
		windowDirectory = new File(resultsDirectory, "Win");
		if (windowDirectory.exists() == false) windowDirectory.mkdir();
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
	 * Also makes gr files for the qValues.
	 * Saves score in last position of the window.scores[].*/
	public void importScores(Window[] windows, File importFile){
		int counter = 0;
		int finalIndex = -1;
		try {
			BufferedReader in = new BufferedReader( new FileReader (importFile));
			String line;
			String[] tokens;

			//skip first three lines
			in.readLine();in.readLine();in.readLine();
			//start
			finalIndex = windows[0].getScores().length -1;
			while ((line= in.readLine()) != null){
				tokens = line.split("\\s+");
				if (tokens.length!=2) continue;
				double score = Double.parseDouble(tokens[1]);
				score = (-10.0) * (Math.log(score)/Num.log10);
				double[] scores= windows[counter].getScores();
				scores[finalIndex] = score;
				windows[counter].setScores(scores);
				counter++;
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		//check counter
		if (counter != windows.length){
			Misc.printExit("\nProblem with QValue R script: the number of lines in the import file ("+counter+
					") do not match the number of Windows ("+windows.length+")!\n");
		}
		if (printPointBars){
			File grDir = new File( resultsDirectory, "QVal");
			makeBarFiles(windows, grDir, finalIndex, false);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}	
		new ScanChromosomes(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		String treatment = null;
		String control = null;
		File oligoPositions = null;
		String chroms = null;
		String resultsDirectoryString = null;
		String range = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': treatment = args[i+1]; i++; break;
					case 'c': control = args[i+1]; i++; break;
					case 'r': resultsDirectoryString = args[i+1]; i++; break;
					case 's': fullPathToPValue = args[i+1]; i++; convertScoresToQValues = true;break;
					case 'q': fullPathToR = args[i+1]; i++; convertScoresToQValues = true; break;
					case 'f': chroms = args[i+1]; i++; break;
					case 'g': range = args[i+1]; i++; break;
					case 'y': usePseudoMedian = true; averageGroups = true; break;
					case 'v': genomeVersion = args[i+1]; i++; break;
					case 'b': split = args[i+1]; i++; break;
					case 'k': strand = args[i+1]; i++; break;
					case 'p': convertScoresToQValues = true; break;
					case 'a': usePseudoMedian = true; break;
					case 'i': printOligoRatios = true; break;
					case 'e': printPointBars = true; break;
					case 'd': printHeatMapBars = true; break;
					case 'j': convertRelDiff2Log2 = true; break;
					case 'm': minOligos = Integer.parseInt(args[i+1]); i++; break;
					case 'z': sizeOfOligoMinusOne =Integer.parseInt(args[i+1])-1;i++; break;
					case 'w': windowSize = Integer.parseInt(args[i+1]); i++; break;
					case 'l': randomizeLabels = true; randomizeIntensities = false; break;
					case 'u': randomizeIntensities = true; randomizeLabels = false; break;
					case 'x': randomizeData = true; break;
					case 'n': numPermutations = Integer.parseInt(args[i+1]); i++; break;
					case 'o': oligoPositions = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for required parameters
		if (treatment == null || resultsDirectoryString == null || oligoPositions == null || genomeVersion == null){
			Misc.printExit("\nPlease complete one or more of the following required parameters: -t, -o, -v, or -r .\n");
		}

		//parse treatments
		treatmentDirectories = IO.extractFiles(treatment);
		if (treatmentDirectories == null) Misc.printExit("\nProblem parsing treatment directories -> "+treatment);

		//parse control
		if (control != null){
			controlDirectories = IO.extractFiles(control);
			if (controlDirectories == null) Misc.printExit("\nProblem parsing control directories -> "+control);
		}
		else {
			pmOnly = false;
			randomizeLabels = false;
			randomizeIntensities = false;
		}

		//particular chromosome processing?
		if (chroms !=null){
			particularChromosomes = chroms.split(",");
			chromosomeOligoPositions = new File[particularChromosomes.length];
			for (int i=0; i<particularChromosomes.length; i++){
				File singleChrome = new File (oligoPositions, particularChromosomes[i]);
				if (singleChrome.exists()) chromosomeOligoPositions[i] = singleChrome;
				else Misc.printExit("\nCannot find your requested chromosome in the oligo positions directory -> "+oligoPositions+" chrom: "+particularChromosomes[i]+"\n");
			}
			//necessary to sort!
			Arrays.sort(chromosomeOligoPositions);
		}
		//parse all oligo positions
		else{
			chromosomeOligoPositions = IO.extractFiles(oligoPositions);
			if (chromosomeOligoPositions == null) Misc.printExit("\nProblem parsing oligo positions -> "+oligoPositions);
		}

		//make and check results directory
		if (split!=null) {
			String[] x = split.split("-");
			if (x.length !=2) Misc.printExit("\nCannot parse your split chromosome(s) processing request?! Check -> "+split+"\n");
			splitBy = Integer.parseInt(x[0]);
			splitUnit = Integer.parseInt(x[1]);
			resultsDirectory = new File(resultsDirectoryString+"_"+splitUnit+"of"+splitBy);
			System.out.println("Splitting chromosome(s) by "+splitBy+", processing part "+splitUnit+".");
		}
		else resultsDirectory = new File(resultsDirectoryString);
		if (resultsDirectory.exists() == false) resultsDirectory.mkdir();

		halfOligoLength = (int)Math.round(((double)sizeOfOligoMinusOne)/2.0);

		if (range !=null){
			String[] minMax = range.split(",");
			minWindowFilter = Double.parseDouble(minMax[0]);
			maxWindowFilter = Double.parseDouble(minMax[1]);
		}

	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Scan Chromosmes: Dec 2007                             **\n" +
				"**************************************************************************************\n" +

				"SC uses a sliding window to calculate to window level statistics on processed cel\n" +
				"file directories (see CelProcessor -r and MakeChromosomeSets) using a Wilcoxon Rank\n" +
				"Sum test and one of the following statistics: a trimmed mean (10%) or pseudo median\n" +
				"on the average replica relative differences or a pseudo median on the all pair\n" +
				"relative differences within each window. To estimate confidence SC wraps Bourgon's\n"+
				"Symmetric Null P-value estimator and John Storey's QValue R package. SC can also use\n"+
				"random label permutation to estimate a p-value and FDR. A variety of bar files can be\n" +
				"saved for direct viewing in Affymetrix' IGB. xxxHM.bar should be viewed as heat maps.\n" +
				"Negative scores/ p-values/ FDRs indicate reduction, positive accumulation.\n\n" +

				"Note, the order of scores associated with each window are:" +
				"\n      [0]= "+SCORE_DESCRIPTION[0]+
				"\n      [1]= "+SCORE_DESCRIPTION[1]+
				"\n      [2]= "+SCORE_DESCRIPTION[2]+
				"\n      [3]= "+SCORE_DESCRIPTION[3]+
				"\n      [4]= "+SCORE_DESCRIPTION[4]+
				"\n\n"+
				"Required:\n"+
				"-o The 'OligoPositions' directory, full path, generated by MakeChromosomeSets.\n" +
				"-r The full path directory text to use in saving the results.\n" +
				"-t Treatment chip set directories, full path, comma delimited, no spaces.\n" +
				"-c Control chip set directories, full path, comma delimited, no spaces.\n" +
				"-v Genome version (ie hg17, dm2, ce2, mm8), get from UCSC Browser,\n" +
				"      http://genome.ucsc.edu/FAQ/FAQreleases, for bar files.\n" +
				"\n"+
				"Optional:\n"+
				"-w Window size, defaults to 675bp.\n"+
				"-m Minimal number of unique oligo positions required in each window. Defaults to\n" +
				"       10. Set to 1 to save all windows.\n" +
				"-z Size of oligo, defaults to 25.\n"+
				"-f Name(s) of chromosomes, comma delimited (e.g. chr21,chr22), to process, others\n" +
				"       will be skipped, defaults to all.\n"+
				"-b Break and process part of each chromosome(s), (e.g. 2-1 to split in 1/2 and\n" +
				"       and process the first half, 3-2 split in thirds, process the 2nd 3rd).\n"+
				"-a Use a layered pseudo median instead of a trimmed mean relative difference.\n" +
				"       This is very slow but much more robust. (For every oligo, all rel diffs are\n"+
				"       calculated between its treatment and control intensities.  These are pooled\n" +
				"       and a median is then calculated on all pairwise means.)\n"+
				"-y Use an average of the treatment and control replicas as input for the relative\n" +
				"       difference instead of making all pairwise relative differences in the\n"+
				"       pseudo median calculation.\n"+
				"-l Use random permutation of the chip labels to estimate a p-value and FDR for the\n" +
				"       window score.\n"+
				"-u Use random permutation of the intensity positions to estimate a p-value and FDR\n" +
				"       for the window score.\n"+
				"-n Number of random permutations, defaults to 10.\n"+
				"-x Randomize within replica intensities.\n"+
				"-p Convert window ratio scores to q-value FDRs, -10log10(multiple test corr pval).\n" +
				"      Requires installing Storey's R Q-Value Pkg and possibly compiling the SymPTest.\n"+
				"-s The full path to the symmetric_p_test, defaults to\n" +
				"      '"+fullPathToPValue+"'\n"+
				"-q The full path to qvalue loaded R, defaults to\n"+
				"      '"+fullPathToR+"'\n"+
				"-i Print an individual oligo log2 ratio (aveT/aveC) bar files.\n"+
				"-e Print point window summary bar files.\n"+
				"-d Print heat map window summary bar files.\n"+
				"-k Strand (either '+' or '-') for stranded data, used when writing bar files.\n"+
				"-j Convert relative difference scores to log2.\n"+
				"-g Exclude windows falling within this range (ie -0.2,0.2) when making heat maps.\n"+


				"\nExample: java -Xmx8000M -jar pathTo/T2/Apps/ScanChromosomesCNV -o /affy/OliPosHWG14 -r\n" +
				"      /affy/res/p53/ -t /Cels/T/B10_ChrNorm,/Cels/T/B11_ChrNorm -c /Cels/C/B10_ChrNorm,\n" +
				"      /Cels/C/B11_ChrNorm,/Cels/C/B12_ChrNorm -m 5 -i -e -l -n 3 -v hg17 \n\n" +

				"\n" +

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

	/**Fetches a float[] from each of the directories given a fileName.  Returns null 
	 * if a file cannot be read or it is not a serialized float[].
	 * @return float[replicas][intensities]*/
	public static float[][] fetchIntensities (File[] directories, String fileName, boolean randomizeData){
		float[][] intensities = new float[directories.length][];
		try {
			for (int i=0; i< directories.length; i++){
				File file = new File (directories[i], fileName);
				if (file.canRead() == false) return null;
				intensities[i] = (float[])IO.fetchObject(file);
				if (randomizeData) Num.randomize(intensities[i], i+System.currentTimeMillis());
			}
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
		return intensities;
	}

	/**Checks to see if the number of values are the same given float[replica#][values].*/
	public static boolean checkSizes(float[][] replicasValues){
		int size = replicasValues[0].length;
		for (int i=1; i< replicasValues.length; i++){
			if (replicasValues[i].length != size) return false;
		}
		return true;
	}

	public float[][] getControlIntensities() {
		return controlIntensities;
	}

	public int getNumPermutations() {
		return numPermutations;
	}

	public int[] getPositions() {
		return positions;
	}

	public File getResultsDirectory() {
		return resultsDirectory;
	}

	public float[][] getTreatmentIntensities() {
		return treatmentIntensities;
	}

	public Window[] getWindows() {
		return windows;
	}

}
