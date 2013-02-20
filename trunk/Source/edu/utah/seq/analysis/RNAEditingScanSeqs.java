package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import edu.utah.seq.useq.apps.Bar2USeq;
import edu.utah.seq.useq.data.Region;
import trans.tpmap.*;

/**Takes PointData and generates sliding window scan statistics. Estimates FDRs using random permutation.
 * @author Nix
 * */
public class RNAEditingScanSeqs {

	//user defined fields
	private File[] pointDataDirs;
	private File saveDirectory;
	private int windowSize = 500;
	private int minimumNumberObservationsInWindow = 5;
	private boolean strandedAnalysis = false;
	private int numberRandomPermutations = 100;
	private double fractionToTrim = 0.2f;
	private float minimumWindowPValue = 0.5f;
	private float minimumBaseFractionEdit = 0.005f;
	private float maximumBaseFractionEdit = 1.0f;
	private String regionsFileName;
	private LinkedHashMap<String,SmoothingWindow[]> chromosomeRegion;

	//internal fields
	private File trMeanDir;
	private File fdrDir;
	private File log2RatioDir;
	private boolean plusStrand = true;
	private boolean minusStrand =  true;
	private HashMap<String,PointData[]> plusPointData;
	private HashMap<String,PointData[]> minusPointData;
	private ArrayList<SmoothingWindowInfo> smiAL = new ArrayList<SmoothingWindowInfo>();

	private WindowMaker windowMaker; 
	private int[][] windows;
	private SmoothingWindow[] smoothingWindow;
	private SmoothingWindowInfo[] smoothingWindowInfo;
	private File swiFile;
	private String[] scoreNames;
	private String[] scoreDescriptions;
	private String[] scoreUnits;
	private String adapterName = "chrAdapter";

	//by chromosome
	private String chromosome;
	private PointData chromPointData = null;
	private PointData chromPlus = null;
	private PointData chromMinus = null;
	private float[][] randomScores;

	//constructor
	public RNAEditingScanSeqs(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);
		
		//make window directories
		trMeanDir = new File(saveDirectory, "TrMean");
		trMeanDir.mkdir();
		fdrDir = new File(saveDirectory, "FDR");
		fdrDir.mkdir();
		log2RatioDir = new File(saveDirectory, "Log2Ratio");
		log2RatioDir.mkdir();

		if (strandedAnalysis){
			plusStrand = true;
			minusStrand =  false;
			scan();
			//clear
			smiAL.clear();
			plusStrand = false;
			minusStrand =  true;
			scan();
		}
		else scan();
		
		//convert window data to useq
		new Bar2USeq(trMeanDir, true);
		new Bar2USeq(fdrDir, true);
		new Bar2USeq(log2RatioDir, true);

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	//methods
	public void scan(){
		//fetch counts
		loadPointDataArrays();

		//for each chromosome
		System.out.print("Scanning Chromosomes");
		
		if (chromosomeRegion == null){
			String[] chromosomes = fetchAllChromosomes();
			for (int i=0; i< chromosomes.length; i++){
				chromosome = chromosomes[i];
				windowScanChromosome();
			}
		}
		else {
			String[] chromosomes = fetchAllChromosomes();
			for (int i=0; i< chromosomes.length; i++){
				chromosome = chromosomes[i];
				windowScanChromosome();
			}
		}
		
		//any data, make arrays
		if (smiAL.size() == 0) Misc.printExit("\nNo windows found. Thresholds too strict?\n");
		smoothingWindowInfo = new SmoothingWindowInfo[smiAL.size()];
		smiAL.toArray(smoothingWindowInfo);
		
		//make EnrichedRegions 
		System.out.println("\nEstimating FDRs...");
		estimateFDRs();

		//save window data 
		String ext = "";
		if (plusStrand == true && minusStrand == true) ext = "";
		else if (plusStrand) ext = "Plus";
		else if (minusStrand) ext = "Minus";

		System.out.println("Writing window objects for the EnrichedRegionMaker...");
		swiFile = new File (saveDirectory, "windowData"+windowSize+"bp"+ext+".swi");
		IO.saveObject(swiFile, smoothingWindowInfo);
		
		System.out.println("Saving graph data...");
		writeBarFileGraphs();

	}

	public void estimateFDRs(){
		//collect all windows
		ArrayList<SmoothingWindow> smAL = new ArrayList<SmoothingWindow>();
		for (int i=0; i< smoothingWindowInfo.length; i++){
			SmoothingWindow[] win = smoothingWindowInfo[i].getSm();
			for (int j = 0; j< win.length; j++){
				smAL.add(win[j]);
			}
		}
		smoothingWindow = new SmoothingWindow[smAL.size()];
		smAL.toArray(smoothingWindow);
		
		//window scores = numberObs, real mean, pval, random mean
		
		//sort windows by real mean score, smallest to largest
		Arrays.sort(smoothingWindow, new ComparatorSmoothingWindowScore(1));
		
		EnrichedRegionMaker erm = new EnrichedRegionMaker(windowSize, smoothingWindowInfo);
		
		//for each window with diff score
		float priorScore = 0;
		float priorFDR = 0;
		for (int i=0; i< smoothingWindow.length; i++){
			//scores = numberObs, real mean, pval, random mean
			float[] currScores = smoothingWindow[i].getScores();
			if (currScores[1] != priorScore){
				//diff score so calc fdr
				double numReal = erm.countEnrichedRegions(false, 1, currScores[1]);
				double numRand = erm.countEnrichedRegions(false, 3, currScores[1]);
				
				if (numRand == 0){
					//no more random regions so assign last FDR * 1%
					priorFDR = priorFDR * 1.01f;
					for (; i< smoothingWindow.length; i++) {
						currScores = smoothingWindow[i].getScores();
						//set fdr
						currScores[2] = priorFDR;
					}
					break;
				}
				float fdr = new Double(-10 * Num.log10(numRand/numReal)).floatValue();
				//only increase the FDR with increasing stringent thresholds
				if (fdr> priorFDR) priorFDR = fdr;
				priorScore = currScores[1];
			}
			//set FDR
			currScores[2] = priorFDR;
		}
		
		//set log2Ratio 
		for (SmoothingWindow w: smoothingWindow){
			float[] currScores = w.getScores();
			currScores[3] = Num.log2(currScores[1]/currScores[3]);
			//if (currScores[2]>13) System.out.println(w);

		}
		
		//final scores = numberObs	realMean	FDR	log2Ratio
	}


	/**Fetches the names of all the chromosomes in the data.*/
	public String[] fetchAllChromosomes(){
		HashSet<String> c = new HashSet<String>();
		Iterator<String> it = plusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = minusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		
		return Misc.hashSetToStringArray(c);
	}

	/**Fetchs the data for a particular chromosome.*/
	public boolean fetchData(){

		//fetch data
		chromPlus = null;
		if (plusPointData.containsKey(chromosome) && plusStrand) chromPlus = PointData.combinePointData(plusPointData.get(chromosome), true);
		chromMinus = null;
		if (minusPointData.containsKey(chromosome) && minusStrand) chromMinus = PointData.combinePointData(minusPointData.get(chromosome), true);
		
		chromPointData = null;
		if (plusStrand && minusStrand == false) chromPointData = chromPlus;
		else if (minusStrand && plusStrand == false) chromPointData = chromMinus;
		else {
			ArrayList<PointData> pdAL = new ArrayList<PointData>();
			if (chromPlus != null ) pdAL.add(chromPlus);
			if (chromMinus != null ) pdAL.add(chromMinus);
			if (pdAL.size() == 2) chromPointData = PointData.combinePointData(pdAL, false);
			else if (pdAL.size() == 1) chromPointData = pdAL.get(0);
		}
		
		//toss low values values
		if (chromPointData != null){
			ArrayList<Integer> goodPos = new ArrayList<Integer>(100000);
			ArrayList<Float> goodVal = new ArrayList<Float>(100000);
			int[] pos = chromPointData.getPositions();
			float[] val = chromPointData.getScores();
			for (int i=0; i< pos.length; i++){
				if (val[i] >= minimumBaseFractionEdit && val[i] < maximumBaseFractionEdit){
					goodPos.add(pos[i]);
					goodVal.add(val[i]);
				}
			}
			pos = Num.arrayListOfIntegerToInts(goodPos);
			val = Num.arrayListOfFloatToArray(goodVal);
			chromPointData.setPositions(pos);
			chromPointData.setScores(val);
		}
		
		if (chromPointData != null) return true;
		return false;
	}

	/**Window scans a chromosome collecting read count data and calculating binomial p-values.*/
	public void windowScanChromosome(){
		//fetch data
		if (fetchData() == false) {
			System.out.println("\n\tSkipping "+chromosome+". Failed to find paired PointData datasets.");
			return;
		}
		//fetch all positions 
		int[] positions = chromPointData.getPositions();
		
		//make windows using all of the reads
		makeWindows(positions);
		
		//any windows?
		if (windows.length == 0){
			System.out.println("\n\tSkipping "+chromosome+". No windows found with minimum observations of "+minimumNumberObservationsInWindow+" within a window size of "+windowSize);
			return;
		}
		System.out.print(".");
		//make SmoothingWindow[] container
		smoothingWindow = new SmoothingWindow[windows.length];
		randomScores = new float [windows.length][numberRandomPermutations];
		
		//scan for real scores
		scanWindowsReal();
		
		//scan permutated base fractions
		scanPermutedWindows();
		
		//calculate p-values
		calculateWindowPValues();
		
		//save results
		addWindowDataToSMIArrayList();
	}

	/**Saves window data to SMI Array.*/
	public void addWindowDataToSMIArrayList(){
		Info info = chromPointData.getInfo();
		HashMap<String,String> notes = new HashMap<String,String>();
		notes.put(BarParser.WINDOW_SIZE, windowSize+"");
		notes.put(BarParser.DESCRIPTION_TAG, Misc.stringArrayToString(scoreNames, ","));
		notes.put(BarParser.UNIT_TAG, Misc.stringArrayToString(scoreUnits, ","));
		info.setNotes(notes);
		info.setStrand(".");		
		smiAL.add(new SmoothingWindowInfo(smoothingWindow, info));
	}

	private void calculateWindowPValues(){
		//for each window 
		ArrayList<SmoothingWindow> goodWindows = new ArrayList<SmoothingWindow>();
		for (int i=0; i< smoothingWindow.length; i++){
			//scores = numberObs, mean, startIndex, stopIndex
			float[] scores = smoothingWindow[i].getScores();
			//get array of permuted mean's
			float[] rs = randomScores[i];
			//count number greater than or equal to
			float badTrials = 0;
			for (int x=0; x< numberRandomPermutations; x++ ){
				if (rs[x] >= scores[1]) badTrials++;
			}
			//pval
			scores[2] = badTrials/(float)numberRandomPermutations;
			if (scores[2] <= minimumWindowPValue) {
				//random mean
				scores[3] = Num.mean(rs);
				goodWindows.add(smoothingWindow[i]);
			}
			//final scores = numberObs, real mean, pval, random mean
			
			
		}
		//replace unfiltered sms
		SmoothingWindow[] sm = new SmoothingWindow[goodWindows.size()];
		goodWindows.toArray(sm);
		smoothingWindow = sm;
	}

	private void scanPermutedWindows(){
		for (int i=0; i<numberRandomPermutations; i++){
			//randomize chromPointData values
			float[] vals = chromPointData.getScores();
			Num.randomize(vals, 0);
			chromPointData.setScores(vals);
			//scan it
			scanWindowsRandom(i);
		}
	}
	
	private void scanWindowsRandom(int randomIndex){
		//get scores
		float[] fractions = chromPointData.getScores();
		//for each window 
		for (int i=0; i< windows.length; i++){
			//get start stop
			//scores = numberObs, mean, startIndex, stopIndex
			float[] winScores = smoothingWindow[i].getScores();
			//num fractions
			int numFrac = (int)(winScores[3]-winScores[2]);
			//calc trimmed mean
			float[] winFrac = new float[numFrac];
			System.arraycopy(fractions, (int)winScores[2], winFrac, 0, numFrac);
			Arrays.sort(winFrac);
			float trimmedMean = (float)Num.trimmedMean(winFrac, fractionToTrim);
			randomScores[i][randomIndex] = trimmedMean;
		}
	}


	private void scanWindowsReal(){
		//get scores
		float[] fractions = chromPointData.getScores();
		//for each window 
		for (int i=0; i< windows.length; i++){
			//fetch indexes
			int[] indexStartStop = chromPointData.findIndexes(windows[i][0], windows[i][1]);
			//get fractions
			int numFrac = indexStartStop[1]-indexStartStop[0];
			float[] winFrac = new float[numFrac];
			System.arraycopy(fractions, indexStartStop[0], winFrac, 0, numFrac);
			Arrays.sort(winFrac);
			float trimmedMean = (float)Num.trimmedMean(winFrac, fractionToTrim);
			//float mean = Num.mean(indexStartStop[0], indexStartStop[1], fractions);
			//scores = numberObs, mean, startIndex, stopIndex
			float[] scores = new float[]{numFrac, trimmedMean ,  (float)indexStartStop[0], (float)indexStartStop[1]};
			//make window
			smoothingWindow[i] = new SmoothingWindow (windows[i][0], windows[i][1], scores);
			/*if (smoothingWindow[i].getStart() == 1011704	 && smoothingWindow[i].getStop() == 1012202) {
				System.out.println(smoothingWindow[i]);
				Misc.printArray(winFrac);
			}*/
		}
	}




	/**Collects and calculates a bunch of stats re the PointData.*/
	private void loadPointDataArrays(){
		//fetch treatment PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (pointDataDirs);
		//remove chrAdapter
		combo[0].remove(adapterName);
		combo[1].remove(adapterName);
		plusPointData = PointData.convertArrayList2Array(combo[0]);
		minusPointData = PointData.convertArrayList2Array(combo[1]);
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
				Info info = smoothingWindowInfo[i].getInfo();
				SmoothingWindow[] sm = smoothingWindowInfo[i].getSm();
				if (plusStrand && minusStrand == false) info.setStrand("+");
				else if (minusStrand && plusStrand == false) info.setStrand("-");
				else info.setStrand(".");
				
				//final scores = numberObs, real mean, FDR, log2Ratio
				saveSmoothedHeatMapData (1, sm, info, trMeanDir, "#FF00FF", false); //magenta
				saveSmoothedHeatMapData (2, sm, info, fdrDir, "#00FF00", false); //green
				saveSmoothedHeatMapData (3, sm, info, log2RatioDir, "#FFFF00", false); //yellow
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
		String fileNames = Misc.stringArrayToString(IO.fetchFileNames(pointDataDirs),",");
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

	/**Sets score names/ descriptions/ units base on whether control data is present.*/
	public void setScoreStrings(){
		//final scores = numberObs, real mean, FDR, log2Ratio
			scoreNames = new String[]{
					"Obs",
					"TrMean",
					"FDR",
					"Log2Ratio"
			};
			scoreDescriptions = new String[]{
					"# non zero base fraction edits in window",
					"Trimmed mean of the non zero base fraction edits in the window",
					"FDR",
					"Log2Ratio(obsMean/randomMean)"
			};
			scoreUnits = new String[]{
					"Count",
					"Fraction",
					"-10Log10(FDR)",
					"Log2Ratio(obsMean/randomMean)"
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
		String strand = null;
		File bedFile = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': pointDataDirs = IO.extractFiles(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'r': saveDirectory = new File(args[++i]); break;
					case 's': strandedAnalysis = true; break;
					case 'w': windowSize = Integer.parseInt(args[++i]); break;
					case 'm': minimumNumberObservationsInWindow = Integer.parseInt(args[++i]); break;
					case 'i': minimumBaseFractionEdit = Float.parseFloat(args[++i]); break;
					case 'x': maximumBaseFractionEdit = Float.parseFloat(args[++i]); break;
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
		if (pointDataDirs == null || pointDataDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your treatment PointData directories(s)!\n");
		//only one directory look deeper
		if (pointDataDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(pointDataDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) pointDataDirs = otherDirs;
		}

		if (windowSize == 0 ) Misc.printExit("\nPlease enter a positive length for the window size if you are going to set the peak shift to 0\n");

		//strand?
		if (strand != null && (strand.equals("+") == false && strand.equals("-") == false)){
			Misc.printExit("\nError: Enter either + or - for a stranded scan.\n");
		}
		plusStrand =  (strand == null || strand.equals("+"));
		minusStrand =  (strand == null || strand.equals("-"));

		//make window maker 
		windowMaker = new WindowMaker(windowSize,minimumNumberObservationsInWindow);

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		if (saveDirectory.exists() == false) saveDirectory.mkdir();


		//set score items
		setScoreStrings();
		
		//load regions
		if (bedFile != null) {
		HashMap<String, Region[]>cr = Region.parseStartStops(bedFile, 0, 0, 0);
		chromosomeRegion = new LinkedHashMap<String, SmoothingWindow[]>();
		for (String chr: cr.keySet()){
			Region[] regions = cr.get(chr);
			SmoothingWindow[] sw = new SmoothingWindow[regions.length];
			for (int i=0; i< regions.length; i++){
				sw[i] = new SmoothingWindow (regions[i].getStart(), regions[i].getStop(), new float[]{0,0,0,0,0,0});
			}
			chromosomeRegion.put(chr, sw);
		}
		regionsFileName = Misc.removeExtension(bedFile.getName());
		}

		

	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            RNA Editing Scan Seqs: Feb 2013                       **\n" +
				"**************************************************************************************\n" +
				"Beta.\n\n" +

				"Options:\n"+
				"-r Results directory, full path.\n"+
				"-p PointData directoryfrom the RNAEditingPileUpParser Base Fraction Edited.\n" +
				"       These should contain stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories. These will be merged when scanning.\n" +
				"-s Perform stranded analysis, defaults to unstranded.\n"+
				"-w Window size, defaults to 500.\n"+
				"-m Minimum number observations in window, defaults to 5. \n" +
				"-i Minimum base fraction editing, defaults to 0.005 \n"+
				"-x Maximum base fraction editing, defaults to 1. Set to > 1 to retain bases with \n" +
				"       100% base fraction editing. These are often SNPs!\n"+
				//"-b A bed file of regions to score (tab delimited: chr start stop ...), defaults to\n" +
				//"       window scanning.\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/ \n\n" +

				"**************************************************************************************\n");

	}

	public SmoothingWindowInfo[] getSmoothingWindowInfo() {
		return smoothingWindowInfo;
	}

	public File getSwiFile() {
		return swiFile;
	}
}
