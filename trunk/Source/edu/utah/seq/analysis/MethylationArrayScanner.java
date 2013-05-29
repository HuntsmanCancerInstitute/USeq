package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import trans.main.WilcoxonRankSumTest;
import trans.main.WilcoxonSignedRankTest;
import trans.tpmap.WindowMaker;
import util.gen.*;
import edu.utah.seq.data.ComparatorPointAscendingScore;
import edu.utah.seq.data.ComparatorPointPosition;
import edu.utah.seq.data.HeatMapMaker;
import edu.utah.seq.data.Info;
import edu.utah.seq.data.MethylationArraySamplePair;
import edu.utah.seq.data.Point;
import edu.utah.seq.data.PointData;
import edu.utah.seq.data.SmoothingWindow;
import edu.utah.seq.data.SmoothingWindowInfo;
import edu.utah.seq.parsers.BarParser;
import edu.utah.seq.useq.apps.Bar2USeq;


/** Application for looking for changes in methylation array data across a genome.  I've hacked in a non paired analysis on top of the original app.
 * @author Nix
 * */
public class MethylationArrayScanner {

	//fields
	private int windowSize = 1000;
	private int minimumNumberObservationsInWindow = 10;
	private float minimumPseudoMedianRatio = 0.2f;
	private int numberRandomTrials = 5;
	private boolean performPairedAnalysis = true;

	private HashMap<String, MethylationArraySamplePair[]> chromMethylationArraySamplePairs;
	private HashMap<String, int[]> chromPositions = new HashMap<String, int[]>();
	private String chromosome;
	private int[] positions;
	private MethylationArraySamplePair[] pairs;
	private WindowMaker windowMaker; 
	private int[][] windows;
	private SmoothingWindow[] smoothingWindow;
	private SmoothingWindow[] allSmoothingWindows;
	private SmoothingWindowInfo[] smoothingWindowInfo;
	private ArrayList<Float> randomPValues = new ArrayList<Float>();
	private ArrayList<SmoothingWindowInfo> smiAL = new ArrayList<SmoothingWindowInfo>();
	private ArrayList<SmoothingWindow> smALAll = new ArrayList<SmoothingWindow>();
	private ArrayList<SmoothingWindow> smALChrom = new ArrayList<SmoothingWindow>();
	private String[] scoreNames;
	private String[] scoreDescriptions;
	private String[] scoreUnits;
	private String versionedGenome;
	private File pseDir;
	private File fdrDir;
	private File baseRatioDir;
	private File saveDirectory;
	private int pseIndex = 0;
	private int pvalueBHFDRIndex = 1;
	private int permFDRIndex = 2;
	private int numObsIndex = 3;
	private File dataDirectory;
	private String treatmentPairedSamples;
	private String controlPairedSamples;
	private int numberMethylationArraySamplePairs;
	private double numberRandomWindows;


	//constructor

	public MethylationArrayScanner(String[] args){	
		long startTime = System.currentTimeMillis();
		//set fields
		processArgs(args);

		//launch
		run();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}

	public void run(){
		//load data hashes
		System.out.println("Loading paired samples...");
		loadMethylationArraySamplePairs();

		//make window maker, compensate for multiple samples
		int minObs = (int)Math.round(((double)minimumNumberObservationsInWindow / (double) numberMethylationArraySamplePairs));
		if (minObs == 1) minObs = 2;
		windowMaker = new WindowMaker(windowSize,(int)minObs);

		//for each chromosome
		System.out.println("\nWindow scanning...");
		for (String chrom: chromMethylationArraySamplePairs.keySet()){
			//fetch data
			chromosome = chrom;		
			
			positions = chromPositions.get(chromosome);
			pairs = chromMethylationArraySamplePairs.get(chromosome);

			//make windows
			windows = windowMaker.makeWindows(positions);
			if (windows.length == 0){
				System.out.println("\t"+chromosome+" Skipping! No windows found with minimum reads of "+minimumNumberObservationsInWindow+" within a window size of "+windowSize);
				continue;
			}
			System.out.println("\t"+chromosome);

			//scan windows
			scanWindows();
			saveWindowDataToSMIArrayList();

			for (int i=0; i< numberRandomTrials; i++){
				//randomize pairs data
				for (MethylationArraySamplePair sp: pairs) sp.randomize();
				//scan for random
				scanWindowsRandom();
			}
		}

		//make arrays
		smoothingWindowInfo = new SmoothingWindowInfo[smiAL.size()];
		smiAL.toArray(smoothingWindowInfo);
		allSmoothingWindows = new SmoothingWindow[smALAll.size()];
		smALAll.toArray(allSmoothingWindows);

		//calculate permutation FDRs, do this first!
		System.out.println("\nCalculating permutation FDRs...");
		calculatePermutationFDRs();

		//multiple test correct the pvalues with B&H, do this second since it modifies the pVal
		System.out.println("Converting wilcoxon pvalues to FDRs with B&H...");
		convertPValuesToFDRs();

		//write out graph data
		System.out.println("Saving graph data...");
		writeBarFileGraphs();

		//convert graph data to useq
		new Bar2USeq(pseDir, true);
		new Bar2USeq(fdrDir, true);
		new Bar2USeq(baseRatioDir, true);

		//save window data 
		System.out.println("Writing window objects for the EnrichedRegionMaker...");
		File swiFile = new File (saveDirectory, "windowData"+windowSize+"bp"+Num.formatNumberOneFraction(minimumPseudoMedianRatio)+"MinPse.swi");
		IO.saveObject(swiFile, smoothingWindowInfo);

	}

	private void calculatePermutationFDRs() {
		//sort RandomScores 
		float[] rs = Num.arrayListOfFloatToArray(randomPValues);
		Arrays.sort(rs);

		//score windows passing pseudoMedian threshold
		for (SmoothingWindow sw: allSmoothingWindows){
			float[] scores = sw.getScores();
			float permFDR = 0;
			if (scores[pseIndex] >= minimumPseudoMedianRatio) {
				permFDR = calculateTransformedPValue(scores[pvalueBHFDRIndex], rs, numberRandomWindows);
			}
			scores[permFDRIndex] = permFDR;
			//if (permFDR !=0) System.out.println("PermFDR "+permFDR+"\tPVal "+scores[pvalueBHFDRIndex]);
		}
	}

	public static float calculateTransformedPValue (float realScore, float[] sortedRandomScores, double totalNumberRandomScores){
		int index = Num.findClosestIndexToValue(sortedRandomScores, realScore);
		double count = sortedRandomScores.length - index;
		if (count == 0.0) count = 1.0;
		return Num.minus10log10Float(count/totalNumberRandomScores);
	}

	private void convertPValuesToFDRs() {
		//collect pvals for B&H corr
		Point[] pvals = new Point[allSmoothingWindows.length];
		for (int i=0; i< allSmoothingWindows.length; i++){
			pvals[i] = new Point(i, allSmoothingWindows[i].getScores()[pvalueBHFDRIndex]);
		}
		//sort by score
		Arrays.sort(pvals, new ComparatorPointAscendingScore());
		//correct
		Point.benjaminiHochbergCorrect(pvals, 0);
		//sort back to original position
		Arrays.sort(pvals, new ComparatorPointPosition());
		//assign FDRs to pVals
		for (int i=0; i< allSmoothingWindows.length; i++){
			float scores[] = allSmoothingWindows[i].getScores();
			float fdr = pvals[i].getScore();
			if (fdr > 0) scores[pvalueBHFDRIndex] = fdr;
			else scores[pvalueBHFDRIndex] = 0.0f;
		}
	}


	/**Writes stair step window bar graph files*/
	public void writeBarFileGraphs(){

		//for each chromosome
		for (int i=0; i< smoothingWindowInfo.length; i++){
			Info in = smoothingWindowInfo[i].getInfo();
			if (in == null) Misc.printErrAndExit("Info is null");
			Info info = smoothingWindowInfo[i].getInfo().copy();
			SmoothingWindow[] sm = smoothingWindowInfo[i].getSm();

			//final scores = mean, PVal, fdr, # obs		
			saveSmoothedHeatMapData (pseIndex, sm, info, pseDir, "#FF00FF"); //magenta
			saveSmoothedHeatMapData (pvalueBHFDRIndex, sm, info, fdrDir, "#00FF00"); //green
			//saveSmoothedHeatMapData (2, sm, info, fdrDir, "#FFFF00"); //yellow
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

	/**Saves window data to SMI Array.*/
	public void saveWindowDataToSMIArrayList(){
		if (smALChrom.size() == 0) return;
		smoothingWindow = new SmoothingWindow[smALChrom.size()];
		smALChrom.toArray(smoothingWindow);
		HashMap<String,String> notes = new HashMap<String,String>();
		notes.put(BarParser.WINDOW_SIZE, windowSize+"");
		notes.put(BarParser.DESCRIPTION_TAG, Misc.stringArrayToString(scoreNames, ","));
		notes.put(BarParser.UNIT_TAG, Misc.stringArrayToString(scoreUnits, ","));
		Info info = new Info(chromosome, versionedGenome, chromosome, ".", 1, notes);		
		smiAL.add(new SmoothingWindowInfo(smoothingWindow, info));
	}


	private void scanWindows() {
		ArrayList<Float> treatmentAL = new ArrayList<Float>();
		ArrayList<Float> controlAL = new ArrayList<Float>();
		smALChrom.clear();

		//for each window (startIndex, stopIndex)
		for (int i=0; i< windows.length; i++){
			treatmentAL.clear();
			controlAL.clear();

			//fetch scores from each sample 
			for (int x=0; x< pairs.length; x++){
				pairs[x].fetchScoresByIndex(windows[i][0], windows[i][1], treatmentAL, controlAL);
			}

			//enough obs?
			if (treatmentAL.size() < minimumNumberObservationsInWindow) continue;

			float[] treatment = Num.arrayListOfFloatToArray(treatmentAL);
			float[] control = Num.arrayListOfFloatToArray(controlAL);
			double pvalue = 0;
			double pse = 0;
			
			//paired analysis
			if (performPairedAnalysis){
				double[] fraction = Num.ratio(treatment, control);
				Num.log2(fraction);
				pse = Num.pseudoMedian(fraction);
				//calculate pval if > 4 pairs
				if (treatment.length > 4){
					WilcoxonSignedRankTest w = new WilcoxonSignedRankTest(treatment, control);
					pvalue = w.getTransformedPValue();
				}
			}
			//non paired analysis
			else{
				//calc ratio
				double pseT = Num.pseudoMedian(treatment);
				double pseC = Num.pseudoMedian(control);
				pse = Num.log2(pseT/pseC);
				//calc pval if more than 9 obs
				if (treatment.length > 9 && control.length > 9){
					WilcoxonRankSumTest w = new WilcoxonRankSumTest();
					pvalue = w.test(treatment, control);
				}
			}
			
			//scores = pse, pvalFDR, permFDR, # obs   
			float[] scores = new float[4];
			scores[pseIndex] = (float) pse;
			scores[pvalueBHFDRIndex] = (float) pvalue;
			scores[numObsIndex] = treatment.length;

			//make window
			SmoothingWindow win = new SmoothingWindow (positions[windows[i][0]], positions[windows[i][1]]+1, scores);
			smALChrom.add(win);

			/*if (win.getStart() == 179456){
				System.out.println(win);
				Misc.printArray(treatment);
				Misc.printArray(control);
				System.exit(0);
			}*/

		}
		smALAll.addAll(smALChrom);

	}

	private void scanWindowsRandom() {
		ArrayList<Float> treatmentAL = new ArrayList<Float>();
		ArrayList<Float> controlAL = new ArrayList<Float>();

		//for each window 
		for (int i=0; i< windows.length; i++){
			treatmentAL.clear();
			controlAL.clear();
			
			//fetch scores from each sample pair, stopIndex is included
			for (int x=0; x< pairs.length; x++){
				pairs[x].fetchScoresByIndex(windows[i][0], windows[i][1], treatmentAL, controlAL);
			}
			//enough obs?
			if (treatmentAL.size() < minimumNumberObservationsInWindow) continue;
			
			
			float[] treatment = Num.arrayListOfFloatToArray(treatmentAL);
			float[] control = Num.arrayListOfFloatToArray(controlAL);
			double pvalue = 0;
			
			//paired analysis
			if (performPairedAnalysis){
				double[] fraction = Num.ratio(treatment, control);
				Num.log2(fraction);
				if (Math.abs(Num.pseudoMedian(fraction)) < minimumPseudoMedianRatio) continue;
				//calculate pval if > 4 pairs
				if (treatment.length > 4){
					WilcoxonSignedRankTest w = new WilcoxonSignedRankTest(treatment, control);
					pvalue = w.getTransformedPValue();
				}
			}
			//non paired analysis
			else{
				//calc ratio
				double pseT = Num.pseudoMedian(treatment);
				double pseC = Num.pseudoMedian(control);
				if (Math.abs(Num.log2(pseT/pseC)) < minimumPseudoMedianRatio) continue;
				//calc pval if more than 9 obs
				if (treatment.length > 9 && control.length > 9){
					WilcoxonRankSumTest w = new WilcoxonRankSumTest();
					pvalue = w.test(treatment, control);
				}
			}
			randomPValues.add((float)pvalue);
			numberRandomWindows++;
			
		}

	}


	public void loadMethylationArraySamplePairs(){
		//parse names
		String[] treatmentNames = treatmentPairedSamples.split(",");
		String[] controlNames = controlPairedSamples.split(",");
		if (treatmentNames.length != controlNames.length) Misc.printErrAndExit("\nError: the number of treatment and control samples differ?");
		numberMethylationArraySamplePairs = treatmentNames.length;

		HashMap<String, ArrayList<MethylationArraySamplePair>> pairs = new HashMap<String, ArrayList<MethylationArraySamplePair>>();

		for (int i=0; i< treatmentNames.length; i++){
			System.out.println("\t"+ treatmentNames[i] +"\t"+ controlNames[i]);
			HashMap<String, PointData> pdT = loadPointData(treatmentNames[i]);
			HashMap<String, PointData> pdC = loadPointData(controlNames[i]);

			//write log2(t/2) graphs
			writeProbeRatioGraph(pdT, pdC, treatmentNames[i]+"_"+controlNames[i]);

			HashSet<String> allChroms = new HashSet<String>();
			allChroms.addAll(pdT.keySet());
			allChroms.addAll(pdC.keySet());
			for (String chrom: allChroms){
				MethylationArraySamplePair sp = new MethylationArraySamplePair(pdT.get(chrom).getScores(), pdC.get(chrom).getScores(), performPairedAnalysis);
				ArrayList<MethylationArraySamplePair> al = pairs.get(chrom);
				if (al == null){
					al = new ArrayList<MethylationArraySamplePair>();
					pairs.put(chrom, al);
				}
				al.add(sp);
			}
		}

		//convert to []s
		chromMethylationArraySamplePairs = new HashMap<String, MethylationArraySamplePair[]>();
		for (String chrom: pairs.keySet()){
			ArrayList<MethylationArraySamplePair> al = pairs.get(chrom);
			if (al.size() != treatmentNames.length) Misc.printErrAndExit("\nError: One or more samples are missing a chromosome PointData set.");
			MethylationArraySamplePair[] sp = new MethylationArraySamplePair[treatmentNames.length];
			al.toArray(sp);
			chromMethylationArraySamplePairs.put(chrom, sp);
		}
	}

	private void writeProbeRatioGraph(HashMap<String, PointData> pdT, HashMap<String, PointData> pdC, String name) {
		File sampleDir = new File (baseRatioDir, name);
		sampleDir.mkdir();
		//for each chromosome
		for (String chrom : pdT.keySet()){
			//make ratios
			PointData t = pdT.get(chrom);
			int[] pos = t.getPositions();
			float[] tScores = t.getScores();
			float[] cScores = pdC.get(chrom).getScores();
			ArrayList<Point> good = new ArrayList<Point>();
			for (int i=0; i< pos.length; i++){
				if (tScores[i] !=0 && cScores[i] !=0) {
					float ratio = Num.log2(tScores[i]/ cScores[i]);
					good.add(new Point(pos[i], ratio));
				}
			}
			//make PointData
			PointData r = Point.extractPositionScores(good);
			Info info = t.getInfo().copy();
			info.setName("Log2Ratio(t/c)");
			//add info to hashmap for writing to bar file
			HashMap<String,String> map = new HashMap<String,String>();		
			//what graph type should be used to display it?
			map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
			//color orange
			map.put(BarParser.GRAPH_TYPE_COLOR_TAG, "#FF6600");
			//description
			map.put(BarParser.DESCRIPTION_TAG, "Log2Ratio of t/c base level contrast");
			//save in info
			info.setNotes(map);
			r.setInfo(info);
			r.writePointData(sampleDir);

		}

	}

	/**Does a variety of checks, also sets the positions for each chrom.*/
	public HashMap<String, PointData> loadPointData(String name){
		File f = new File (dataDirectory, name);
		if (f.exists() == false) Misc.printErrAndExit("\nError: failed to find sample -> "+f);
		if (f.isDirectory() == false) Misc.printErrAndExit("\nError: sample doesn't appear to be a directory -> "+f);
		HashMap<String, PointData> pd = PointData.fetchPointData (f, null, false);
		if (pd == null || pd.keySet().size() == 0) Misc.printErrAndExit("\nError: failed to load PointData from -> "+f);
		
		//check and or set positions
		for (String chrom: pd.keySet()){
			int[] existingPos = chromPositions.get(chrom);
			PointData p = pd.get(chrom);
			int[] pdPos = p.getPositions();
			if (versionedGenome == null) versionedGenome = p.getInfo().getVersionedGenome();
			if (existingPos == null) chromPositions.put(chrom, pdPos);
			else {
				if (existingPos.length != pdPos.length) Misc.printErrAndExit("\nError: lengths of PointData don't match.");
				for (int i=0; i< existingPos.length; i++){
					if (existingPos[i] != pdPos[i]) Misc.printErrAndExit("\nError: positions in PointData don't match.");
				}
			}
		}
		return pd;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MethylationArrayScanner(args);
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
					case 'd': dataDirectory = new File (args[++i]); break;
					case 's': saveDirectory = new File (args[++i]); break;
					case 't': treatmentPairedSamples = args[++i]; break;
					case 'c': controlPairedSamples = args[++i]; break;
					case 'w': windowSize = Integer.parseInt(args[++i]); break;
					case 'o': minimumNumberObservationsInWindow = Integer.parseInt(args[++i]); break;
					case 'r': numberRandomTrials = Integer.parseInt(args[++i]); break;
					case 'p': minimumPseudoMedianRatio = Float.parseFloat(args[++i]); break;
					case 'n': performPairedAnalysis = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}	
		}
		//look for data dir
		if (dataDirectory == null || dataDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: cannot find your data directory containing sample bar file directories.\n");

		//samples
		if (treatmentPairedSamples == null || controlPairedSamples == null) Misc.printErrAndExit("\nError: please enter at least one paired treatment control sample set to contrast.\n");

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		if (saveDirectory.exists() == false) saveDirectory.mkdir();

		pseDir = new File(saveDirectory, "WindowLog2Rto");
		pseDir.mkdirs();
		fdrDir = new File(saveDirectory, "WindowWilcoxFDR");
		fdrDir.mkdirs();
		baseRatioDir = new File(saveDirectory, "BPLog2Ratio");
		baseRatioDir.mkdirs();

		scoreNames = new String[]{
				"Log2Rto",
				"WilcoxFDR",
				"PermFDR",
				"#Obs",
		};
		scoreDescriptions = new String[]{
				"Pseudomedian of log2(treatment/control) ratios",
				"B&H Corrected Wilcoxon Signed Rank test p-value",
				"Permutation extimated FDR",
				"Number paired observations",
		};
		scoreUnits = new String[]{
				"pseudomedian(log2(treatment/control) ratios)",	
				"-10Log10(WilcoxFDR)",
				"-10Log10(PermFDR)",
				"count",
		};
		//change for non paired
		if (performPairedAnalysis == false){
			scoreDescriptions[0] = "log2( pseudomedian treatments/ pseudomedian controls) ratio";
			scoreUnits[0] = "log2(pseT/pseC)";
		}
	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Methylation Array Scanner: May 2013                      **\n" +
				"**************************************************************************************\n" +
				"MAS takes paired or non-paired sample PointData representing beta values (0-1) from\n" +
				"arrays and scores regions with enriched/ reduced signal using a sliding window\n" +
				"approach. A B&H corrected Wilcoxon signed rank (or rank sum test for non-paired),\n" +
				"pseudo median of the log2(treat/control) ratios (or log2(pseT/pseC) for non-paired),\n" +
				"and permutation test FDR is calculated for each window. Use the EnrichedRegionMaker\n" +
				"to identify enriched and reduced regions by picking thresholds (e.g. -i 0,1 -s 0.2,13).\n" +
				"MAS generates several data tracks for visualization in IGB including paired sample bp\n" +
				"log2 ratios, window level Wilcoxon FDRs, and window level pseudomedian log2 ratios. \n" +
				"Note, non-paired analysis are very underpowered and require > 30 obs/ window to see\n" +
				"any significant FDRs.\n"+

				"Required Options:\n"+
				"-s Path to a directory for saving the results.\n"+
				"-d Path to a directory containing individual sample PointData directories, each of\n"+
				"      which should contain chromosome split bar files (e.g. chr1.bar, chr2.bar, ...)\n"+
				"-t Names of the treatment sample directories in -d, comma delimited, no spaces.\n"+
				"-c Ditto but for the control samples, the ordering is critical and describes how to\n"+
				"      pair the samples for a paired analysis.\n"+

				"\nAdvanced Options:\n"+
				"-n Run a non-paired analysis where t and c are treated as groups and pooled.\n"+
				"-w Window size, defaults to 1000.\n"+
				"-o Minimum number observations in window, defaults to 10.\n" +
				"-p Minimum pseudomedian log2 ratio for estimating the permutation FDR, defaults to 0.2\n" +
				"-r Number permutations, defaults to 5\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/MethylationArrayScanner -s ~/MAS/Res\n" +
				"     -v H_sapiens_Feb_2009 -d ~/MAS/Bar/ -t Early1,Early2,Early3 -c Late1,Late2,Late3\n" +
				"     -w 1500\n\n" +

		"**************************************************************************************\n");
	}
}
