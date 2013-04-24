package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import trans.main.WilcoxonSignedRankTest;
import trans.tpmap.WindowMaker;
import util.bio.annotation.Bed;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.gen.*;
import edu.utah.seq.analysis.multi.Condition;
import edu.utah.seq.analysis.multi.GeneCount;
import edu.utah.seq.analysis.multi.Replica;
import edu.utah.seq.data.ComparatorPointAscendingScore;
import edu.utah.seq.data.ComparatorPointPosition;
import edu.utah.seq.data.HeatMapMaker;
import edu.utah.seq.data.Info;
import edu.utah.seq.data.Point;
import edu.utah.seq.data.PointData;
import edu.utah.seq.data.SmoothingWindow;
import edu.utah.seq.data.SmoothingWindowInfo;
import edu.utah.seq.parsers.BarParser;
import edu.utah.seq.useq.apps.Bar2USeq;
import edu.utah.seq.useq.apps.Text2USeq;


/**
 * @author Nix
 * */
public class MethylationArrayScanSeqs {

	//fields
	private int windowSize = 100000;
	private int minimumNumberObservationsInWindow = 20;
	private HashMap<String, SamplePair[]> chromSamplePairs;
	private HashMap<String, int[]> chromPositions;
	private String chromosome;
	private int[] positions;
	private SamplePair[] pairs;
	private WindowMaker windowMaker; 
	private int[][] windows;
	private SmoothingWindow[] smoothingWindow;
	private SmoothingWindowInfo[] smoothingWindowInfo;
	private ArrayList<Float> randomPValues = new ArrayList<Float>();
	private ArrayList<SmoothingWindowInfo> smiAL = new ArrayList<SmoothingWindowInfo>();
	private ArrayList<SmoothingWindow> smAL = new ArrayList<SmoothingWindow>();
	private String[] scoreNames;
	private String[] scoreDescriptions;
	private String[] scoreUnits;
	private String versionedGenome;
	private File pseDir;
	private File fdrDir;
	private File saveDirectory;
	private int pseIndex = 0;
	private int pvalueFDRIndex = 1;
	private int permPValIndex = 2;
	private int numObsIndex = 3;
	private File dataDirectory;
	private String pairedSamples;
	private float badValue = -666.0f;
	private String treatmentPairedSamples;
	private String controlPairedSamples;
	

	//constructor

	public MethylationArrayScanSeqs(String[] args){	
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
		makeSamplePairs();

		//make window maker 
		windowMaker = new WindowMaker(windowSize,minimumNumberObservationsInWindow);

		//for each chromosome
		for (String chrom: chromSamplePairs.keySet()){
			//fetch data
			chromosome = chrom;
			positions = chromPositions.get(chromosome);
			pairs = chromSamplePairs.get(chromosome);

			//make windows
			windows = windowMaker.makeWindows(positions);
			if (windows.length == 0){
				System.out.println("\t"+chromosome+" Skipping! No windows found with minimum reads of "+minimumNumberObservationsInWindow+" within a window size of "+windowSize);
				return;
			}
			System.out.println("\t"+chromosome);
			//make SmoothingWindow[] container
			smoothingWindow = new SmoothingWindow[windows.length];

			//scan windows
			scanWindows();
			saveWindowDataToSMIArrayList();

			//randomize pairs
//Misc.randomize(pairs, System.currentTimeMillis());
			
			//scan for random
//scanWindowsRandom();
		}
		
		smoothingWindowInfo = new SmoothingWindowInfo[smiAL.size()];
		smiAL.toArray(smoothingWindowInfo);
		
		//multiple test correct the pvalues with B&H
		System.out.println("Converting pvalues to FDRs with B&H...");
		convertPValuesToFDRs();

		//write out graph data
		System.out.println("\nSaving graph data...");
		writeBarFileGraphs();
		
		//convert graph data to useq
		new Bar2USeq(pseDir, true);
		new Bar2USeq(fdrDir, true);
		//new Bar2USeq(fdrDir, true);

		//save window data 
		System.out.println("Writing window objects for the EnrichedRegionMaker...");
		File swiFile = new File (saveDirectory, "windowData"+windowSize+"bp.swi");
		IO.saveObject(swiFile, smoothingWindowInfo);

	}
	
	public void makeSamplePairs(){
		//private HashMap<String, SamplePair[]> chromSamplePairs;
		//private HashMap<String, int[]> chromPositions;
	}
	
	private void convertPValuesToFDRs() {
		SmoothingWindow[] all = new SmoothingWindow[smAL.size()];
		smAL.toArray(all);
		
		//collect pvals for B&H corr
		Point[] pvals = new Point[all.length];
		for (int i=0; i< all.length; i++){
			pvals[i] = new Point(i, all[i].getScores()[pvalueFDRIndex]);
		}
		//sort by score
		Arrays.sort(pvals, new ComparatorPointAscendingScore());
		//correct
		Point.benjaminiHochbergCorrect(pvals, 0);
		//sort back to original position
		Arrays.sort(pvals, new ComparatorPointPosition());
		//assign FDRs to pVals
		for (int i=0; i< all.length; i++){
			float scores[] = all[i].getScores();
			scores[pvalueFDRIndex] = pvals[i].getScore();
		}
	}
	
	/*
	private void calcucontrolWindowFDRs(){
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
				float transPVal = randomScores[windowLength].calcucontrolTransformedPValue(scores[0]);
				
				//final scores = numberObs, pseMedian, pval
				scores = new float[]{windowLength, mean, transPVal};
				smoothingWindow[j].setScores(scores);	
			}
		}
	}*/

	/**Writes stair step window bar graph files*/
	public void writeBarFileGraphs(){

			//for each chromosome
			for (int i=0; i< smoothingWindowInfo.length; i++){
				Info info = smoothingWindowInfo[i].getInfo().copy();
				SmoothingWindow[] sm = smoothingWindowInfo[i].getSm();
				
				//final scores = mean, PVal, fdr, # obs		
				saveSmoothedHeatMapData (pseIndex, sm, info, pseDir, "#FF00FF"); //magenta
				saveSmoothedHeatMapData (pvalueFDRIndex, sm, info, fdrDir, "#00FF00"); //green
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

		//for each window 
		for (int i=0; i< windows.length; i++){
			treatmentAL.clear();
			controlAL.clear();

			//fetch scores from each sample pair
			for (int x=0; x< pairs.length; x++){
				pairs[x].fetchScores(windows[i][0], windows[i][1], treatmentAL, controlAL);
			}

			float[] treatment = Num.arrayListOfFloatToArray(treatmentAL);
			float[] control = Num.arrayListOfFloatToArray(controlAL);
			double[] fraction = Num.ratio(treatment, control);
			float pse = (float)Num.log2(Num.pseudoMedian(fraction));
			WilcoxonSignedRankTest w = new WilcoxonSignedRankTest(treatment, control);
			float pvalue = (float)w.getTransformedPValue();

			//scores = pse, pvalFDR, permPVal, # obs             
			float[] scores = new float[4];
			scores[pseIndex] = pse;
			scores[pvalueFDRIndex] = pvalue;
			scores[numObsIndex] = treatment.length;
			                           
			//make window
			smoothingWindow[i] = new SmoothingWindow (positions[windows[i][0]], positions[windows[i][1]]+1, scores);
			smAL.add(smoothingWindow[i]);
		}

	}

	private void scanWindowsRandom() {
		ArrayList<Float> treatmentAL = new ArrayList<Float>();
		ArrayList<Float> controlAL = new ArrayList<Float>();

		//for each window 
		for (int i=0; i< windows.length; i++){
			treatmentAL.clear();
			controlAL.clear();
			//fetch scores
			for (int x=windows[i][0]; x< windows[i][1]; x++){
				pairs[x].fetchScores(windows[i][0], windows[i][1], treatmentAL, controlAL);
			}
			float[] treatment = Num.arrayListOfFloatToArray(treatmentAL);
			float[] control = Num.arrayListOfFloatToArray(controlAL);
			WilcoxonSignedRankTest w = new WilcoxonSignedRankTest(treatment, control);
			float pvalue = (float)w.getTransformedPValue();
			randomPValues.add(pvalue);
		}

	}


	public void loadSamplePairs(){
		//parse names
		String[] treatmentNames = treatmentPairedSamples.split(",");
		String[] controlNames = controlPairedSamples.split(",");
		if (treatmentNames.length != controlNames.length) Misc.printErrAndExit("\nError: the number of treatment and control samples differ?");
		
		HashMap<String, ArrayList<SamplePair>> pairs = new HashMap<String, ArrayList<SamplePair>>();
		
		for (int i=0; i< treatmentNames.length; i++){
			HashMap<String, PointData> pdT = loadPointData(treatmentNames[i]);
			HashMap<String, PointData> pdC = loadPointData(controlNames[i]);
			
			//write log2(t/2) graphs
			writeProbeRatioGraph(pdT, pdC);
			
			HashSet<String> allChroms = new HashSet<String>();
			allChroms.addAll(pdT.keySet());
			allChroms.addAll(pdC.keySet());
			for (String chrom: allChroms){
				SamplePair sp = new SamplePair(pdT.get(chrom).getScores(), pdC.get(chrom).getScores());
				ArrayList<SamplePair> al = pairs.get(chrom);
				if (al == null){
					al = new ArrayList<SamplePair>();
					pairs.put(chrom, al);
				}
				al.add(sp);
			}
		}
		
		//convert to []s
		chromSamplePairs = new HashMap<String, SamplePair[]>();
		for (String chrom: pairs.keySet()){
			ArrayList<SamplePair> al = pairs.get(chrom);
			if (al.size() != treatmentNames.length) Misc.printErrAndExit("\nError: One or more samples are missing a chromosome PointData set.");
			SamplePair[] sp = new SamplePair[treatmentNames.length];
			al.toArray(sp);
			chromSamplePairs.put(chrom, sp);
		}
	}
	
	private void writeProbeRatioGraph(HashMap<String, PointData> pdT, HashMap<String, PointData> pdC) {
		// TODO Auto-generated method stub
		
	}

	class SamplePair{

		float[] treatment;
		float[] control;
		
		public SamplePair(float[] treatment, float[] control) {
			this.treatment = treatment;
			this.control = control;
		}

		public void fetchScores(int startIndex, int stopIndex, ArrayList<Float> treatmentAL, ArrayList<Float> controlAL){
			for (int i=startIndex; i< stopIndex; i++){
				if (treatment[i] !=0.0f && control[i] != 0.0f){
					treatmentAL.add(treatment[i]);
					controlAL.add(control[i]);
				}
			}
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
			int[] pdPos = pd.get(chrom).getPositions();
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
		new MethylationArrayScanSeqs(args);
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
		if (pairedSamples == null) Misc.printErrAndExit("\nError: please enter at least one paired sample set to contrast.\n");
		
		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		if (saveDirectory.exists() == false) saveDirectory.mkdir();
		
		pseDir = new File(saveDirectory, "Pse");
		pseDir.mkdir();
		fdrDir = new File(saveDirectory, "FDR");
		fdrDir.mkdir();

		scoreNames = new String[]{
				"Pse",
				"FDR",
				"permPVal",
				"#Obs",
		};
		scoreDescriptions = new String[]{
				"Pseudomedian of treatment/control ratios",
				"B&H Corrected Wilcoxon Signed Rank test p-value",
				"Permutation adjusted p-value",
				"Number paired observations",
		};
		scoreUnits = new String[]{
				"log2(pse(treatment/control ratios))",	
				"-10Log10(FDR)",
				"-10Log10(permPVal)",
				"count",
		};
	}	
	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Defined Region Differential Seq: Sept 2012                **\n" +
				"**************************************************************************************\n" +
				"DRDS takes bam files, one per replica, minimum one per condition, minimum two\n" +
				"conditions (e.g. treatment and control or a time course/ multiple conditions) and\n" +
				"identifies differentially expressed genes under any pairwise comparison using DESeq.\n" +
				"DESeq's variance corrected count data is used to heirachically cluster the \n" +
				"differentially expressed genes as well as the samples. See the DESeq manual for\n" +
				"details. Alternative splicing is estimated using a chi-square test of independence.\n" +
				"In addition to the cluster plots, a spread sheet is created with the pValue,\n" +
				"FDR, and variance corrected log2Ratios for each of the pairwise comparisons as well as\n" +
				"the raw, FPKM, and log2 variance corrected alignment counts.  Use the controlr for\n" +
				"subsequent clustering and distance estimations.\n"+

				"\nOptions:\n"+
				"-s Save directory.\n"+
				"-c Conditions directory containing one directory for each condition with one xxx.bam\n" +
				"       file per biological replica and their xxx.bai indexs. 3-4 reps recommended per\n" +
				"       condition. The BAM files should be sorted by coordinate using Picard's SortSam.\n" +
				"       All spice junction coordinates should be converted to genomic coordinates, see\n" +
				"       USeq's SamTranscriptomeParser.\n" +
				"-r Full path to R loaded with DESeq library, defaults to '/usr/bin/R' file, see\n" +
				"       http://www-huber.embl.de/users/anders/DESeq/ . Type 'library(DESeq)' in\n" +
				"       an R terminal to see if it is installed. \n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. Tab delimited, see RefSeq Genes\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (uniqueName1 name2(optional) chrom\n" +
				"       strand txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts\n" +
				"       (commaDelimited)exonEnds). Example: ENSG00000183888 C1orf64 chr1 + 16203317\n" +
				"       16207889 16203385 16205428 2 16203317,16205000 16203467,16207889 . NOTE:\n" +
				"       this table should contain only ONE composite transcript per gene (e.g. use\n" +
				"       Ensembl genes NOT transcripts). Use the MergeUCSCGeneTable app to collapse\n" +
				"       transcripts. See http://useq.sourceforge.net/usageRNASeq.html for details.\n"+
				"-b (Or) a bed file (chr, start, stop,...), full path, See,\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1\n"+
				"-g Genome Version  (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +

				"\nAdvanced Options:\n"+
				"-v Filter for variance outliers in DESeq, defaults to not filtering.\n"+
				"-m Mask overlapping gene annotations, recommended for well annotated genomes.\n"+
				"-x Max per base alignment depth, defaults to 50000. Genes containing such high\n"+
				"       density coverage are ignored. Warnings are thrown.\n"+
				"-n Max number repeat alignments. Defaults to all.  Assumes 'IH' tags have been set by\n" +
				"       processing raw alignments with the SamTranscriptomeProcessor.\n"+
				"-f Minimum FDR threshold for sorting, defaults to 10 (-10Log10(FDR=0.1)).\n"+
				"-l Minimum absolute varCorLog2Rto threshold for sorting, defaults to 1 (2x).\n"+
				"-e Minimum number alignments per gene/ region, defaults to 20.\n"+
				"-i Score introns instead of exons.\n"+
				"-p Perform a stranded analysis. Only collect reads from the same strand as the\n" +
				"      annotation.\n" +
				"-j Reverse stranded analysis.  Only collect reads from the opposite strand of the\n" +
				"      annotation.  This setting should be used for the Illumina's strand-specific dUTP protocol.\n" +
				"-k Second read flipped. This setting can be used to flip the strand of the second read in a pair.\n" +
				"      This setting makes it easier to view in IGB, but can break other downstream applications.\n" +
				"-a Perform a permutation based chi-square test for differential exon usage.\n"+
				"      Needs 4 or more replicas per condition.\n"+
				"-t Don't delete temp files (R script, R results, Rout, etc..).\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/DefinedRegionDifferentialSeq -c\n" +
				"      /Data/TimeCourse/ESCells/ -s /Data/TimeCourse/DRDS -g H_sapiens_Feb_2009\n" +
				"     -u /Anno/mergedHg19EnsemblGenes.ucsc.gz\n\n" +

		"**************************************************************************************\n");

	}
}
