package edu.utah.seq.analysis;
import java.io.*;

import edu.utah.seq.data.FilterPointData;
import edu.utah.seq.data.PointData;
import edu.utah.seq.data.PointDataManipulator;
import edu.utah.seq.data.ReadCoverage;
import edu.utah.seq.data.SmoothingWindowInfo;
import edu.utah.seq.data.SubSamplePointData;
import edu.utah.seq.parsers.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.*;

import edu.utah.seq.useq.*;
import edu.utah.seq.useq.apps.Bar2USeq;
import util.bio.annotation.Coordinate;
import util.bio.parsers.UCSCGeneModelTableReader;
import util.gen.*;


/**Application for chaining together the many steps in a ChIP-Seq analysis.*/
public class ChIPSeq {

	//user defined fields
	private File rApplication = new File ("/usr/bin/R");
	private File resultsDirectory;
	private File windowFilterFile = null;
	private boolean convert2USeq = false;
	private float minimumFDR = 0.9f;
	private float maximumAlignmentScore = 60;
	private float minimumMappingQualityScore = 13;
	private String genomeVersion;
	private int minimumNumberReadsInWindow = 10;
	/*currently only sam, eland, novoalign, bed are supported*/
	private String alignmentType;
	private String alignmentTypeSam = "sam";
	private String alignmentTypeEland = "eland";
	private String alignmentTypeNovoalign = "novoalign";
	private String alignmentTypeBed = "bed";
	private File[] treatmentReplicaDirectories;
	private File[] controlReplicaDirectories;
	private boolean runMultipleReplicaAnalysis = true;
	private boolean findReducedRegions = true;
	private boolean verbose = false;
	private boolean bypassVarianceOutlierFiltering = false;

	//internal fields
	private File readCoverageTracks;
	private boolean skipParsingPointData = false;
	private File[] treatmentPointData;
	private File[] controlPointData;

	//ScanSeqs
	private int peakShift = -1;
	private int windowSize = -1;
	private int defaultPeakShift = 150;
	private int defaultWindowSize = 250;
	private SmoothingWindowInfo[] smoothingWindowInfo;
	private File swiFile;
	private File scanSeqs;

	//for MultipleReplicaScanSeqs, FDR, Log2Ratio
	private int[] scoreIndexes = new int[]{0,1};
	private float[] scoreThresholds = new float[]{10,1};


	public ChIPSeq (String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);

		//parse PointData?
		parsePointData();

		//filter duplicate reads
		filterDuplicates();

		//make ReadCoverage tracks, one for each replica
		readCoverage();

		//attempt to find peakShift
		peakShift();

		//run MultipleReplicaScanSeqs
		scanSeqs();

		//run EnrichedRegionMaker with various blocked exons
		enrichedRegionMaker();

		//Clean up
		cleanUp();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}



	public void peakShift(){
		System.out.println("\n*******************************************************************************");

		if (peakShift !=-1 && windowSize !=-1) {
			System.out.println("Skipping peak shift estimation since you explicity defined the peak shift and window size.");
			return;
		}
		String name = "ScanSeqs";
		if (runMultipleReplicaAnalysis) name = "MultipleReplicaScanSeqs";

		scanSeqs = new File (resultsDirectory, name);
		if (scanSeqs.exists() && scanSeqs.isDirectory()){
			System.out.println("\t"+name+" folder exists, Peak Shift estimation not needed.  Delete it to reprocess.");
			return;
		}
		System.out.println("Using the PeakShiftFinder to estimate the peak shift in your data...");

		File results = new File (resultsDirectory, "PeakShiftFinder");
		results.mkdir();
		PeakShiftFinder p = new PeakShiftFinder(treatmentPointData, controlPointData, results, true);
		int[] peakShiftAndWindowSize = p.fetchPeakShiftAndWindowSize();
		if (peakShiftAndWindowSize == null) {
			if (peakShift == -1) peakShift = defaultPeakShift;
			if (windowSize == -1) windowSize = defaultWindowSize;
			System.out.println("\nCould not robustly identify the peak shift in your data, using the default "+peakShift+"bp peak shift and "+windowSize+"bp window size.");
		}
		else {
			if (peakShift == -1) peakShift = peakShiftAndWindowSize[0];
			if (windowSize == -1) windowSize = peakShiftAndWindowSize[1];
			System.out.println("\nUsing a "+peakShift+"bp peak shift and "+windowSize+"bp window size.");
		}
	}

	public void cleanUp(){
		System.out.println("\n*******************************************************************************");
		System.out.println("Converting data...");
		if (convert2USeq){
			//convert ReadCoverageTracks
			if (readCoverageTracks!= null) new Bar2USeq(readCoverageTracks, true);
			//convert ScanSeqs
			if (scanSeqs!= null) new Bar2USeq(scanSeqs,true);
		}
	}

	public void enrichedRegionMaker(){
		System.out.println("\n*******************************************************************************");
		System.out.println("Building lists of chIP-seq peaks (enriched and reduced) from the window data using the EnrichedRegionMaker...");
		if (swiFile == null){
			if (runMultipleReplicaAnalysis) System.out.print("\tWARNING: MultipleReplicaScanSeqs ");
			else {
				System.out.print("\tWARNING: ScanSeqs ");
			}
			System.out.println("was not run or produced no significant windows, skipping EnrichedRegionMaker.  Delete it to reprocess.");
			return;
		}
		if (runMultipleReplicaAnalysis == false){
			//reset thresholds and indexes
			scoreIndexes = new int[]{1,4};
			scoreThresholds = new float[]{30,1};
		}
		new EnrichedRegionMaker(smoothingWindowInfo, swiFile, windowSize, treatmentPointData, controlPointData, false, windowFilterFile, 0, rApplication, scoreIndexes, scoreThresholds, verbose);
		if (findReducedRegions) new EnrichedRegionMaker(smoothingWindowInfo, swiFile, windowSize, treatmentPointData, controlPointData, true, windowFilterFile, 0, rApplication, scoreIndexes, scoreThresholds, verbose);
	}

	public void scanSeqs(){
		System.out.println("\n*******************************************************************************");

		if (runMultipleReplicaAnalysis){
			System.out.println("Scanning chromosomes for differential enrichment/reduction with MultipleReplicaScanSeqs...");
			scanSeqs = new File (resultsDirectory, "MultipleReplicaScanSeqs");
			if (scanSeqs.exists() && scanSeqs.isDirectory()){
				System.out.println("\tWARNING: MultipleReplicaScanSeqs folder exists, skipping analysis.  Delete it to reprocess.");
				return;
			}
			scanSeqs.mkdir();
			MultipleReplicaScanSeqs ss = new MultipleReplicaScanSeqs(treatmentPointData, controlPointData, scanSeqs, rApplication, windowSize, peakShift, minimumNumberReadsInWindow, minimumFDR, bypassVarianceOutlierFiltering, verbose);
			swiFile = ss.getSwiFile();
		}
		else {
			System.out.println("Scanning chromosomes for differential enrichment/reduction with ScanSeqs...");
			scanSeqs = new File (resultsDirectory, "ScanSeqs");
			if (scanSeqs.exists() && scanSeqs.isDirectory()){
				System.out.println("\tWARNING: ScanSeqs folder exists, skipping analysis.  Delete it to reprocess.");
				return;
			}
			scanSeqs.mkdir();
			ScanSeqs ss = new ScanSeqs(treatmentPointData, controlPointData, scanSeqs, rApplication, windowSize, peakShift, minimumNumberReadsInWindow, findReducedRegions, verbose);
			swiFile = ss.getSwiFile();	
		}

		//must reload since object is modified when making graphs
		if (swiFile != null && swiFile.exists()) smoothingWindowInfo = (SmoothingWindowInfo[])IO.fetchObject(swiFile);
	}

	public void readCoverage(){
		System.out.println("\n*******************************************************************************");
		System.out.println("Making relative ReadCoverage tracks (the number of reads per million mapped that intersect each bp) from the PointData...");
		readCoverageTracks = new File (resultsDirectory, "ReadCoverageTracks");
		if (readCoverageTracks.exists() && readCoverageTracks.isDirectory()){
			System.out.println("\tWARNING: ReadCoverage folder exists, skipping track generation.  Delete it to reprocess.");
			return;
		}
		readCoverageTracks.mkdir();
		for (int i=0; i< treatmentPointData.length; i++){
			File rc = new File(readCoverageTracks, "T_"+treatmentPointData[i].getName());
			rc.mkdir();
			new ReadCoverage(rc, new File[]{treatmentPointData[i]}, false);
		}
		if (treatmentPointData.length !=1) {
			File all = new File(readCoverageTracks, "T_Combine");
			all.mkdir();
			new ReadCoverage(all, treatmentPointData, false);
		}
		for (int i=0; i< controlPointData.length; i++){
			File rc = new File(readCoverageTracks, "C_"+controlPointData[i].getName());
			rc.mkdir();
			new ReadCoverage(rc, new File[]{controlPointData[i]}, false);
		}
		if (controlPointData.length !=1) {
			File all = new File(readCoverageTracks, "C_Combine");
			all.mkdir();
			new ReadCoverage(all, controlPointData, false);
		}
	}

	public void parsePointData(){
		System.out.println("\n*******************************************************************************");
		System.out.println("Parsing PointData from raw alignments...");
		File pointDataDirectory = new File(resultsDirectory, "PointData");
		if (pointDataDirectory.exists()){
			skipParsingPointData = true;
			System.out.println("\tWARNING: PointData directory exists, skipping parsing and using files within.  Delete it to reprocess.");
		}
		pointDataDirectory.mkdir();
		File t = new File (pointDataDirectory, "TreatmentPointData");
		if (skipParsingPointData == false) System.out.println("\tParsing treatment PointData...");
		treatmentPointData = parsePointData(t, treatmentReplicaDirectories);
		File c = new File (pointDataDirectory, "ControlPointData");
		if (skipParsingPointData == false) System.out.println("\tParsing control PointData...");
		controlPointData = parsePointData(c, controlReplicaDirectories);
	}

	public void filterDuplicates(){
		System.out.println("\n*******************************************************************************");
		System.out.println("Filtering PointData for duplicate reads (same strand and start position)...");
		File dupFilt = new File(resultsDirectory, "DupFiltPointData");
		if (dupFilt.exists()){
			System.out.println("\tWARNING: DupFiltPointData directory exists, skipping filtering.  Delete it to reprocess.");
			//add files
			treatmentPointData = IO.extractOnlyDirectories(new File (dupFilt, "TreatmentPointData"));
			controlPointData = IO.extractOnlyDirectories(new File (dupFilt, "ControlPointData"));
			if (treatmentPointData.length !=0 && controlPointData.length !=0) return;
			else System.out.println("\tNo filtered data found, attempting to re filter PointData?!");
		}
		dupFilt.mkdir();
		File t = new File (dupFilt, "TreatmentPointData");
		t.mkdir();
		for (int i=0; i< treatmentPointData.length; i++){			
			File repDir = new File(t,"Rep"+i);
			repDir.mkdir();
			PointDataManipulator pdm = new PointDataManipulator(new File[]{treatmentPointData[i]}, repDir);
			System.out.println("\t\tT "+repDir.getName()+"\t"+Num.formatPercentOneFraction(pdm.fetchFractionUniqueReads())+" Unique");
			//reassign
			treatmentPointData[i] = repDir;
		}
		//for each unfiltered treatment pd dir
		File c = new File (dupFilt, "ControlPointData");
		c.mkdir();
		for (int i=0; i< controlPointData.length; i++){
			File repDir = new File(c,"Rep"+i);
			repDir.mkdir();
			PointDataManipulator pdm = new PointDataManipulator(new File[]{controlPointData[i]}, repDir);
			System.out.println("\t\tC "+repDir.getName()+"\t"+Num.formatPercentOneFraction(pdm.fetchFractionUniqueReads())+" Unique");
			//reassign
			controlPointData[i] = repDir;
		}
	}

	public File[] parsePointData (File saveDir, File[] replicaDirectories){

		saveDir.mkdir();
		//make PointData
		File[] pd = new File[replicaDirectories.length];
		//for each replica directory
		for (int i=0; i< replicaDirectories.length; i++){
			//make a dir to same pd
			File repDir = new File (saveDir,"Rep"+i);
			repDir.mkdir();
			pd [i] = new File(repDir,"PointData");
			if (skipParsingPointData) continue;
			System.out.print("\t"+repDir+"\t");
			//parse data
			if (alignmentType.equals(alignmentTypeEland)){
				pd[i].mkdir();
				ElandParser p = new ElandParser(pd[i], IO.extractFiles(replicaDirectories[i]), maximumAlignmentScore,genomeVersion);
				System.out.println(p.getTotalNumMatch());
			}
			else if (alignmentType.equals(alignmentTypeSam)){
				SamParser p = new SamParser(repDir, IO.extractFiles(replicaDirectories[i]), minimumMappingQualityScore, maximumAlignmentScore, genomeVersion);
				System.out.println(p.getNumberPassingAlignments());
			}
			else if (alignmentType.equals(alignmentTypeNovoalign)){
				NovoalignParser p = new NovoalignParser(repDir, IO.extractFiles(replicaDirectories[i]), minimumMappingQualityScore, maximumAlignmentScore, genomeVersion, null);
				System.out.println(p.getNumberPassingAlignments());
			}
			else if (alignmentType.equals(alignmentTypeBed)){
				Tag2Point p = new Tag2Point(pd[i], IO.extractFiles(replicaDirectories[i]), genomeVersion);
				System.out.println(p.getNumberParsedRegions());
			}
		}
		return pd;
	}

	public void checkArgs(){
		System.out.println("\nChecking parameters...");
		boolean passed = true;
		StringBuilder notes = new StringBuilder();

		//find results directory
		if (resultsDirectory != null){
			if (resultsDirectory.exists()) notes.append("Your save results directory exits, may overwrite the files within.\n");
			else if (resultsDirectory.mkdirs() == false){
				notes.append("\tCannot create your results directory? Does the parent directory exist?\n");
				passed = false;
			}
		}
		else {
			notes.append("Please enter a directory in which to save your results.\n");
			passed = false;
		}

		//check for R and required libraries
		if (rApplication == null || rApplication.canExecute()== false) {
			notes.append("Cannot find or execute the R application -> "+rApplication+"\n");
			passed = false;
		}
		else {
			//DESeq?
			if (this.runMultipleReplicaAnalysis){
				String errors = IO.runRCommandLookForError("library(DESeq)", rApplication, resultsDirectory);
				if (errors == null || errors.length() !=0){
					passed = false;
					notes.append("\nError: Cannot find the required R library.  Did you install DESeq " +
							"(http://www-huber.embl.de/users/anders/DESeq/)?  See the author's websites for installation instructions. Once installed, " +
							"launch an R terminal and type 'library(DESeq)' to see if it is present. Error message:\n\t\t"+errors+"\n\n");
				}
			}
			else {
				String errors = IO.runRCommandLookForError("library(qvalue)", rApplication, resultsDirectory);
				if (errors == null || errors.length() !=0){
					Misc.printExit("\nError: Cannot find the required R library.  Did you install qvalue " +
							"(http://genomics.princeton.edu/storeylab/qvalue/)?  See the author's websites for installation instructions. Once installed, " +
							"launch an R terminal and type 'library(qvalue)' to see if it is present. R error message:\n\t\t"+errors+"\n\n");
				}
			}
		}

		//alignment type
		if (alignmentType != null){
			alignmentType = alignmentType.toLowerCase();
			if (alignmentType.startsWith("s")) alignmentType = alignmentTypeSam;
			else if  (alignmentType.startsWith("e")) alignmentType = alignmentTypeEland;
			else if  (alignmentType.startsWith("n")) alignmentType = alignmentTypeNovoalign;
			else if  (alignmentType.startsWith("b")) alignmentType = alignmentTypeBed;
			else {
				notes.append("Your alignment type does not appear to match any of the recognized formats? \n");
				passed = false;
			}
		}
		else {
			notes.append("Please indicate what type of alignments you are providing. Either "+alignmentTypeSam+", "+alignmentTypeEland+", or "+alignmentTypeNovoalign+" .\n");
			passed = false;
		}

		//genome version
		if (genomeVersion != null){
			//check to see if verisonedGenome follows form 
			if (ArchiveInfo.DAS2_VERSIONED_GENOME_FORM.matcher(genomeVersion).matches() == false) {
				notes.append("WARNING! Your versioned genome does not follow recommended form (e.g. H_sapiens_Mar_2006) -> "+genomeVersion+"\n");
			}
		}
		else {
			notes.append("\tPlease provide a versioned genome.\n");
			passed = false;
		}

		//check filter regions file
		if (windowFilterFile !=null){
			if (windowFilterFile.canRead() == false) {
				notes.append("\tCannot find your filter regions bed file? .\n");
				passed = false;
			}

		}

		//check alignment directories 
		if (treatmentReplicaDirectories == null || controlReplicaDirectories == null){
			notes.append("Please enter, at minimum, one or more treatment and control alignment file directories.\n");
			passed = false;
		}
		else {
			//only one directory look deeper
			if (treatmentReplicaDirectories.length == 1){
				File[] otherDirs = IO.extractOnlyDirectories(treatmentReplicaDirectories[0]);
				if (otherDirs != null && otherDirs.length > 0) treatmentReplicaDirectories = otherDirs;
			}
			if (controlReplicaDirectories.length == 1){
				File[] otherDirs = IO.extractOnlyDirectories(controlReplicaDirectories[0]);
				if (otherDirs != null && otherDirs.length > 0) controlReplicaDirectories = otherDirs;
			}
			for (int i=0; i<treatmentReplicaDirectories.length; i++){
				if (treatmentReplicaDirectories[i].isDirectory() == false){
					notes.append("This treatment directory can't be found or isn't a directory? -> "+treatmentReplicaDirectories[i]+"\n");
					passed = false;
				}
			}
			for (int i=0; i<controlReplicaDirectories.length; i++){
				if (controlReplicaDirectories[i].isDirectory() == false){
					notes.append("This control directory can't be found or isn't a directory? -> "+controlReplicaDirectories[i]+"\n");
					passed = false;
				}
			}
		}

		//print out statement
		if (passed){
			printParameters();
			if (notes.length()!=0) System.out.print("\nNotes: "+notes);
		}
		else {
			Misc.printErrAndExit("\nThe following problems were encountered when processing your parameter file.  Correct and restart. -> \n"+notes);
		}

	}

	public void printParameters(){
		System.out.println("\tApp Version = "+IO.fetchUSeqVersion());
		System.out.println("\tResults directory = "+resultsDirectory);
		System.out.println("\tWindow filter file = "+windowFilterFile);
		System.out.println("\tAlignment type = "+alignmentType);
		System.out.println("\tMaximum alignment score = "+maximumAlignmentScore);
		System.out.println("\tMinimum mapping quality score = "+minimumMappingQualityScore);
		System.out.println("\tGenome version = "+genomeVersion);
		System.out.println("\tTreatment replica directories:\n\t\t"+IO.concatinateFileFullPathNames(treatmentReplicaDirectories, "\n\t\t"));
		System.out.println("\tControl replica directories:\n\t\t"+IO.concatinateFileFullPathNames(controlReplicaDirectories, "\n\t\t"));
		System.out.println("\tConvert bar graph files to useq format = "+convert2USeq);
		System.out.println("\tRun multiple replica DESeq analysis = "+runMultipleReplicaAnalysis);
		System.out.println("\tFind reduced regions = "+findReducedRegions);
		System.out.println("\tMinimum number reads in each window = "+minimumNumberReadsInWindow);
	}

	public static void main(String[] args) {
		if (args.length<1){
			printDocs();
			System.exit(0);	
		}
		new ChIPSeq(args);
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
					case 's': resultsDirectory = new File(args[++i]); break;
					case 't': treatmentReplicaDirectories = IO.extractFiles(args[++i]); break;
					case 'c': controlReplicaDirectories = IO.extractFiles(args[++i]); break;
					case 'y': alignmentType = args[++i]; break;
					case 'v': genomeVersion = args[++i]; break;
					case 'f': windowFilterFile = new File(args[++i]); break;
					case 'r': rApplication = new File(args[++i]); break;
					case 'a': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'd': minimumFDR = Float.parseFloat(args[++i]); break;
					case 'q': minimumMappingQualityScore = Float.parseFloat(args[++i]); break;
					case 'u': convert2USeq = true; break;
					case 'g': verbose = true; break;
					case 'e': findReducedRegions = false; break;
					case 'm': runMultipleReplicaAnalysis = false; break;
					case 'b': bypassVarianceOutlierFiltering = true; break;
					case 'i': minimumNumberReadsInWindow = Integer.parseInt(args[++i]); break;
					case 'w': windowSize = Integer.parseInt(args[++i]);break;
					case 'p': peakShift = Integer.parseInt(args[++i]);break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		checkArgs();
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                  ChIPSeq: April 2012                             **\n" +
				"**************************************************************************************\n" +
				"The ChIPSeq application is a wrapper for processing ChIP-Seq data through a variety of\n" +
				"USeq applications. It:\n" +
				"   1) Parses raw alignments (sam, eland, bed, or novoalign) into binary PointData\n" +
				"   2) Filters PointData for duplicate alignments\n"+
				"   3) Makes relative ReadCoverage tracks from the PointData (reads per million mapped)\n" +
				"   4) Runs the PeakShiftFinder to estimate the peak shift and optimal window size\n" +
				"   5) Runs the MultipleReplicaScanSeqs to window scan the genome generating enrichment\n" +
				"        tracks using DESeq's negative binomial pvalues and B&H's FDRs\n" +
				"   6) Runs the EnrichedRegionMaker to identify likely chIP peaks (FDR < 1%, >2x).\n" +

				"\nNote, given R's poor memory management, this app requires 64bit R and >6-8G RAM.\n"+

				"\nOptions:\n"+
				"-s Save directory, full path.\n"+
				"-t Treatment alignment file directories, full path, comma delimited, no spaces, one\n" +
				"       for each biological replica. These should each contain one or more text\n" +
				"       alignment files (gz/zip OK) for a particular replica. Alternatively, provide\n" +
				"       one directory that contains multiple alignment file directories.\n" +
				"-c Control alignment file directories, ditto. \n" +
				"-y Type of alignments, either novoalign, sam, bed, or eland (sorted or export).\n"+
				"-v Genome version (e.g. H_sapiens_Feb_2009, M_musculus_Jul_2007), see UCSC FAQ,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-r Full path to 64bit R, defaults to '/usr/bin/R'. Be sure to install Ander's DESeq\n" +
				"       (http://www-huber.embl.de/users/anders/DESeq/) R library.\n"+

				"\nAdvanced Options:\n"+
				"-m Combine any replicas and run single replica analysis (ScanSeqs), defaults to\n" +
				"      using DESeq.\n"+
				"-d Minimum FDR threshold for filtering windows, defaults to 0.9\n"+
				"-a Maximum alignment score. Defaults to 60, smaller numbers are more stringent.\n"+
				"-q Minimum mapping quality score. Defaults to 13, bigger numbers are more stringent.\n" +
				"      This is a phred-scaled posterior probability that the mapping position of read\n" +
				"      is incorrect. Set to 0 for RNASeq data.\n" +
				"-p Peak shift, defaults to the PeakShiftFinder peak shift or 150bp. Set to 0 for\n" +
				"      RNASeq data.\n"+
				"-w Window size, defaults to the PeakShiftFinder peak shift + stnd dev or 250bp.\n"+
				"-i Minimum number reads in window, defaults to 10.\n"+
				"-u Convert bar graph folders to xxx.useq format.\n"+
				"-f Filter bed file (tab delimited: chr start stop) to use in excluding intersecting\n" +
				"      windows while making peaks, e.g. satelliteRepeats.bed .\n"+
				"-g Print verbose output from each application.\n"+
				"-e Don't look for reduced regions.\n"+
				"-b Bypass DESeq's variance outlier filtering. Recommended for first pass.\n"+
				

				"\n"+
				"Example: java -Xmx2G -jar pathTo/USeq/Apps/ChIPSeq -y eland -v D_rerio_Dec_2008 -t \n" +
				"      /Data/PolIIRep1/,/Data/PolIIRep2/ -c /Data/PolIINRep1/,/Data/PolIINRep2/ -s\n" +
				"      /Data/Results/WtVsNull -f /Anno/satelliteRepeats.bed\n\n" +

		"**************************************************************************************\n");		
	}
}
