package edu.utah.seq.analysis;
import java.io.*;
import edu.utah.seq.data.Sam2USeq;
import edu.utah.seq.data.SmoothingWindowInfo;
import edu.utah.seq.parsers.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.*;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

import edu.utah.seq.useq.*;
import edu.utah.seq.useq.apps.Bar2USeq;
import util.bio.annotation.Coordinate;
import util.bio.parsers.UCSCGeneModelTableReader;
import util.gen.*;


/**Application for chaining together the many steps in an RNA-Seq analysis.*/
public class RNASeq {

	//user defined fields
	private File rApplication = new File ("/usr/bin/R");
	private File resultsDirectory;
	private boolean stranded = false;
	private boolean performReverseStrandedAnalysis = false;
	private boolean secondStrandFlipped = false;
	private boolean bothFlipped = false;
	private File geneTableFile;
	private boolean filterGeneTable = true;
	private float maximumAlignmentScore = 120;
	private String genomeVersion;
	private File treatmentDirectory;
	private File controlDirectory;
	private File[] treatmentReplicaDirectories;
	private File[] controlReplicaDirectories;
	private float minimumFDR = 0.5f;
	private int maxMatches = 1;
	private boolean useDESeq = true;
	private boolean verbose = false;
	private int maxAlignmentsDepth = 50000;

	//internal fields
	private File filteredGeneTableFile;
	private File geneRegionFile;
	private File geneExonFile;
	private File geneIntronFile;
	private File readCoverageTracks;
	private boolean bamFilesWereGenerated = true;
	private File[] treatmentPointData;
	private File[] controlPointData;
	private File treatmentBamDirectory;
	private File controlBamDirectory;
	private File[] treatmentBamFiles;
	private File[] controlBamFiles;

	//ScanSeqs
	private int scanSeqsWindowSize = 150;
	private SmoothingWindowInfo[] smoothingWindowInfo;
	private File swiFile;
	private File scanSeqs;

	//EnrichedRegionMaker, QValFDR, Log2Ratio
	private int[] bpBuffers = new int[]{0,500, 5000};



	public RNASeq (String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);

		//convert to BAM?
		convertSamFiles();
		
		//parse gene table
		parseGeneTable();
		
		//parse PointData?
		parsePointData();
		
		//known gene analysis with exons and introns
		differentialExpression();
		
		//make ReadCoverage tracks, one for each replica
		readCoverage();

		//run window scanning for novels
		scanForNovels();

		//run EnrichedRegionMaker with various blocked exons
		enrichedRegionMaker();

		//Clean up
		cleanUp();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}

	public void convertSamFiles(){
		System.out.println("\n*******************************************************************************");
		System.out.println("Parsing raw sam files for unique alignments with the SamTranscriptomeParser....");
		try {
			boolean convertSam = true;
			boolean tBamMade = true;
			//Any prexisting treatment bam files?
			treatmentBamFiles = IO.extractFiles(treatmentBamDirectory, ".bam");
			if (treatmentBamFiles != null){
				//clear the existing bam files?
				if (treatmentBamFiles.length != treatmentReplicaDirectories.length) {
					IO.deleteFiles(treatmentBamFiles);
					IO.deleteFiles(IO.extractFiles(treatmentBamDirectory, ".bai"));
				}
				//nope ok so don't make em
				else {
					System.out.println("\tFound treatment bam files. Delete to reprocess.");
					convertSam = false;
					tBamMade = false;
				}
			}
			if (convertSam){
				treatmentBamFiles = new File[treatmentReplicaDirectories.length];
				//run STP on each treatmentReplicaDirectory
				for (int i=0; i< treatmentReplicaDirectories.length; i++){
					System.out.println("\n"+treatmentReplicaDirectories[i]);
					treatmentBamFiles[i] = new File (treatmentBamDirectory, treatmentReplicaDirectories[i].getName()+".bam");
					File[] samFiles = fetchSamFiles(treatmentReplicaDirectories[i]);
					SamTranscriptomeParser stp = new SamTranscriptomeParser(samFiles, treatmentBamFiles[i], maximumAlignmentScore, genomeVersion, maxMatches, verbose, secondStrandFlipped, bothFlipped);
					if (stp.getNumberPassingAlignments() < 1000) Misc.printErrAndExit("\nError: too few passing alignments?  Are these in sam format?\n");
				}
			}
			System.out.println();
			convertSam = true;
			boolean cBamMade = true;
			//Any prexisting control bam files?
			controlBamFiles = IO.extractFiles(controlBamDirectory, ".bam");
			if (controlBamFiles != null){
				//clear the existing bam files?
				if (controlBamFiles.length != controlReplicaDirectories.length) {
					IO.deleteFiles(controlBamFiles);
					IO.deleteFiles(IO.extractFiles(controlBamDirectory, ".bai"));
				}
				//nope ok so don't make em
				else {
					System.out.println("\tFound control bam files. Delete to reprocess.");
					convertSam = false;
					cBamMade = false;
				}
			}
			if (convertSam){
				controlBamFiles = new File[controlReplicaDirectories.length];
				//run STP on each controlReplicaDirectory
				for (int i=0; i< controlReplicaDirectories.length; i++){
					System.out.println(controlReplicaDirectories[i]);
					controlBamFiles[i] = new File (controlBamDirectory, controlReplicaDirectories[i].getName()+".bam");
					File[] samFiles = fetchSamFiles(controlReplicaDirectories[i]);
					SamTranscriptomeParser stp = new SamTranscriptomeParser(samFiles, controlBamFiles[i], maximumAlignmentScore, genomeVersion, maxMatches, verbose, secondStrandFlipped, bothFlipped);
					if (stp.getNumberPassingAlignments() < 1000) Misc.printErrAndExit("\nError: too few passing alignments?  Are these in sam format?\n");
				}
			}
			if (tBamMade == false && cBamMade == false) bamFilesWereGenerated = false;


		} catch (IOException e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing sam files!\n");
		}

	}

	/**Returns the sam files or an empty array.*/
	public static File[] fetchSamFiles (File dir){
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(dir,".sam");
		tot[1] = IO.extractFiles(dir,".sam.gz");
		tot[2] = IO.extractFiles(dir,".sam.zip");
		File[] f = IO.collapseFileArray(tot);
		
		if (f.length == 0) f = IO.extractFiles(dir,".gz");
		if (f.length == 0) Misc.printErrAndExit("\nError: failed to find any sam alignment files in -> "+dir);
		return f;
	}

	public void cleanUp(){
		System.out.println("\n*******************************************************************************");
		System.out.println("Compressing and possibly converting data...");
		//compress swi file
		if (swiFile !=null && swiFile.exists()) IO.zipAndDelete(swiFile);
	}

	public void differentialExpression(){
		System.out.println("\n*******************************************************************************");
		String name;
		if (useDESeq) {
			name = "DefinedRegionDifferentialSeq";
		}
		else {
			name = "DefinedRegionScanSeqs";
		}
		System.out.println("Scanning genes for differential expression with "+name+"...");
		File mrdrss = new File (resultsDirectory, name);
		
		//were bam files parsed if so then delete it
		if (bamFilesWereGenerated) {
			IO.deleteDirectory(mrdrss);
			mrdrss = new File(resultsDirectory, name);
			mrdrss.mkdirs();
		}
		else if (mrdrss.exists() && mrdrss.list().length !=0){
			System.out.println("\t"+name+" folder exists, skipping analysis.  Delete to reprocess.");
			return;
		}
		mrdrss.mkdir();
		File exons = new File (mrdrss, "Exons");
		exons.mkdir();
		File introns = new File (mrdrss, "Introns");
		introns.mkdir();

		System.out.println("\nExons....");
		
		if (useDESeq) {
			//score exons
			DefinedRegionDifferentialSeq o = new DefinedRegionDifferentialSeq(treatmentBamDirectory, controlBamDirectory, genomeVersion, exons, rApplication, filteredGeneTableFile, false, stranded, verbose, this.performReverseStrandedAnalysis, this.secondStrandFlipped, this.maxAlignmentsDepth);
			System.out.println("Introns....");
			//score introns
			o= new DefinedRegionDifferentialSeq(treatmentBamDirectory, controlBamDirectory, genomeVersion, introns, rApplication, filteredGeneTableFile, true, stranded, verbose, this.performReverseStrandedAnalysis, this.secondStrandFlipped, this.maxAlignmentsDepth);
	
		}
		else {
			//score exons
			DefinedRegionScanSeqs m = new DefinedRegionScanSeqs(treatmentPointData, controlPointData, exons, rApplication, filteredGeneTableFile, false, verbose);
			System.out.println("Introns....\n");
			//score introns
			m = new DefinedRegionScanSeqs(treatmentPointData, controlPointData, introns, rApplication, filteredGeneTableFile, true, verbose);
		}
	}

	public String fetchFractionExonicIntronicIntergenic(double exonic, double intronic, double total){
		double exonicFraction = exonic/ total;
		double intronicFraction = intronic/ total;
		double intergenic = total - (exonic + intronic);
		double intergenicFraction = intergenic / total;
		StringBuilder sb = new StringBuilder();
		sb.append(Num.formatNumber(exonicFraction, 3)+" ("+(int)exonic+"/"+(int)total+")\t");
		sb.append(Num.formatNumber(intronicFraction, 3)+" ("+(int)intronic+"/"+(int)total+")\t");
		sb.append(Num.formatNumber(intergenicFraction, 3)+" ("+(int)intergenic+"/"+(int)total+")");
		return sb.toString();
	}

	public void enrichedRegionMaker(){
		System.out.println("\n*******************************************************************************");
		System.out.println("Building lists of differentially expressed transfrags from the window data using the EnrichedRegionMaker...");
		int[] scoreIndexes;
		float[] scoreThresholds;
		String name;
		if (useDESeq) {
			name = "MultipleReplicaScanSeqs";
			scoreIndexes = new int[]{0,1};
			scoreThresholds = new float[]{20,1};
		}
		else {
			name = "ScanSeqs";
			scoreIndexes = new int[]{1,4};
			scoreThresholds = new float[]{20,1};
		}

		if (swiFile == null){
			System.out.println("\t"+name+" was not run, skipping EnrichedRegionMaker.  Delete to reprocess. ");
			return;
		}
		for (int i=0; i<bpBuffers.length; i++){
			File dummy = new File(swiFile.getParentFile(),bpBuffers[i]+"BPBuffer.swi");
			new EnrichedRegionMaker(smoothingWindowInfo, dummy, scanSeqsWindowSize, treatmentPointData, controlPointData, false, geneExonFile, bpBuffers[i], rApplication, scoreIndexes, scoreThresholds, verbose);
			new EnrichedRegionMaker(smoothingWindowInfo, dummy, scanSeqsWindowSize, treatmentPointData, controlPointData, true, geneExonFile, bpBuffers[i], rApplication, scoreIndexes, scoreThresholds, verbose);
		}
	}

	public void scanForNovels(){
		System.out.println("\n*******************************************************************************");
		String name;
		if (useDESeq) name = "MultipleReplicaScanSeqs";
		else name = "ScanSeqs";
		System.out.println("Scanning chromosomes for differential transcription with "+name+"...");
		scanSeqs = new File (resultsDirectory, name);
		
		//were bam files parsed if so then delete it
		if (bamFilesWereGenerated) {
			IO.deleteDirectory(scanSeqs);
			scanSeqs = new File(resultsDirectory, name);
			scanSeqs.mkdirs();
		}
		else if (scanSeqs.exists() && scanSeqs.isDirectory()){
			System.out.println("\t"+name+" folder exists, skipping analysis.  Delete to reprocess.");
			return;
		}
		scanSeqs.mkdir();
		
		if (useDESeq) {
			System.out.println();
			MultipleReplicaScanSeqs ss = new MultipleReplicaScanSeqs(treatmentPointData, controlPointData, scanSeqs, rApplication, scanSeqsWindowSize, 0, 10, minimumFDR, true, verbose);
			swiFile = ss.getSwiFile();
		}
		else {
			ScanSeqs ss = new ScanSeqs(treatmentPointData, controlPointData, scanSeqs, rApplication, scanSeqsWindowSize, 0, 10, true, verbose);
			swiFile = ss.getSwiFile();
		}
		//must reload since object is modified when making graphs
		smoothingWindowInfo = null;
		if (swiFile != null && swiFile.exists()) smoothingWindowInfo = (SmoothingWindowInfo[])IO.fetchObject(swiFile);
		
		//convert internal window bar dirs to useq format 
		if (scanSeqs!= null) new Bar2USeq(scanSeqs,true);

	}

	public void readCoverage(){
		System.out.println("\n*******************************************************************************");
		System.out.println("Making relative ReadCoverage tracks (the number of reads per million mapped that intersect each bp) from the bam alignments using the Sam2USeq app...");
		
		readCoverageTracks = new File (resultsDirectory, "RelativeReadCoverageTracks");
		
		//were bam files parsed if so then delete it
		if (bamFilesWereGenerated) {
			IO.deleteDirectory(readCoverageTracks);
			readCoverageTracks = new File(resultsDirectory, "RelativeReadCoverageTracks");
			readCoverageTracks.mkdirs();
		}
		
		File tDir = new File (readCoverageTracks, "Treatment");
		tDir.mkdirs();
		
		//for each treatment bam file 
		for (int i=0; i< treatmentBamFiles.length; i++) {
			String name = treatmentBamFiles[i].getName();
			name = name.substring(0, name.length()-4);
			File rc = new File (tDir, name+".useq");
			System.out.print("\t"+treatmentBamFiles[i]);
			if (rc.exists() == false) {
				//make track
				Sam2USeq s2u = new Sam2USeq (new File[]{treatmentBamFiles[i]}, rc, genomeVersion, true, stranded, true, verbose);
				System.out.println("\t"+s2u.getNumberPassingAlignments());
			}
			else System.out.println("\tAlready processed. Delete to reparse.");
		}
		
		//for all treatments
		File combineTreatment = new File (tDir, "treatment.useq");
		System.out.print("\t"+combineTreatment);
		if (combineTreatment.exists() == false){
			//make track
			Sam2USeq s2u = new Sam2USeq (treatmentBamFiles, combineTreatment, genomeVersion, true, stranded, true, verbose);
			System.out.println("\t"+s2u.getNumberPassingAlignments());
		}
		else System.out.println("\tAlready processed. Delete reparse.");
		
		System.out.println();
		
		//control
		File cDir = new File (readCoverageTracks, "Control");
		cDir.mkdirs();
		
		//for each control bam file 
		for (int i=0; i< controlBamFiles.length; i++) {
			String name = controlBamFiles[i].getName();
			name = name.substring(0, name.length()-4);
			File rc = new File (cDir, name+".useq");
			System.out.print("\t"+controlBamFiles[i]);
			if (rc.exists() == false) {
				//make track
				Sam2USeq s2u = new Sam2USeq (new File[]{controlBamFiles[i]}, rc, genomeVersion, true, stranded, true, verbose);
				System.out.println("\t"+s2u.getNumberPassingAlignments());
			}
			else System.out.println("\tAlready processed. Delete to reparse.");
		}
		
		//for all controls
		File combineControl = new File (cDir, "control.useq");
		System.out.print("\t"+combineControl);
		if (combineControl.exists() == false){
			//make track
			Sam2USeq s2u = new Sam2USeq (controlBamFiles, combineControl, genomeVersion, true, stranded, true, verbose);
			System.out.println("\t"+s2u.getNumberPassingAlignments());
		}
		else System.out.println("\tAlready processed. Delete to reparse.");
		
		
		
	}

	public void parseGeneTable(){
		System.out.println("\n*******************************************************************************");
		System.out.println("Parsing gene table for gene regions and exons...");
		//make directory to hold annotation
		File annotation = new File(resultsDirectory, "Annotation");
		annotation.mkdir();
		filteredGeneTableFile = new File(annotation, "geneModels.ucsc");
		geneRegionFile = new File (annotation, "geneRegions.bed");
		geneExonFile = new File (annotation, "geneExons.bed");
		geneIntronFile = new File (annotation, "geneIntrons.bed");
		if (filteredGeneTableFile.exists() && geneRegionFile.exists() && geneExonFile.exists()){
			System.out.println("\tAnnotation folder exists.  Delete to reprocess.");
			return;
		}

		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(geneTableFile, 0);
		
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your regions's coordinates are reversed. Check that each start is less than the stop.\n");
		//check gene name is unique
		if (reader.uniqueGeneNames() == false) Misc.printExit("\nDuplicate gene names were found in your gene / bed file, these must be unique.\n");
		//check that genes are stranded
		if (stranded && reader.checkStrand() == false) Misc.printExit("\nError: you have indicated to perform a stranded analysis yet one or more of your genes/ regions is unstranded, aborting.\n");


		if (filterGeneTable){
			System.out.print("\tRemoving overlapping exons");
			String deletedGenes = reader.removeOverlappingExons();
			if (deletedGenes.length() !=0) System.out.println("\t\tWARNING: the following genes had more than 1/2 of their exonic bps removed -> "+deletedGenes);
		}

		System.out.println("\tWriting gene regions, exons, introns, and unique gene models...");
		//write out gene table sans overlaps
		reader.writeGeneTableToFile(filteredGeneTableFile);
		//make gene region file
		Coordinate[] geneRegions = reader.fetchGeneRegions();
		Coordinate.writeToFile(geneRegions, geneRegionFile);
		//make exon file
		Coordinate[] exons = reader.fetchExons();
		Coordinate.writeToFile(exons, geneExonFile);
		//make intron file
		reader.swapIntronsForExons();
		Coordinate[] introns = reader.fetchExons();
		Coordinate.writeToFile(introns, geneIntronFile);
	}

	public void parsePointData(){
		System.out.println("\n*******************************************************************************");
		System.out.println("Parsing PointData from bam alignments...");
		File pointDataDirectory = new File(resultsDirectory, "PointData");

		//were bam files parsed if so then clear PointData dir
		if (bamFilesWereGenerated) {
			IO.deleteDirectory(pointDataDirectory);
			pointDataDirectory = new File(resultsDirectory, "PointData");
			pointDataDirectory.mkdirs();
		}

		File tDir = new File (pointDataDirectory, "Treatment");
		tDir.mkdirs();
		treatmentPointData = new File[treatmentBamFiles.length];
		//for each treatment bam file 
		for (int i=0; i< treatmentBamFiles.length; i++) {
			String name = treatmentBamFiles[i].getName();
			name = name.substring(0, name.length()-4);
			treatmentPointData[i] = new File (tDir, name);
			System.out.print("\t"+treatmentBamFiles[i]);
			if (treatmentPointData[i].exists() == false || IO.extractFiles(treatmentPointData[i], ".bar.zip").length == 0) {
				//remove any old files
				IO.deleteFiles(treatmentPointData[i].listFiles());
				//make pointData
				SamParser p = new SamParser(treatmentPointData[i], treatmentBamFiles[i], 0, maximumAlignmentScore, genomeVersion, verbose);
				System.out.println("\t"+p.getNumberPassingAlignments());
			}
			else System.out.println("\tAlready parsed. Delete to reparse.");
		}

		System.out.println();
		
		File cDir = new File (pointDataDirectory, "Control");
		cDir.mkdirs();
		controlPointData = new File[controlBamFiles.length];
		//for each control bam file 
		for (int i=0; i< controlBamFiles.length; i++) {
			String name = controlBamFiles[i].getName();
			name = name.substring(0, name.length()-4);
			controlPointData[i] = new File (cDir, name);
			System.out.print("\t"+controlBamFiles[i]);
			if (controlPointData[i].exists() == false || IO.extractFiles(controlPointData[i], ".bar.zip").length == 0) {
				//remove any old files
				IO.deleteFiles(controlPointData[i].listFiles());
				//make pointData
				SamParser p = new SamParser(controlPointData[i], controlBamFiles[i], 0, maximumAlignmentScore, genomeVersion, verbose);
				System.out.println("\t"+p.getNumberPassingAlignments());
			}
			else System.out.println("\tAlready parsed. Delete to reparse.");
		}

	}

	

	public void checkArgs(){
		System.out.println("\nChecking parameters...");
		boolean passed = true;
		StringBuilder notes = new StringBuilder();

		//check java
		/*
		if (IO.checkJava() == false){
			notes.append("\tYour java application is not >= 1.6 (type 'java -version' on the cmd line). Install the most recent java from http://www.java.com/en/download/ .\n");
			passed = false;
		}*/

		//find results directory
		if (resultsDirectory != null){
			if (resultsDirectory.exists()) notes.append("\tYour save results directory exits, may overwrite files within.\n");
			else if (resultsDirectory.mkdirs() == false){
				notes.append("\tCannot create your results directory? Does the parent directory exist?\n");
				passed = false;
			}
		}
		else {
			notes.append("\tPlease enter a directory in which to save your results.\n");
			passed = false;
		}

		//check for R and required libraries
		if (rApplication == null || rApplication.canExecute()== false) {
			notes.append("\tCannot find or execute the R application -> "+rApplication+"\n");
			passed = false;
		}
		else {
			String errors = IO.runRCommandLookForError("library(DESeq)", rApplication, resultsDirectory);
			if (errors == null || errors.length() !=0){
				passed = false;
				notes.append("\nError: Cannot find the required R library.  Did you install DESeq " +
						"(http://www-huber.embl.de/users/anders/DESeq/)?  See the author's websites for installation instructions. Once installed, " +
						"launch an R terminal and type 'library(DESeq)' to see if it is present. Error message:\n\t\t"+errors+"\n\n");
			}
			if (useDESeq == false){
				//look for qvalue library
				errors = IO.runRCommandLookForError("library(qvalue)", rApplication, resultsDirectory);
				if (errors == null || errors.length() !=0){
					passed = false;
					notes.append("\nError: Cannot find the required R library.  Did you install qvalue " +
							"(http://genomics.princeton.edu/storeylab/qvalue/)?  See the author's websites for installation instructions. Once installed, " +
							"launch an R terminal and type 'library(qvalue)' to see if it is present. R error message:\n\t\t"+errors+"\n\n");
				}
				
				//look for estimateDispersions() function
				boolean useEstimateDispersions = DefinedRegionDifferentialSeq.estimateDispersions(rApplication, resultsDirectory);
				if (useEstimateDispersions == false){
					Misc.printErrAndExit("\nError: Please upgrade DESeq to the latest version, see http://www-huber.embl.de/users/anders/DESeq/ \n");
				}
			}
		}

		//parse gene table
		if (geneTableFile == null || geneTableFile.canRead()== false) {
			notes.append("\tCannot find your geneTable -> "+geneTableFile+"\n");
			passed = false;
		}

		//genome version
		if (genomeVersion != null){
			//check to see if verisonedGenome follows form 
			if (ArchiveInfo.DAS2_VERSIONED_GENOME_FORM.matcher(genomeVersion).matches() == false) {
				notes.append("\tWARNING! Your versioned genome does not follow recommended form (e.g. H_sapiens_Mar_2006) -> "+genomeVersion+"\n");
			}
		}
		else {
			notes.append("\tPlease provide a versioned genome.\n");
			passed = false;
		}

		//check alignment directories 
		if (treatmentDirectory == null || treatmentDirectory.isDirectory() == false || controlDirectory == null || controlDirectory.isDirectory() == false){
			notes.append("\tPlease enter both a treatment and control directory containing replica directories of alignment files.\n");
			passed = false;
		}
		else {
			//parse treatment replica dirs
			treatmentReplicaDirectories = IO.extractOnlyDirectories(treatmentDirectory);
			if (treatmentReplicaDirectories == null || treatmentReplicaDirectories.length == 0) treatmentReplicaDirectories = new File[] {treatmentDirectory};

			//parse control replica dirs
			controlReplicaDirectories = IO.extractOnlyDirectories(controlDirectory);
			if (controlReplicaDirectories == null || controlReplicaDirectories.length == 0) treatmentReplicaDirectories = new File[] {controlDirectory};

			//make dirs to hold bam files
			File bamFiles = new File (resultsDirectory,"BamAlignments");
			bamFiles.mkdirs();
			treatmentBamDirectory = new File (bamFiles, "Treatment");
			treatmentBamDirectory.mkdirs();
			controlBamDirectory = new File (bamFiles, "Control");
			controlBamDirectory.mkdirs();
			
			//too few reps?
			int numReps = treatmentReplicaDirectories.length + controlReplicaDirectories.length;
			if (numReps < 6) notes.append("\tWarning: Too few biological replicas. DESeq tends to be very conservative when analyzing data with < 4 replicas.\n");
			
		}
		
		if (this.performReverseStrandedAnalysis && this.stranded) {
			notes.append("\tStranded and reverse-stranded are both set to true.  Please select one and try again\n");
			passed = false;
		}

		//print out statement
		if (passed){
			printParameters();
			if (notes.length()!=0) System.out.print("\nNotes:\n"+notes);
		}
		else {
			Misc.printErrAndExit("\nThe following problems were encountered when processing your parameters.  Correct and restart. -> \n"+notes);
		}

	}

	/**Does some minimal error checking on a bam alignment file.
	 * @return null if no problems, otherwise an error.*/
	public static String checkBamFile(File bamFile) {
		String message = null;
		SAMFileReader reader = null;
		Pattern oneTwoDigit = Pattern.compile("\\w{1,2}");
		try {
			reader = new SAMFileReader(bamFile);
			//index present?
			if (reader.hasBrowseableIndex() == false) return bamFile+ " does not have an associated xxx.bai index? Use Picard's SortSam to sort and index your parsed sam files.";
			//check sort order
			SAMFileHeader h = reader.getFileHeader();
			if (h.getSortOrder().compareTo(SAMFileHeader.SortOrder.coordinate) !=0) return  bamFile+" does not appear to be sorted by coordinate. Did you sort with SAMTools?  Use Picard's SortSam instead."; 
			//check that their chromosomes aren't 1,2,3, should be chr1, chr2, chr3
			List<SAMSequenceRecord> chroms = h.getSequenceDictionary().getSequences();
			StringBuilder badChroms = new StringBuilder();
			boolean badMito = false;
			for (SAMSequenceRecord r: chroms){
				if (oneTwoDigit.matcher(r.getSequenceName()).matches()) badChroms.append(r.getSequenceName()+" ");
				if (r.getSequenceName().equals("chrMT")) badMito = true;
			}
			//just warn
			if (badChroms.length() !=0) System.err.println("\t\tWarning: "+bamFile+" bam file contains chromosomes that are 1-2 letters/ numbers long. For DAS compatibility they should start with 'chr' for chromosomes and something longish for contigs/ unassembled segments, see -> "+badChroms+"\n");
			if (badMito) System.err.println("\t\tWarning: "+bamFile+" bam file contains alignments with a chrMT chromosome. For DAS compatibility convert it to chrM.");
			//read an alignment
			SAMRecordIterator it = reader.iterator();
			if (it.hasNext()) it.next();
			//clean up
			reader.close();
		} catch (Exception e){
			message = e.getMessage();
		} finally {
			if (reader != null) reader.close();
		}
		return message;
	}

	public void printParameters(){
		System.out.println("\tApp Version = "+IO.fetchUSeqVersion());
		System.out.println("\tResults directory = "+resultsDirectory);
		System.out.println("\tStranded = "+stranded);
		System.out.println("\tGene table = "+geneTableFile);
		System.out.println("\tRemove overlapping exons = "+filterGeneTable);
		System.out.println("\tMaximum alignment score = "+maximumAlignmentScore);
		System.out.println("\tMinimum window FDR = "+minimumFDR);
		System.out.println("\tGenome version = "+genomeVersion);
		System.out.println("\tRun multiple replica DESeq analysis = "+useDESeq);
		System.out.println("\tTreatment replicas:\n\t\t"+IO.concatinateFileFullPathNames(treatmentReplicaDirectories, "\n\t\t"));
		System.out.println("\tControl replicas:\n\t\t"+IO.concatinateFileFullPathNames(controlReplicaDirectories, "\n\t\t"));
	}

	public static void main(String[] args) {
		if (args.length<1){
			printDocs();
			System.exit(0);	
		}
		new RNASeq(args);
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
					case 's': resultsDirectory = new File(args[++i]); break;
					case 't': treatmentDirectory = new File(args[++i]); break;
					case 'c': controlDirectory = new File(args[++i]); break;
					case 'v': genomeVersion = args[++i]; break;
					case 'g': geneTableFile = new File(args[++i]); break;
					case 'r': rApplication = new File(args[++i]); break;
					case 'a': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'o': filterGeneTable = false; break;
					case 'd': minimumFDR = Float.parseFloat(args[++i]); break;
					case 'n': stranded = true; break;
					case 'e': verbose = true; break;
					case 'm': useDESeq = false; break;
					case 'j': performReverseStrandedAnalysis = true; break;
					case 'k': secondStrandFlipped = true; break;
					case 'b': bothFlipped = true; break;
					case 'x': maxAlignmentsDepth = Integer.parseInt(args[++i]); break;
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
				"**                                   RNASeq: July 2012                              **\n" +
				"**************************************************************************************\n" +
				"The RNASeq application is a wrapper for processing RNA-Seq data through a variety of\n" +
				"USeq applications. It uses the DESeq package for calling significant differential\n" +
				"expression.  3-4 biological replicas per condition are strongly recommended. See \n" +
				"http://useq.sourceforge.net/usageRNASeq.html for details constructing splice indexes,\n" +
				"aligning your reads, and building a proper gene (NOT transcript) table.\n\n" +

				"The pipeline:\n" +
				"   1) Converts raw sam alignments containing splice junction coordinates into genome\n "+
				"         coordinates outputting sorted bam alignemnts.\n"+
				"   2) Makes relative read depth coverage tracks.\n"+
				"   3) Scores known genes for differential exonic and intronic expression using DESeq\n" +
				"         and alternative splicing with a chi-square test.\n" +
				"   4) Identifies unannotated differentially expressed transfrags using a window\n" +
				"         scan and DESeq.\n" +

				"\nUse this application as a starting point in your transcriptome analysis.\n\n" +

				"Options:\n"+
				"-s Save directory, full path.\n"+
				"-t Treatment alignment file directory, full path.  Contained within should be one\n" +
				"       directory per biological replica, each containing one or more raw\n" +
				"       SAM (.gz/.zip OK) files.\n" +
				"-c Control alignment file directory, ditto.  \n" +
				"-n Data is stranded. Only analyze reads from the same strand as the annotation.\n"+
				"-j Reverse stranded analysis.  Only count reads from the opposite strand of the\n" +
				"       annotation.  This setting should be used for the Illumina's strand-specific dUTP protocol.\n" +
				"-k Second read flipped. This setting can be used to flip the strand of the second read in a pair.\n" +
				"       This setting makes it easier to view in IGB, but can break other downstream applications.\n" +
				"-b Reverse the strand of both pairs.  Use this option if you would like the orientation\n" +
				"      of the alignments to match the orientation of the annotation in Illumina stranded \n" +
				"      UTP sequencing.\n" +
				"-x Max per base alignment depth, defaults to 50000. Genes containing such high\n"+
				"       density coverage are ignored. Warnings are thrown.\n"+
				"-v Genome version (e.g. H_sapiens_Feb_2009, M_musculus_Jul_2007), see UCSC FAQ,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-g UCSC RefFlat or RefSeq gene table file, full path. Tab delimited, see RefSeq Genes\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (uniqueName1 name2(optional) chrom\n" +
				"       strand txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts\n" +
				"       (commaDelimited)exonEnds). Example: ENSG00000183888 C1orf64 chr1 + 16203317\n" +
				"       16207889 16203385 16205428 2 16203317,16205000 16203467,16207889 . NOTE:\n" +
				"       this table should contain only ONE composite transcript per gene (e.g. use\n" +
				"       Ensembl genes NOT transcripts). Use the MergeUCSCGeneTable app to collapse\n" +
				"       transcripts to genes. See the RNASeq usage guide for details.\n"+
				"-r Full path to R, defaults to '/usr/bin/R'. Be sure to install Ander's DESeq\n" +
				"       (http://www-huber.embl.de/users/anders/DESeq/) R library.\n"+
				
				
				"\nAdvanced Options:\n"+
				"-m Combine replicas and run single replica analysis using binomial based statistics,\n" +
				"       defaults to DESeq and a negative binomial test.\n"+
				"-a Maximum alignment score. Defaults to 120, smaller numbers are more stringent.\n"+
				"-d Minimum FDR threshold for filtering windows, defaults to 0.5\n"+
				"-o Don't delete overlapping exons from the gene table.\n"+
				"-e Print verbose output from each application.\n"+				

				"\n"+
				"Example: java -Xmx2G -jar pathTo/USeq/Apps/RNASeq -v D_rerio_Dec_2008 -t \n" +
				"      /Data/PolIIMut/ -c /Data/PolIIWT/ -s\n" +
				"      /Data/Results/MutVsWT -g /Anno/zv8Genes.ucsc \n\n" +

		"**************************************************************************************\n");		
	}

	public static LinkedHashMap <String,String> loadKeyValueFile(File file){
		LinkedHashMap<String,String> map = new LinkedHashMap<String,String>();
		String line = null;
		try{
			BufferedReader in = new BufferedReader(new FileReader(file));

			Pattern equal = Pattern.compile("(.+)\\s*=\\s*(.+)");
			while ((line = in.readLine())!=null){
				line = line.trim();
				if (line.startsWith("#") || line.length() ==0) continue;
				Matcher mat = equal.matcher(line);
				if (mat.matches()) map.put(mat.group(1), mat.group(2));
				else Misc.printErrAndExit("\tProblem parsing parameter file line -> "+line);
			}
		}catch(Exception e){
			System.out.println("Prob with loadin parameter file.\n\t"+line);
			e.printStackTrace();
		}
		return map;
	}
}
