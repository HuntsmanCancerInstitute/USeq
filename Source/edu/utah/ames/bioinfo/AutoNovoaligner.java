package edu.utah.ames.bioinfo;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.util.MultidimensionalCounter.Iterator;

import trans.tpmap.WindowMaker;

/**
 * This class contains methods to identify and align newly deposited fastq sequences 
 * against the appropriate genome build. It grabs a file containing full path to 
 * fastq files and genome build, creates a cmd.txt file in the user's Tomato jobs
 * directory (with my email address), soft links the fastq files there, runs the 
 * alignments using appropriate params with Tomato, and then populates an analysis report 
 * in GNomEx. Bisulphite files are split into smaller chunks, aligned individually, 
 * and then the SAM files are merged. Completion success is reported back to sender 
 * (me) via email. Uses Tomato's email warning options to report back if there are 
 * warnings and/or errors with alignments.
 * 
 * @author darren.ames@hci.utah.edu
 *
 */

public class AutoNovoaligner {
	
	//fields
	 private static String requester; 
	 private static String lab; 
	 //private static String labNum; //number used by Tomato to create analysis report for appropriate lab
	 private static String sequencingApplicationCode; 
	 private static String organism;
	 private static String requestNum; 
	 private static String sampleID; 
	 private static String singleOrPairedEnd; //RunType
	 private static String myEmail = "darren.ames@hci.utah.edu"; //change this? Who else needs notifications?
	 private static String freshDataReport = "/Users/darren/Desktop/novoalignerTestDir/freshDataReport.txt";
	 private static String lane;
	 private static String numberOfCycles;
	 private static String requestYear;
	 private static String genome;
	 private static String novoindex;
	 private static String buildCode;
	 private static String tomatoDataDir = "/Users/darren/Desktop/novoalignerTestDir/";
	 //private String fastqFileName;
	 private String smtpHostName = "mail.inscc.utah.edu"; //default
	 static String EMAILREGEX = "^#e\\s+(.+@.+)"; //email address pattern
	 static String LABNAMEREGEX = "\\w+\\s\\w+(?=\\WLab)"; //matches first and last name of lab using positive look
	 private String cmdFilePath;
	 private static String jobsFolder = "/Users/darren/Desktop/novoalignerTestDir/";
	 private static String experimentFolder = "/Repository/MicrarrayData/" + requestYear + "/" + requestNum + "/Fastq/";
	 private static File paramsFile;
	 //private long fastqReadyDate = 0;
	 private static File logFile;
	 private HashSet<String> emails = new HashSet<String>();
	 private static boolean doNotAlign = false;
	 private static boolean isGenomicDNA = false;
	 private static boolean isSingleEnd = false;
	 private static boolean isPairedEnd = false;
	 private static boolean isBisSeq = false;
	 private static boolean isRNASeq = false;
	 private static boolean isSmallRNA = false;
	 private static boolean sequenceFilesAreNew = false;
	 private static boolean makeCmdFile = false;
	 private static boolean softLinkSuccess = false;
	 private static boolean bisulphiteSplitSuccess = false;
	 private static boolean mergeSamFilesSuccess = false;
	 private static boolean hasNovoindex = true;
	 private static boolean makeNewNovoindex = false;
	 private static boolean alignmentSuccess = false;
	 private static boolean populateAnalysisReport = false;
	
	 /**
	  * ***Coerce this thing skip creation of cmd.txt when lacking novoindex or build info
	  * ***or other necessary info. Right now, it's throwing nulls into the cmdFile string
	  * ***rather than skipping.
	  */
	 
	//constructor
	public AutoNovoaligner(String[] args) {	
		processArgs(args);
	}
	
	public static void main(String[] args) throws Exception {
		
		//check for args
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
			AutoNovoaligner an = new AutoNovoaligner(args);
			an.parseFreshDataReport();
	}
	
	/**
	 * Defines alignment parameters based on fields in freshDataReport.
	 * Fetches mandatory args for cmd.txt file.
	 * @throws Exception 
	 */
	public void parseFreshDataReport() throws Exception {
		//create buffered reader to read fresh data report
		BufferedReader br = new BufferedReader(new FileReader(freshDataReport));
		String line;
		//skip first line in file
		br.readLine();
		
		//loop the read block until all lines in the file are read
		while ((line = br.readLine()) != null) {

			//split contents of tab-delimited file 
			String dataValue[] = line.split("\t");

			//instantiate sample object
			Sample s = new Sample(dataValue);
			
			//parse out "Lab" from lab name string
			Pattern p = Pattern.compile(".+(?=Lab)");
			Matcher m = p.matcher(s.getLab());
			while (m.find()) {
				s.setLab(m.group());
			}
			//assign contents to fields
			if (s.getOrganism() == "Unknown" || s.getOrganism() == "Other") {
				doNotAlign = true;
				continue;
			}
			String genome = dataValue[s.genomeIndex];
			s.setGenome(genome);
			//exit (continue) if genome index unspecified
			if (s.getGenome() == "None") {
				//System.exit(0);
				hasNovoindex = false;
				doNotAlign = true;
				continue;
			}
			//set sequencing application codes if other redundant names are used
			if (s.getSequencingApplicationCode() == "APP2" || 
					s.getSequencingApplicationCode() == "APP3") {
				s.setSequencingApplicationCode("MRNASEQ");
			}
			//TODO make sure to remove the extra "_" characters in these two names below when
			//building the appropriate index name
			else if (s.getSequencingApplicationCode() == "APP1" || 
					s.getSequencingApplicationCode() == "TDNASEQ" || 
					s.getSequencingApplicationCode() == "EXCAPNIM") {
				s.setSequencingApplicationCode("");
			}
			else if (s.getSequencingApplicationCode() == "APP4" || 
					s.getSequencingApplicationCode() == "DNASEQ" || 
					s.getSequencingApplicationCode() == "APP5" || 
					s.getSequencingApplicationCode() == "CHIPSEQ") {
				s.setSequencingApplicationCode("");
			}
			
			novoindex = parseGenomeIndex(s);
		}
		br.close();
	}
	
	/**
	 * 
	 * @param sourceStr
	 * @return
	 */
	public String manipulateStr(String sourceStr) {
		StringBuffer sb = new StringBuffer(sourceStr);
		
		return sourceStr;
	}
	
	/**
	 * Parses sample object fields to set appropriate novoindex name 
	 * @param sample
	 * @return
	 * @throws Exception 
	 */
	public static String parseGenomeIndex(Sample s) throws Exception {
		
		String BUILDCODE = "\\w+(?=;\\s\\w+;)"; //match only build code in sample genome string
		
		//parse build code from genome string
		if (s.getGenome() != null) {
			Pattern p = Pattern.compile(BUILDCODE);
			Matcher m = p.matcher(s.getGenome());
			while (m.find()) {
				s.setBuildCode(m.group());
				String build = s.getBuildCode();
			}
		}
		
		//cast numberOfCycles String to int so I can do some math on it
		//int i = Integer.parseInt(s.getNumberOfCycles());
		//int radius = (i - 4); //splice junction radius to search for in index name
		//cast int radius back to String
		//String parsedRadius = String.valueOf(radius);
		//pattern to match a normal (non-bisSeq) index
		//String NORMINDEX = build + "\\w+" + parsedRadius + "\\w+(\\.nov\\.illumina\\.nix)";
		//pattern to match all bisulphite indices (not exclusive)
		//String BISINDEXALL = build + "\\w+Lamb\\w+(\\.nov\\.bisulphite\\.nix)";
		
		//instantiate Directory object so I can use its methods
		//Directory d = new Directory(tomatoDataDir);
		//do {
			//if not bisulfite, do this
			//if (!s.getSequencingApplicationCode().equalsIgnoreCase("BISSEQ")) {
				//d.getFiles(s, NORMINDEX);
				//System.out.println(s.getNovoindex());
			//}
			//else {
				//if bisulfite, do that
				//Pattern p1 = Pattern.compile(BISINDEXALL, Pattern.CASE_INSENSITIVE);
			//} 
			//load index names into hash for downstream use
			loadNovoindex(s);
		return novoindex;
		//}
		//while (tomatoDataDir != null);
	}
	
	/**
	 * Loads a list of novoindices into a hash. First column (build_sequencingApplicationCode_radius) 
	 * is the key, second (shortened index name) is the value.
	 * @throws IOException 
	 */
	public static HashMap<String,String> loadNovoindex(Sample s) throws IOException{
		
		//**TODO**Check hashmap implementation in test.java to see if it would be more concise.***
		
		//create new hashmap, key and value are string
		HashMap<String,String> indices = new HashMap<String,String>();
		//pattern to match novoindices
		String SHORTINDEXNAME = "([a-z0-9]+_[A-Z]+)_{0,1}([0-9]+){0,1}";
		String path = tomatoDataDir + "/indices/";
		//create new directory object
		Directory d = new Directory(path);
		try{
			//create new file object
			File fs[] = null;
			File f = new File(path);
			if (f.isDirectory()) {
				fs = d.getFiles();
				//make sure directory isn't empty
				if (!(fs.length > 0)) {
					System.out.println("Error, no indices found");
					//System.exit(1);
				}
			}
			//TODO FIX THIS GARBAGE!!!
			//TODO this doesn't correctly populate the hash so a match isn't taking place.
			//TODO fix problem where it prints empty cmd files w/o an index match
			
			//populate the hash
			for (int i = 0; i < fs.length; i++) {
				String filename = fs[i].getName();
				//System.out.println(filename);
				Pattern p = Pattern.compile(SHORTINDEXNAME);
				Matcher m = p.matcher(filename);
				boolean found = m.find();
				if (found) {
					//System.out.println(m.group(0));
					String keyName = m.group(0);
					indices.put(keyName, m.group(0));
				}
				//TODO This throws array out of bounds error if fresh data report 
				//contains unsupported genome name
				if (found == false) {
					System.out.println("Problem finding appropriate index file in the directory.");
				}
			}
		}catch(Exception e){
			System.out.println("Problem loading the Hash");
			e.printStackTrace();
		}
		//check to see if no build code is specified in the report
		if (s.getBuildCode() == null) {
			//if empty, set key to null
			String key = null;
			hasNovoindex = false;
			doNotAlign = true;
		}
		else {
			//build short index name to use in hashmap
			String key = (s.getBuildCode() + "_" + s.getSequencingApplicationCode() 
					+ "_" + s.getNumberOfCycles()); 
			System.out.println(key);
			//check hashmap for key (short index name), set index code
			if (indices.containsKey(key)) {
				s.setIndexCode(indices.get(s.getIndexCode()));
			}
		}
		/*
		//iterate through the has and print key/value pairs
		java.util.Iterator<String> iterator = indices.keySet().iterator();
		while (iterator.hasNext()) {
			String k = iterator.next().toString();
			String v = indices.get(k).toString();
			System.out.println(k + "\t" + v);
		}
		*/
		createCmdFile(s);
		return indices;
	}

	/**
	 * Makes cmd.txt file to execute the alignment.
	 */
	public static File createCmdFile(Sample s) throws IOException {
		
		//create an executable cmd.txt file for submitting new job
		File cmdFile = new File(AutoNovoaligner.getCmdFileMessageGeneral(s));
		String jobDirPath = jobsFolder + s.getRequestNumber();
		//TODO Create one folder for all samples of same requestNumber, then subfolders for each sample
		//TODO Create new analysis number for all samples from same requestNumber, then put all samples
		//TODO in the same analysis number
		try {
			//first make sure an index exists before continuing
			if ((s.getGenome().isEmpty() || s.getOrganism().isEmpty()) || 
					(s.getBuildCode().isEmpty())) {
				//System.out.println("Missing genome or index info for this sample. Skipping sample...\n");
			}
			else {
				//check to see if job dir exists for requestNum
				
				File dir = new File(jobDirPath);
				//directory exists
				if (dir.exists()) {
					//create cmd.txt file
					FileWriter fw = new FileWriter(jobDirPath + "/" + "cmd.txt");
					BufferedWriter out = new BufferedWriter(fw);
					out.write(AutoNovoaligner.getCmdFileMessageGeneral(s));

					//set checkpoint
					makeCmdFile = true;
					
					//close output stream
					out.close();
				}
				else {
					//if dir does not exist
					if (!dir.exists()) {
						//make new directory
						dir.mkdir();
						//make new cmd.txt file
						FileWriter fw = new FileWriter(jobDirPath + "/" + "cmd.txt"); 
						BufferedWriter out = new BufferedWriter(fw);
						out.write(AutoNovoaligner.getCmdFileMessageGeneral(s));
						//set checkpoint
						makeCmdFile = true;
						//close output stream
						out.close();
					}
					else {
						makeCmdFile = false;
						doNotAlign = true;
						throw new Exception("Couldn't make cmd.txt or directory");
					}
				}
			}
		}
		catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
			System.exit(0);
		}
		return cmdFile;
	}
	
	/**
	 * Generates the body of the cmd.txt file using sample-specific params.
	 * Calls appropriate methods to do the actual populating based on whether
	 * sequencing application is bisulfite, miRNA/small RNA, genomic DNA, or other (mRNA...)
	 * 
	 * @param sample
	 * @return
	 * @throws IOException
	 */
	public static String getCmdFileMessageGeneral(Sample s) throws IOException {
		//set string for non-variable part of cmd.txt params
		String nonVariablePartOfCmdFile = "#e " + myEmail + "\n#l " + s.getLab() + "\n## Novoalignments of " 
						+ s.getRequestNumber() + " " + s.getProjectName() + "\n## Aligning to " 
						+ s.getBuildCode() + "\n\n@align -novoalign " 
						+ s.getAlignParams() + " -g " + s.getNovoindex();
		
		//instantiate new StringBuffer object for holding cmd.txt file's body
		StringBuffer sb = new StringBuffer();
		
		//don't generate cmd.txt file if appropriate novoindex isn't available
		while (s.getIndexCode() != null) {
			//bisulfite? 
			if (s.getSequencingApplicationCode().equalsIgnoreCase("APP6")) {
				//get the bisulfite alignment params
				AutoNovoaligner.getCmdFileMessageBisulfite(s);
				//build the command
				sb.append(nonVariablePartOfCmdFile);
			}
			//maybe mRNA or exome?
			else if (s.getSequencingApplicationCode().equalsIgnoreCase("APP2") || 
					s.getSequencingApplicationCode().equalsIgnoreCase("APP3") || 
					s.getSequencingApplicationCode().equalsIgnoreCase("MRNASEQ") || 
					s.getSequencingApplicationCode().equalsIgnoreCase("APP1") || 
					s.getSequencingApplicationCode().equalsIgnoreCase("TDNASEQ") || 
					s.getSequencingApplicationCode().equalsIgnoreCase("EXCAPNIM")) {
				//get standard alignment params
				AutoNovoaligner.getCmdFileMessageStdParams(s);
				//build the command
				sb.append(nonVariablePartOfCmdFile);
				}
			//how about genomic DNA or ChIP-Seq then?
			else if (s.getSequencingApplicationCode().equalsIgnoreCase("APP4") || 
					s.getSequencingApplicationCode().equalsIgnoreCase("DNASEQ") || 
					s.getSequencingApplicationCode().equalsIgnoreCase("APP5") || 
					s.getSequencingApplicationCode().equalsIgnoreCase("CHIPSEQ")) {
				//get standard alignment params
				AutoNovoaligner.getCmdFileMessageStdParams(s);
				//build the command
				sb.append(nonVariablePartOfCmdFile);
			}
			else if (s.getSequencingApplicationCode().equalsIgnoreCase("SMRNASEQ")) {
				//get small RNA alignment params
				AutoNovoaligner.getCmdFileMessageSmallRNA(s);
				//build the command
				sb.append(nonVariablePartOfCmdFile);
			}
			else {
				//it's something else not currently supported in AutoNovoaligner (exit and run custom)
				doNotAlign = true;
				continue;	
			}
		}
		return sb.toString();
	}
	
	//bisulfite alignment params
	public static String getCmdFileMessageBisulfite(Sample s) {
		s.setAlignParams(" [-o SAM -r Random -t 240 -h 120 -b 2] -i *.gz -p bisulphite -gzip\n");
		return null;	
	}
	
	//small RNA alignment params
	public static String getCmdFileMessageSmallRNA(Sample s) {
		//include 3' adapter stripping prior to alignment
		//this is NOT the standard Illumina Gex Adapter 2 used as default by novoalign
		s.setAlignParams(" [-o SAM -r All 50 -m -a \"ATCTCGTATGCCGTCTTCTGCTTG\" -l 15 -t 30] -i *.gz -gzip");
		return null;
	}
	
	//standard mRNA, DNA or exome alignment params
	public static String getCmdFileMessageStdParams(Sample s) {
		s.setAlignParams(" [-o SAM -r All 50] -i *.gz -gzip\n");
		return null; 
	}
	
	/**
	 * Parses the indexed column of a tab delimited file, all lines included.
	 * Will return an empty String if index column doesn't exist.
	 * @param inputFile
	 * @param index
	 * @return
	 */
	/*
	public static String[] parseColumn (String inputFile, int index){
		ArrayList al = new ArrayList();
		inputFile = "/Users/darren/Desktop/novoalignerTestDir/labNamesAndNumbers.txt";
		try{
			BufferedReader in = new BufferedReader (new FileReader(inputFile));
			String line;
			String[] tokens;
			while ((line=in.readLine()) != null){
				tokens = line.split("\\t");
				if (tokens.length <= index) al.add("");
				else al.add(tokens[index].trim());
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
			Misc.printExit("\nError: problem parsing matcher file, aborting.\n");
		}
		String[] col = new String[al.size()];
		al.toArray(col);
		return col;
	}
	*/
	
	/**
	 * Soft-link the new fastq files to the appropriate Tomato jobs dir.
	 * TODO: maybe unnecessary. I could just use fastq file path instead for simplicity.
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void linkFastqFiles() throws IOException, InterruptedException {
		String oldName = null;
		String newName = null;
		Process link = Runtime.getRuntime().exec(new String[] {"ln", "-s", oldName, newName});
		link.waitFor();
		link.destroy();
	}
	
	/**
	 * Break BisSeq files into smaller manageable bits.
	 */
	public static void chunkBisulphiteFiles() {
		
		FileSplitter split = new FileSplitter(null);
		if (isBisSeq == true) {
			//TODO: add stuff here so it breaks BisSeq fastq file into manageable chunks.
			//find filename and assign to Sample.setFastqFileName
			//use sampleID + "_(last two chars of requestYear)*.fastq.gz"
			
			/**
			 * modify FileSplitter so I can call the methods from here and specify params
			 * Do like this (from ScanSeqs)
			 *
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
			 */
		}
	}
	
	/**
	 * This method creates a new empty "b" file if the file with the same name does 
	 * not already exist. Returns true if the file with the same name did not 
	 * exist and it was created successfully, false otherwise. This starts the 
	 * Tomato job. 
	 */
	public static void runTomatoAlignments() {
		
		//create File object
		File file = new File(jobsFolder + requestNum + "/" + "b");
		boolean bFileCreated = false;
		try {
			bFileCreated = file.createNewFile();
			bFileCreated = true;
		} catch (IOException ioe) {
			bFileCreated = false;
			doNotAlign = true;
			System.out.println("Error while creating the new b file: " + ioe); 
		}
	}
	
	/**
	 * Report back to me when alignments successfully finish w/o errors.
	 * Send me emails if warnings/errors are reported from Tomato.
	 */
	public static void emailAlignmentUpdates() {

	}
	
	/**
	 * Use Tony's CreateAnalysis app for this.
	 */
	public static void populateGNomExAnalysisReport() {
		
	}
	
	/**
	 * This method will process each argument and assign new variables.
	 * @param args
	 */
	public void processArgs(String[] args) {
		Pattern pat = Pattern.compile("-[a-z]");
		String programArgs = Misc.stringArrayToString(args, ",");
		boolean verbose = false;
		if (verbose) System.out.println("\nArguments: " + programArgs + "\n");
		for (int i = 0; i < args.length; i++) {
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()) {
				char test = args[i].charAt(1);
				try {
					switch (test) {
					//case 'f': fastqFileName = new String(args[++i]); break;
					//case 'c': cmdFilePath = new String(args[++i]); break;
					//case 'g': genomeBuild = args[i+1]; i++; break;
					case 'p': paramsFile = new File(args[++i]); break;
					//case 'u': submitter = new String(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem--unknown option used!" + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -" + test + "\n");
				}
			}
		}
		//if (versionedGenome == null) Misc.printErrAndExit
			//("\nPlease provide a versioned genome (e.g. H_sapiens_Mar_2006).\n");
		//if (cmdFilePath == null) Misc.printErrAndExit("\nPlease provide full path to cmd.txt file.\n");
		//if (fastqFileName == null) Misc.printErrAndExit("\nPlease provide full path to fastq file(s) to align.\n");
		//if (paramsFile == null) Misc.printErrAndExit("\nPlease provide full path to alignment params file.\n");
		//if (userName == null) Misc.printErrAndExit("\nPlease provide name of user whose sequences you're aligning.\n");
		
	}
	
	public static void printDocs() {
		System.out.println("\n" +
				"**********************************************************************************\n" +
				"**                        MakeNovoalignments: Oct 2012                          **\n" +
				"**********************************************************************************\n" +
				"AutoNovoaligner does somewhat automatic alignments using novoalign\n" + 
				"**********************************************************************************\n");
	}
}