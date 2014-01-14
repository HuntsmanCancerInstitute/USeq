package edu.utah.ames.bioinfo;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;
import util.gen.Misc;


/**
 * This class contains methods to automatically align fastq sequences after they come out of the sequencing pipeline.  
 * All the necessary input data is contained in the fresh data report that comes out of the pipeline.
 * Samples are run individually. A new analysis report is created in GNomEx (linked to experiment numbers) for all 
 * samples of the same request number. Bisulphite files are first split into smaller chunks and aligned separately. 
 * Tomato will email the user only if there are exceptions and/or job failures (-ef option, i.e. #e email@nobody.com -ef). 
 * All input fastq files are first run through fastqc app. A log is created for each run, with data reported on how 
 * each sample is handled. Run parameters supplied through the autoalign_run_configs.txt file.
 * 
 * @author darren.ames@hci.utah.edu
 *
 */

public class Autoaligner {

	//fields
	private String analysisNumber;
	private String myEmail = "darren.ames@hci.utah.edu"; 
	private String autoalignReport = null;
	private String novoindexNames = null;
	private String parsedReports = null;
	private String tomatoJobDir = null;
	private File configFile = null;
	private String logDir = null;
	private String createAnalysisMain = null;
	static String EMAILREGEX = "^#e\\s+(.+@.+)"; //email address pattern
	static String LABNAMEREGEX = "\\w+\\s\\w+(?=\\WLab)"; //matches first and last name of lab using positive look
	

	private HashMap<String, String> genomeIndex;
	boolean isSmallRNA = false;
	private Logger logFile;

	//constructor
	public Autoaligner(String[] args) {	
		processArgs(args);
		this.parseConfig();

		//write info to log file
		Date date = new Date();
		logFile = new Logger(logDir,"autoalign","INFO");
		logFile.writeInfoMessage("Running Autoaligner: " + date);
	}

	public static void main(String[] args) throws Exception {
		//check for args
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		Autoaligner a = new Autoaligner(args);
		a.parseFreshDataReport();
		a.closeLogger();
	}

	public void closeLogger() { 
		logFile.closeLogger(); 
	}

	/**
	 * prints details for each sample 
	 * @param s
	 */
	public void printLogs(Sample s) {
		logFile.writeInfoMessage("\n");
		logFile.writeInfoMessage("Request number:\t" + s.getRequestNumber());
		logFile.writeInfoMessage("Sample ID:\t" + s.getSampleID());
		logFile.writeInfoMessage("Lab:\t" + s.getLab());
		logFile.writeInfoMessage("Requester:\t" + s.getRequester());
		logFile.writeInfoMessage("Requester email:\t" + s.getRequesterEmail());
		logFile.writeInfoMessage("Request year:\t" + s.getRequestYear());
		logFile.writeInfoMessage("Project name:\t" + s.getProjectName());
		logFile.writeInfoMessage("Organism:\t" + s.getOrganism());
		logFile.writeInfoMessage("Build code:\t" + s.getBuildCode());
		logFile.writeInfoMessage("Number of cycles:\t" + s.getNumberOfCycles());
		logFile.writeInfoMessage("Lane number:\t" + s.getLane());
		logFile.writeInfoMessage("Sequencing application code:\t" + s.getSequencingApplicationCode());

		if (s.toAlign() == true) {
			logFile.writeInfoMessage("Known sequencing application code:\t" + s.knownSeqAppCode());
			if (s.isBisSeq() == true) logFile.writeInfoMessage("Sample type:\tbisulfite"); 
			else if (s.isChIPSeq() == true) logFile.writeInfoMessage("Sample type:\tChIP-Seq");
			//RNA-Seq
			else if (s.isRNASeq() == true) { 
				//small RNA
				if (s.isSmallRNA() == true) {
					logFile.writeInfoMessage("Sample type:\tsmall RNA");
				}
				//other RNA
				else {
					logFile.writeInfoMessage("Sample type:\tRNA");
				}
			}
			//DNA
			else if (s.getSequencingApplication().toString().equals("_DNASEQ")) {
				logFile.writeInfoMessage("Sample type:\tDNA");
			}
			//unknown sample type
			else {
				logFile.writeInfoMessage("Sample type:\tunknown");
			}
			logFile.writeInfoMessage("Single or paired-end reads:\t" + s.getSingleOrPairedEnd());
			//matching novoindex
			if (s.hasNovoindex() == true) {
				logFile.writeInfoMessage("Novoindex:\t" + s.getNovoindex());
			}
			//no matching novoindex
			else {
				logFile.writeInfoMessage("Novoindex:\t" + s.hasNovoindex());
			}
			//paired-end reads
			if (s.isPairedEnd() == true) {
				logFile.writeInfoMessage("Fastq file 1:\t" + s.getFastqFile1());
				logFile.writeInfoMessage("Fastq file 2:\t" + s.getFastqFile2());
			}
			//single-end reads
			else {
				logFile.writeInfoMessage("Fastq file:\t" + s.getFastqFileName());
			}
			//adapters present
			if (s.adaptersIncluded() == true) {
				logFile.writeInfoMessage("Read adapter 1:\t" + s.getRead1Adapter());
				logFile.writeInfoMessage("Read adapter 2:\t" + s.getRead2Adapter());
			}
			//adapters not included
			else {
				logFile.writeInfoMessage("Read adapters:\tnone");
			}
			//ok to align
			if (s.toAlign() == true) {
				logFile.writeInfoMessage("Analysis number:\t" + s.getAnalysisNumber());
				logFile.writeInfoMessage("Aligned:\t" + s.toAlign());
			}
		}
		else {
			if (s.isBisSeq() == true) {
				logFile.writeInfoMessage("**Sample skipped...Bisulfite**");
			}
			else {
				logFile.writeInfoMessage("**Sample skipped...missing required info**");
			}
		}
	}

	/**
	 * Loads a list of novoindices into a hash. First column (build_sequencingApplication[_radius?]) 
	 * is the key, abbreviated index name used by Tomato in second column is the value.
	 * @param s
	 * @return
	 * @throws IOException
	 */
	public HashMap<String,String> loadNovoindexList() {
		//make new hashmap
		HashMap<String,String> indices = new HashMap<String,String>(10000);
		try{
			//read in the file of index codes and names
			BufferedReader br = new BufferedReader(new FileReader(novoindexNames));
			String line;
			String[] keyValue;
			//split file on tabs
			while ((line = br.readLine()) != null){
				keyValue = line.split("\t");
				//put the key/value pairs in the hashmap
				indices.put(keyValue[0].trim(), keyValue[1].trim());
			}
		} 
		catch(Exception e) {
			System.out.println("Problem loading the Hash");
			e.printStackTrace();
		}
		return indices;
	}

	/**
	 * Parses the fresh data report and populates all the fields of a Sample object.
	 * @throws Exception 
	 */
	public void parseFreshDataReport() throws Exception {
		//load novoindex hash
		genomeIndex = loadNovoindexList();
		//create buffered reader to read fresh data report
		BufferedReader br = new BufferedReader(new FileReader(autoalignReport));
		String line;
		//skip first line in file
		br.readLine();
		//make hashmap
		HashMap<String,String> hm = new HashMap<String,String>();

		//loop the read block until all lines in the file are read
		while ((line = br.readLine()) != null) {

			//split contents of tab-delimited file 
			String dataValue[] = line.split("\t");

			//check for correct number of columns in fresh data report
			if (dataValue.length < 18) {
				//continue;
			}
			else {
				Sample s = new Sample(dataValue);

				//catch B37 genome if specified
				//TODO catch B37 only with RNA-Seq, not exome stuff
				//TODO make sure we have B37 index in data dir
				if (s.getGenome().toString().equals("H_sapiens_Jun_2003")) {
					s.setGenome("hg19; H_sapiens_Feb_2009; GRCh37");
				}
				
				//if paired-end
				if (s.getSingleOrPairedEnd().toString().equals("Paired-end")) {

					s.setPairedEnd(true);

					//split on comma
					String[] pairedFiles = s.getFastqFilePath().split(",");

					//set paths
					s.setFastqFilePath1(pairedFiles[0]);
					String f1 = pairedFiles[0].toString();
					s.setFastqFile1(f1);
					s.setFastqFilePath2(pairedFiles[1]);
					String f2 = pairedFiles[1].toString();
					s.setFastqFile2(f2);
				}
				//single-end
				else {
					//set fastq file name
					String path = s.getFastqFilePath();
					String fastqFileName = new File(path).getName();
					s.setFastqFileName(fastqFileName);
				}
				//parse out "Lab" from lab name string
				Pattern p = Pattern.compile(".+(?=Lab)");
				Matcher m = p.matcher(s.getLab());
				if (m.find()) {
					s.setLab(m.group());
				}

				String genome = dataValue[s.genomeIndex];
				//skip samples without organism or index info
				if (s.getOrganism().toString().equals("Unknown") || s.getOrganism().toString().equals("Other") || 
						s.getOrganism().toString().equals(null) || s.getGenome().toString().equals("None") || 
						s.getGenome().toString().equals(null)) {

					//do not align
					s.setAlign(false);
					//continue;
				}
				else {
					s.setGenome(genome);
					//call downstream methods to match samples with the correct novoindex
					this.setCondSeqAppCode(s);
				}

				//check for matching novoindex for the sample before creating new analysis report 
				if (!(s.getNovoindex() == null)) {
					//get new analysis report number for samples that have a matching novoindex
					analysisNumber = hm.get(s.getRequestNumber());
					if (analysisNumber == null) {
						//populate with request numbers
						analysisNumber = this.getAnalysisNumber(s);
						s.setAnalysisNumber(analysisNumber);
						hm.put(s.getRequestNumber(), analysisNumber);
					}
					else if (hm.containsKey(s.getRequestNumber())) {
						s.setAnalysisNumber(hm.get(s.getRequestNumber()));
					}
					this.createCmdFile(s);
					this.sequencingAdapterState(s);
				}
				else {
					s.setAlign(false);
					//continue;
				}
				//print logs for each sample
				this.printLogs(s);
			}
		}
		//close the reader
		br.close();

		//move fresh data report to processedReports dir
		File ar = new File(autoalignReport);
		boolean success = ar.renameTo(new File(parsedReports, ar.getName()));
		if (!success) {
			//write to log file
			logFile.writeInfoMessage("Problem moving parsed fresh data report "
					+ autoalignReport.toString());
		}
	}

	/**
	 * Checks if Read1Adapter and Read2Adapter fields are populated or contain "None"
	 * @param s
	 * @return
	 */
	public boolean sequencingAdapterState(Sample s) {
		//check if adapter strings are "None"
		if (s.getRead1Adapter().toString().equalsIgnoreCase("None") && 
				s.getRead2Adapter().toString().equalsIgnoreCase("None")) {
			s.setAdaptersIncluded(false);
		}
		//adapters are there
		else {
			s.setAdaptersIncluded(true);
		}
		return s.adaptersIncluded();
	}

	/**
	 * Runs an external shell script, writing params to stdin and collecting output in stderr/stdout. 
	 * The method grabs a new analysis number in GNomEx for each unique request number that has 
	 * a matching novoindex in the tomato/data directory. New analysis reports and their corresponding
	 * experiment numbers are linked. Makes the new directory after fetching the new analysis number.
	 * @param s
	 * @throws IOException
	 */
	public String getAnalysisNumber(Sample s) throws IOException {
		String line;
		InputStream stdout = null;
		InputStream stderr = null;
		String analysisNum = null;
		String analysisPath = null;

		//feed input params into string array
		String cmd[] = {"./httpclient_create_analysis.sh", "-properties", "gnomex_httpclient.properties", 
				"-serverURL", "https://bioserver.hci.utah.edu", "-name", s.getProjectName(), 
				"-folderName", "novoalignments", "-lab", s.getLab(), "-organism", s.getOrganism(), 
				"-genomeBuild", s.getBuildCode(), "-analysisType", "Alignment", "-seqLane", s.getSequenceLaneNumber()};

		//launch the script and grab stdin/stdout and stderr
		Process process = Runtime.getRuntime().exec(cmd, null, new File(createAnalysisMain));

		stderr = process.getErrorStream();
		stdout = process.getInputStream();

		//clean up if any output in stdout
		BufferedReader brCleanup = new BufferedReader(new InputStreamReader(stdout));

		//fetch output analysis number
		while (((line = brCleanup.readLine()) != null)) {

			System.out.println("[stdout] " + line);
			if (line.indexOf("idAnalysis=") > 0) {
				//grab the new analysis number in output
				String matchANum = "([A][0-9]+)";
				Pattern p = Pattern.compile(matchANum);
				Matcher m = p.matcher(line);
				if (m.find()) {
					//assign analysis number to var
					analysisNum = m.group();
					s.setAnalysisNumber(analysisNum);
				}
				else {
					s.setAlign(false);
				}
			}

			//capture path for new analysis report
			String matchPath = "[/\\w]+[A][0-9]+";
			Pattern p = Pattern.compile(matchPath);
			Matcher m = p.matcher(line);
			if (m.find()) {
				//set the path in sample object
				analysisPath = m.group() + "/";
				s.setAnalysisNumberPath(analysisPath);

				//write to log file
				logFile.writeInfoMessage("Path for " + s.getAnalysisNumber() + ": "
						+ s.getAnalysisNumberPath());
			}
		}

		//set path
		String path = s.getAnalysisNumberPath();

		//check for new Analysis Report directory
		if (path.isEmpty() || s.getAnalysisNumber() == null) {
			System.out.println("Problem creating new Analysis Report " + s.getAnalysisNumber());

			s.setAlign(false);
		}
		//close reader
		brCleanup.close();

		return analysisNum;
	}

	/**
	 * Sets condensed sequencing application codes for matching sample to appropriate novoindex
	 * @param s
	 * @throws Exception 
	 */
	public void setCondSeqAppCode(Sample s) throws Exception {
		
		//set isPaired flag
		if (s.getSingleOrPairedEnd().toString().equals("Paired-end")) {
			s.setPairedEnd(true);
		}
		//set reverse strand flag for paired-end stranded (m)RNA sequencing
		if ((s.getSequencingApplicationCode().toString().equals("APP2") || 
				s.getSequencingApplicationCode().toString().equals("APP3")) && 
				s.isPairedEnd() == true) {
			s.setReverseStrand(true);
		}
		else {
			s.setReverseStrand(false);
		}

		//check to see if small RNA, and if so, set the boolean flag
		if (s.getSequencingApplicationCode().toString().equals("SMRNASEQ")) {
			isSmallRNA = true;
			s.setKnownSeqAppCode(true);
		}
		else {
			isSmallRNA = false;
		}

		//set condensed sequencing application codes 
		if (s.getSequencingApplicationCode().toString().equals("APP2") || 
				s.getSequencingApplicationCode().toString().equals("APP3") || 
				s.getSequencingApplicationCode().toString().equals("MRNASEQ") || 
				s.getSequencingApplicationCode().toString().equals("SMRNASEQ") || 
				s.getSequencingApplicationCode().toString().equals("APP27") ||
				s.getSequencingApplicationCode().toString().equals("APP9") || 
				s.getSequencingApplicationCode().toString().equals("DMRNASEQ")) {
			s.setSequencingApplication("_MRNASEQ_");
			s.setKnownSeqAppCode(true);
			s.setRNASeq(true);
		}
		else if (s.getSequencingApplicationCode().toString().equals("APP1") ||  
				s.getSequencingApplicationCode().toString().equals("APP4") || 
				s.getSequencingApplicationCode().toString().equals("APP8") || 
				s.getSequencingApplicationCode().toString().equals("APP11") || 
				s.getSequencingApplicationCode().toString().equals("APP12") ||
				s.getSequencingApplicationCode().toString().equals("APP13") || 
				s.getSequencingApplicationCode().toString().equals("APP20") || 
				s.getSequencingApplicationCode().toString().equals("APP21") || 
				s.getSequencingApplicationCode().toString().equals("APP22") || 
				s.getSequencingApplicationCode().toString().equals("APP23") || 
				s.getSequencingApplicationCode().toString().equals("APP25") || 
				s.getSequencingApplicationCode().toString().equals("APP26") || 
				s.getSequencingApplicationCode().toString().equals("APP29") || 
				s.getSequencingApplicationCode().toString().equals("APP30") ||
				s.getSequencingApplicationCode().toString().equals("APP33") ||
				s.getSequencingApplicationCode().toString().equals("APP42") ||
				s.getSequencingApplicationCode().toString().equals("CHIPSEQ") ||
				s.getSequencingApplicationCode().toString().equals("EXCAPSSC") ||
				s.getSequencingApplicationCode().toString().equals("TDNASEQ")) {
			s.setSequencingApplication("_DNASEQ");
			s.setKnownSeqAppCode(true);
		}
		else if (s.getSequencingApplicationCode().toString().equals("APP6") || 
				s.getSequencingApplicationCode().toString().equals("APP24")) {
			s.setSequencingApplication("_BISSEQ");
			s.setBisSeq(true);
			s.setKnownSeqAppCode(true);
			s.setAlign(false);
			emailExceptions(s);
		}
		//unknown sequencing application code
		else {
			s.setKnownSeqAppCode(false);
			s.setAlign(false);
			emailExceptions(s);
		}
		if (s.knownSeqAppCode() == true) {
			parseGenomeIndex(s);
		}
	}

	/**
	 * Sends me emails for the two most common errors
	 * @param s
	 * @throws IOException
	 */
	public void emailExceptions(Sample s) throws IOException {
		if (s.isBisSeq() == true) {
			//email myself so I know to check out the bisulfite sequences that are going to be split and aligned
			try {
				this.postMail(myEmail, "Bisulfite alignments", ("Hello kind Sir,\n\nBisulfite sequences have come out of pipeline. " +
						"Please check this fresh data report for details:\t" + this.autoalignReport.toString() +
						"\n\nLive long and prosper.\n\nDetails:\nRequest number:\t" + s.getRequestNumber() + 
						"\nSample ID:\t" + s.getSampleID() + "\nGenome:\t" + s.getGenome() + "\nFile:\t" +
						autoalignReport.toString()), "autoalign@nobody.com");
			} catch (MessagingException e) {
				e.printStackTrace();
			}
		}
		else if (s.knownSeqAppCode() == false) {
			try {
				this.postMail(myEmail, "unknown seqAppCode", ("Hello kind Sir,\n\nAlignments for experiment " + s.getRequestNumber() 
						+ " were skipped due to unknown sequencing application code: " + s.getSequencingApplicationCode() + "."
						+ " Please add this code and the appropriate alignment params to the list of knowns in the autoaligner." +
						"\n\n" + this.autoalignReport.toString()), 
						"autoalign@nobody.com");
			} catch (MessagingException e) {
				e.printStackTrace();
			}
		}
	}

	/**
	 * Parses abbreviated genome name (hg19, mm10...) from larger build string
	 * @param sample
	 * @return
	 * @throws IOException 
	 * @throws Exception 
	 */
	public void parseGenomeIndex(Sample s) throws IOException {
		//match build code (mm10, zv9, hg19, etc)
		String BUILDCODE = "([A-Za-z]+[0-9]+)(?=[;,:])"; 
		//match versioned genome (D_rerio_Jul_2010, etc)
		String VERSIONEDGENOME = "\\w_\\w+_\\w{3}_[0-9]{4}";
		
		if (s.getGenome() != null) {
			Pattern p = Pattern.compile(BUILDCODE);
			Matcher m = p.matcher(s.getGenome());
			if (m.find()) {
				s.setBuildCode(m.group().toLowerCase());
			}
			else {
				s.setAlign(false);
			}
			//parse versioned genome from genome string
			Pattern p2 = Pattern.compile(VERSIONEDGENOME);
			Matcher m2 = p2.matcher(s.getGenome());
			if (m2.find()) {
				s.setVersionedGenome(m2.group().toLowerCase());
			}
			//missing versioned genome?
			else {
				s.setAlign(false);
			}
		}
		buildShortIndexName(s);
	}

	/**
	 * Builds the shortened novoindex name from data in the fresh data report.
	 * Uses this shortened name as key for matching sample to the associated matching novoindex
	 * name that's recognized by Tomato
	 * @param s
	 * @return
	 * @throws IOException
	 */
	public void buildShortIndexName(Sample s) throws IOException {
		//check to see if no build code is specified in the report
		if (s.getBuildCode() == null) {
			s.setAlign(false);
		}
		else {
			//if RNA: build short index name for key in hashmap
			if (s.getSequencingApplication().toString().equals("_MRNASEQ_")) {
				//build the key
				String key = (s.getBuildCode() + s.getSequencingApplication() + s.getNumberOfCycles()); 
				//check hashmap for key (short index name), set index code
				if (genomeIndex.containsKey(key)) {
					s.setIndexCode(key);
					s.setNovoindex(genomeIndex.get(key));
					s.setHasNovoindex(true);
					s.setAlign(true);
				}
				//missing index
				else {
					s.setAlign(false);
					s.setHasNovoindex(false);
				}
			}
			//if DNA or bisulfite, don't include #cycles in short index name key in hashmap
			else if (s.getSequencingApplication().toString().equals("_DNASEQ")) {
				//build the key
				String key = (s.getBuildCode() + s.getSequencingApplication());
				if (genomeIndex.containsKey(key)) {
					s.setIndexCode(key);
					s.setNovoindex(genomeIndex.get(key));
					s.setAlign(true);
					s.setHasNovoindex(true);
				}
				//missing novoindex
				else {
					s.setAlign(false);
					s.setHasNovoindex(false);
				}
			}
			else if (s.getSequencingApplication().toString().equals("_BISSEQ")) {
				String key = (s.getBuildCode() + s.getSequencingApplication());
				if (genomeIndex.containsKey(key)) {
					s.setIndexCode(key);
					s.setNovoindex(genomeIndex.get(key));
					s.setHasNovoindex(true);
					//TODO remove this line below to align bisulfite seqs
					s.setAlign(false);
				}
				else {
					s.setAlign(false);
					s.setHasNovoindex(false);
				}
			}
		}
		createCmdFile(s);
	}

	/**
	 * Checks to see if a directory exists for the request number, and makes a new directory if it doesn't. 
	 * Subdirectories are made for each separate sample ID and populated with cmd.txt file. The cmd.txt file
	 * for bisulfite sequences first contains file splitting parameters, and upon completion, an intermediary
	 * file (bisAlignWait.txt) containing alignment params is moved to each new chunked bisulfite directory, 
	 * renamed cmd.txt and executed.
	 * @param s
	 * @return
	 * @throws IOException
	 */
	public void createCmdFile(Sample s) throws IOException {
		String jobDirPath = tomatoJobDir + s.getRequestNumber() + "/" + s.getSampleID();
		try {
			//first make sure an index exists before continuing
			if (!(s.getIndexCode() == null)) {
				//check to see if job dir exists for requestNum
				File dir = new File(jobDirPath);
				//make appropriate new directory and subdirectories
				dir.mkdirs();
				//change dir permissions so Tomato can execute 
				Runtime.getRuntime().exec("chmod 777 " + dir);
				//create cmd.txt file
				FileWriter fw = new FileWriter(jobDirPath + "/" + "cmd.txt");
				BufferedWriter out = new BufferedWriter(fw);

				//make bisAlignWait file for bisulfite sequences
				if (s.getSequencingApplication().toString().equals("_BISSEQ")) {
					//make bisAlignWait.txt file that will later be the @align cmd.txt file after file splitting
					FileWriter bisAlignWait = new FileWriter(jobDirPath + "/" + "bisAlignWait.txt");
					BufferedWriter bout = new BufferedWriter(bisAlignWait);
					bout.write(this.getCmdFileMessageBisulfite(s));
					bout.close();
				}
				//write the messages
				out.write(this.getCmdFileMsgGen(s));
				//close output stream
				out.close();
			}
		}
		catch (Exception e) {
			s.setAlign(false);
			System.exit(0);
		}
	}

	/**
	 * This method soft links the appropriate input fastq file in Repository 
	 * to its respective Tomato job directory.
	 * @param s
	 * @param dirPath
	 * @throws IOException
	 */
	public void softLinkFiles(Sample s, String dirPath) throws IOException {
		//job folder
		String jobDirPath = tomatoJobDir + s.getRequestNumber() + "/" + s.getSampleID();

		//paired-end reads?
		if (s.isPairedEnd() == true) {
			Process p1 = Runtime.getRuntime().exec(new String[] {"ln", "-s", s.getFastqFilePath1(), jobDirPath});
			Process p2 = Runtime.getRuntime().exec(new String[] {"ln", "-s", s.getFastqFilePath2(), jobDirPath});
		}
		//single-end reads
		else {
			Process p = Runtime.getRuntime().exec(new String[] {"ln", "-s", s.getFastqFilePath(), jobDirPath});
		}
	}

	/**
	 * Generates the general body of the cmd.txt file using sample-specific params and calls appropriate methods to do 
	 * the actual populating of application-specific params based on whether sequencing application is 
	 * bisulfite, miRNA/small RNA, exome, genomic DNA, or other (mRNA...)
	 * 
	 * @param sample
	 * @return
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public String getCmdFileMsgGen(Sample s) throws IOException, InterruptedException {
		//set string for first non-variable part of cmd.txt params
		String msg = "#e " + myEmail + "\n#a " + s.getAnalysisNumber() + "\n## Novoalignments of " 
				+ s.getRequestNumber() + " " + s.getProjectName() + "\n## Aligning to " 
				+ s.getBuildCode() + " for " + s.getRequester() + " in the " + s.getLab() + "Lab "
				+ "\n\nfastqc *.gz --noextract" 
				+ "\n\n@align -novoalign ";
		//set string for first non-variable part of cmd.txt params
		String msg2 = " -g " + s.getNovoindex() + " -i *.gz" + " -gzip\n\n";

		//command string for splitting bisulfite files and preparing dirs for alignment
		//splits into files with 2 million reads, which with our current hardware at CHPC, takes ~1 hr to align/file
		String bisulfiteMsg = "#e " + myEmail + " -ef" + "\n#a " + s.getAnalysisNumber() + "\n## Splitting bisulfite fastq files of " 
				+ s.getRequestNumber() + " " + s.getProjectName() + " for " + s.getRequester() + " in the "
				+ s.getLab() + "Lab " + "\n\nfastqc " + s.getSampleID() + "*.gz" + " --noextract" 
				+ "\n\nFileSplitter.jar" + " -f " + s.getSampleID() + "*_1.txt.gz" + " -n 8000000 " + "-g"
				+ "\n\nFileSplitter.jar" + " -f " + s.getSampleID() + "*_2.txt.gz" + " -n 8000000 " + "-g"
				+ "\n\nfor i in *_" + s.getSampleID() + "*_1.txt.gz; do o=${i%*_1.txt.gz}; n=${i%_" + s.getSampleID() + "*}; mkdir $n; " 
				+ "mv ${o}* $n; cp bisAlignWait.txt $n; mv $n/bisAlignWait.txt $n/cmd.txt; touch $n/b; done\n";

		//instantiate new StringBuffer object for holding cmd.txt file's body
		StringBuffer sb = new StringBuffer();

		//don't generate cmd.txt file if appropriate novoindex isn't available
		if (!(s.getIndexCode() == null)) {
			//is it be bisulfite? 
			if (s.getSequencingApplication().equalsIgnoreCase("_BISSEQ")) {
				//get the bisulfite alignment params
				sb.append(bisulfiteMsg);
				//sb.append(this.getCmdFileMessageBisulfite(s));
			}
			//all right then, how about small RNA?
			else if (s.getSequencingApplication().equalsIgnoreCase("_MRNASEQ_") && 
					(isSmallRNA == true)) {
				//build the cmd.txt contents
				sb.append(msg);
				sb.append(this.getCmdFileMessageSmallRNA(s));
				sb.append(msg2);
			}
			//hmmm, could it be regular mRNA or RNA?
			else if (s.getSequencingApplication().equalsIgnoreCase("_MRNASEQ_") && 
					(isSmallRNA == false)) {
				sb.append(msg);
				sb.append(this.getCmdFileMessageStdParams(s));
				sb.append(msg2);
			}
			//ok then, it's got to be genomic DNA, mononucleosome or ChIP-Seq
			else if (s.getSequencingApplication().equalsIgnoreCase("_DNASEQ")) {
				sb.append(msg);
				sb.append(this.getCmdFileMessageGenomic(s));
				sb.append(msg2);
			}
			//nope, it's something freaky that isn't currently supported in Autoaligner 
			else {
				s.setAlign(false);
			}
		}
		//soft link input files
		String dirPath = "/Repository/MicroarrayData/" + s.getRequestYear() + "/" + s.getRequestNumber() + "/Fastq/"; 
		this.softLinkFiles(s, dirPath);
		//start the job by creating the b file
		this.startTomatoJob(s);
		return sb.toString();
	}

	/**
	 * Exome, genomic DNA, mononucleosome, ChIP alignment params
	 * @param s
	 * @return
	 */
	public String getCmdFileMessageGenomic(Sample s) {
		//no adapter sequences available?
		if (s.adaptersIncluded() == false) {
			s.setParams("[-o SAM -r None -H -k]");
		}
		//adapter sequences available
		else {
			s.setParams("[-o SAM -r None -a " + s.getRead1Adapter() + " " + s.getRead2Adapter() + " -H -k]");
		}
		return s.getParams();
	}

	/**
	 * Bisulfite alignment params
	 * @param s
	 * @return
	 */
	public String getCmdFileMessageBisulfite(Sample s) {
		s.setParams("#e " + myEmail + "\n#a " + s.getAnalysisNumber() + "\n## Bisulfite novoalignments of " + 
				s.getRequestNumber() + " " + s.getProjectName() + "\n## Aligning to " + s.getBuildCode() +
				" for " + s.getRequester() + " in the " + s.getLab() + "Lab " + "\n\n@align -novoalign " +
				"[-o SAM -r Random -t 240 -h 120 -b 2] -i *_" + s.getSampleID() + "_*.txt.gz -g " 
				+ s.getNovoindex() + " -p bisulphite -gzip\n");
		return s.getParams();	
	}

	/**
	 * Small RNA alignment params
	 * @param s
	 * @return
	 */
	public String getCmdFileMessageSmallRNA(Sample s) {
		//includes 3' adapter stripping prior to alignment
		//***this is NOT the standard Illumina Gex Adapter 2 used as default by novoalign
		s.setParams("[-o SAM -r All 50 -m -a " + s.getRead1Adapter() + " -l 18 -h 60]");
		return s.getParams();
	}

	/**
	 * Standard mRNA alignment params. 
	 * @param s
	 * @return
	 */
	public String getCmdFileMessageStdParams(Sample s) {
		//paired-end data
		if (s.isPairedEnd() == true) {
			//adapters included
			if (s.adaptersIncluded() == true) {
				s.setParams("[-o SAM -r All 50 -a " + s.getRead1Adapter() + " " + s.getRead2Adapter() + "]");
			}
			//adapters NOT included
			else {
				s.setParams("[-o SAM -r All 50]");
			}
		}
		//single-end data
		else if (s.isPairedEnd() == false) {
			//adapter included
			if (s.adaptersIncluded() == true) {
				s.setParams("[-o SAM -r All 50 -a " + s.getRead1Adapter() + "]");
			}
			//adapter NOT included
			else {
				s.setParams("[-o SAM -r All 50]");
			}
		}
		return s.getParams(); 
	}

	/**
	 * This method creates a new empty "b" file if the file with the same name does 
	 * not already exist. Returns true if the file with the same name did not 
	 * exist and it was created successfully, false otherwise. This starts the 
	 * Tomato job. 
	 */
	public void startTomatoJob(Sample s) {
		if (s.toAlign() == true) {
			//create b file to start the jobs
			File f = new File(tomatoJobDir + s.getRequestNumber() + "/" + s.getSampleID() + "/" + "b");
			try {
				f.createNewFile();
			}
			catch (IOException ioe) {
				s.setAlign(false);
			}
		}
	}

	/**
	 * Sends an email when bisulfite sequences are being split
	 * @param recipients
	 * @param subject
	 * @param message
	 * @param from
	 * @throws MessagingException
	 */
	public void postMail(String recipients, String subject, String message, String from) throws MessagingException {
		//set the host smtp address
		Properties props = new Properties();
		props.put("mail.smtp.host", "hci-mail.hci.utah.edu");

		//create some properties and get the default Session
		javax.mail.Session session = javax.mail.Session.getDefaultInstance(props, null);

		//create message
		Message msg = new MimeMessage(session);

		//set the from and to address
		InternetAddress addressFrom = new InternetAddress(from);
		msg.setFrom(addressFrom);
		msg.setRecipients(Message.RecipientType.TO, InternetAddress.parse(recipients, false));

		//optional: can also set custom headers here in the email if wanted
		//msg.addHeader("MyHeaderName", "myHeaderValue");

		//setting the Subject and Content type
		msg.setSubject(subject);
		msg.setContent(message, "text/plain");
		Transport.send(msg);
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
					case 'f': autoalignReport = new String(args[++i]); break; 
					case 'c': configFile = new File(args[++i]); break;
					default: Misc.printErrAndExit("Problem--unknown option used: " + mat.group() 
							+ "\nUsage: java -jar Autoaligner.jar -f autoalign_2013-01-01.txt "
							+ "-c configFile.txt");
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -" + test + "\n");
				}
			}
		}
		if (configFile == null) {
			System.out.println("Please set config file path");
			System.exit(1);
		}
		if (!configFile.exists()) {
			System.out.println("config file is missing!!!");
			System.exit(1);
		}
	}

	/**
	 * Checks for necessary variables in the config file
	 * @param name
	 * @param hm
	 * @return
	 */
	private String checkExistance(String name, HashMap<String,String> hm) {
		String var = "";
		if (hm.containsKey(name)) {
			var = hm.get(name);
		} else {
			System.out.println("Config file doesn't contain variable " + name);
			System.exit(1);
		}
		return var;
	}
	
	/**
	 * Parses the config run file to populate necessary parameters
	 */
	private void parseConfig() {
		HashMap<String,String> hm = new HashMap<String,String>();

		try {
			BufferedReader br = new BufferedReader(new FileReader(this.configFile));
			String tmp = null;

			while ((tmp = br.readLine()) != null) {
				String[] items = tmp.split("=");
				hm.put(items[0], items[1]);
			}
			this.novoindexNames = this.checkExistance("novoindexNames", hm);
			this.parsedReports = this.checkExistance("parsedReports", hm);
			this.tomatoJobDir = this.checkExistance("tomatoJobDir", hm);
			this.logDir = this.checkExistance("logDir", hm);
			this.createAnalysisMain = this.checkExistance("createAnalysisMain", hm);

		} catch (FileNotFoundException e) {
			System.out.println("Please set config file path");
			System.exit(1);
		} catch (IOException e) {
			System.out.println("Error reading file: " + this.configFile.getAbsolutePath());
			System.exit(1);
		}
	}

	public static void printDocs() {
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Autoaligner: Jan 2014                               **\n" +
				"**************************************************************************************\n" + 
				"This class contains methods to automatically align fastq sequences using novoalign\n" +
				"via Tomato after they come out of the HiSeq Pipeline. All the necessary input data\n" +
				"is contained in the autoalign report that comes out of the HiSeq Pipeline. Files\n" +
				"are aligned individually. A new analysis report is created in GNomEx (linked to\n" +
				"experiment numbers) for all samples of the same request number. The requester is\n" +
				"notified via email when jobs start and finish (except bisulfite at the moment).\n" +
				"Bisulphite files are split into smaller chunks and aligned individually. All input\n" +
				"fastq files are first run through fastqc app. New novoindices must be added to the\n" +
				"list of acceptable indices in novoindexNamesTable.txt, along with its associated\n" +
				"short index name. Each run logs what happens to each sample in the report.\n" + 
				"\nParameters:\n\n" +
				"-f filename for autoalign report to process\n" +
				"-c config file with run parameters\n" +
				"\nUsage:\n" +
				"java -jar pathTo/Autoaligner.jar -f autoalign_2014-01-01.txt -c configFile.txt\n" +
				"**************************************************************************************\n");
	}
}