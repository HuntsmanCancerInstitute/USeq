package util.bio.wrappers;

import java.io.*;
import java.net.*;
import java.util.regex.*;
import java.util.zip.GZIPOutputStream;
import java.util.*;

import util.gen.*;

public class CHPCAligner {

	//fields
	//On CHPC
	private File genomeIndex;
	private File resultsDirectory;
	private File splitFastqDataDirectory1;
	private ArrayList<File> splitFastqFiles1 = new ArrayList<File>();
	private File splitFastqDataDirectory2;
	private ArrayList<File> splitFastqFiles2 = new ArrayList<File>();
	private File shellScriptDirectory;
	private File logDirectory;
	private File alignmentDirectory;
	private String name = null;
	private boolean launch = false;
	private String smtpHostName = "smtp.utah.edu";
	private boolean filterForChrLines = false;
	private boolean relaunchBadJobs = true;
	private String chpcAccount = null;
	private double numberOfJobs = 0;
	private boolean stripSAMSQHeaders = false;

	//aligner
	private File alignerApp = new File("/uufs/chpc.utah.edu/common/home/u0028003/BioApps/Novocraft/novocraft/novoalign");
	private String alignerName;

	//fastq file(s)
	private String fastqFile1;
	private File actualFastqFile1;
	private String fastqFile2;
	private File actualFastqFile2;
	private String fastqFileUserNameAndServer = "hci_u0028003@hci-bio3.hci.utah.edu";
	private String splitDataFileExtension = ".gz";

	//Final results archive, should contain the actual directory in which to place the compressed alignments
	private String archiveDirectory;
	private String archiveDirectoryUserNameAndServer = "hci_u0028003@hci-bio3.hci.utah.edu";

	//misc
	private String alignerParams = " ";
	private int numberReadsPerSlice = 500000;
	private String adminEmail = "david.nix@hci.utah.edu";
	private int hrsWallTime = 24;
	private int numberCPUs = 12;
	private ClusterJob job;
	private ClusterJob[] jobs;
	private ArrayList<String> log = new ArrayList<String>();

	//email
	private String[] recipientEmailAddresses;
	private Pattern analysisPattern = Pattern.compile(".+/(A\\d+)/*");

	public CHPCAligner (String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process and check args
		processArgs(args);

		//get and split data
		getAndSplitFastqData();

		//build and launch jobs
		buildAndLaunchJobs();

		if (launch)  {
			//monitor
			monitor();

			//look for bad jobs and clean up good ones
			if (relaunchBadJobs) {
				for (int i=0; i< 3; i++){
					int numBadJobs = cleanUpJobs().size();
					if (numBadJobs > 0) {
						relaunchBadJobs(); 
						monitor();
					}
					else break;
				}
			}

			//calculate run time stats
			int[] runTimes = fetchRunTimes();
			if (runTimes.length !=0){
				Arrays.sort(runTimes);
				double mean = Num.mean(runTimes);
				double median = Num.median(runTimes);
				int[] minMax = Num.findMinMaxIntValues(runTimes);
				print("\n\tJob run time stats (in minutes):");
				print("\t\t"+jobs.length+"\t# Jobs");
				print("\t\t"+minMax[0]+"\tMinimum");
				print("\t\t"+minMax[1]+"\tMaximum");
				print("\t\t"+(int)mean+"\tMean");
				print("\t\t"+(int)median+"\tMedian");
				print("\t\tRun times: "+ Num.intArrayToString(runTimes, " "));
			}

			//clean up jobs
			ArrayList<ClusterJob> badJobs = cleanUpJobs();
			boolean success = false;
			if (badJobs.size() == 0){
				//write log file to log dir
				IO.writeArrayList(log, new File (logDirectory, "runLog.txt"));
				//compress and send
				print("\tCompressing  and sending data to archive...");
				success = compressAndSend();
			}
			else {
				print("\n\tWARNING: Bad jobs found! Aborting! Check:\n\tName\tID\tLogFile\tAlignmentFile\tSplitFastqFile\tShellFile");
				for (int i=0; i< badJobs.size(); i++){
					ClusterJob job = badJobs.get(i);
					print("\t"+job);
				}
				//write log file to log dir
				IO.writeArrayList(log, new File (logDirectory, "runLog.txt"));
			}

			//email results!
			emailResults(success);

			//delete files? this often doesn't work due to locks by the system
			if (success) IO.deleteDirectory(resultsDirectory);

		}
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		print("\nDone! "+Math.round(diffTime)+" minutes\n");
	}

	private void relaunchBadJobs(){
		ArrayList<ClusterJob> badJobs = cleanUpJobs();
		int numBadJobs = badJobs.size();
		if (numBadJobs !=0){
			print("\nRelaunching bad jobs:");
			for (int i=0; i< badJobs.size(); i++){
				ClusterJob job = badJobs.get(i);
				//set new log file
				job.setLog(new File(logDirectory, job.getName()+ "_ReDo.log"));
				print("\t"+job);
				//launch
				String[] qsubOut = IO.executeCommandLineReturnAll(new String[]{"qsub",job.getShellScript().toString()});
				//set job id
				job.setJobID(qsubOut[0]);

			}
		}
	}

	public void emailResults(boolean success){
		String subject;
		String message;
		String[] toEmail;
		if (success) {
			print("\n\tRun completed, no issues.");
			subject = "Alignment "+name+" completed.";
			Matcher mat = analysisPattern.matcher(archiveDirectory);
			if (mat.matches()){
				message = "Your alignments may be downloaded from https://bioserver.hci.utah.edu/gnomex/gnomexFlex.jsp?launchWindow=AnalysisDetail&analysisNumber="+mat.group(1)+"\n\n";
			}
			else message = "Fetch your alignments from "+archiveDirectory+"\n";
			if (recipientEmailAddresses == null) toEmail = new String[]{adminEmail};
			else {
				toEmail = new String[recipientEmailAddresses.length+1];
				System.arraycopy(recipientEmailAddresses, 0, toEmail, 0, recipientEmailAddresses.length);
				toEmail[recipientEmailAddresses.length] = adminEmail;
			}
		}
		else {
			print("\n\tWARNING: Problems found with run.");
			subject = "Alignment of "+name+" failed!";
			message = "";
			toEmail = new String[]{adminEmail};
		}
		message = message +"Log:\n"+ Misc.stringArrayListToString(log, "\n");
		/*if(smtpPassword != null){
			Email e = new Email (smtpHostName, smtpUser, smtpPassword, "text/plain");
			boolean sent = e.postMail(toEmail, subject, message, adminEmail);
			print("Email sent? "+sent);
		}*/
		boolean sent = Email.postMailNoAuthentication(toEmail, subject, message, adminEmail, "text/plain", smtpHostName);
		print("\n\tEmail sent? "+sent);
	}

	public int[] fetchRunTimes(){
		ArrayList<Integer> timesAL = new ArrayList<Integer>();
		for (int i=0; i< jobs.length; i++){
			int rt = jobs[i].getRunTime();
			if (rt > 0) timesAL.add(new Integer(rt));
		}
		return Num.arrayListOfIntegerToInts(timesAL);
	}

	public void print(String s){
		System.out.println(s);
		log.add(s);
	}

	public void printError(String s){
		System.err.println(s);
		log.add(s);
	}

	public boolean extractTarFile(File tarFile, File dir, ArrayList<File> dataFiles){

		try {
			String cmd = "tar -C "+dir.getCanonicalPath()+" -xf "+ tarFile.getCanonicalPath();
			String[] results = IO.executeShellScript(cmd, dir);
			//shouldn't return anything
			if (results == null || results.length !=0) {
				print("Could not extract files from your tar archive?!");
				print("Command: "+ cmd);
				print("Results: "+Misc.stringArrayToString(results, " "));
				throw new Exception();
			}
			//fetch all files, move if needed, and add those gzipped to the dataFiles ArrayList
			ArrayList<File> allFiles = IO.fetchAllFilesRecursively(dir);
			ArrayList<File> toDelete = new ArrayList<File>();
			for (int i=0; i< allFiles.size(); i++){
				File f = allFiles.get(i);
				if (f.getName().endsWith(splitDataFileExtension)) {
					if (f.getParentFile().equals(dir) == false) {
						File m = new File (dir, f.getName());
						f.renameTo(m);
						dataFiles.add(m);
					}
					else dataFiles.add(f);
				}
				else toDelete.add(f);
			}
			//delete the non data and tar files
			for (int i=0; i< toDelete.size(); i++) {
				toDelete.get(i).delete();
			}
			//delete any enclosed directories from tar extraction
			File[] enclosedDirs = IO.extractOnlyDirectories(dir);
			for (int i=0; i< enclosedDirs.length; i++) enclosedDirs[i].delete();
			tarFile.delete();
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		return true;
	}

	public void getAndSplitFastqData(){
		print("\tFetching and splitting fastq data...");
		//first fastq file
		int numLinesInFirst = -1;
		int numFilesInFirst = 0;
		if (splitFastqDataDirectory1.exists() == false){
			//download file?
			BufferedReader fq = fetchFastqDataReader(fastqFile1, true);
			//make directory to hold split files
			splitFastqDataDirectory1.mkdir();
			//split file?
			if (actualFastqFile1.getName().endsWith(".tar")) {
				boolean extracted = extractTarFile(actualFastqFile1, splitFastqDataDirectory1, splitFastqFiles1);
				if (extracted == false) Misc.printErrAndExit("\t\tError: Failed to extract your tar fastq archive!");
				numFilesInFirst = IO.numberFilesExist(splitFastqDataDirectory1, splitDataFileExtension);
				if (numFilesInFirst ==0) Misc.printErrAndExit("\t\tError: No split (xxx"+splitDataFileExtension+") fastq files found in your tar archive?!" );
			}
			//just one job so move over?
			else if (numberOfJobs == 1){
				File indir1 = new File(splitFastqDataDirectory1, actualFastqFile1.getName());
				actualFastqFile1.renameTo(indir1);
				splitFastqFiles1.add(indir1);
				numFilesInFirst = 1;
			}
			else {
				//break file into set number of chunks?
				if (numberOfJobs > 1){
					double numLines = this.countLinesInFile(actualFastqFile1);
					numberReadsPerSlice = (int)Math.round(numLines/(numberOfJobs * 4.0));
					int r = (int)numLines % (int)(numberOfJobs * 4.0);
					if (r !=0) numberReadsPerSlice++;
					System.out.println("\t\tSplitting reads into "+numberReadsPerSlice + " read chunks to launch "+numberOfJobs+" jobs. Total read count = "+(int)(numLines/4));
				}
				numLinesInFirst = splitFastqDataGZipOutput(fq, true);
			}
			//delete original file?
			if (actualFastqFile1 != null && numberOfJobs !=1) actualFastqFile1.delete();
		}
		else {
			print("\tSplitFastqData1 directory exists. Delete it to download and split...");
			File[] files = IO.extractFiles(splitFastqDataDirectory1, splitDataFileExtension);
			for (int i=0; i< files.length; i++) splitFastqFiles1.add(files[i]);
			numFilesInFirst = files.length;
			if (numFilesInFirst ==0) Misc.printErrAndExit("\t\tError: No split (xxx"+splitDataFileExtension+") fastq files found in your tar archive?!" );

		}


		//second fastq file?
		if (fastqFile2 != null){
			if (splitFastqDataDirectory2.exists() == false){
				//download file?
				BufferedReader fq = fetchFastqDataReader(fastqFile2, false);
				//make directory to hold split files
				splitFastqDataDirectory2.mkdir();
				//split file?
				int numLinesInSecond = 0;
				if (actualFastqFile2.getName().endsWith(".tar")) {
					boolean extracted = extractTarFile(actualFastqFile2, splitFastqDataDirectory2, splitFastqFiles2);
					if (extracted == false) Misc.printErrAndExit("\t\tError: Failed to extract your tar fastq archive!");
				}
				//just one job so move over?
				else if (numberOfJobs == 1){
					File indir2 = new File(splitFastqDataDirectory2, actualFastqFile2.getName());
					actualFastqFile2.renameTo(indir2);
					splitFastqFiles2.add(indir2);
				}
				//nope split it
				else numLinesInSecond = splitFastqDataGZipOutput(fq, false);
				//check number of lines?
				if (numLinesInFirst !=-1 && numLinesInFirst != numLinesInSecond) Misc.printErrAndExit("\t\tError: The number of sequences in your paired fastq data files differs! ");
				if (actualFastqFile2 != null && numberOfJobs !=1) actualFastqFile2.delete();
			}
			else {
				print("\tSplitFastqData2 directory exists. Delete it to download and split...");
				File[] files = IO.extractFiles(splitFastqDataDirectory2, splitDataFileExtension);
				for (int i=0; i< files.length; i++) splitFastqFiles2.add(files[i]);
			}
			//check that number of files is the same
			int numFilesInSecond = IO.extractFiles(splitFastqDataDirectory2, splitDataFileExtension).length;
			if (numFilesInSecond ==0) Misc.printErrAndExit("\t\tError: No split (xxx"+splitDataFileExtension+") fastq files found in your tar archive?!" );
			if (numFilesInFirst != numFilesInSecond) Misc.printErrAndExit("\t\tProblem: The number of split files in your paired fastq data differs! "+ numFilesInFirst +" vs "+ numFilesInSecond );
		}
	}

	public boolean compressAndSend(){
		StringBuilder sb = new StringBuilder();
		sb.append("cd "+resultsDirectory.getParent());
		sb.append("\n");
		sb.append("tar -cf "+name+".tar "+name);
		sb.append("\n");
		sb.append("rsync -L -e ssh "+name+".tar "+archiveDirectoryUserNameAndServer+ ":"+ archiveDirectory);
		sb.append("\n");
		String[] results = IO.executeShellScript(sb.toString(), logDirectory);
		//shouldn't return anything
		if (results == null || results.length !=0) {
			print("\t\tFailed to compress and or transfer your alignments?!");
			print("\t\tCommand: "+ sb.toString());
			print("\t\tResults: "+Misc.stringArrayToString(results, " "));
			return false;
		}
		//delete tar file
		File tar = new File (resultsDirectory+".tar");
		tar.delete();
		return true;
	}

	public ArrayList<ClusterJob> cleanUpJobs(){
		ArrayList<ClusterJob> badJobs = new ArrayList<ClusterJob>();
		for (int i=0; i< jobs.length; i++){
			if (jobs[i].isGzippedResultFileOK()) {
				jobs[i].getData1().delete();
				if (fastqFile2 != null) jobs[i].getData2().delete();
			}
			else badJobs.add(jobs[i]);
		}
		if (badJobs.size() ==0) {
			IO.deleteDirectory(splitFastqDataDirectory1);
			IO.deleteDirectory(splitFastqDataDirectory2);
		}
		return badJobs;
	}

	public ArrayList<ClusterJob> fetchIncompleteJobs(){
		ArrayList<ClusterJob> j = new ArrayList<ClusterJob>();
		for (int i=0; i< jobs.length; i++){
			if (jobs[i].isLogFilePresent() == false) j.add(jobs[i]);
		}
		return j;
	}

	public boolean monitor(){
		print("\tMonitoring job progress...");
		double numJobs = jobs.length;
		long startTime = System.currentTimeMillis();
		double oldFractionComplete = 0;
		while (true){
			try {
				//check every 5 min
				Thread.sleep(1000 * 60 * 5);
				double numIncom = fetchIncompleteJobs().size();
				double fractionComplete = 1- numIncom/numJobs;
				if (oldFractionComplete != fractionComplete){
					oldFractionComplete = fractionComplete;
					String percentComplete = Num.formatPercentOneFraction(fractionComplete);
					double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
					print("\t\t"+percentComplete +" Complete at "+ Math.round(diffTime)+" min");
				}
				if (fractionComplete == 1) break;
			} catch (InterruptedException e) {
				e.printStackTrace();
				break;
			}
		}
		return true;
	}

	public void buildAndLaunchJobs(){
		print("\tBuilding jobs for cluster...");
		jobs = new ClusterJob[splitFastqFiles1.size()];

		//fetch files and sort
		File[] splitFiles1 = new File[splitFastqFiles1.size()];
		splitFastqFiles1.toArray(splitFiles1);
		Arrays.sort(splitFiles1);
		File[] splitFiles2 = null;
		if (fastqFile2 != null){
			splitFiles2 = new File[splitFastqFiles2.size()];
			splitFastqFiles2.toArray(splitFiles2);
			Arrays.sort(splitFiles2);
		}
		//for each split data file
		for (int i=0; i< splitFiles1.length; i++){
			String index = (i+1)+"";
			File scriptFile = new File (shellScriptDirectory, index+".sh");
			String temp;
			if (splitFiles1.length == 1) temp = name;
			else temp = (i+1)+"_"+name;
			//make job
			job = new ClusterJob();
			job.setName(temp);
			job.setResults(new File(alignmentDirectory, temp+"."+ alignerName));
			job.setData1(splitFiles1[i]);
			if (fastqFile2 != null) job.setData2(splitFiles2[i]);
			job.setShellScript(scriptFile);
			job.setChpcAligner(this);
			job.setFilterForChrLines(filterForChrLines);
			job.setLog(new File(logDirectory, job.getName()+ ".log"));
			job.writeShellScript();
			if (launch){
				//launch
				String[] qsubOut = IO.executeCommandLineReturnAll(new String[]{"qsub",job.getShellScript().toString()});
				//set job id
				job.setJobID(qsubOut[0]);
				//add
				jobs[i] = job;
			}
		}
		if (launch) print("\tLaunched "+jobs.length+" jobs...");
		else print("\tRerun with the -l option to launch the jobs on the cluster.");
	}
	
	

	public int splitFastqDataGZipOutput(BufferedReader in, boolean firstFile){
		int totalLines = 0;
		try {
			//make split files
			String line = null;
			int fileNumber = 1;
			int counter = 0;
			int numLinesPerFile = numberReadsPerSlice * 4;
			File splitDir;
			ArrayList<File> al;
			if (firstFile) {
				splitDir = this.splitFastqDataDirectory1;
				al = this.splitFastqFiles1;
			}
			else {
				splitDir = this.splitFastqDataDirectory2;
				al = this.splitFastqFiles2;
			}
			//make first writer
			File f = new File (splitDir, fileNumber +splitDataFileExtension);
			GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(f));
			fileNumber++;
			al.add(f);
			byte[] cr = "\n".getBytes();
			byte[] b = null;
			while ((line = in.readLine()) !=null) {
				//skip blank lines
				if (line.length()==0) continue;
				//less than counter
				if (counter < numLinesPerFile) {
					counter++;
				}
				else {
					counter =1;
					out.finish();
					out.close();
					f = new File (splitDir, fileNumber +splitDataFileExtension);
					fileNumber++;
					al.add(f);
					out = new GZIPOutputStream(new FileOutputStream(f));
				}
				totalLines++;
				b=line.getBytes();
				out.write(b);
				out.write(cr);
			}
			//close final PrintWriter
			out.finish();
			out.close();
			in.close();
			//check number of lines
			if ((totalLines % 4) != 0) Misc.printErrAndExit("Error: the number of text lines in your fastq file is not divisible by 4!");
			else print("\t\t"+(totalLines/4)+" reads found");

		} catch (Exception e){
			print("Problem splitting your fastq data.");
			e.printStackTrace();
		}
		return totalLines;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CHPCAligner(args);
	}		

	public BufferedReader fetchFastqDataReader(String fastqFile, boolean firstFile){
		//fetch buffered reader
		BufferedReader in = null;
		try {
			//is it a web address?
			if (fastqFile.startsWith("http") || fastqFile.startsWith("ftp")){
				URL u = null;
				u = new URL (fastqFile);
				in = IO.fetchBufferedReader(u);
			}	
			else {
				//nope it's a file or directory of split files so download it
				File actualFastqFile = new File (resultsDirectory, fastqFile.substring(fastqFile.lastIndexOf("/")+1));				
				String cmd = "rsync -L -e ssh  "+ fastqFileUserNameAndServer+":"+fastqFile +" "+ actualFastqFile.getCanonicalPath();
				//String cmd = "scp "+ fastqFileUserNameAndServer+":"+fastqFile +" "+ actualFastqFile.getCanonicalPath();
				String[] results = IO.executeShellScript(cmd, resultsDirectory);

				//shouldn't return anything
				if (results == null || results.length !=0) {
					//print("Could not fetch your fastq file?!");
					print("Command: "+ cmd);
					print("Results: "+Misc.stringArrayToString(results, " "));
					throw new Exception();
				}

				//make a reader
				in = IO.fetchBufferedReader(actualFastqFile);
				//set file
				if (firstFile) actualFastqFile1 = actualFastqFile;
				else actualFastqFile2 = actualFastqFile;
			}
		} catch (Exception e) {
			print("Failed to aquire a BufferedReader on your fastq data.");
			e.printStackTrace();
		} 
		return in;
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String s = Misc.stringArrayToString(args, " ");
		print("\nCHPCAligner Arguments: "+s+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'i': genomeIndex = new File(args[++i]); break;
					case 'r': resultsDirectory = new File(args[++i]); break;
					case 'j': alignerApp = new File(args[++i]); break;
					case 'f': fastqFile1 = args[++i]; break;
					case 's': fastqFile2 = args[++i]; break;
					case 'd': fastqFileUserNameAndServer = args[++i]; break;
					case 'g': archiveDirectoryUserNameAndServer = args[++i]; break;
					case 'a': archiveDirectory = args[++i]; break;
					case 'p': alignerParams = args[++i]; break;
					case 'e': adminEmail = args[++i]; break;
					case 'o': chpcAccount = args[++i]; break;
					case 'c': recipientEmailAddresses = args[++i].split(","); break;
					case 'n': numberReadsPerSlice = Integer.parseInt(args[++i]); break;
					case 'k': numberOfJobs = Integer.parseInt(args[++i]); break;
					case 'l': launch = true; break;
					case 'b': relaunchBadJobs = false; break;
					case 't': filterForChrLines = true; break;
					case 'q': stripSAMSQHeaders = true; break;
					case 'w': hrsWallTime = Integer.parseInt(args[++i]); break;
					case 'x': numberCPUs = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}		

		//check params
		if (genomeIndex == null || genomeIndex.canRead()== false) Misc.printErrAndExit("Error: Cannot find or read your novoalign index file -> "+genomeIndex);
		if (resultsDirectory == null) Misc.printErrAndExit("Error: Please enter a results directory on CHPC.");
		if (fastqFile1 == null) Misc.printErrAndExit("Error: Please enter the fastq file to align.");
		if (archiveDirectory == null) Misc.printErrAndExit("Error: Please enter a results archive directory to save the alignment files.");

		//make results directory on scratch disk
		if (resultsDirectory.exists()) print("\tScratch results directory exists, data will be overwritten...");
		else resultsDirectory.mkdir();
		name = resultsDirectory.getName();

		//instantiate directory object to hold split fastq files, don't actually make it though
		splitFastqDataDirectory1 = new File(resultsDirectory, "SplitFastqData1");
		if (fastqFile2 != null) splitFastqDataDirectory2 = new File(resultsDirectory, "SplitFastqData2");

		//check aligner
		if (alignerApp == null) Misc.printErrAndExit("Error: Please enter some novoalign arguments.");
		alignerName = alignerApp.getName();
		if (alignerName.contains("novoalign") == false) Misc.printErrAndExit("Error: Only novoalign is supported at this time.");
		if (alignerParams.contains("SAM")) alignerName = alignerName + ".sam";

		//test access to fastq file(s)
		print("\tChecking access to fastq data...");
		if (splitFastqDataDirectory1.exists() == false && cannotAccessFastQFile(fastqFile1)) Misc.printErrAndExit("Error: cannot access your fastq file -> "+fastqFile1 + " on "+fastqFileUserNameAndServer);
		if (fastqFile2 != null){
			if (splitFastqDataDirectory2.exists() == false && cannotAccessFastQFile(fastqFile2)) Misc.printErrAndExit("Error: cannot access your fastq file -> "+fastqFile2 + " on "+fastqFileUserNameAndServer);
		}

		//test connection to final results directory on hci-ma where the tar file will be placed
		if (cannotAccessArchiveDirectory()) Misc.printErrAndExit("Error: cannot access your final archive directory -> "+archiveDirectory +" on "+archiveDirectoryUserNameAndServer);

		//make directory to hold shell scripts
		shellScriptDirectory = new File (resultsDirectory, "ShellScripts");
		shellScriptDirectory.mkdir();

		//make directory to hold logs results
		logDirectory = new File (resultsDirectory, "Logs");
		logDirectory.mkdir();

		//make directory to hold alignment results
		alignmentDirectory = new File (resultsDirectory, "Alignments");
		alignmentDirectory.mkdir();

	}

	public boolean cannotAccessFastQFile(String fastQFile){
		//is it a web address?
		if (fastQFile.startsWith("http")){
			URL u = null;
			try {
				u = new URL (fastQFile);
				BufferedReader in = IO.fetchBufferedReader(u);
				//read 10 lines
				for (int i=0; i< 10; i++) in.readLine();
				in.close();
				return false;
			} catch (Exception e) {
				e.printStackTrace();
				return true;
			}
		}
		//nope it's a file
		String command = "ssh " + fastqFileUserNameAndServer+ " 'ls "+fastQFile+"'";
		String[] results = IO.executeShellScript(command, resultsDirectory);
		if (results== null || results[0].equals(fastQFile) == false) {
			print("Command: "+command);
			print("Results: "+Misc.stringArrayToString(results, " "));
			return true;
		}
		return false;
	}

	private boolean cannotAccessArchiveDirectory(){
		print("\tChecking access to final archive directory...");
		String command = "ssh "+ archiveDirectoryUserNameAndServer+ " 'test -d "+archiveDirectory+" && echo true'";
		String[] results = IO.executeShellScript(command, resultsDirectory);		
		if (results== null || results[0].equals("true") == false) {
			print("Command: "+command);
			print("Results: "+Misc.stringArrayToString(results, " "));
			return true;
		}
		return false;
	}
	
	public int countLinesInFile(File file){
		String command;
		if (file.getName().endsWith(".gz")) command = "gunzip -c "+ file+ " | wc -l ";
		else command = "wc -l "+ file;
		String[] results = IO.executeShellScript(command, resultsDirectory);
		if (results== null || results[0].equals("")) {
			print("Problem counting the lines in the file -> "+file);
			print("Command: "+command);
			print("Results: "+Misc.stringArrayToString(results, " "));
			return -1;
		}
		results = results[0].trim().split("\\s+");
		return Integer.parseInt(results[0]);
	}


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                CHPC Aligner: Sept 2011                           **\n" +
				"**************************************************************************************\n" +
				"Wrapper for running novoalign on the CHPC clusters. You will need to configure ssh\n" +
				"keys from CHPC to your other servers. See http://linuxproblem.org/art_9.html . Run\n" +
				"this app at the CHPC.\n\n" +

				"Required Options:\n"+
				"-i Genome index file on CHPC\n"+
				"-r Results directory on CHPC, this also defines the name of the final data archive\n" +
				"-f First fastq file on the raw data server\n"+
				"-s (Optional) Second paired end read fastq file on the raw data server\n"+
				"-a Archive directory on the analysis server for saving the final alignments\n"+

				"\nDefault Options:\n"+
				"-l Launch jobs, defaults to not launching jobs, inspect and test the shell scripts\n"+
				"     before committing.\n"+
				"-w Wall time in hours, defaults to 24. (Max on ember general is 24.)\n"+
				"-x Number CPUs, defaults to 12 for ember, set to 4 for sanddune.\n"+
				"-e Administrator email address, defaults to david.nix@hci.utah.edu\n"+
				"-c (Optional) Client email addresses, comma delimited, no spaces.\n"+
				"-b Don't relaunch bad jobs, defaults to making 3 attempts before aborting.\n"+
				"-o CHPC account to draw hours from (e.g. kaplan, kaplan-em, cairns, etc)\n"+
				"-d Raw data user name and server, defaults to hci_u0028003@hci-bio3.hci.utah.edu\n"+
				"-g Final alignment data user name and server, defaults to\n" +
				"     hci_u0028003@hci-bio3.hci.utah.edu\n"+
				"-j Aligner application, defaults to \n" +
				"     '/uufs/chpc.utah.edu/common/home/u0028003/BioApps/Novocraft/novocraft/novoalign'\n" +
				"-p Aligner cmd line options\n"+
				"-n Number of reads to process per job, defaults to 500000\n"+
				"-k Number of jobs to run, defaults to number of reads per job setting\n"+
				"-t Filter results for lines containing a 'chr' string, defaults to all.\n"+
				"-q Strip @SQ: lines from SAM alignment results, recommended for transcriptomes.\n"+


				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/CHPCAlign  -p '-F ILMFQ -t60 -rRandom' \n"+ 
				"    -i ~/Genomes/hg19Splices34bpAdaptersNovo.index\n" +
				"    -r /scratch/serial/u0028003/7317X1_100602 \n"+
				"    -f /mnt/hci-ma/MicroarrayData/2010/7317R/GAII/100602_7317X1_s_7_1_sequence.txt.gz \n"+
				"    -s /mnt/hci-ma/MicroarrayData/2010/7317R/GAII/100602_7317X1_s_7_2_sequence.txt.gz \n"+
				"    -a /mnt/hci-ma/AnalysisData/2010/A115  -w 6 -e nix@gmail.com -t -b \n\n" +

		"**************************************************************************************\n");

	}

	public File getGenomeIndex() {
		return genomeIndex;
	}

	public File getAlignerApp() {
		return alignerApp;
	}

	public String getArchiveDirectory() {
		return archiveDirectory;
	}

	public String getArchiveDirectoryUserNameAndServer() {
		return archiveDirectoryUserNameAndServer;
	}

	public String getEmail() {
		return adminEmail;
	}

	public int getHrsWallTime() {
		return hrsWallTime;
	}

	public void setFastQFileUserNameAndServer(String fastQFileUserNameAndServer) {
		this.fastqFileUserNameAndServer = fastQFileUserNameAndServer;
	}

	public String getAlignerName() {
		return alignerName;
	}

	public String getAlignerParams() {
		return alignerParams;
	}

	public File getLogDirectory() {
		return logDirectory;
	}

	public String getChpcAccount() {
		return chpcAccount;
	}

	public int getNumberCPUs() {
		return numberCPUs;
	}

	public boolean isStripSAMSQHeaders() {
		return stripSAMSQHeaders;
	}


}
