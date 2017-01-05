package edu.utah.pysano;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.LinkedHashMap;

import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;
import javax.mail.Session;

import edu.utah.pysano.commands.Command;
import edu.utah.pysano.commands.CommandAlign;
import edu.utah.pysano.commands.CommandHaplotype;
import edu.utah.pysano.commands.CommandMerge;
import edu.utah.pysano.commands.CommandMetrics;
import edu.utah.pysano.commands.CommandPostProcess;
import edu.utah.pysano.commands.CommandVariantMerge;
import edu.utah.pysano.daemon.CThread;
import edu.utah.pysano.daemon.CThreadDaemon;
import edu.utah.pysano.utils.Constants;
import edu.utah.pysano.utils.Logger;
import util.gen.IO;
import util.gen.Misc;



public class TomatoFarmer {
	private Logger logFile;	
	private String username;
	
	private File projectDirectory = null;
	private String email = null;
	private boolean fullGenome = false;
	private boolean suppressEmail = true;
	private boolean use1KGenomes = false;
	private boolean validateFastq = false;
	private String studyName = "STUDY";
	private String cluster = null;
	private Integer heartbeat = 120;
	private Integer wallTime = null;
	private File targetFile = null;
	private int jobNumber = 5;
	
	private HashMap<String,String> properties = new HashMap<String,String>();
	private LinkedHashMap<String,Boolean> analysisSteps = new LinkedHashMap<String,Boolean>();
	
	
	private boolean tomatoFarmerFailed = true;
	private boolean useNovo = false;
	private boolean useVqsr = false;
	private boolean metricsFinished = false;
	
	
	private CThreadDaemon daemon = null;

	//Email validation
	public static final Pattern VALID_EMAIL_ADDRESS_REGEX = 
		    Pattern.compile("^[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,6}$", Pattern.CASE_INSENSITIVE);

	public static boolean validateEmail(String emailStr) {
	        Matcher matcher = VALID_EMAIL_ADDRESS_REGEX .matcher(emailStr);
	        return matcher.find();
	}

	
	public TomatoFarmer(String args[]) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		processArgs(args);
		
		
		
		for (String step: analysisSteps.keySet()) {
			//Check to see if we are going to run the command
			if (!analysisSteps.get(step)) {
				continue;
			}
			
			//Create command depending on step
			ArrayList<ArrayList<Command>> commandList = null;
			int maxFailures = 5;
			if (step.equals("align")) {
				commandList = CommandAlign.checkForInputs(projectDirectory, properties, logFile, studyName, useNovo, validateFastq, targetFile, analysisSteps, 
						cluster, email, wallTime, suppressEmail);
			} else if (step.equals("merge")) {
				commandList = CommandMerge.checkForInputs(projectDirectory, properties, logFile, studyName, targetFile, analysisSteps, 
						cluster, email, wallTime, suppressEmail);
			} else if (step.equals("postprocess")) {
				commandList = CommandPostProcess.checkForInputs(projectDirectory, properties, logFile, analysisSteps, 
						cluster, email, wallTime, suppressEmail);
			} else if (step.equals("metrics")) {
				commandList = CommandMetrics.checkForInputs(projectDirectory, properties, logFile, studyName, targetFile, analysisSteps, 
						cluster, email, wallTime, suppressEmail);
			} else if (step.equals("haplotype")) {
				commandList = CommandHaplotype.checkForInputs(projectDirectory, properties, logFile, targetFile, analysisSteps, 
						cluster, email, wallTime, suppressEmail);
			} else if (step.equals("genotype")) {
				commandList = CommandVariantMerge.checkForInputs(projectDirectory, properties, logFile, targetFile, studyName, useVqsr, use1KGenomes, true, analysisSteps, 
						cluster, email, wallTime, suppressEmail);
			}
			
			
			//Skip step if empty, probably shouldn't happen
			if (commandList == null) {
				logFile.writeWarningMessage("[TF] Nothing to process for step: " + step + ", moving to next phase.");
				continue;
			}
			
			//Process jobs
			ArrayList<Command> toRun = commandList.get(0);
			ArrayList<Command> toCleanup = commandList.get(1);
			if (toRun.size() > 0) {
				logFile.writeInfoMessage("[TF] Found " + toRun.size() + " jobs to run for step " + step + ".");
				
				//Start deamon
				daemon = new CThreadDaemon(logFile,step,toRun.size(),jobNumber);
				this.daemon.start();
				int threadCount = 1;
				for (Command c: toRun) {
					CThread thread = c.prepareJobDirectory(threadCount++, heartbeat, maxFailures);
					daemon.addJob(thread);	
				}
				
				//Wait for command to finish
				try {
					this.daemon.join();
					Thread.sleep(5000);
					if (this.daemon.getFailed()) {
						System.exit(1);
					}
				} catch (InterruptedException ie) {
					logFile.writeErrorMessage("[TF] Daemon interrupted",true);
					System.exit(1);
				}
			}
			if (toCleanup.size() > 0) {
				logFile.writeInfoMessage("[TF] Found " + toCleanup.size() + " jobs to cleanup for step " + step + ".");
				
				for (Command c: toCleanup) {
					c.postProcess();
				}
			}
			
			//Non-pysano command bits
			if (step.equals("metrics") || (analysisSteps.containsKey("metrics") && isCurrentStepAfterStep(step,"metrics") && !metricsFinished)) {
				CommandMetrics.mergeMetrics(projectDirectory, properties, studyName, logFile);
				metricsFinished = true;
			}
		}
		
		try {
			CThread.postMail(this.email, "TF finished successfully", "You can find your results here: " + projectDirectory.getAbsolutePath());
		} catch (MessagingException e) {
			logFile.writeErrorMessage("[TF] Failed to send ending email", true);
			e.printStackTrace();
		} catch (NoClassDefFoundError ncdfe) {
 			logFile.writeErrorMessage("[TF] Failed to send ending email", true);
 		}
		
		this.tomatoFarmerFailed = false;
		logFile.writeInfoMessage("[TF] TomatoFarmer finished successfully!");
	}

	private boolean isCurrentStepAfterStep(String currentStep, String testStep) {
		int pos1 = 0;
		int pos2 = 0;
		int counter = 0;
		for (String step: analysisSteps.keySet()) {
			if (step.equals(currentStep)) {
				pos1 = counter;
			} else if (step.equals(testStep)) {
				pos2 = counter;
			}
			counter++;
		}
		if (pos1 > pos2 ) {
			return true;
		} else {
			return false;
		}
	}
	
	private void validatePropertiesFile(File propertiesFile) {
		ArrayList<String> properties = new ArrayList<String>();
		properties.add("DATA_PATH");
		properties.add("TEMPLATE_PATH");
		properties.add("BACKGROUND_PATH");
		properties.add("BWA_PATH");
		properties.add("PICARD_PATH");
		properties.add("GATK_PATH");
		properties.add("SAM_PATH");
		properties.add("R_PATH");
		properties.add("USEQ_PATH");
		properties.add("NOVOALIGN_PATH");
		properties.add("TABIX_PATH");
		properties.add("SAM_PATH_LOCAL");
		properties.add("VCF_PATH_LOCAL");
		properties.add("TABIX_PATH_LOCAL");
		properties.add("USEQ_PATH_LOCAL");
		properties.add("TARGET_DEFAULT");
		properties.add("SAMBAMBA_PATH");
		properties.add("SAMBLAST_PATH");
		properties.add("SCRATCH_PATH");
		properties.add("JAVA_PATH");
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(propertiesFile));
			String line = null;
			
			Pattern p = Pattern.compile("(.+?)=(.+)");
			
			while ((line = br.readLine()) != null) {
				Matcher m = p.matcher(line);
				if (m.matches()) {
					this.properties.put(m.group(1), m.group(2));
				}
			}
			
			br.close();
				
		} catch (IOException ioex) {
			System.out.println("Error reading properties file: " + ioex.getMessage());
			System.exit(1);
		}
		
		boolean passed = true;
		for (String p: properties) {
			if (!this.properties.containsKey(p)) {
				System.out.println("Properties file does not contain required key: " + p);
				passed = false;
			}
		}
		
		if (!passed) {
			System.exit(1);
		}
	}
	
	
	
	
	private void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		
		String targetType = null;
		File directory = null;
		File propertiesFile = null;
		String email = null;
		
		String logLevel = "INFO";
		String analysisType = null;
		String cluster = null;
		
	
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': directory = new File(args[++i]); break;
					case 'e': email = args[++i]; break;
					case 'w': wallTime = Integer.parseInt(args[++i]); break;
					case 'y': analysisType = args[++i]; break;
					case 't': targetType = args[++i]; break;
					case 'l': logLevel = args[++i]; break;
					case 'a': heartbeat = Integer.parseInt(args[++i]); break;
					case 'j': jobNumber = Integer.parseInt(args[++i]); break;
					case 's': studyName = args[++i]; break;
					case 'x': suppressEmail = false; break; 
					case 'n': useNovo = true; break;
					case 'r': useVqsr = true; break;
					case 'p': propertiesFile = new File(args[++i]); break; 
					case 'h': printDocs(); System.exit(0);
					case 'v': validateFastq = true; break;
					case 'b': use1KGenomes = true; break;
					case 'g': fullGenome = true; break;
					case 'c': cluster = args[++i]; break;
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//read in properties file
		if (propertiesFile == null) {
			System.out.println("The properties ( -p ) file was not specified, exiting");
			System.exit(1);
		}
		if (!propertiesFile.exists()) {
			System.out.println("Specified properties file does not exist, exiting");
			System.exit(1);
		}
		validatePropertiesFile(propertiesFile);
		
		//***********************************************
		// Set up logging file
		//***********************************************
		
		//Get username
		this.username = System.getProperty("user.name");
				
		//Initialize log files
		logFile = new Logger(new File(this.properties.get("TEMPLATE_PATH")),projectDirectory,this.username,logLevel);
		
		//Write Username to file
		logFile.writeSystemMessage("[TF] Username: " + this.username);
		
		//************************************************
		// Parse command line arguments
		//************************************************
	
		//Write out system arguments
		logFile.writeInfoMessage("[TF] " + IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " "));
		
		if (directory == null) {
			logFile.writeErrorMessage("[TF] You must set directory name",false);
			System.exit(1);
		}
		
	    projectDirectory = directory;
		
		if (email == null) {
			logFile.writeErrorMessage("[TF] You must set email",false);
			System.exit(1);
		}
		
		this.email = email;
		
		if (!directory.canWrite() || !directory.canRead() || !directory.canExecute()) {
			logFile.writeErrorMessage("[TF] Your directory is not accessable to TF, please edit permissions",false);
			System.exit(1);
		}
		
		if (!logLevel.equals("INFO") && !logLevel.equals("WARNING") && !logLevel.equals("ERROR")) {
			logFile.writeWarningMessage("[TF] Logging level must be INFO, WARNING or ERROR.  Reverting to INFO");
		}
		
		if (jobNumber > 20 || jobNumber < 1) {
			logFile.writeErrorMessage("[TF] Job number must be between 1 and 20",false);
		}
		
		logFile.writeInfoMessage("[TF] Using " + jobNumber + " threads for this project");
		
		if (!validateEmail(email)) {
			logFile.writeErrorMessage("[TF] Email does not pass our format check, please check if it was entered properly.  Offending email: " + email,false);
			System.exit(1);
		}
	
		//Check to make sure the run type is valid
		if (analysisType == null) {
			logFile.writeErrorMessage("[TF] You must set runtype",false);
			System.exit(1);
		}
		
		//Make sure analysisType is a valid menu option
		if (!Constants.validMenu.contains(analysisType)) {
			logFile.writeErrorMessage("[TF] " + analysisType + " is not a valid analysis type. Please see help menu for options.", false);
			System.exit(1);
		}
		
		
		//Make sure the target capture is valid
		
		
		if (fullGenome) {
			targetFile = new File(this.properties.get("TARGET_GENOME"));
		} else {
			if (targetType != null) {
				File targetDirectory = new File(new File(this.properties.get("TEMPLATE_PATH")),"captureRegions");
				if (Constants.validTargets.contains(targetType)) {
					targetFile = new File(targetDirectory,targetType + ".bed");
					if (!targetFile.exists()) {
						logFile.writeErrorMessage("[TF] The target capture bed file is missing: " + targetFile.getAbsolutePath(), true);
						System.exit(1);
					}
				} else {
					logFile.writeInfoMessage("[TF] Non-standard target capture, looking for custom target file");
					targetFile = new File(targetType);
					if (!targetFile.exists()) {
						logFile.writeErrorMessage("[TF] Specified target region file does not exist, exiting", false);
						System.exit(1);
					} else {
						targetFile = this.checkTarget(targetFile);
					}
				}
			} else {
				targetFile = new File(this.properties.get("TARGET_DEFAULT"));
			}
		}
		
		//Check cluster
		if (cluster != null) {
			if (!Constants.availClusters.contains(cluster)) {
				System.out.println("[TF] The specified cluster is not available: " + cluster);
				System.exit(1);
			} else {
				this.cluster = cluster;
			}
		}
		
		
		//Build up the command list
		if (analysisType.equals("ugp_full")) {
			analysisSteps.put("align", true);
			analysisSteps.put("merge", true);
			analysisSteps.put("postprocess", true);
			analysisSteps.put("metrics", true);
			analysisSteps.put("haplotype", true);
			analysisSteps.put("genotype", true);
		} else if (analysisType.equals("ugp_align")) {
			analysisSteps.put("align", true);
			analysisSteps.put("merge", true);
			analysisSteps.put("postprocess", true);
		} else if (analysisType.equals("metrics")) {
			analysisSteps.put("metrics", true);
		} else if (analysisType.equals("ugp_variant")) {
			analysisSteps.put("haplotype", true);
			analysisSteps.put("genotype", true);
		} else if (analysisType.equals("ugp_align_metrics")) {
			analysisSteps.put("align", true);
			analysisSteps.put("merge", true);
			analysisSteps.put("postprocess", true);
			analysisSteps.put("metrics", true);
		} else if (analysisType.equals("ugp_haplotype_only")) {
			analysisSteps.put("haplotype",true);
		} else {
			logFile.writeErrorMessage("Analysis " + analysisType + " is not a valid analysis type.", false);
			System.exit(1);
		}
		
		//************************************
		// Attach shutdown
		//***********************************
		attachShutDownHook();
		
	}	
	
	private File checkTarget(File origTarget) {
		
		String name = origTarget.getName();
		File directory = origTarget.getParentFile();
		
		File usedFile = origTarget;
		
		//Unzip target if zipped
		if (origTarget.getName().endsWith(".gz")) {
			this.logFile.writeInfoMessage("[TF]  Found gzipped target file, decompressing");
			byte[] buffer = new byte[1024];
			
			try {
				GZIPInputStream gzis = new GZIPInputStream(new FileInputStream(origTarget));
				File unzipped = new File(directory,name.substring(0, name.length()-3));
			    FileOutputStream fos = new FileOutputStream(unzipped);
			 
		        int len;
		        while ((len = gzis.read(buffer)) > 0) {
		        	fos.write(buffer, 0, len);
		        }
		 
		        gzis.close();
		    	fos.close();
		    	
		    	usedFile = unzipped;
			} catch (IOException ioex) {
				this.logFile.writeErrorMessage("[TF]  Error decompressing target file, exiting", true);
				System.exit(1);
			}
		} else if (origTarget.getName().endsWith(".zip")) {
			this.logFile.writeInfoMessage("[TF]  Found zipped target file, decompressing");
			byte[] buffer = new byte[1024];
			
			try {
				ZipInputStream zis = new ZipInputStream(new FileInputStream(origTarget));
				ZipEntry ze = zis.getNextEntry();
				
				if (ze == null) {
					this.logFile.writeErrorMessage("[TF]  Zip file appears corruped, exiting", true);
					System.exit(1);
				}
				
				String entryName = ze.getName();
				File unzipped = new File(directory,entryName);
				FileOutputStream fos = new FileOutputStream(unzipped);
				
				int len;
				while ((len = zis.read(buffer)) > 0) {
	                fos.write(buffer, 0, len);
	            }
			
				fos.close();
				zis.close();
				
				usedFile = unzipped;
				
			} catch (IOException ioex) {
				this.logFile.writeErrorMessage("[TF]  Error decompressing target file, exiting", true);
				System.exit(1);
			}	
		}
		
		
		//Go through the file, strip "chr"
//		File cleanedFile = new File(directory,usedFile.getName() + "_Cleaned.bed");
		
		
//		try {
//			BufferedReader br = new BufferedReader(new FileReader(usedFile));
//			BufferedWriter bw = new BufferedWriter(new FileWriter(cleanedFile));
//			
//			String temp = null;
//			boolean editMade = false;
//			while ( (temp = br.readLine() ) != null ) {
//				if (temp.startsWith("chr")) {
//					editMade = true;
//					String fixed = temp.substring(3);
//					bw.write(fixed + "\n");
//				} else {
//					bw.write(temp + "\n");
//				}
//			}
//			
//			br.close();
//			bw.close();
//			
//			if (editMade) {
//				this.logFile.writeInfoMessage("[TF] Stripped 'chr' from target chromosomes");
//				usedFile = cleanedFile;
//			} else {
//				cleanedFile.delete();
//			}
//			
//		} catch (IOException ioex) {
//			this.logFile.writeErrorMessage("[TF] Error writing cleaned target file", true);
//			System.exit(1);
//		}
		
		
		return usedFile;
	}
	
	private static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                TomatoFarmer: June 2015                           **\n" +
				"**************************************************************************************\n" +
				"TomatoFarmer controls an exome analysis from start to finish.  It creates alignment \n" +
				"jobs for each of the samples in your directory, waits for all jobs to finish and then \n" +
				"launches metrics and variant calling jobs.  Jobs will be resubmitted up to a set \n" +
				"number of times to combat spurious CHPC erors.  Job directories are left behind so \n" +
				"you can save your log files.\n\n" +
				
				"\nRequired Arguments:\n\n"+
				"-d Job directory.  This directory must be a subdirectory of /tomato/version/job. Can \n" +
				"       can be several directory levels below /tomato/version/job. \n" +
				"       Example: '-d /tomato/dev/job/krustofsky/demo'. \n" +
				"-e Email address.  TomatoFarmer emails you once the job completes/fails. You can also\n" +
				"       opt to get all tomato emails as individual jobs start/end (see option -x).\n" +
				"       Example: '-e hershel.krustofsky@hci.utah.edu'.\n" +
				"-y Analysis pipeline. The analysis pipeline or step to run.  Current options are: \n" +
				"          Full pipline:\n" +
				"          1) ugp_full - UGP v1.3.0 pipeline, includes: Alignment, metrics and variant\n" +
				"             calls. Defaults to bwa and raw variant filtering\n" +
				"          Indivdual Steps\n"+
				"          2) ugp_align - Alignment/recalibration only. Defaults to bwa.\n" +
				"          6) metrics - Sample QC metrics only. Requires sample level  *final.bam \n " +
				"             files from aligment steps.\n" +
				"          7) ugp_variant - Variant detection and filtering. Requires 1 or more \n" + 
				"             *final.bam files from alignment step. Defaults to raw variant filtering.\n" +
				"      Example: '-y ugp_full'.\n" +  
				"-p Properties file.  This file contains a list of cluster-specific paths and options \n" +
				"      this file doesn't need to be changed by the user. Example: '-p properties.txt' \n" +
				"\nOptional Arguments:\n\n"+
				"-n Use novoalign. This option will change the aligner from bwa to novoalign.\n" +
				"-r Use variant recalibration. This option will change variant filtration from raw \n" +
				"   to vqsr.  This should only be used if there are enough samples in the study or \n" +
				"   1K genome background files are used.\n" +
				"-g Full genome.  Use this option if you want to detect variants genome-wide\n" +
				"-t Target regions. Setting this argument will restrict coverage metrics and variant \n" +
				"      detection to targeted regions.  This speeds up the variation detection process\n" +
				"      and reduces noise. Options are:\n" +
				"          1) AgilentAllExonV4\n" +
				"          2) AgilentAllExonV5\n" +
				"          3) AgilentAllExonV5UTR\n" +
				"          4) AgilentAllExon50MB\n" +
				"          5) NimbleGenEZCapV2\n" +
				"          6) NimbleGenEZCapV3\n" +
				"          7) TruSeq\n" +
				"          8) path to custom targed bed file.\n" +
				"      If nothing is specifed for this argument, gatk exome boundaries are used. \n" + 
				"      Example: '-t truseq'.\n" +
				"-b 1K Genome samples.  Use this option if you want to spike in 200 1K genome samples \n" +
				"      as the background sample set.  This should improve VQSR variant calling and \n" +
				"      VAAST, but it will take a lot more time to process. BETA, only works for core \n" +
				"      users!\n" +
				"-s Study name.  Set this if you want your VCF files to have a prefix other than \n" +
				"      'STUDY'. Example: '-s DEMO'.\n" +
				"-v Validate fastq files.  TomatoFarmer will validate your fastq files before running \n" +
				"      This is required if any of your samples are ASCII-64\n" +
				"-c Cluster. Specify cluster, defaults to all available clusters\n" +
			
				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/TomatoFarmer -d /tomato/version/job/demo/\n" +
				"      -e herschel.krustofsky@hci.utah.edu -y ugp_bwa -r -b -s DEMO \n" +
				"      -t AgilentAllExon50MB\n\n" +

				"**************************************************************************************\n");

	}

	
	private void attachShutDownHook(){
		  Runtime.getRuntime().addShutdownHook(new Thread() {
		   @Override
		   public void run() {
			 if (daemon != null && daemon.isAlive()) {
		    		logFile.writeWarningMessage("[TF] Recieved termination signal, interrupting daemon");
		    		daemon.interrupt();
		    		while(daemon.isAlive()){}
			 }
		      
		     if (tomatoFarmerFailed) {
		    	if (logFile.getFailed()) {
		    		try {
			 			CThread.postMail(email, "TF failed, exiting.",logFile.getLastErrorMessage());
			 		} catch (MessagingException e) {
			 			logFile.writeErrorMessage("[TF] Failed to send ending email", true);
			 		} catch (NoClassDefFoundError ncdfe) {
			 			logFile.writeErrorMessage("[TF] Failed to send ending email", true);
			 		}
		    	} else {
		    		try {
			 			CThread.postMail(email, "TF was stopped prematurely.","If you did not halt TF yourself, please contact the core.");
			 		} catch (MessagingException e) {
			 			logFile.writeErrorMessage("[TF] Failed to send ending email", true);
			 		} catch (NoClassDefFoundError ncdfe) {
			 			logFile.writeErrorMessage("[TF] Failed to send ending email", true);
			 		}
		    	}
		     }
		     
		     try {
		    	 logFile.closeLogger();
		     } catch (Exception ex) {}
		     
		   }
		  });
	}
	
	public static void main(String[] args) {
		new TomatoFarmer(args);
		

	}
	
}