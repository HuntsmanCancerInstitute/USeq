package edu.utah.tomato;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;
import javax.mail.Session;

import edu.utah.tomato.TFConstants;


import util.gen.IO;
import util.gen.Misc;



public class TomatoFarmer {
	private TFLogger logFile;	
	private String username;
	private File rootDirectory = null;
	private String email = null;
	
	private boolean isFull = false;
	private HashMap<String,String> properties = new HashMap<String,String>();
	
	
	private boolean tomatoFarmerFailed = true;
	
	//Sample list
	private ArrayList<TFSampleInfo> sampleList= null;
	
	//Analysis step management
	private ArrayList<TFCommand> commandList = null;
	
	 
	
	
	
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
		
		try {
			BufferedWriter br = new BufferedWriter(new FileWriter(new File(this.rootDirectory,"run.conf")));
			
			
			for (TFCommand command: this.commandList) {
				br.write(command.commandString);
				command.run(sampleList);
				
			}
			
			br.close();
		} catch (IOException ioex) {
			logFile.writeErrorMessage("Could not write to run configuration file",true);
			System.exit(1);
		}
		
		
		
		
		
		try {
			this.postMail(this.email, "TomatoFarmer finished successfully", "", "Farmer@tomatosareawesome.org");
		} catch (MessagingException e) {
			logFile.writeErrorMessage("Failed to send ending email", true);
			e.printStackTrace();
		} catch (NoClassDefFoundError ncdfe) {
 			logFile.writeErrorMessage("Failed to send ending email", true);
 		}
		
		this.tomatoFarmerFailed = false;
		logFile.writeInfoMessage("TomatoFarmer finished successfully!");
	}

	
	private void validatePropertiesFile(File propertiesFile) {
		ArrayList<String> properties = new ArrayList<String>();
		properties.add("DATA_PATH");
		properties.add("DATA2_PATH");
		properties.add("BWA_PATH");
		properties.add("PICARD_PATH");
		properties.add("GATK_PATH");
		properties.add("SAM_PATH");
		properties.add("R_PATH");
		properties.add("USEQ_PATH");
		properties.add("TEMPLATE_PATH");
		properties.add("SCRATCH_PATH");
		properties.add("JAVA_MEM");
		properties.add("JAVA_PATH");
		properties.add("NOVOALIGN_PATH");

		
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


	private void validateFileSet(TFSampleInfo fi) {
		logFile.writeInfoMessage("Validating fastq files for sample: " + fi.getSampleName());
		logFile.writeInfoMessage("Validating paired-end file 1");
		int lines1 = validateFastq(fi,TFConstants.FILE_FASTQ1);
		if (lines1 == 0) {
			logFile.writeErrorMessage("Paired-end fastq file 1 is empty",false);
			System.exit(1);
		}
		logFile.writeInfoMessage("Validated " + lines1 + " records");
		logFile.writeInfoMessage("Validating paired-end file 2");
		int lines2 = validateFastq(fi,TFConstants.FILE_FASTQ2);
		if (lines2 == 0) {
			logFile.writeErrorMessage("Paired-end fastq file 2 is empty",false);
			System.exit(1);
		}
		logFile.writeInfoMessage("Validated " + lines2 + " records");
		if (lines2 != lines1) {
			logFile.writeErrorMessage("Paired-end fastq files have unequal numbers of lines: " + lines1 + " " + lines2,false);
			System.exit(1);
		}
		fi.setRecordCount(lines1);
	}
	
	private int validateFastq(TFSampleInfo fi, String type) {
		int retval = 0;
		try {
			GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(fi.getFile(type)));
			BufferedReader br = new BufferedReader(new InputStreamReader(gzip));
			int counter = 0;
			int readLength = 0;
			Pattern p = Pattern.compile("^[ACTGN]+$");
			String line = null;
			int min = 100000000;
			int max = 0;
			int rl = 0;
			
		
			while((line = br.readLine()) != null) {
				if (counter % 4 == 0) {
					if (!line.startsWith("@")) {
						logFile.writeErrorMessage("Fastq line 1 does not start with @,  Error at line#: " + (counter + 1) + ". Contents: " + line,false);
						System.exit(1);
					}
				} else if (counter % 4 == 1) {
					Matcher m = p.matcher(line);
					readLength = line.length();
					if (readLength > rl) {
						rl = readLength;
					}
					if (!m.matches()) {
						logFile.writeErrorMessage("Fastq line 2 has characters other than A,C,G,T or N.  Error at line#: " + (counter + 1) + ". Contents: " + line,false);
						System.exit(1);
					}
				} else if (counter % 4 == 2) {
					if (!line.startsWith("+")) {
						logFile.writeErrorMessage("Fastq line 3 does not start with +.  Error at line#: " + (counter + 1) + ". Contents: " + line,false);
						System.exit(1);
					}
				} else if (counter % 4 == 3) {
					if (line.length() != readLength) {
						logFile.writeErrorMessage("Fastq line 4 length does not match fastq line 2,  Error at line#: " + (counter + 1) + ". Contents: " + line,false);
						System.exit(1);
					}
					//Quality check
					for (int i = 0; i < line.length(); i++){
					    char c = line.charAt(i);        
					  
					    int val = (int)c;
					    if (val > max) {
					    	max = val;
					    } else if (val < min) {
					    	min = val;
					    }
					    
					}
				} else {
					logFile.writeErrorMessage("Should not reach this line",true);
					System.exit(1);
				}
				counter++;
			}
			
			if (min-33 >= 31) {
				logFile.writeInfoMessage("Fastq detected as ASCII-64.  Min: " + (min-64) + ". Max: " + (max-64));
				fi.setQual64();
			} else {
				logFile.writeInfoMessage("Fastq detected as ASCII-33.  Min: " + (min-33) + ". Max: " + (max-33));
				
			}
			
			
			if (counter % 4 != 0) {
				logFile.writeErrorMessage("Fastq file number not divible by 4, exiting",false);
				System.exit(1);
			}
			
			br.close();
			retval = counter / 4;
			
			fi.setReadLength(rl);
		} catch (Exception ioex) {
			logFile.writeErrorMessage("Error reading fastq file, exiting: " + ioex.getMessage(),true) ;
			ioex.printStackTrace();
			System.exit(1);
		}
		return retval;
	}
	
	private void postMail(String recipients, String subject, String message, String from) throws MessagingException {

		//set the host smtp address
		Properties props = new Properties();
		props.put("mail.smtp.host", "hci-mail.hci.utah.edu");

		//create some properties and get the default Session
		Session session = Session.getDefaultInstance(props, null);

		//create message
		Message msg = new MimeMessage(session);

		//set the from and to address
		InternetAddress addressFrom = new InternetAddress(from);
		msg.setFrom(addressFrom);
		msg.setRecipients(Message.RecipientType.TO, InternetAddress.parse(recipients, false));

		//optional: can also set custom headers here in the Email if wanted
		//msg.addHeader("MyHeaderName", "myHeaderValue");

		//setting the Subject and Content type
		msg.setSubject(subject);
		msg.setContent(message, "text/plain");
		Transport.send(msg);
	}
	
	
	private void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		
		String targetType = null;
		File directory = null;
		File configFile = null;
		File propertiesFile = null;
		String email = null;
		Integer wallTime = null;
		String logLevel = "INFO";
		String studyName  = "STUDY";
		String analysisType = null;
		
		int heartbeat = 30;
		int jobNumber = 5;
		boolean manualFail = false;
		boolean splitChroms = false;
		boolean noSplit = false;
		boolean suppressEmail = true;
		boolean use1KGenomes = false;
		boolean deleteMetricsBam = true;
		boolean deleteReducedBam = true;
		
		
		
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
					case 'c': splitChroms = true; break;
					case 'n': noSplit = true; break;
					case 'y': analysisType = args[++i]; break;
					case 't': targetType = args[++i]; break;
					case 'f': manualFail = true; break;
					case 'l': logLevel = args[++i]; break;
					case 'b': heartbeat = Integer.parseInt(args[++i]); break;
					case 'j': jobNumber = Integer.parseInt(args[++i]); break;
					case 's': studyName = args[++i]; break;
					case 'x': suppressEmail = false; break; 
					case 'p': propertiesFile = new File(args[++i]); break; 
					case 'h': printDocs(); System.exit(0);
					case 'u': deleteMetricsBam = false; break;
					case 'v': deleteReducedBam = false; break;
					case 'g': use1KGenomes = true; break;
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
		logFile = new TFLogger(new File(this.properties.get("TEMPLATE_PATH")),this.rootDirectory,this.username,logLevel);
		
		//Write Username to file
		logFile.writeSystemMessage("Username: " + this.username);
		logFile.writeInfoMessage("");
		
		//************************************************
		// Parse command line arguments
		//************************************************
		
		String splitType = "chunk";
		
		//Determine run type
		if (splitChroms == true) {
			if (noSplit == true) {
				logFile.writeErrorMessage("Argument '-c' and '-n' can't be used together.", true);
				System.exit(1);
			} else {
				splitType = "chrom";
			}
		} else if (noSplit == true) {
			splitType = "none";
		}
		
		
		
		//Write out system arguments
		logFile.writeInfoMessage(IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " "));
		
		if (directory == null) {
			logFile.writeErrorMessage("You must set directory name",false);
			System.exit(1);
		}
		
		this.rootDirectory = directory;
		
		if (email == null) {
			logFile.writeErrorMessage("You must set email",false);
			System.exit(1);
		}
		
		this.email = email;
		

//		if (!directory.getAbsolutePath().startsWith("/tomato/version/job")) {
//			logFile.writeErrorMessage("Your job directory must be contained within /tomato/version/job.  Offending directory: " + directory.getAbsolutePath(),false);
//			System.exit(1);
//		}
		
		if (!directory.canWrite() || !directory.canRead() || !directory.canExecute()) {
			logFile.writeErrorMessage("Your directory is not accessable to tomatoFarmer, please edit permissions",false);
			System.exit(1);
		}
		
		
		if (!logLevel.equals("INFO") && !logLevel.equals("WARNING") && !logLevel.equals("ERROR")) {
			logFile.writeWarningMessage("Logging level must be INFO, WARNING or ERROR.  Reverting to INFO");
		}
		
		if (jobNumber > 20 || jobNumber < 1) {
			logFile.writeErrorMessage("Job number must be between 1 and 20",false);
		}
		
		logFile.writeInfoMessage("Using " + jobNumber + " threads for this project");
		
		if (!validateEmail(email)) {
			logFile.writeErrorMessage("Email does not pass our format check, please check if it was entered properly.  Offending email: " + email,false);
			System.exit(1);
		}
	
		//Check to make sure the run type is valid
		if (analysisType == null) {
			logFile.writeErrorMessage("You must set runtype",false);
			System.exit(1);
		}
	
		//Specify the configuration file
		if (TFConstants.validTypes.contains(analysisType)) {
			configFile = new File(this.properties.get("TEMPLATE_PATH"),"defaults/" + analysisType + ".default.txt");
			if (!configFile.exists()) {
				logFile.writeErrorMessage("Default configuration file for runtype " + analysisType + " does not exist", true);
				System.exit(1);
			}
			
			if (analysisType.equals("exome_bw_raw") || analysisType.equals("exome_novo_raw") || analysisType.equals("exome_novo_vqsr") ||
					analysisType.equals("exome_bw_vqsr")) {
				this.isFull = true;
			}
		} else {
			configFile = new File(analysisType);
			logFile.writeInfoMessage("Non-standard analysis type, looking for custom configuration file");
			if (!configFile.exists()) {
				logFile.writeErrorMessage("Specified configuration file does not exist, exiting",false);
				System.exit(1);
			}
		}
		
		
		
		//Make sure the target capture is valid
		File targetFile = null;
		
		if (targetType != null) {
			File targetDirectory = new File(new File(this.properties.get("TEMPLATE_PATH")),"captureRegions");
			if (TFConstants.validTargets.contains(targetType)) {
				targetFile = new File(targetDirectory,targetType + ".bed");
				if (!targetFile.exists()) {
					logFile.writeErrorMessage("The target capture bed file is missing: " + targetFile.getAbsolutePath(), true);
					System.exit(1);
				}
			} else {
				logFile.writeInfoMessage("Non-standard target capture, looking for custom target file");
				targetFile = new File(targetType);
				if (!targetFile.exists()) {
					logFile.writeErrorMessage("Specified target region file does not exist, exiting", false);
					System.exit(1);
				}
			}
		} else {
			targetFile = new File(this.properties.get("TARGET_DEFAULT"));
		}
		

		//Validate configuration file
		commandList = new ArrayList<TFCommand>();
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(configFile));
			String commandString = null;
			while ((commandString = br.readLine()) != null) {
				Integer failmax = null;
				if (manualFail) {
					Scanner scanner = new Scanner(System.in);
					while(true) {
						logFile.writeInputMessage("User selected manual failmax.  Please enter your selection for " + commandString + ". [Leave blank for default]");
						String input = scanner.nextLine();
						try {
							if (input.equals("")) {
								break;
							} else {
								failmax = Integer.parseInt(input);
								if (failmax < 1 || failmax > 20) {
									logFile.writeWarningMessage("Failmax must be between 1 and 20, try again");
								} else {
									break;
								}
							}
						} catch(NumberFormatException nfe) {
							logFile.writeWarningMessage(input + " is not a number, try again");
						}
					}
					scanner.close();
				} 
				TFCommand cmd = TFCommand.getCommandObject(directory, commandString, logFile, email, wallTime, heartbeat, 
						failmax, jobNumber, suppressEmail, deleteMetricsBam, deleteReducedBam, isFull, use1KGenomes, studyName, 
						splitType, targetFile, this.properties);
				commandList.add(cmd);
			}
			br.close();
		} catch (IOException ioex) {
			logFile.writeErrorMessage("Could not read configuration file",true);
			System.exit(1);
		}
		
		
		//**************************************
		// Look through the files in the directory and assign filetypes
		//**************************************
		
		//Create file lists
		sampleList = TFSampleInfo.identifyAndValidateSamples(directory, commandList.get(0).getCommandType(), this.logFile);
		
		logFile.writeInfoMessage("Found " + this.sampleList.size() + " samples in the study " + studyName);
		
		if (sampleList.size() == 0) {
			logFile.writeErrorMessage("No valid input file found",false);
			System.exit(1);
		}
		
		
		//************************************
		// Validate fastq files, if needed
		//************************************
		
		
		if (commandList.get(0).getCommandType().equals(TFConstants.ANALYSIS_EXOME_ALIGN_BWA) || 
				commandList.get(0).getCommandType().equals(TFConstants.ANALYSIS_EXOME_ALIGN_NOVO)) {
			for (TFSampleInfo fi: sampleList) {
				validateFileSet(fi);
			}
		}

		
		//************************************
		// Attach shutdown
		//***********************************
		attachShutDownHook();
		
	}	
	
	private static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                TomatoFarmer: June 2013                           **\n" +
				"**************************************************************************************\n" +
				"TomatoFarmer controls an exome analysis from start to finish.  It creates alignment \n" +
				"jobs for each of the samples in your directory, waits for all jobs to finish and then \n" +
				"launches metrics and variant calling jobs.  Jobs will be resubmitted up to a set \n" +
				"number of times to combat spurious CHPC erors.  Job directories are left behind so \n" +
				"you can save your log files.\n\n" +
				
				"\nRequired Arguments:\n\n"+
				"-d Job directory.  This directory must be a subdirectory of /tomato/version/job. Can \n" +
				"       can be several directory levels below /tomato/version/job. \n" +
				"       Example: '-d /tomato/version/job/krustofsky/demo'. \n" +
				"-e Email address.  TomatoFarmer emails you once the job completes/fails. You can also\n" +
				"       opt to get all tomato emails as individual jobs start/end (see option -x).\n" +
				"       Example: '-e hershel.krustofsky@hci.utah.edu'.\n" +
				"-y Analysis pipeline. The analysis pipeline or step to run.  Current options are: \n" +
				"          Full Pipeline \n" +
				"          1) exome_bwa_raw - Full exome analysis, using bwa and raw filtering.\n" +
				"          2) exome_bwa_vqsr - Full exome analysis, using bwa and vqsr filtering\n" +
				"          3) exome_novo_raw - Full exome analysis, using novoalign and raw filtering\n" +
				"          4) exome_novo_vqsr - Full exome analysis, using novoalign and vqsr filtering\n\n"+
				"          A la carte\n"+
				"          5) exome_align_bwa - Alignment/recalibration only (bwa).\n" +
				"          6) exome_align_novo - Align/recalibration only (novoalign).\n" +
				"          7) exome_metrics - Sample QC metrics only. Requires *mate.bam and \n" +
				"             *split.lane.bam from 1, 2, 3, 4, 5 or 6 in the launch directory.\n" +
				"          8) exome_variant_raw - Variant detection and filtering (raw settings). \n" + 
				"             Requires *reduced.bam from 1, 2, 3, 4, 5 or 6 in the launch directory.\n" + 
				"          9) exome_variant_vqsr - Variant detection and filtering (vqsr) \n" +
				"             Requires *reduced.bam from 1, 2, 3, 4, 5 or 6 in the launch directory.\n" + 
				"          10) path to configuration file.  This allows you to use older versions of \n" +
				"             the pipeline.\n" +		
				"      If you want run older versions of the pipeline, you must provide a configuration\n" +
				"      file, otherwise the most recent versions of each command are used. \n" +
				"      Example: '-y exome_bwa'.\n" +  
				"-p Properties file.  This file contains a list of cluster-specific paths and options \n" +
				"      this file doesn't need to be changed by the user. Example: '-p properties.txt' \n" +
				"\nOptional Arguments:\n\n"+
				"-t Target regions. Setting this argument will restrict coverage metrics and variant \n" +
				"      detection to targeted regions.  This speeds up the variation detection process\n" +
				"      and reduces noise. Options are:\n" +
				"          1) AgilentAllExonV4\n" +
				"          2) AgilentAllExonV5\n" +
				"          3) AgilentAllExon50MB\n" +
				"          4) NimbleGenEZCapV2\n" +
				"          5) NimbleGenEZCapV3\n" +
				"          6) TruSeq\n" +
				"          7) path to custom targed bed file.\n" +
				"      If nothing is specifed for this argument, the full genome will be queried for \n" +
				"      variants and ccds exomes will be used for capture metrics. Example: '-t truseq'.\n" +
				"-g 1K Genome samples.  Use this option if you want to spike in 200 1K genome samples \n" +
				"      as the background sample set.  This should improve VQSR variant calling and \n" +
				"      VAAST, but it will take a lot more time to process. BETA, ONLY WORKS FOR CORE \n" +
				"      USERS!!!\n" +
				"-w Wall time.  Use this option followed by a new wall time, in hours, if you want less\n" +
				"      wall time than the default 240 hours. Useful when there is upcoming CHPC \n" +
				"      downtime. Example: '-w 40'. \n" +
				"-s Study name.  Set this if you want your VCF files to have a prefix other than \n" +
				"      'STUDY'. Example: '-s DEMO'.\n" +
				"-c Split chromsomes.  Set this option if you want to run variant calling on each \n" +
				"      chromosome separately. By default, the genome is split by callable region. \n" +
				"      Example: '-c'.\n" +
				"-n No splitting.  Set this option if you want to run variant calling on the entire\n" +
				"      genome at one time.  This is only suggested when you have a very small capture\n" +
				"      region.  By default, the genome is split by callable region,  Example '-n'. \n" +
				"-x Unsuppress tomato emails.  Receive both tomato and TomatoFarmer emails. \n" +
				"      Example: '-x'.\n" +
				"-u Don't delete metrics bams (full pipeline specific).  Use this option if you want \n "+
				"      to keep the bams used to generate the metrics output.  The two bam files are \n" +
				"      the unprocessed lane-level sample bam and the processed lane-level sample bam \n" +
				"      files.  Example '-u' \n" +
				"-v Don't delete reduced bams (full pipeline specific.) Use this option if you want \n " +
				"      to keep the reduced bams used to generate variant calls.  Reduced bams are \n" +
				"      efficient in GATK pipelines but don't necessarily work in non-GATK software. \n" +
				"      We recommend using the files labeled SAMPLE_NAME.bam for downstream apps\n. " +
				"      Example '-v'\n" +
				
				"\nAdmin Arguments:\n\n" +
				"-f Manually set failure level.  If this variable is set, the user will be prompted to \n" + 
				"      override default allowed failures for each analysis step before the analysis \n" + 
				"      begins.  If exome_align is set to 3, once four failures are reached across all \n" +
				"      spawned threads, everything will be shut down.  You might think about changing \n" +
				"      the default settings if you're running lots of samples. Note that outright \n" + 
				"      tomato failures (job never actually starts) don't count against the cap.  Must \n" +
				"      be between 1 and 20. \n" + 	
				"-l Logging level.  Level of logging you want to see.  Options INFO, WARNING, ERROR.\n" +
				"      ERROR just displays error messages, WARNING displays warning and error messages.\n" +
				"      INFO shows all three levels. Default: INFO.\n" +
				"-b Heartbeat frequency.  How often you want a thread heartbeat message in minutes. \n" + 
				"      Default: 30 mins.\n" +
				"-j Number of jobs at a time. Don't abuse this, big brother is watching you. Default:\n" +
				"      5 jobs.\n" +
			
				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/TomatoFarmer -d /tomato/version/job/demo/\n" +
				"      -e herschel.krustofsky@hci.utah.edu -y exome_bwa -s DEMO -c -t AgilentAllExon50MB\n\n" +

				"**************************************************************************************\n");

	}
	
	private void attachShutDownHook(){
		  Runtime.getRuntime().addShutdownHook(new Thread() {
		   @Override
		   public void run() {
		     for (TFCommand command: commandList) {
		    	 command.shutdown();
		     }
		      
		     if (tomatoFarmerFailed) {
		    	try {
		    		
		 			postMail(email, "TomatoFarmer failed, check logs ","", "Farmer@tomatosareawesome.org");
		 		} catch (MessagingException e) {
		 			logFile.writeErrorMessage("Failed to send ending email", true);
		 		} catch (NoClassDefFoundError ncdfe) {
		 			logFile.writeErrorMessage("Failed to send ending email", true);
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

