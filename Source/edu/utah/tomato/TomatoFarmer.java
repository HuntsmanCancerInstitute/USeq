package edu.utah.tomato;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
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
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;
import javax.mail.Session;

import edu.utah.tomato.model.TFCommand;
import edu.utah.tomato.model.TFSampleInfo;
import edu.utah.tomato.util.TFConstants;
import edu.utah.tomato.util.TFLogger;


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
		
		
		ArrayList<TFSampleInfo> sampleList = new ArrayList<TFSampleInfo>();
		
		for (TFCommand command: this.commandList) {
			sampleList = command.run(sampleList);
		}
			
		
		try {
			this.postMail(this.email, "TomatoFarmer finished successfully", "You can find your results here: " + this.rootDirectory.getAbsolutePath(), "Farmer@tomatosareawesome.org");
		} catch (MessagingException e) {
			logFile.writeErrorMessage("[TomatoFarmer] Failed to send ending email", true);
			e.printStackTrace();
		} catch (NoClassDefFoundError ncdfe) {
 			logFile.writeErrorMessage("[TomatoFarmer] Failed to send ending email", true);
 		}
		
		this.tomatoFarmerFailed = false;
		logFile.writeInfoMessage("[TomatoFarmer] TomatoFarmer finished successfully!");
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
		
		properties.add("SCRATCH_PATH");
		//properties.add("JAVA_MEM");
		properties.add("JAVA_PATH");
		//properties.add("NCTHREAD");
		//properties.add("NTHREAD");
		
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
		boolean suppressEmail = true;
		boolean use1KGenomes = false;
		boolean validateFastq = false;
		boolean fullGenome = false;
		boolean useNovoalign = false;
		boolean useVariantRecal = false;
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
					case 'f': manualFail = true; break;
					case 'l': logLevel = args[++i]; break;
					case 'a': heartbeat = Integer.parseInt(args[++i]); break;
					case 'j': jobNumber = Integer.parseInt(args[++i]); break;
					case 's': studyName = args[++i]; break;
					case 'x': suppressEmail = false; break; 
					case 'n': useNovoalign = true; break;
					case 'r': useVariantRecal = true; break;
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
		logFile = new TFLogger(new File(this.properties.get("TEMPLATE_PATH")),this.rootDirectory,this.username,logLevel);
		
		//Write Username to file
		logFile.writeSystemMessage("[TomatoFarmer] Username: " + this.username);
		
		//************************************************
		// Parse command line arguments
		//************************************************
	
		//Write out system arguments
		logFile.writeInfoMessage("[TomatoFarmer] " + IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " "));
		
		if (directory == null) {
			logFile.writeErrorMessage("[TomatoFarmer] You must set directory name",false);
			System.exit(1);
		}
		
		this.rootDirectory = directory;
		
		if (email == null) {
			logFile.writeErrorMessage("[TomatoFarmer] You must set email",false);
			System.exit(1);
		}
		
		this.email = email;
		

//		if (!directory.getAbsolutePath().startsWith("/tomato/version/job")) {
//			logFile.writeErrorMessage("Your job directory must be contained within /tomato/version/job.  Offending directory: " + directory.getAbsolutePath(),false);
//			System.exit(1);
//		}
		
		if (!directory.canWrite() || !directory.canRead() || !directory.canExecute()) {
			logFile.writeErrorMessage("[TomatoFarmer] Your directory is not accessable to tomatoFarmer, please edit permissions",false);
			System.exit(1);
		}
		
		
		if (!logLevel.equals("INFO") && !logLevel.equals("WARNING") && !logLevel.equals("ERROR")) {
			logFile.writeWarningMessage("[TomatoFarmer] Logging level must be INFO, WARNING or ERROR.  Reverting to INFO");
		}
		
		if (jobNumber > 20 || jobNumber < 1) {
			logFile.writeErrorMessage("[TomatoFarmer] Job number must be between 1 and 20",false);
		}
		
		logFile.writeInfoMessage("[TomatoFarmer] Using " + jobNumber + " threads for this project");
		
		if (!validateEmail(email)) {
			logFile.writeErrorMessage("[TomatoFarmer] Email does not pass our format check, please check if it was entered properly.  Offending email: " + email,false);
			System.exit(1);
		}
	
		//Check to make sure the run type is valid
		if (analysisType == null) {
			logFile.writeErrorMessage("[TomatoFarmer] You must set runtype",false);
			System.exit(1);
		}
		
		//Make sure analysisType is a valid menu option
		if (!TFConstants.validMenu.contains(analysisType)) {
			logFile.writeErrorMessage("[TomatoFarmer] " + analysisType + " is not a valid analysis type. Please see help menu for options.", false);
			System.exit(1);
		}
		
		//Add aligner to command (if necessary)
		
		String origAnalysisType = analysisType;
		if (TFConstants.alignOptions.contains(origAnalysisType)) {
			String aligner = "bwa";
			if (useNovoalign) {
				aligner = "nov";
			}
			analysisType += "_" + aligner;
		} else if (useNovoalign) {
			logFile.writeWarningMessage("[TomatoFarmer] The option -n is not valid for your analysis type, ignoring");
		}
		
		
		//Add variant command option (if necessary)
		if (TFConstants.variantOptions.contains(origAnalysisType)) {
			String variantMethod = "raw";
			if (useVariantRecal) {
				variantMethod = "vqsr";
			}
			analysisType += "_" + variantMethod;
		} else if (useVariantRecal) {
			logFile.writeWarningMessage("[TomatoFarmer] The option -r is not valid for your analysis type, ignoring");
		}
		

		if (TFConstants.validTypes.contains(analysisType)) {
			configFile = new File(this.properties.get("TEMPLATE_PATH"),"defaults/" + analysisType + ".default.txt");
			if (!configFile.exists()) {
				logFile.writeErrorMessage("[TomatoFarmer] Default configuration file for runtype " + analysisType + " does not exist", true);
				System.exit(1);
			}
			
			if (analysisType.equals("exome_bwa_raw") || analysisType.equals("exome_novo_raw") || analysisType.equals("exome_novo_vqsr") ||
					analysisType.equals("exome_bwa_vqsr") || analysisType.equals("exome_best")) {
				this.isFull = true;
			}
		} else {
			//I used support a-la carte analyses, but no-one used it.
//			configFile = new File(analysisType);
//			logFile.writeInfoMessage("[TomatoFarmer] Non-standard analysis type, looking for custom configuration file");
//			if (!configFile.exists()) {
//				logFile.writeErrorMessage("[TomatoFarmer] Specified configuration file does not exist, exiting",false);
//				System.exit(1);
//			}
			logFile.writeErrorMessage("[TomatoFarmer] " + analysisType + " is not a valid analysis type. Please see help menu for options.", false);
			System.exit(1);
		}
		
		
		
		//Make sure the target capture is valid
		File targetFile = null;
		
		if (fullGenome) {
			targetFile = new File(this.properties.get("TARGET_GENOME"));
		} else {
			if (targetType != null) {
				File targetDirectory = new File(new File(this.properties.get("TEMPLATE_PATH")),"captureRegions");
				if (TFConstants.validTargets.contains(targetType)) {
					targetFile = new File(targetDirectory,targetType + ".bed");
					if (!targetFile.exists()) {
						logFile.writeErrorMessage("[TomatoFarmer] The target capture bed file is missing: " + targetFile.getAbsolutePath(), true);
						System.exit(1);
					}
				} else {
					logFile.writeInfoMessage("[TomatoFarmer] Non-standard target capture, looking for custom target file");
					targetFile = new File(targetType);
					if (!targetFile.exists()) {
						logFile.writeErrorMessage("[TomatoFarmer] Specified target region file does not exist, exiting", false);
						System.exit(1);
					} else {
						targetFile = this.checkTarget(targetFile);
					}
				}
			} else {
				targetFile = new File(this.properties.get("TARGET_DEFAULT"));
			}
		}
		
		
		//Validate configuration file
		commandList = new ArrayList<TFCommand>();
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(configFile));
			String commandString = null;
			Scanner scanner = new Scanner(System.in);
			while ((commandString = br.readLine()) != null) {
				Integer failmax = null;
				if (manualFail) {
					while(true) {
						logFile.writeInputMessage("[TomatoFarmer] User selected manual failmax.  Please enter your selection for " + commandString + ". [Leave blank for default]");
						String input = scanner.nextLine();
						try {
							if (input.equals("")) {
								break;
							} else {
								failmax = Integer.parseInt(input);
								if (failmax < 0 || failmax > 20) {
									logFile.writeWarningMessage("[TomatoFarmer] Failmax must be between 0 and 20, try again");
								} else {
									break;
								}
							}
						} catch(NumberFormatException nfe) {
							logFile.writeWarningMessage("[TomatoFarmer] " + input + " is not a number, try again");
						}
					}
					
				} 
				TFCommand cmd = TFCommand.getCommandObject(directory, commandString, logFile, email, wallTime, heartbeat, 
						failmax, jobNumber, suppressEmail, false, false, isFull, use1KGenomes, validateFastq, studyName, 
						targetFile, this.properties, cluster);
				commandList.add(cmd);
			}
			scanner.close();
			br.close();
		} catch (IOException ioex) {
			logFile.writeErrorMessage("[TomatoFarmer] Could not read configuration file",true);
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
			this.logFile.writeInfoMessage("[TomatoFarmer]  Found gzipped target file, decompressing");
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
				this.logFile.writeErrorMessage("[TomatoFarmer]  Error decompressing target file, exiting", true);
				System.exit(1);
			}
		} else if (origTarget.getName().endsWith(".zip")) {
			this.logFile.writeInfoMessage("[TomatoFarmer]  Found zipped target file, decompressing");
			byte[] buffer = new byte[1024];
			
			try {
				ZipInputStream zis = new ZipInputStream(new FileInputStream(origTarget));
				ZipEntry ze = zis.getNextEntry();
				
				if (ze == null) {
					this.logFile.writeErrorMessage("[TomatoFarmer]  Zip file appears corruped, exiting", true);
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
				this.logFile.writeErrorMessage("[TomatoFarmer]  Error decompressing target file, exiting", true);
				System.exit(1);
			}	
		}
		
		
		//Go through the file, strip "chr"
		File cleanedFile = new File(directory,usedFile.getName() + "_Cleaned.bed");
		
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(usedFile));
			BufferedWriter bw = new BufferedWriter(new FileWriter(cleanedFile));
			
			String temp = null;
			boolean editMade = false;
			while ( (temp = br.readLine() ) != null ) {
				if (temp.startsWith("chr")) {
					editMade = true;
					String fixed = temp.substring(3);
					bw.write(fixed + "\n");
				} else {
					bw.write(temp + "\n");
				}
			}
			
			br.close();
			bw.close();
			
			if (editMade) {
				this.logFile.writeInfoMessage("[TomatoFarmer] Stripped 'chr' from target chromosomes");
				usedFile = cleanedFile;
			} else {
				cleanedFile.delete();
			}
			
		} catch (IOException ioex) {
			this.logFile.writeErrorMessage("[TomatoFarmer] Error writing cleaned target file", true);
			System.exit(1);
		}
		
		
		return usedFile;
	}
	
	private static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                TomatoFarmer: September 2014                      **\n" +
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
				"          Full pipline:\n" +
				"          1) ugp_full - UGP v1.0.5 pipeline, includes: Alignment, metrics and variant\n" +
				"             calls. Defaults to bwa and raw variant filtering\n" +
				"          2) core_full - Core best practices, includes: Alignment, metrics and variant\n" +
				"             calls.  Defaults to bwa, only raw variant filtering allowed\n" +
				"          Indivdual Steps\n"+
				"          4) ugp_align - Alignment/recalibration only. Defaults to bwa.\n" +
				"          5) core_align - Alignment/dedup only. Defalts to bwa.\n" +
				"          6) metrics - Sample QC metrics only. Requires sample level *raw.bam and \n" +
				"             *final.bam files from aligment steps.\n" +
				"          7) ugp_variant - Variant detection and filtering. Requires 1 or more \n" + 
				"             *final.bam files from alignment step. Defaults to raw variant filtering.\n" +
				"          8) core_variant - Variant detection Requires 1 or more *final.bam from \n" +
				"             alignment step.\n" +
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
		     for (TFCommand command: commandList) {
		    	 command.shutdown();
		     }
		      
		     if (tomatoFarmerFailed) {
		    	if (logFile.getFailed()) {
		    		try {
			 			postMail(email, "TomatoFarmer failed, exiting.",logFile.getLastErrorMessage(), "noreply@tomatofarmer.org");
			 		} catch (MessagingException e) {
			 			logFile.writeErrorMessage("[TomatoFarmer] Failed to send ending email", true);
			 		} catch (NoClassDefFoundError ncdfe) {
			 			logFile.writeErrorMessage("[TomatoFarmer] Failed to send ending email", true);
			 		}
		    	} else {
		    		try {
			 			postMail(email, "TomatoFarmer was stopped prematurely.","If you did not halt TF yourself, please contact the core.", "noreply@tomatofarmer.org");
			 		} catch (MessagingException e) {
			 			logFile.writeErrorMessage("[TomatoFarmer] Failed to send ending email", true);
			 		} catch (NoClassDefFoundError ncdfe) {
			 			logFile.writeErrorMessage("[TomatoFarmer] Failed to send ending email", true);
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

