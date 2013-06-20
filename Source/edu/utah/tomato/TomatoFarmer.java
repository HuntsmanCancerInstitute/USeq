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


import util.gen.IO;
import util.gen.Misc;



public class TomatoFarmer {
	private TFLogger logFile;	
	private String username;
	private File rootDirectory = null;
	private String email = null;
	
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
		}
		this.tomatoFarmerFailed = false;
		logFile.writeInfoMessage("TomatoFarmer finished successfully!");
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
		
			while((line = br.readLine()) != null) {
				if (counter % 4 == 0) {
					if (!line.startsWith("@")) {
						logFile.writeErrorMessage("Fastq line 1 does not start with @,  Error at line#: " + (counter + 1) + ". Contents: " + line,false);
						System.exit(1);
					}
				} else if (counter % 4 == 1) {
					Matcher m = p.matcher(line);
					readLength = line.length();
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
		javax.mail.Session session = javax.mail.Session.getDefaultInstance(props, null);

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
		
		File configFile = null;
		String analysisType = null;
		File directory = null;
		String email = null;
		Integer wallTime = null;
		String logLevel = "INFO";
		String studyName  = "STUDY";
		
		int heartbeat = 30;
		int jobNumber = 5;
		boolean manualFail = false;
		boolean splitChroms = false;
		boolean suppressEmail = true;
		
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
					case 'y': analysisType = args[++i].toLowerCase(); break;
					case 'g': configFile = new File(args[++i]); break;
					case 'f': manualFail = true; break;
					case 'l': logLevel = args[++i]; break;
					case 'b': heartbeat = Integer.parseInt(args[++i]); break;
					case 'j': jobNumber = Integer.parseInt(args[++i]); break;
					case 's': studyName = args[++i]; break;
					case 'x': suppressEmail = false; break; 
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//***********************************************
		// Set up logging file
		//***********************************************
		
		//Get username
		this.username = System.getProperty("user.name");
				
		//Initialize log files
		logFile = new TFLogger(TFConstants.templateDir,this.username,logLevel);
		
		//Write Username to file
		logFile.writeSystemMessage("Username: " + this.username);
		logFile.writeInfoMessage("");
		
		//************************************************
		// Parse command line arguments
		//************************************************
		
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
		

		if (!directory.getAbsolutePath().startsWith("/tomato/version/job")) {
			logFile.writeErrorMessage("Your job directory must be contained within /tomato/version/job.  Offending directory: " + directory.getAbsolutePath(),false);
			System.exit(1);
		}
		
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
		} else if (!TFConstants.validTypes.contains(analysisType)) {
			logFile.writeErrorMessage("Specifed run type is not valid: " + analysisType + ". See help for options.",false);
			System.exit(1);
		}
		
		//Specify the configuration file
		if (analysisType.equals(TFConstants.ANALYSIS_CUSTOM)) {
			if (configFile == null) {
				logFile.writeErrorMessage("Configuration file (-g) not specified",false);
			} else if (!configFile.exists()) {
				logFile.writeErrorMessage("Specified configuration file does not exist, exiting",false);
			} 
		} else {
			configFile = new File(TFConstants.templateDir,"defaults/" + analysisType + ".default.txt");
			if (!configFile.exists()) {
				logFile.writeErrorMessage("Default configuration file does not exist", true);
				System.exit(1);
			}
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
				TFCommand cmd = TFCommand.getCommandObject(TFConstants.templateDir, directory, commandString, logFile, email, wallTime, heartbeat, failmax, jobNumber, suppressEmail, studyName, splitChroms);
				commandList.add(cmd);
			}
			br.close();
		} catch (IOException ioex) {
			logFile.writeErrorMessage("Could not read configuation file",true);
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
				commandList.get(0).getCommandType().equals(TFConstants.ANALYSIS_EXOME_ALIGN_NOVOALIGN)) {
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
				"**                               TomatoFarmer: June 2013                            **\n" +
				"**************************************************************************************\n" +
				"TomatoFarmer controls an exome analysis from start to finish.  It creates alignment jobs\n" +
				"for each of the samples in your directory, waits for all jobs to finish, then launches\n" +
				"metrics and variant calling jobs.  Jobs will be resubmitted up to a set number of times \n" +
				"if they fail to account for spurious CHPC errors.\n" +
				
				"\n\nRequired Arguments:\n\n"+
				"-d Job directory.  This directory must be a subdirectory of /tomato/version/job, but can\n"+
				"       be several levels below.\n"+
				"-e Email address.  You will be emailed as tomato jobs start and finish. This will happen\n" +
				"       several times during the run,\n" + 
				"-y Analysis pipeline. The analysis pipeline you want to run.  Current options are: \n" +
				"       1) exome_bwa - Full exome analysis (bwa)\n" +
				"       2) exome_novoalign - Full exome analysis (novoalign) \n" +
				"       3) exome_align_bwa - Alignment/recalibration only (bwa)\n" +
				"       4) exome_align_novoalign - Align/recalibration only (novoalign)\n" +
				"       5) exome_metrics - Sample QC metrics only\n" +
				"       6) exome_variant_raw - Variant detection and filtering (raw)\n " + 
				"       7) custom - Specify your own configuration file.\n" +		
				"   If you specify 'custom', you must use the option -g and specify an exiting configuration file. \n" +
				"   'Custom' allows you to to run older versions of the pipeline.  Non-custom pipelines always use \n " + 
				"   the most recent versions of each command.\n" +

				"\n\nOptional Arguments:\n\n"+
				"-g Configuration file.  If 'custom' is specified as the Analysis Pipeline, you must specify \n" +
				"      an existing configuration file to work off of. \n" +
				"-w Wall time.  Set this if you want less wall time than the default 240 hours. Useful\n" +
				"      when there is upcoming CHPC downtime\n" +
				"-s Study name.  Set this if you want your VCF files to have a prefix other than 'STUDY'\n" +
				"-c Split chromsome.  Set this option if you want to run variant calling on each \n" +
				"      chromosome separately.  This is good for large projects (>15 exomes). They  \n" +
				"      will be merged once they are all finished.\n" +
				"-x Unsuppress tomato emails.  Recieve both tomato and TomatoFarmer emails\n" +
				
				
				"\n\nAdmin Arguments:\n\n" +
				"-f Manually set failure level.  If this variable is set, the user will be prompted to set the number \n" +
				"      of allowed failures for each analysis step before the analysis begins.  If exome_align is set to \n" +
				"      3, once four failures are reached across all spawned threads, everything will be shut down.  You might think \n" +
				"      about changing the default settings if you're running lots of samples. Note that outright tomato failures \n" +
				"      (job never actually starts) don't count against the cap.  Must be between 1 and 20. \n" + 	
				"-l Logging level.  Level of logging you want to see.  Options INFO, WARNING, ERROR.\n" +
				"      ERROR just displays error messages, WARNING displays warning and error messages. \n" +
				"      INFO shows all three levels. Default value is INFO\n" +
				"-b Heatbeat frequency.  How often you want a thread heartbeat message in minutes, default 30\n" +
				"-j Number of jobs at a time, don't abuse this.  Default 5 jobs.\n" +
			
				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/TomatoFarmer -d /tomato/version/test/ -e \n" +
				"              bob.bobson@hci.utah.edu -s DEMO -c\n\n" +

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
		 			postMail(email, "TomatoFarmer failed, check logs", "", "Farmer@tomatosareawesome.org");
		 		} catch (MessagingException e) {
		 			logFile.writeErrorMessage("Failed to send ending email", true);
		 			e.printStackTrace();
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
