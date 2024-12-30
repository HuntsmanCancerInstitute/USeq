package edu.utah.seq.run.caris;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.hci.bioinfo.smm.Subject;
import edu.utah.hci.bioinfo.smm.SubjectMatchMaker;
import util.gen.IO;
import util.gen.Misc;

/**Queries the Caris bucket, finds old files, downloads and deletes the objects, and builds TNRunner dirs.
 * 
 * Test setup
 for x in TN22-121409_2022-03-25_10_45.xml TN22-121409_220321_A00747_0565_BHY55JDSX2_2022-03-25_10_45.vcf  
 do
 echo $x
 aws s3 cp s3://hci-caris/$x s3://hcibioinfo-nix-test/$x
 done
 aws s3 cp s3://hcibioinfo-nix-test/TN22-121409_220321_A00747_0565_BHY55JDSX2_2022-03-25_10_45.vcf s3://hcibioinfo-nix-test/DNA_TN22-121409_S55.R1.fastq.gz
 aws s3 cp s3://hcibioinfo-nix-test/TN22-121409_220321_A00747_0565_BHY55JDSX2_2022-03-25_10_45.vcf s3://hcibioinfo-nix-test/DNA_TN22-121409_S55.R2.fastq.gz
 aws s3 cp s3://hcibioinfo-nix-test/TN22-121409_220321_A00747_0565_BHY55JDSX2_2022-03-25_10_45.vcf s3://hcibioinfo-nix-test/RNA_TN22-121409_S54.R1.fastq.gz
 aws s3 cp s3://hcibioinfo-nix-test/TN22-121409_220321_A00747_0565_BHY55JDSX2_2022-03-25_10_45.vcf s3://hcibioinfo-nix-test/RNA_TN22-121409_S54.R2.fastq.gz
 * */
public class CarisDataWrangler {

	//user defined fields
	private String bucketURI = null;
	private File jobsDirectory = null;
	private File phiDirectory = null;
	private File smmDirectory = null;
	private File smmRegistryDirectory = null;
	private int minimumHrsOld = 12;
	private String awsPath = "aws";
	private String profile = "default";
	private boolean deleteAfterDownload = false;
	private boolean verbose = false;
	private String[] filePrefixesToIgnore = null;
	
	//internal
	private HashMap<String, CarisPatient> namePatient = new HashMap<String, CarisPatient>();
	private HashMap<String, String> envPropToAdd = new HashMap <String, String>();
	private ArrayList<String> errorMessages = new ArrayList<String>();
	private ArrayList<String> s3Objects2Delete = new ArrayList<String>();
	private ArrayList<String> partialPatientTestIDs = new ArrayList<String>();
	private ArrayList<CarisPatient> workingPatients = new ArrayList<CarisPatient>();
	private Subject[] queries = null;

	public CarisDataWrangler (String[] args) {
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			checkAwsCli();
			
			//pull bucket objects and load patient lines
			loadPatientLines();
			
			//age, group, download, delete
			checkPatients();
			
			//any errors?
			if (errorMessages.size()!=0) {
				Misc.printErrAndExit("\nERRORS were found, fix the following and restart:\n"+ Misc.stringArrayListToString(errorMessages, "\n"));
			}
			
			//download files building TNRunner Analysis dirs
			if (workingPatients.size() !=0) {
				pl("\n"+workingPatients.size()+" patient datasets ready for processing from "+namePatient.size()+" tests...");
				downloadXml();
				downloadDatasets();
				deleteDownloadedDatasets();
			}
			else pl("\nNo patient datasets ready for processing.");
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			pl("\nDone! "+Math.round(diffTime)+" Min\n");
			
		} catch (Exception e) {
			IO.el("\nERROR running the CarisDataWrangler, aborting");
			e.printStackTrace();
			System.exit(1);
		}
	}


	private void deleteDownloadedDatasets() throws Exception {
		if (deleteAfterDownload) {
			pl("\nDeleting downloaded datasets from S3 "+bucketURI+"...");
			File log = new File (jobsDirectory.getParentFile(), "downloadedDeleted"+Misc.getDateNoSpaces()+".txt");
			PrintWriter out = new PrintWriter( new FileWriter(log));
			out.println("Deleting aws s3 objects from "+bucketURI+ " on "+Misc.getDateTime());

			for (String s3Uri: s3Objects2Delete) {
				String[] cmd = {awsPath, "s3", "rm", s3Uri, "--profile", profile};
				int exitCode = executeReturnExitCode(cmd, true, null);
				if (exitCode != 0) {
					out.println("Failed to delete "+s3Uri);
					out.close();
					throw new IOException("Error: failed to delete "+ s3Uri);
				}
				else out.println("Deleted\t"+s3Uri);
			}
			out.close();
			
			//upload the log
			pl("\nUploading the deletion dataset log to "+bucketURI+ log.getName());
			String[] cmd = {awsPath, "s3", "cp", log.getCanonicalPath(), bucketURI+ log.getName(), "--profile", profile};
			int exitCode = executeReturnExitCode(cmd, true, null);
			if (exitCode != 0) throw new IOException("Error: failed to upload the deletion log "+ log+" to "+bucketURI);
			log.delete();
		}
		else pl("\nSkipping the deletion of downloaded datasets per option -d from "+bucketURI);
	}

	private void downloadXml() throws Exception {
		pl("\nLoading xml files for ID matching...");
		
		//make IO for writing out queries
		File queryFile = new File (smmDirectory, "carisPatientQueries_PHI.txt");
		File queryResDir = new File (smmDirectory, "SMMQueryResults_PHI");
		queryResDir.mkdirs();
		
		PrintWriter out = new PrintWriter (new FileWriter( queryFile));
		
		//write out each for lookup, check some xmls are broken
		ArrayList<CarisPatient> readyPatients = new ArrayList<CarisPatient>();
		for (CarisPatient cp : workingPatients) {
			cp.fetchXmlAndLoad();
			if (cp.isReady()) {
				String att = cp.getCarisXml().fetchSubjectMatchMakerLine();
				out.println(att);
				readyPatients.add(cp);
			}
		}
		out.close();
		workingPatients = readyPatients;
		
		//run the SMM
		String[] args = {
				"-r", smmRegistryDirectory.getCanonicalPath(),
				"-q", queryFile.getCanonicalPath(),
				"-o", queryResDir.getCanonicalPath(),
				"-v", //turn off verbosity
				"-a", //adding queries not found to registry
				"-c", //making name matching case-insensitive for EDW
				"-u"  //fill in info not found to registry entries from queries, e.g. missing mrns, otherIds
		};
		SubjectMatchMaker smm = new SubjectMatchMaker(args);
		queries = smm.getQuerySubjects();
	}
	
	private void downloadDatasets() throws Exception {
		pl("\nDownloading files and building TNRunner run folders...");
		for (int i=0; i< workingPatients.size(); i++) {
			CarisPatient cp = workingPatients.get(i);
			if (cp.isReady()) {
				Subject s = queries[i];
				String coreId = s.getCoreIdNewOrMatch();
				pl("\t"+ coreId+"/"+cp.getTestID());
				cp.makeJobDirsMoveXml(coreId);
				cp.downloadDatasets();
			}
		}
	}


	private void checkPatients() throws Exception {
		pl("\nChecking patient datasets...");
		for (String testID : namePatient.keySet()) {
			pl(testID);
			CarisPatient cp = namePatient.get(testID);
			cp.parseFileLines(minimumHrsOld);
			if (cp.isReady()) workingPatients.add(cp);
			else if (cp.isTooYoung() == false) partialPatientTestIDs.add(testID);
		}
		pl(partialPatientTestIDs.size()+" PartialTestDatasets: "+Misc.stringArrayListToString(partialPatientTestIDs, " "));
	}


	private void loadPatientLines() throws IOException {
		pl("\nPulling bucket object list...");
		String[] cmd = {awsPath, "s3", "ls", bucketURI, "--profile", profile};
		String[] out = executeViaProcessBuilder(cmd, false, null);
		
		Pattern pat = Pattern.compile(".+[\\s_]+TN\\d\\d-.+"); // looks for '_TN##-' or ' TN##-'

		for (String line: out) {
			line = line.trim();
			Matcher mat = pat.matcher(line);
			if (mat.matches() && line.endsWith(".log") == false) {
				//pull the testID
				//2022-03-16 12:00:11 11912006509 DNA_TN22-116244_S32.R1.fastq.gz
				//    0          1         2                  3
				String[] tokens = Misc.WHITESPACE.split(line);
				if (tokens.length != 4) continue;
				String fileName = tokens[3];
				String[] splitName = Misc.UNDERSCORE.split(tokens[3]);
				String testID = null;
				if (fileName.startsWith("DNA_TN") || fileName.startsWith("RNA_TN"))  testID = splitName[1];
				else if (fileName.startsWith("TN")) testID = splitName[0];
				else {
					boolean errorIt = true;
					//one of the prefixes to ignore?
					if (filePrefixesToIgnore != null) {
						for (String pre: filePrefixesToIgnore) {
							if (fileName.startsWith(pre)) {
								errorIt = false;
								break;
							}
						}
					}
					if (errorIt) {
						errorMessages.add("Failed to parse the test ID from "+fileName+ " in "+line);
						pl("\tERROR\t"+line);
					}
				}

				//fetch or create the patient
				if (testID != null) {
					CarisPatient cp = namePatient.get(testID);
					if (cp == null) {
						cp = new CarisPatient(testID, this);
						namePatient.put(testID, cp);
					}
					pl("\tAdding\t"+line);
					cp.addObjectLine(tokens);
				}
			}
			else pl("\tSkipping\t"+line);
		}
		
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CarisDataWrangler(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{

			IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
			Pattern pat = Pattern.compile("-[a-zA-Z]");
			File tmpDir = null;
			for (int i = 0; i<args.length; i++){
				Matcher mat = pat.matcher(args[i]);
				if (mat.matches()){
					char test = args[i].charAt(1);
					try{
						switch (test){
						case 'b': bucketURI = args[++i]; break;
						case 'j': jobsDirectory = new File(args[++i]); break;
						case 'r': phiDirectory = new File(args[++i]); break;
						case 's': tmpDir = new File(args[++i]); break;
						case 'c': smmRegistryDirectory = new File(args[++i]); break;
						case 't': minimumHrsOld = Integer.parseInt(args[++i]); break;
						case 'p': profile =args[++i]; break;
						case 'd': deleteAfterDownload = true; break;
						case 'f': filePrefixesToIgnore = Misc.COMMA.split(args[++i]); break;
						default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
						}
					}
					catch (Exception e){
						e.printStackTrace();
						Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
					}
				}
			}
			
			//check bucket
			if (bucketURI == null) Misc.printErrAndExit("Error: please provide an S3 URI to the folder containing the patient datasets, e.g. s3://hci-caris/");
			if (bucketURI.startsWith("s3://")==false) Misc.printErrAndExit("Error: your bucket URS does not start with s3:// ? Aborting. "+bucketURI);
			if (bucketURI.endsWith("/")==false) bucketURI = bucketURI+"/";
			
			//check credentials
			File awsCredentialsDir = new File (System.getProperty("user.home")+"/.aws");
			if (awsCredentialsDir.exists() == false)  Misc.printErrAndExit("Error: failed to find the aws credentials directory "+awsCredentialsDir);
			File credentialsFile = new File (awsCredentialsDir, "credentials").getCanonicalFile();
			if (awsCredentialsDir.exists() == false)  Misc.printErrAndExit("Error: failed to find your aws credentials file "+credentialsFile);
			
			//check jobDirectory
			if (jobsDirectory !=null) jobsDirectory.mkdirs();
			if (jobsDirectory == null || jobsDirectory.exists() == false) Misc.printErrAndExit("Error: cannot find or make your jobs directory "+jobsDirectory);
			
			//check phiDirectory
			if (phiDirectory !=null) phiDirectory.mkdirs();
			if (phiDirectory == null || phiDirectory.exists() == false) Misc.printErrAndExit("Error: cannot find or make your xml report directory "+phiDirectory);
			
			//check smm registry directory
			if (smmRegistryDirectory == null || smmRegistryDirectory.exists() == false || smmRegistryDirectory.isDirectory()==false) Misc.printErrAndExit("Error: cannot find your Subject Match Maker registry directory, "+smmRegistryDirectory);
			
			//check tmpDir for the subject match maker
			if (tmpDir == null) Misc.printErrAndExit("Error: cannot find your temp subject ID directory "+tmpDir);
			tmpDir.mkdirs();
			if (tmpDir.exists() == false) Misc.printErrAndExit("Error: cannot make your temp subject ID directory "+tmpDir);
			smmDirectory = new File (tmpDir, System.currentTimeMillis()+"_DeleteMe");
			smmDirectory.mkdir();
			
			//set env var for the aws cli 
			envPropToAdd.put("AWS_SHARED_CREDENTIALS_FILE", credentialsFile.getCanonicalPath());

			//Only needed in Eclipse
			if (false) {
				envPropToAdd.put("PATH", "/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin");
				awsPath="/usr/local/bin/aws";
			}
			
	}

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                             Caris Data Wrangler : Oct 2024                       **\n" +
				"**************************************************************************************\n" +
				"The Caris Data Wrangler downloads complete patient datasets from an AWS bucket, parses\n"+
				"the xml test file for patient info, fetches/makes coreIds using the SubjectMatchMaker\n"+
				"app and assembles TNRunner compatible analysis directories. If indicated, it will\n"+
				"delete AWS objects that were successfully downloaded. A deidentified version of the\n"+
				"xml report is saved. Do look at the log output at the failing patient datasets and\n"+
				"resolve the issues (e.g. multiple xmls, multiple vcf, for the same Caris ID).\n"+
				
				"\nOptions:\n"+
				"-b S3 URI where Caris files are placed\n"+
				"-j Directory to build patient job folders\n" +
				"-r Directory to place PHI Caris xml reports\n"+
				"-s Directory to place temp files with PHI for Subject ID matching\n"+
				"-c Directory containing the SubjectMatchMaker 'currentRegistry_' file\n"+
				"-t Minimum hours old before downloading, defaults to 12\n"+
				"-p Credentials profile, defaults to 'default'\n"+
				"-d Delete S3 objects after successful download\n"+
				"-f S3 file prefixes to ignore, comma delimited, no spaces\n"+
				

				"\nExample: java -jar pathToUSeq/Apps/CarisDataWrangler -b s3://hci-caris \n"+
				"     -j ~/Scratch/Caris/CJobs/ -r ~/Scratch/Caris/XmlReports_PHI/ -c \n"+
				"     ~/Scratch/SMM/Registry_PHI -s ~/Scratch/Caris/SMM_Tmp_PHI\n"+

				"\n**************************************************************************************\n");
	}
	
	/**Uses ProcessBuilder to execute a cmd, combines standard error and standard out into one and returns their output.
	 * @throws IOException */
	public String[] executeViaProcessBuilder(String[] command, boolean printStandardOut, File workingDirectory) throws IOException{
		if (verbose) pl ("Executing: '"+Misc.stringArrayToString(command, " ")+"'");
		ArrayList<String> al = new ArrayList<String>();
		ProcessBuilder pb = new ProcessBuilder(command);
		//add enviro props?
		pb.environment().putAll(envPropToAdd);
		if (workingDirectory !=null) pb.directory(workingDirectory);
		
		pb.redirectErrorStream(true);
		Process proc = pb.start();
		BufferedReader data = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		String line;
		while ((line = data.readLine()) != null){			
			al.add(line);
			if (printStandardOut) pl(line);
		}
		data.close();
		String[] res = new String[al.size()];
		al.toArray(res);
		return res;
	}

	/**Uses ProcessBuilder to execute a cmd, combines standard error and standard out into one and printsToLogs if indicated.
	 * Returns exit code, 0=OK, >0 a problem
	 * @throws IOException */
	public int executeReturnExitCode(String[] command, boolean printStandardOut, File workingDirectory) throws Exception{
		if (verbose) pl ("Executing: "+Misc.stringArrayToString(command, " "));
		ProcessBuilder pb = new ProcessBuilder(command);
		//add enviro props?
		pb.environment().putAll(envPropToAdd);
		if (workingDirectory !=null) pb.directory(workingDirectory);
		pb.redirectErrorStream(true);
		Process proc = pb.start();
		BufferedReader data = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		StringBuilder sb = new StringBuilder();
		String line;
		while ((line = data.readLine()) != null) {
			if (printStandardOut) pl(line);
			sb.append(line);
			sb.append("\n");
		}
		int exitCode = proc.waitFor();
		if (exitCode !=0) pl(sb.toString());
		return exitCode;
	}
	
	private void checkAwsCli() throws Exception {
		pl("\nChecking the aws cli...");
		
		String[] cmd = {awsPath, "--version"};
		int exitCode = executeReturnExitCode(cmd, false, null);
		if ( exitCode != 0) throw new IOException("Error: 'aws --version' failed to return an appropriate response.");
	}

	private void pl(Object message) {
		System.out.println(message.toString());
	}
	private void el(Object message) {
		System.err.println(message.toString());
	}

	private ArrayList<String> fetchFilesInS3Dir (String s3DirUri) throws IOException{
		String[] cmd = {awsPath, "s3", "ls", s3DirUri, "--profile", profile};
		String[] out = executeViaProcessBuilder(cmd, false, null);
		ArrayList<String> names = new ArrayList<String>();
		for (String s: out) {
			s= s.trim();
			//empty
			if (s.length() == 0) continue;
			//folder?
			if (s.startsWith("PRE ")) continue;
			String[] fields = Misc.WHITESPACE.split(s);
			String pathName = fields[fields.length-1];
			String fileName = "";
			if (pathName.contains("/")) {
				int lastSlash = pathName.lastIndexOf("/") +1;
				fileName = pathName.substring(lastSlash);
			}
			else fileName = pathName;
			names.add(fileName);
		}
		return names;
	}


	public ArrayList<String> getErrorMessages() {
		return errorMessages;
	}


	public void cp(String awsObjectName, File localPath, boolean addToDelete) throws Exception {
		String fullS3Uri = bucketURI+awsObjectName;
		if (localPath.exists()) IO.pl("\tSkipping, exists "+localPath);
		else {
			String[] cmd = {awsPath, "s3", "cp", fullS3Uri, localPath.getCanonicalPath(), "--profile", profile, "--no-progress"};
			int exitCode = executeReturnExitCode(cmd, true, null);
			if ( exitCode != 0) {
				localPath.delete();
				throw new IOException("Error: failed to cp "+ fullS3Uri + " to "+localPath); 
			}
		}
		s3Objects2Delete.add(fullS3Uri);
	}
	
	public void addObjectToDelete(String awsObjectName) {
		String fullS3Uri = bucketURI+awsObjectName;
		s3Objects2Delete.add(fullS3Uri);
	}


	public File getJobsDirectory() {
		return jobsDirectory;
	}
	public File getPhiDirectory() {
		return phiDirectory;
	}


}
