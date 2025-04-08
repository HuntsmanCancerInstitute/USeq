package edu.utah.seq.run.tempus;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.hci.bioinfo.smm.Subject;
import edu.utah.hci.bioinfo.smm.SubjectMatchMaker;
import util.gen.IO;
import util.gen.Misc;

/**Queries the Tempus bucket, finds old files, downloads the objects, and builds TNRunner dirs.
 * 
 Seeing multiple test jsons, one contained new rnaFindings but none of the older dna seq info.  Thus need to merge these.
 This makes processing these difficult since we don't know if and when new json files will come over. Wait a week? a month?
 How trigger a re analysis.* */
public class TempusDataWrangler {

	//user defined fields
	private String bucketURI = null;
	private File jobsDirectory = null;
	private File phiDirectory = null;
	private File smmDirectory = null;
	private File smmRegistryDirectory = null;
	private int minimumHrsOld = 48;
	private String awsPath = "aws";
	private String profile = "default";
	private boolean verbose = false;
	private File testIDsToSkip = null;
	
	//internal
	private HashMap<String, TempusPatient> namePatient = new HashMap<String, TempusPatient>();
	private ArrayList<String[]> jsonFilePaths = new ArrayList<String[]>();
	private HashMap<String, String> envPropToAdd = new HashMap <String, String>();
	private ArrayList<String> errorMessages = new ArrayList<String>();
	private HashSet<String> testIds2Skip = new HashSet<String>();
	private ArrayList<String> partialPatientTestIDs = new ArrayList<String>();
	private ArrayList<TempusPatient> workingPatients = new ArrayList<TempusPatient>();
	
	// 2018-11-16 13:47:22 5839004176 TL-18-41018B/RNA/FastQ/TL-18-41018B_-RNA-fastq.tar.gz
	// 2018-11-16 14:08:20 4770633766 TL-18-C174D3/TL-18-C174D3/DNA/FastQ/TL-18-C174D3_TL-18-C174D3-DNA-fastq.tar.gz
	private static final Pattern tarPattern = Pattern.compile(".+(TL-\\d\\d-.+)/.NA/.+tar.gz");
	
	// 2025-02-27 09:27:16      32556 TL-25-1CGRW7DEBN/DNA/TL-25-1CGRW7DEBN_20250227.soma.freebayes.vcf
	// 2025-02-27 09:27:28      96070 TL-25-1CGRW7DEBN/DNA/TL-25-1CGRW7DEBN_20250227.germ.pindel.vcf
	private static final Pattern vcfPattern = Pattern.compile(".+(TL-\\d\\d-.+)/.NA/.+vcf");
	
	// 2020-01-16 10:43:23      12520 result_json/result_69d44328-7c5e-4295-a6af-0da8097acc27.json
	private static final Pattern jsonPattern = Pattern.compile(".+\\sresult_json/.+\\.json");
	public static final SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	
	private Subject[] queries = null;

	public TempusDataWrangler (String[] args) {
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			checkAwsCli();
			
			//load test IDs to skip, from prior processing, can't modify the tempus tar files in their bucket
			if (testIDsToSkip != null) testIds2Skip = IO.loadFileIntoHashSet(testIDsToSkip);
			
			//pull bucket objects and load patient tar and json lines
			loadPatientLines();

			//download the json test files and parse the test id as well as phi for the subject match maker
			downloadJsonFiles();
			
			//identify tempus datasets ready for processing
			checkPatients();
			
			//any errors?
			if (errorMessages.size()!=0) {
				Misc.printErrAndExit("\nERRORS were found, fix the following and restart:\n"+ Misc.stringArrayListToString(errorMessages, "\n"));
			}
			
			//download files building TNRunner Analysis dirs
			if (workingPatients.size() !=0) {
				pl("\n"+workingPatients.size()+" patient datasets ready for processing from "+namePatient.size()+" tests...");
				findMakeJobDir();
				downloadDatasets();
				writeOutProcessedTestIDs();
				
			}
			else pl("\nNo patient datasets ready for processing.");
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			pl("\nDone! "+Math.round(diffTime)+" Min\n");
			
		} catch (Exception e) {
			IO.el("\nERROR running the TempusDataWrangler, aborting");
			e.printStackTrace();
			System.exit(1);
		}
	}


	private void downloadJsonFiles() throws Exception {
		pl("\nDownloading and parsing json files...");
		
		//download and parse them
		File[] jsonFiles = new File[jsonFilePaths.size()];
		for (int i=0; i< jsonFiles.length; i++) {
			//2020-01-16   10:43:23  12520 result_json/result_69d44328-7c5e-4295-a6af-0da8097acc27.json
			//     0          1         2                         3
			String[] tokens = jsonFilePaths.get(i);
			String[] t = Misc.FORWARD_SLASH.split(tokens[3]);
			jsonFiles[i] = new File (phiDirectory, t[1]);
			cp(tokens[3], jsonFiles[i], false);
			
			TempusJsonParser jd = new TempusJsonParser(jsonFiles[i], tokens[3]);
			// did it parse? some are broken
			if (jd.isParsed()) {
				//in skip file?
				if (testIds2Skip.contains(jd.getTestId()) == false) {
					//pull patient and if it exists, add all of the json files
					TempusPatient tp = namePatient.get(jd.getTestId());
					if (tp != null) tp.getJsonDatasets().add(jd);
				}
			}
		}
	}

	private void writeOutProcessedTestIDs() throws Exception {
		pl("\nWriting out processed test TL-xxxx IDs...");
		
		//add working dataset ids to hash
		for (TempusPatient cp : workingPatients) testIds2Skip.add(cp.getTestID());
		
		File tmp = new File (testIDsToSkip.getCanonicalPath()+".tmp.delme");
		
		PrintWriter out = new PrintWriter (new FileWriter(tmp));
		for (String id: testIds2Skip) out.println(id);	
		out.close();
		
		//now replace original
		tmp.renameTo(testIDsToSkip);
	}

	private void findMakeJobDir() throws Exception {
		pl("\nMatching patient PHI for job directory generation...");
		
		//make IO for writing out queries
		File queryFile = new File (smmDirectory, "tempusPatientQueries_PHI.txt");
		File queryResDir = new File (smmDirectory, "SMMQueryResults_PHI");
		queryResDir.mkdirs();
		
		PrintWriter out = new PrintWriter (new FileWriter( queryFile));
		
		//write out each for lookup
		for (TempusPatient cp : workingPatients) {
			String att = cp.getJsonDatasets().get(0).fetchSubjectMatchMakerLine();
			out.println(att);	
		}
		out.close();
		
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
		//in same order as that which was printed out above from workingPatients
		queries = smm.getQuerySubjects();
	}
	
	private void downloadDatasets() throws Exception {
		pl("\nDownloading files and building TNRunner run folders...");
		for (int i=0; i< workingPatients.size(); i++) {
			TempusPatient cp = workingPatients.get(i);
			Subject s = queries[i];
			String coreId = s.getCoreIdNewOrMatch();
			pl("\t"+ coreId+"/"+cp.getTestID());
			for (TempusJsonParser tjp: cp.getJsonDatasets()) tjp.printDeIDDocument();
			cp.makeJobDirsMoveJson(coreId);
			cp.downloadDatasets();
		}
	}


	private void checkPatients() throws Exception {
		
		pl("\nChecking patient datasets...");
		for (String testID : namePatient.keySet()) {
			TempusPatient cp = namePatient.get(testID);
			cp.parseFileLines(minimumHrsOld);
			if (cp.isReady()) workingPatients.add(cp);
			else partialPatientTestIDs.add(testID);
		}
		pl(workingPatients.size()+" ready for processing.");
		pl(partialPatientTestIDs.size()+" not supported or missing files. ");
	}


	private void loadPatientLines() throws Exception {
		pl("\nPulling bucket object list...");
		String[] cmd = {awsPath, "s3", "ls", bucketURI, "--profile", profile, "--recursive"};
		String[] out = executeViaProcessBuilder(cmd, false, null);

		//pull tar.gz and json files
		for (String line: out) {
			line = line.trim();
			//tar.gz?
			Matcher matTar = tarPattern.matcher(line);
			if (matTar.matches()) {
				//pull the testID
				String testID = matTar.group(1);

				//2018-11-16 13:47:32 6877541013 TL-18-29F99A/TL-18-29F99A/DNA/FastQ/TL-18-29F99A_TL-18-29F99A-DNA-fastq.tar.gz
				//    0          1         2                  3
				String[] tokens = Misc.WHITESPACE.split(line);
				if (tokens.length != 4) throw new IOException ("\nERROR: not seeing 4 parts of the .tar.gz line: "+line);

				//skip?
				if (testID != null && testIds2Skip.contains(testID) == true) pl("\tIn skip file\t"+line);
				//check the age
				else if (checkAge(tokens) == false) pl("\tToo new\t"+line);
				else {
					//fetch or create the patient
					if (testID != null) {
						TempusPatient cp = namePatient.get(testID);
						if (cp == null) {
							cp = new TempusPatient(testID, this);
							namePatient.put(testID, cp);
						}
						pl("\tAdding\t"+line);
						cp.addObjectLine(tokens);
					}
				}
			}
			else {
				//vcf?
				Matcher matvcf = vcfPattern.matcher(line);
				if (matvcf.matches()) {
					//soma or germ?
					if (line.contains(".soma") || line.contains(".germ")) {
						//pull the testID
						String testID = matvcf.group(1);

						//2018-11-16 13:47:32 6877541013 TL-25-1CGRW7DEBN/DNA/TL-25-1CGRW7DEBN_20250227.germ.pindel.vcf
						//    0          1         2                  3
						String[] tokens = Misc.WHITESPACE.split(line);
						if (tokens.length != 4) throw new IOException ("\nERROR: not seeing 4 parts of the .vcf line: "+line);

						//skip?
						if (testID != null && testIds2Skip.contains(testID) == true) pl("\tIn skip file\t"+line);
						//check the age
						else if (checkAge(tokens) == false) pl("\tToo new\t"+line);
						else {
							//fetch or create the patient
							if (testID != null) {
								TempusPatient cp = namePatient.get(testID);
								if (cp == null) {
									cp = new TempusPatient(testID, this);
									namePatient.put(testID, cp);
								}
								pl("\tAdding\t"+line);
								cp.getVcfPaths().add(tokens[3]);
							}
						}
					}
				}
				//json?
				else {
					//2020-01-16 10:43:23      12520 result_json/result_69d44328-7c5e-4295-a6af-0da8097acc27.json
					matTar = jsonPattern.matcher(line);
					if (matTar.matches()) {
						String[] tokens = Misc.WHITESPACE.split(line);
						if (tokens.length != 4) throw new IOException ("\nERROR: not seeing 4 parts of the .json line: "+line);
						
						//check the age
						if (checkAge(tokens) == false) pl("\tToo new\t"+line);
						else jsonFilePaths.add(tokens);
					}
					//else pl("\tSkipping\t"+line);
				}
			}
		}
	}

	public boolean checkAge(String[] dataLine) throws ParseException {
		//2018-11-16 13:47:32 6877541013 TL-18-29F99A/TL-18-29F99A/DNA/FastQ/TL-18-29F99A_TL-18-29F99A-DNA-fastq.tar.gz
		//    0          1         2                  3
		Date objectDate = sdf.parse(dataLine[0]+" "+dataLine[1]);
		long diff = System.currentTimeMillis() - objectDate.getTime();
		long diffHours = TimeUnit.HOURS.convert(diff, TimeUnit.MILLISECONDS);
		if (diffHours< minimumHrsOld) return false;
		return true;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new TempusDataWrangler(args);
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
						case 'i': testIDsToSkip = new File(args[++i]); break;
						default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
						}
					}
					catch (Exception e){
						e.printStackTrace();
						Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
					}
				}
			}
			
			//check prior proc files
			if (testIDsToSkip == null || testIDsToSkip.exists() == false) {
				Misc.printErrAndExit("Error: please provide a file containing Tempus TL-xxx IDs to skip, may be empty. Will be updated with the new tests that were successfully downloaded.");
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
			//if (1==2) {
			//	envPropToAdd.put("PATH", "/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin");
			//	awsPath="/usr/local/bin/aws";
			//}
			
	}

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                          Tempus Data Wrangler : March 2025                       **\n" +
				"**************************************************************************************\n" +
				"The Tempus Data Wrangler downloads complete patient datasets from an AWS bucket, parses\n"+
				"the json test file for patient info, fetches/makes coreIds using the SubjectMatchMaker\n"+
				"app and assembles TNRunner compatible analysis directories. If indicated, it will\n"+
				"delete AWS objects that were successfully downloaded. Deidentified versions of the\n"+
				"json reports are saved. Do look at the log output at the failing patient datasets and\n"+
				"resolve the issues (e.g. broken/missing jsons, multiple DNA/RNA tars, etc.).\n"+
				
				"\nOptions:\n"+
				"-b S3 URI where Tempus files are placed\n"+
				"-j Directory to build patient job folders\n" +
				"-r Directory to place PHI Tempus json reports\n"+
				"-s Directory to place temp files with PHI for Subject ID matching\n"+
				"-c Directory containing the SubjectMatchMaker 'currentRegistry_' file\n"+
				"-t Minimum hours old before downloading, defaults to 48\n"+
				"-p Credentials profile, defaults to 'default'\n"+
				"-i File containing Tempus TL-xxx IDs to skip, one per line.  This will be read and then\n"+
				"      appended with datasets downloaded. Required. Can be empty.\n"+
				

				"\nExample: java -jar pathToUSeq/Apps/TempusDataWrangler -b s3://tm-huntsman/ \n"+
				"     -j ~/Scratch/Tempus/TJobs/ -r ~/Scratch/Tempus/JsonReports_PHI/ -c \n"+
				"     ~/Scratch/SMM/Registry_PHI -s ~/Scratch/Tempus/SMM_Tmp_PHI -i ~/Scratch/Tempus/\n"+
				"     priorProcessedTestIds.txt\n"+

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
		pl("Checking the aws cli...");
		String[] cmd = {awsPath, "--version"};
		int exitCode = executeReturnExitCode(cmd, false, null);
		if ( exitCode != 0) throw new IOException("Error: 'aws --version' failed to return an appropriate response.");
	}

	private void pl(Object message) {
		System.out.println(message.toString());
	}

	public ArrayList<String> getErrorMessages() {
		return errorMessages;
	}


	public void cp(String awsObjectName, File localPath, boolean addToDeleteList) throws Exception {
		String fullS3Uri = bucketURI+awsObjectName;
		if (localPath.exists() == false) {
			String[] cmd = {awsPath, "s3", "cp", fullS3Uri, localPath.getCanonicalPath(), "--profile", profile, "--no-progress"};
			int exitCode = executeReturnExitCode(cmd, true, null);
			if ( exitCode != 0) {
				localPath.delete();
				throw new IOException("Error: failed to cp "+ fullS3Uri + " to "+localPath); 
			}
		}
	}
	public File getJobsDirectory() {
		return jobsDirectory;
	}
	public File getPhiDirectory() {
		return phiDirectory;
	}
	public String getBucketURI() {
		return bucketURI;
	}


}
