package edu.utah.seq.run.tempus.v3;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.seq.impl.NewAssembledSymbolList;
import org.json.JSONException;
import edu.utah.hci.bioinfo.smm.Subject;
import edu.utah.hci.bioinfo.smm.SubjectMatchMaker;
import edu.utah.seq.pmr.PMRSearch;
import edu.utah.seq.run.tempus.TempusJsonParser;
import edu.utah.seq.run.tempus.TempusPatient;
import edu.utah.seq.vcf.json.tempusv3.TempusV3Json2Vcf;
import edu.utah.seq.vcf.json.tempusv3.TempusV3JsonCollection;
import edu.utah.seq.vcf.json.tempusv3.TempusV3JsonSummary;
import edu.utah.seq.vcf.json.tempusv3.TempusV3Order;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Queries the Tempus bucket, finds old files, downloads the objects, and builds TNRunner dirs.
 * 
 Works with v3.3+ Tempus json test results, otherwise they are ignored, use TempusDataWrangler for legacy analysis.* */
public class TempusDataWranglerV3 {

	//user defined fields
	private String bucketURI = "s3://tm-huntsman";
	private File jobsDirectory = null;
	private File phiDirectory = null;
	private File smmDirectory = null;
	private File smmRegistryDirectory = null;
	private int minimumDaysOld = 7;
	private String numProcThreads = "10";
	private String awsPath = "aws";
	private String profile = "tempus";
	private boolean verbose = false;
	private File pmrParsedForTempus = null;
	private File tempusDataBucketList = null;
	
	//from the Tempus data bucket
	private ArrayList<String> jsonObjectLines = new ArrayList<String>();
	private HashMap<String, String[]> jsonFileNameDataLine = new HashMap<String, String[]>();
	private File[] allJsonFiles = null;
	private String vcfKeys = null;
	private String tarKeys = null;
	
	//internal
	private String processDate = null;
	private HashMap<String, String> envPropToAdd = new HashMap <String, String>();
	private HashSet<String> orderIds2Process = new HashSet<String>();
	private ArrayList<TempusDataWranglerTumorV3> finalOrdersToProcess  = null; //final set of orders to process
	private HashSet<String> processedIdSets = null;
	private HashMap<String, TempusV3JsonCollection> orderCollections = null;

	public static final SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
	
	//for pulling processed accessionIds
	//2025-04-07 11:26:34      13110 Patients/AAK6HT3yad/Tempus/25tnlyzo_20250527/ClinicalReport/25tnlyzo_20250527_TL-25-6KOAZ9D9H1_TL-25-Q4TTL72L1T_TL-25-T3L8TUDC61_TL-25-UYR2M4OVJR_TL-25-VA94DV6NPJ.manifest.txt
	
	private static final Pattern accessionPattern = Pattern.compile(".+_\\d+_(TL-.+).manifest.txt$");

	public TempusDataWranglerV3 (String[] args) {
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			checkAwsCli();
			
			//parse Tempus data bucket for json, vcf, and tar files from >= 2025 
			parseTempusDataBucketFile();
			
			//download all of the json files not present locally 
			downloadJsonFiles();
			
			parseManifestFile();
			
			parseJsonFiles();
			
			fetchRequiredResourceFiles();
			
			//anything to process?
			if (finalOrdersToProcess.size()!=0) {
				fetchPmrIds();
				downloadDatasetsV3();
			}
			else IO.pl("\nNothing to do, exiting...");
			

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			pl("\nDone! "+Math.round(diffTime)+" Min\n");
			
		} catch (Exception e) {
			IO.el("\nERROR running the TempusDataWrangler, aborting");
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void parseJsonFiles() throws Exception {
		
		IO.pl("\nParsing json files...");
		ArrayList<TempusV3JsonSummary> sums = new ArrayList<TempusV3JsonSummary>();
		for (File json: allJsonFiles) {
			if (json == null) continue;
			if (verbose) IO.pl("\t"+json.getName());
			TempusV3JsonSummary sum = TempusV3Json2Vcf.parseJsonNoVariants(json);
			if (sum != null) sums.add(sum);
		}
		IO.pl("\t"+sums.size()+" v3 schema jsons");

		IO.pl("\nGrouping by tempusOrderId...");
		orderCollections = new HashMap<String, TempusV3JsonCollection>();
		for (TempusV3JsonSummary s: sums) {
			String tempusOrderId = s.getTempusOrder().getTempusOrderId();
			TempusV3JsonCollection c = orderCollections.get(tempusOrderId);
			if (c==null) {
				c = new TempusV3JsonCollection(tempusOrderId);
				orderCollections.put(tempusOrderId, c);
			}
			c.getJsonSummaries().add(s);

		}
		IO.pl("\t"+orderCollections.size()+"\tOrders");


		//Check if all of the json reports in an order meet age requirements and check if they have already been processed and uploaded to the PMR
		IO.pl("\nChecking file age and whether tumor json collection was already processed...");
		for (String tempusOrderId: orderCollections.keySet()) {
			if (verbose) IO.pl("TempusOrderId: "+tempusOrderId);
			//get the collection
			TempusV3JsonCollection collection = orderCollections.get(tempusOrderId);

				ArrayList<TempusV3JsonSummary> al = collection.getJsonSummaries();

				//look at the reports associated with the tumor,
				boolean allOk = true;
				String[] accessionIds = new String[al.size()];
				for (int i=0; i< accessionIds.length; i++) {
					TempusV3JsonSummary sum = al.get(i);
					TempusV3Order order = sum.getTempusOrder();
					String jsonFileName = sum.getTempusReport().getJsonFile().getName();
					String[] dataLine = jsonFileNameDataLine.get(jsonFileName);
					boolean ok = checkAge(dataLine);
					if (ok==false) allOk = false;
					if (verbose) IO.pl("\tReportId: "+order.getAccessionId()+" "+order.getTestCode()+" "+order.getTempusOrderId()+" "+jsonFileName);
					accessionIds[i] = order.getAccessionId();
				}
				//are the ages OK? and thus ready for processing?
				if (allOk) {
					collection.setPassMinAge(true);
					Arrays.sort(accessionIds);
					String colName = Misc.stringArrayToString(accessionIds, "_");
					if (verbose) IO.pl("\tAge passes");
							
					//has it been processed before and uploaded into the PMR?
					//with new proc method this will exist:
					//Patients/AAK6HT3yad/Tempus/25tnlyzo_20250527/ClinicalVars/25tnlyzo_20250527_TL-25-6KOAZ9D9H1_TL-25-Q4TTL72L1T_TL-25-T3L8TUDC61_TL-25-UYR2M4OVJR_TL-25-VA94DV6NPJ.manifest.txt
					//List of TL reports included
					if (processedIdSets.contains(colName) == false) {
						if (verbose) IO.pl("\tNot in PMR, process: "+colName);
						orderIds2Process.add(collection.getTempusOrderId());
						collection.setAlreadyProcessed(false);
					}
					else {
						if (verbose) IO.pl("\tIn PMR, skip: "+colName);
						collection.setAlreadyProcessed(true);
					}
				}
				else if (verbose) IO.pl("\tAge fails");
			}
		IO.pl("\t"+orderIds2Process.size()+" Patients with new data to process.");
	}

	private void downloadJsonFiles() throws Exception {
		pl("\nDownloading and parsing json files year 2025 and later...");
		
		//download and parse them
		allJsonFiles = new File[jsonObjectLines.size()];
		for (int i=0; i< allJsonFiles.length; i++) {
			//2020-01-16   10:43:23  12520 result_json/result_69d44328-7c5e-4295-a6af-0da8097acc27.json
			//     0          1         2                         3
			String[] tokens = Misc.WHITESPACE.split(jsonObjectLines.get(i));
			String[] t = Misc.FORWARD_SLASH.split(tokens[3]);
			allJsonFiles[i] = new File (phiDirectory, t[t.length-1]);
			cp(tokens[3], Long.parseLong(tokens[2]), allJsonFiles[i]);
			jsonFileNameDataLine.put(allJsonFiles[i].getName(), tokens);
		}
		IO.pl("\t"+jsonFileNameDataLine.size()+" json files");
	}
	
	/**Returns provided year - 2025,  2020-01-16 = -5 */
	private static int diffFrom2025(String ymd) {
		String[] split = Misc.DASH.split(ymd);
		int year = Integer.parseInt(split[0]);
		return year - 2025;
	}

	private void fetchPmrIds() throws Exception {
		pl("\nRunning the Subject Match Maker to find or generate new PMR IDs using patient PHI...");
		
		//make IO for writing out queries
		File queryFile = new File (smmDirectory, "tempusPatientQueries_PHI.txt");
		File queryResDir = new File (smmDirectory, "SMMQueryResults_PHI");
		queryResDir.mkdirs();
		
		//Hash on patient subject match maker line
		HashMap<String, ArrayList<TempusDataWranglerTumorV3>> smm2Tumor = new HashMap<String, ArrayList<TempusDataWranglerTumorV3>>();
		for (TempusDataWranglerTumorV3 o : finalOrdersToProcess) {
			String smm = o.getCollection().getJsonSummaries().get(0).getTempusPatient().fetchSubjectMatchMakerLine();
			ArrayList<TempusDataWranglerTumorV3> al = smm2Tumor.get(smm);
			if (al == null) {
				al = new ArrayList<TempusDataWranglerTumorV3>();
				smm2Tumor.put(smm, al);
			}
			al.add(o);
		}
		
		//write out each for lookup
		PrintWriter out = new PrintWriter (new FileWriter( queryFile));
		for (String smm : smm2Tumor.keySet()) out.println(smm);	
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
		Subject[] queries = smm.getQuerySubjects();
		
		//add pmrId to the first json summary's TempusPatient
		int i=0;
		for (String key : smm2Tumor.keySet()) {
			Subject s = queries[i++];
			String coreId = s.getCoreIdNewOrMatch();
			ArrayList<TempusDataWranglerTumorV3> al = smm2Tumor.get(key);
			for (TempusDataWranglerTumorV3 o: al) {
				o.getCollection().getJsonSummaries().get(0).getTempusPatient().setPmrId(coreId);
			}
		}
	}
	
	private void downloadDatasetsV3() throws Exception {
		pl("\nDownloading files and building TNRunner run folders...");

		ArrayList<String> cmdsToExecute = new ArrayList<String>();
		ArrayList<File> tarFiles = new ArrayList<File>();
		String awsCmdProfBucket = "aws --profile "+profile+" s3 cp "+bucketURI;
				
		for (TempusDataWranglerTumorV3 o : finalOrdersToProcess) {
			//make TJob ClinicalReport and Fastq dirs with PMR Id
			o.makeJobDirectories(jobsDirectory, processDate);
			//write deidentified jsons to ClinicalReport
			o.writeDeIdentifiedJsons();
			//collect download cmds
			o.addVcfTarDownloadCmds(cmdsToExecute, awsCmdProfBucket, verbose, tarFiles);
		}		
		
		//anything to download
		if (cmdsToExecute.size()!=0) executeInParallel(cmdsToExecute);
		
		//anything to untar?
		if (tarFiles.size()!=0) {
			IO.pl("Untaring fastq datasets...");
			unTarArchives(tarFiles);
			
			//check RNA and move DNA to correct folders (normal or tumor)
			for (TempusDataWranglerTumorV3 o : finalOrdersToProcess) {
				o.checkRnaFastqDir();
				o.moveDNAFastq();
			}
		}
		
		//write the manifest.txt files
		for (TempusDataWranglerTumorV3 o : finalOrdersToProcess)  o.writeManifest();
	}
	
	private void unTarArchives(ArrayList<File> tarFiles) throws Exception {
		ArrayList<String> cmdsToExecute = new ArrayList<String>();
		for (File t: tarFiles) {
			String tarPath = t.getCanonicalPath();
			String cmd = "tar -xf "+ tarPath+ " -C "+ t.getCanonicalFile().getParent()+"/";	
			cmdsToExecute.add(cmd);
		}
		executeInParallel(cmdsToExecute);
	}
	
	private void executeInParallel(ArrayList<String> cmdsToExecute) throws Exception {
		//write out the cmds
		File tmpExec = new File(System.getProperty("java.io.tmpdir")+"/tempusDataWranglerCmdsForParallel.txt");
		IO.writeArrayList(cmdsToExecute, tmpExec);
		String[] cmd = {
				"parallel", "--will-cite", "--jobs", numProcThreads, "--halt", "soon,fail=1", "--arg-file", tmpExec.getCanonicalPath()
		};
		if (verbose)IO.pl("Executing: "+Misc.stringArrayToString(cmd, " "));
		int exitCode = executeReturnExitCode(cmd, verbose, null);
		if (exitCode !=0) throw new IOException("Parallel execution failed, for "+Misc.stringArrayToString(cmd, " "));
	}

	private void parseTempusDataBucketFile() throws IOException {
		BufferedReader in = IO.fetchBufferedReader(tempusDataBucketList);
		String line;
		
		HashSet<String> vcfObjectLines = new HashSet<String>();
		HashSet<String> tarObjectLines = new HashSet<String>();
		
		while ((line = in.readLine())!=null) {
			line = line.trim();
			//2020-01-16   10:43:23  12520 result_json/result_69d44328-7c5e-4295-a6af-0da8097acc27.json
			//     0          1         2                         3
			String[] tokens = Misc.WHITESPACE.split(line);
			int yearDiff = diffFrom2025(tokens[0]);
			if (yearDiff < 0) continue;
			
			if (line.endsWith(".json")) jsonObjectLines.add(line);
			else if (line.endsWith(".vcf")) vcfObjectLines.add(tokens[2]+"\t"+tokens[3]);
			else if (line.endsWith(".tar.gz")) tarObjectLines.add(tokens[2]+"\t"+tokens[3]);
		}
		
		in.close();
		
		// for searching
		vcfKeys = PMRSearch.fetchKeysSearchString(vcfObjectLines);
		tarKeys = PMRSearch.fetchKeysSearchString(tarObjectLines);
	}
	
	private void parseManifestFile() throws IOException {
		processedIdSets = new HashSet<String>();
		BufferedReader in = IO.fetchBufferedReader(pmrParsedForTempus);
		String line;
		Matcher mat;
		while ((line = in.readLine())!=null) {
			line = line.trim();
			IO.pl("Line: "+line);
			mat = accessionPattern.matcher(line);
			if (mat.matches()) processedIdSets.add(mat.group(1));
		}
		in.close();
	}
	
	private void fetchRequiredResourceFiles() {
		IO.pl("\nAttempting to fetch vcf and tar files for each tumor order...");
		finalOrdersToProcess = new ArrayList<TempusDataWranglerTumorV3>();
		ArrayList<TempusDataWranglerTumorV3> missingFiles = new ArrayList<TempusDataWranglerTumorV3>();
		for (String tempusOrderId: orderIds2Process) {
			//get the collection
			TempusV3JsonCollection collection = orderCollections.get(tempusOrderId);
			TempusDataWranglerTumorV3 tumor = new TempusDataWranglerTumorV3(this, collection);
			if (tumor.isReady()) finalOrdersToProcess.add(tumor);
			else missingFiles.add(tumor);
		}
		IO.pl("\t"+finalOrdersToProcess.size()+"\tOrders ready for processing");
		IO.pl("\t"+missingFiles.size()+"\tOrders skipped due to missing files");
		
		if (missingFiles.size()!=0) {
			IO.pl("\nThe following are missing one or more tar.gz and or vcf files given the test codes. Contact Tempus!\n");
			for (TempusDataWranglerTumorV3 bad: missingFiles) IO.pl(bad+"\n");
		}

	}
	
	public boolean checkAge(String[] dataLine) throws ParseException {
		//2018-11-16 13:47:32 6877541013 TL-18-29F99A/TL-18-29F99A/DNA/FastQ/TL-18-29F99A_TL-18-29F99A-DNA-fastq.tar.gz
		//    0          1         2                  3
		Date objectDate = sdf.parse(dataLine[0]+" "+dataLine[1]);
		long diff = System.currentTimeMillis() - objectDate.getTime();
		long diffDays = TimeUnit.DAYS.convert(diff, TimeUnit.MILLISECONDS);
		if (diffDays< minimumDaysOld) return false;
		return true;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new TempusDataWranglerV3(args);
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
						case 't': tmpDir = new File(args[++i]); break;
						case 'c': smmRegistryDirectory = new File(args[++i]); break;
						case 'o': minimumDaysOld = Integer.parseInt(args[++i]); break;
						case 'p': profile =args[++i]; break;
						case 'm': pmrParsedForTempus = new File(args[++i]); break;
						case 'd': tempusDataBucketList = new File(args[++i]); break;
						case 'v': verbose = true; break;
						default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
						}
					}
					catch (Exception e){
						e.printStackTrace();
						Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
					}
				}
			}
			
			//check for the pmr list file, can use the entire thing best to pull just /Tempus/ containing lines
			if (pmrParsedForTempus == null || pmrParsedForTempus.exists() == false) {
				Misc.printErrAndExit("\nError: please provide a file containing the s3 object list of the PMR bucket parsed for /Tempus/ lines.\n");
			}
			
			//check bucket
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
			
			processDate = Misc.getDateInNumbersNoSpaces();
			
			//set env var for the aws cli 
			envPropToAdd.put("AWS_SHARED_CREDENTIALS_FILE", credentialsFile.getCanonicalPath());

			//Only needed in Eclipse
			if (1==2) {
				envPropToAdd.put("PATH", "/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin");
				awsPath="/usr/local/bin/aws";
			}
			
	}

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                          Tempus Data Wrangler V3 : May 28 2025                   **\n" +
				"**************************************************************************************\n" +
				"The Tempus Data Wrangler downloads complete patient datasets from an AWS bucket, parses\n"+
				"the json test files for patient info, fetches/makes coreIds using the SubjectMatchMaker\n"+
				"app and assembles TNRunner compatible analysis directories. Deidentified versions of the\n"+
				"json reports are saved. Do look at the log output at the failing patient datasets and\n"+
				"resolve the issues (e.g. broken/missing jsons, multiple DNA/RNA tars, etc.). This tool\n"+
				"works with json schema's 3.x+ and ignores others. For each patient tumor sample, all\n"+
				"of the associated json reports (RS, cfDNA, xT, PDL, MMR, etc. ) and associated vcf\n"+
				"files in the ClinicalReport folder for processing with the USeq TempusV3Json2Vcf app.\n"+
				"Tar.gz files are uncompressed and the fastq.gz files placed in the Fastq/ folders.\n"+
				"\nRequired system apps: AWS CLI (https://aws.amazon.com/cli) and parallel \n"+
				"(https://www.gnu.org/software/parallel)\n"+
				
				"\nOptions:\n"+
				"-j Directory to build patient job folders\n" +
				"-r Directory to download Tempus json reports with PHI\n"+
				"-c Directory containing the SubjectMatchMaker 'currentRegistry_' file\n"+
				"-t Directory for writing tmp files with PHI for SubjectMatchMaker.\n"+
				"-o Minimum days old before processing any patient tumor data collection, defaults to 7\n"+
				"      Tempus releases reports as they become available and does not hold reports before\n"+
				"      the entire set is ready.  Immediately processing reports will lead to partial\n"+
				"      datasets. Each order results in multiple reports.\n"+
				"-b S3 bucket URI where test results files are placed by Tempus, defaults to \n"+
				"      's3://tm-huntsman'\n"+
				"-p AWS credentials profile for accessing your Tempus bucket, defaults to 'tempus'\n"+
				"-d File containing the list of files in your Tempus bucket, e.g. \n"+
				"      'aws --profile tempus s3 ls s3://tm-huntsman --recursive > tm-huntsman.list.txt'\n"+
				"-m File containing the list of Tempus ClinicalReport 'xxx.manifest.txt' files in your\n"+
				"      patient molecular repository bucket, e.g. 'aws s3 ls \n"+
				"      s3://hcibioinfo-patient-molecular-repo --recursive | grep _manifest.txt | grep\n"+
				"      ClinicalReport > manifestFiles.txt'\n"+
				"-v Verbose output.\n"+
				
				"\nExample: java -jar pathToUSeq/Apps/TempusDataWranglerV3 -o 5 -j ~/Tempus/TJobs/ \n"+
				"      -r ~/Tempus/JsonReports_PHI/ -c ~/SMM/Registry_PHI -t ~/SMM/SMM_Tmp_PHI \n"+
				"      -d tm-huntsman.list.txt -m manifests.txt\n"+

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

	public void cp(String awsObjectName, long awsObjectSize, File localPath) throws Exception {
		String fullS3Uri = bucketURI+awsObjectName;
		if (localPath.exists() == false || localPath.length() != awsObjectSize) {
			String[] cmd = null;
			if (verbose) {
				IO.pl("\tcp "+fullS3Uri+" exists? "+localPath.exists()+" diff size? "+(localPath.length() != awsObjectSize));
				cmd = new String[] {awsPath, "s3", "cp", fullS3Uri, localPath.getCanonicalPath(), "--profile", profile, "--no-progress"};
			}
			else cmd = new String[] {awsPath, "s3", "cp", "--quiet", fullS3Uri, localPath.getCanonicalPath(), "--profile", profile, "--no-progress"};
			int exitCode = executeReturnExitCode(cmd, true, null);
			if (exitCode != 0) {
				localPath.delete();
				throw new IOException("Error: failed to cp "+ fullS3Uri + " to "+localPath); 
			}
		}
		else if (verbose) IO.pl("\tNot copying, already exists with correct size -> "+localPath);
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

	public String getVcfKeys() {
		return vcfKeys;
	}

	public String getTarKeys() {
		return tarKeys;
	}

	public boolean isVerbose() {
		return verbose;
	}

	public String getProcessDate() {
		return processDate;
	}


}
