package edu.utah.seq.pmr;

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
import util.gen.IO;
import util.gen.Misc;

/**Use to identify older redundant Avatar analysis that should be deleted.  E.g. when M2Gen sends along new or updated samples that complete a prior partial analysis 
 * or provide a matched sample when a prior run used a non platform matched sample.  This won't catch everything but should reduce the duplicates.*/
public class PMRAvatarDedup {

	//user defined fields
	private File clinicalReportDir = null;
	private String awsPatientDirUri = "s3://hcibioinfo-patient-molecular-repo/Patients/"; //must end with /
	private String awsPath = "aws";
	private String profile = "default";
	private HashMap<String, String> envPropToAdd = new HashMap <String, String>();
	private File awsRepoList = null;
	private HashMap<String, Patient> idPatient = new HashMap<String, Patient>();

	//names
	private String clinicalReportDirName = "ClinicalReport/";
	private String avatarSourceName = "Avatar";

	public PMRAvatarDedup (String[] args) {
		long startTime = System.currentTimeMillis();
		try {
			processArgs(args);

			downloadAwsFileList();

			//only Avatar patients are created
			loadPatientFiles();

			fetchAvatarClinicalJsons();

			loadAvatarClinicalJsons();

			deDupAvatarDatasets();
			
			//hmm, can't delete since the fastq's aren't all copied into the most recent version, must copy them over
			
			//look for ADD_SEQ_DATA ?
			
			//see /Ce8yq7PT/Avatar/

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
		} catch (Exception e) {
			IO.el("\nERROR running the SearchPMR...");
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void deDupAvatarDatasets() throws IOException {
		IO.pl("Examining Avatar analysis for legacy partial datasets and those that have been updated...");
		ArrayList<AvatarDataset> complete = new ArrayList<AvatarDataset>();
		ArrayList<AvatarDataset> incomplete = new ArrayList<AvatarDataset>();
		ArrayList<String> toDelete = new ArrayList<String>();

		for (Patient p: idPatient.values()) {
			//all of these are Avatar
			//just one dataset? then skip
			if (p.getIdDataSets().size() < 2) {
				IO.pl("Skipping "+p.getCoreId()+", just one Avatar analysis, nothing to dedup.");
			}
			else {
				IO.pl("Checking "+p.getCoreId()+"...");

				//for each patient split their datasets by tumor sample id, some will be NA
				HashMap <String, ArrayList<AvatarDataset>> tumorDatasets = splitByTumorId(p);

				//for each tumor id
				ArrayList<AvatarDataset> naDatasets = tumorDatasets.get("NA");
				for (String tumorId: tumorDatasets.keySet()) {
					//Skip NA
					if (tumorId.equals("NA")) continue;

					ArrayList<AvatarDataset> ads = tumorDatasets.get(tumorId);
					IO.pl("\tTumorId: "+tumorId);

					//just one, skip
					if (ads.size()== 1) {
						IO.pl("\t\t"+ads.get(0).getDataset().getDatasetId()+" just one, so skip");
						ArrayList<AvatarDataset> naToRemove = findMatchingNasToRemove(ads.get(0), naDatasets);
						for (AvatarDataset na: naToRemove) {
							IO.pl("\t\t\t"+na.getDataset().getDatasetId()+" is a matched NA, so DELETE");
							toDelete.add(p.getCoreId()+"/Avatar/"+na.getDataset().getDatasetId());
						}
					}
					//more than one so split into complete and incomplete
					else {
						complete.clear();
						incomplete.clear();
						for (AvatarDataset ad: ads) {
							//has all three datasets and isn't mixed
							if (ad.isComplete()) complete.add(ad);
							else incomplete.add(ad);
						}

						//best case one complete
						if (complete.size() == 1) {
							IO.pl("\t\t"+complete.get(0).getDataset().getDatasetId()+" is complete, so skip");
							for (AvatarDataset ad: incomplete) {
								IO.pl("\t\t\t"+ad.getDataset().getDatasetId()+" is incomplete, so DELETE");
								toDelete.add(p.getCoreId()+"/Avatar/"+ad.getDataset().getDatasetId());
							}
							if (naDatasets != null) {
								for (AvatarDataset na: naDatasets) {
									IO.pl("\t\t\t"+na.getDataset().getDatasetId()+" is NA, so DELETE");
									toDelete.add(p.getCoreId()+"/Avatar/"+na.getDataset().getDatasetId());
								}
							}
						}
						//more than one complete
						else if (complete.size() > 1) {
							IO.pl("\t\tMore than one complete, REVIEW:");
							for (AvatarDataset ad: complete) IO.pl("\t\t\t"+ad.getDataset().getDatasetId()+" is complete");
							for (AvatarDataset ad: incomplete) IO.pl("\t\t\t"+ad.getDataset().getDatasetId()+" is incomplete");
							if (naDatasets != null) {
								for (AvatarDataset na: naDatasets) {
									IO.pl("\t\t\t"+na.getDataset().getDatasetId()+" is NA, so DELETE");
									toDelete.add(p.getCoreId()+"/Avatar/"+na.getDataset().getDatasetId());
								}
							}
						}
						//no complete
						else {
							IO.pl("\t\tNo complete, leave:");
							for (AvatarDataset ad: incomplete) IO.pl("\t\t\t"+ad.getDataset().getDatasetId()+" is incomplete");	
						}
					}

				}
				//see lots of completes that just differ with the transcriptome ID?  anyway to prioritize these?
				//if same ORIENSpecimenID, take the more recent?  Some don't have analysis, eg Zz2Yb9Aw, do need to check if fastq?
			}
		}
		
		if (toDelete.size() == 0) IO.pl("No redundant datasets to delete!");
		else {
			IO.pl(toDelete.size()+"\tRedundant datasets to delete, consider executing:");
			for (String ad: toDelete) IO.pl("aws s3 rm --recursive "+awsPatientDirUri+ ad);
		}
		
//Hmm, many need analysis, just have fastq and the json
	}



	private ArrayList<AvatarDataset> findMatchingNasToRemove(AvatarDataset avatarDataset, ArrayList<AvatarDataset> naDatasets) {
		ArrayList<AvatarDataset> toDelete = new ArrayList<AvatarDataset>();
		if (naDatasets != null) {
			String completeNormalId = avatarDataset.getNormalExomeId();
			String completeTumorTrans = avatarDataset.getTumorTranscriptomeId();

			for (AvatarDataset na: naDatasets) {
				if (na.getNormalExomeId().equals(completeNormalId) || na.getTumorTranscriptomeId().equals(completeTumorTrans)) toDelete.add(na);
			}
		}
		return toDelete;
	}

	private HashMap<String, ArrayList<AvatarDataset>> splitByTumorId(Patient p) throws IOException {
		HashMap <String, ArrayList<AvatarDataset>> tumorDatasets = new HashMap <String, ArrayList<AvatarDataset>>();
		for (Dataset d: p.getIdDataSets().values()) {
			// A032049_SL419345_SL419548_SL420681
			AvatarDataset ad = new AvatarDataset(d);
			ArrayList<AvatarDataset> al = tumorDatasets.get(ad.getTumorExomeId());
			if (al == null) {
				al = new ArrayList<AvatarDataset>();
				tumorDatasets.put(ad.getTumorExomeId(), al);
			}
			al.add(ad);
		}
		return tumorDatasets;
	}

	private void loadAvatarClinicalJsons() {
		IO.pl("Loading Avatar Clinical Json files...");
		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				if (d.getSource().equals(avatarSourceName)) {
					if (d.getClinicalInfoFiles() == null) IO.el("WARNING! "+p.getCoreId()+" "+d.getDatasetId()+" Avatar dataset is missing a json!");
					else if (d.getClinicalInfoFiles().size() !=1) IO.el("WARNING! "+p.getCoreId()+" "+d.getDatasetId()+" Avatar dataset has more than one json!");
					else {
						AvatarClinicalInfo aci =  new AvatarClinicalInfo(d.getClinicalInfoFiles().get(0));
						d.setAvatarClinicalInfo(aci);
					}
				}
			}
		}
	}

	private void fetchAvatarClinicalJsons() throws Exception {
		IO.pl("\nDownloading Avatar Clinical Json files...");
		int numAvatarJsons = 0;
		for (Patient p: idPatient.values()) {
			for (Dataset d: p.getIdDataSets().values()) {
				//IO.pl(d.getSource()+" -> "+d.getDatasetId());
				if (d.getSource().equals(avatarSourceName)) {
					for (String partPath: d.getPartialPaths()) {
						//IO.pl(partPath);
						if (partPath.contains(clinicalReportDirName) && partPath.endsWith("json")) {
							numAvatarJsons++;
							String relativePath = p.getCoreId()+"/"+d.getSource()+"/"+d.getDatasetId()+"/"+partPath;
							File j = new File (clinicalReportDir, relativePath);
							if (j.exists() == false) {
								if (download(awsPatientDirUri+relativePath, j) == false) throw new IOException("Failed to download "+j);
							}
							d.setClinicalInfoFile(j);
						}
					}
				}
			}
		}
		IO.pl("\t"+numAvatarJsons);
	}

	private boolean download(String s3Uri, File j) throws Exception {
		IO.pl("Fetching file from AWS "+s3Uri+ " -> "+j);
		String[] cmd = {awsPath, "s3", "cp", s3Uri, j.getCanonicalPath(), "--profile", profile};
		int exitCode = executeReturnExitCode(cmd);
		if (exitCode != 0 || j.exists()==false) return false;
		return true;

	}

	private void loadPatientFiles() throws IOException {
		IO.pl("Loading Patients...");
		BufferedReader in = IO.fetchBufferedReader(awsRepoList);
		String line = null;
		String[] fields = null;
		String[] keys = null;
		String toFind = "/"+ avatarSourceName + "/";
		while ((line = in.readLine())!=null) {
			//2023-03-01   09:04:43   2938    Patients/AA2mF6Vy/Avatar/A032049_SL419345_SL419548_SL420681/ClinicalReport/A032049_SL419345_SL419548_SL420681_IDTv1_SAR_F.json
			//  date         time     size       key
			//    0            1        2         3
			// is it an Avatar json file line?
			if (line.contains(toFind) == false || line.endsWith(".json") == false) continue;

			fields = Misc.TAB.split(line);
			if (fields.length != 4) throw new IOException("Failed to find 4 fields in "+line);

			//Patients   AA2mF6Vy   Avatar   A032049_SL419345_SL419548_SL420681     ClinicalReport/A032049_SL419345_SL419548_SL420681_IDTv1_SAR_F.json
			//   0          1         2                 3                                  4
			keys = Misc.FORWARD_SLASH.split(fields[3]);
			if (keys.length < 5) throw new IOException("Failed to find more then 4 dirs in "+line);

			Patient p = fetchPatient(keys[1]);
			p.addDataFile(keys);
		}
		in.close();
	}

	private Patient fetchPatient(String id) {
		Patient p = idPatient.get(id);
		if (p == null) {
			p = new Patient(id);
			idPatient.put(id, p);
		}
		return p;
	}

	/*Looks for the repo file list file, downloads it if not found for today and saves it.*/
	private void downloadAwsFileList() throws IOException {
		String date = Misc.getDateNoSpaces();
		awsRepoList = new File (clinicalReportDir, "awsRepoList"+date+".txt");
		if (awsRepoList.exists() == false) fetchFilesInRepo();
	}

	/**Uses ProcessBuilder to execute a cmd
	 * Returns exit code, 0=OK, >0 a problem
	 * @throws IOException */
	public int executeReturnExitCode(String[] command) throws Exception{
		IO.pl ("Executing: "+Misc.stringArrayToString(command, " "));
		ProcessBuilder pb = new ProcessBuilder(command);
		//add enviro props?
		pb.environment().putAll(envPropToAdd);
		pb.redirectErrorStream(true);
		Process proc = pb.start();
		return proc.waitFor();
	}


	/**Uses ProcessBuilder to execute a cmd, combines standard error and standard out into one and returns their output.
	 * @throws IOException */
	public String[] executeViaProcessBuilder(String[] command, boolean printStandardOut, File workingDirectory) throws IOException{
		IO.pl ("Executing: '"+Misc.stringArrayToString(command, " ")+"'");
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
			if (printStandardOut) IO.pl(line);
		}
		data.close();
		String[] res = new String[al.size()];
		al.toArray(res);
		return res;
	}

	private void fetchFilesInRepo() throws IOException{
		IO.pl("Fetching file list from AWS for '"+awsPatientDirUri+ "'...");
		String[] cmd = {awsPath, "s3", "ls", "--recursive", awsPatientDirUri, "--profile", profile};
		String[] res = executeViaProcessBuilder(cmd, false, null);
		PrintWriter out = new PrintWriter(new FileWriter(awsRepoList));
		for (String s: res) {
			s= s.trim();
			//empty
			if (s.length() == 0) continue;
			//folder?
			if (s.startsWith("PRE ")) continue;
			String[] fields = Misc.WHITESPACE.split(s);
			//2023-03-01   09:04:43   2938    Patients/AA2mF6Vy/Avatar/A032049_SL419345_SL419548_SL420681/ClinicalReport/A032049_SL419345_SL419548_SL420681_IDTv1_SAR_F.json
			//   date         time     size       key
			out.println(Misc.stringArrayToString(fields, "\t"));
		}
		out.close();
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new PMRAvatarDedup(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException {
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': clinicalReportDir = new File(args[++i]).getCanonicalFile(); break;
					case 'x': awsPatientDirUri = args[++i]; break;
					case 'p': profile =args[++i]; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//root patient dirs?
		if (clinicalReportDir == null ) {
			Misc.printErrAndExit("Error: failed to find your directory for saving ClinicalReport json and xml files, see -c ? "+clinicalReportDir);
		}
		clinicalReportDir.mkdirs();

		//check credentials
		File awsCredentialsDir = new File (System.getProperty("user.home")+"/.aws");
		if (awsCredentialsDir.exists() == false)  Misc.printErrAndExit("Error: failed to find the aws credentials directory "+awsCredentialsDir);
		File credentialsFile = new File (awsCredentialsDir, "credentials").getCanonicalFile();
		if (awsCredentialsDir.exists() == false)  Misc.printErrAndExit("Error: failed to find your aws credentials file "+credentialsFile);

		//set env var for the aws cli 
		envPropToAdd.put("AWS_SHARED_CREDENTIALS_FILE", credentialsFile.getCanonicalPath());

		//Only needed in Eclipse

		if (true) {
			envPropToAdd.put("PATH", "/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin");
			awsPath="/usr/local/bin/aws";
		}
	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                     Patient Molecular Repo Avatar Dedup : May 2023               **\n" +
				"**************************************************************************************\n" +
				"Interactive searching of the clinical and sample attribute information in the json/xml\n"+
				"reports in the HCI PMR /ClinicalReport/ folders to identify datasets for analysis.\n"+
				"Press return after loading to see the search menu. Assumes the AWS CLI:\n"+
				"https://docs.aws.amazon.com/cli/ is installed in your path and you have read access\n"+
				"to the PMR. For searching for particular mutations, fusions, or cnvs, use GQuery:\n"+
				"https://github.com/HuntsmanCancerInstitute/GQuery\n"+

				"\nOptions:\n"+
				"-d  Directory to save the PHI redacted clinical xml and json reports.\n"+
				"-u  S3 URI containing the patient molecular repo, defaults to\n"+
				"      s3://hcibioinfo-patient-molecular-repo/Patients/ \n"+
				"-p  AWS credential profile, defaults to 'default'\n"+
				"-v  Verbose output.\n"+


				"\nExample: java -jar pathToUSeq/Apps/PMRSearch -d ~/PMRFiles/\n"+

				"**************************************************************************************\n");
	}
}
