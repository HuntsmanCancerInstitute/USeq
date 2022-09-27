package edu.utah.seq.run.avproj.adw;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.hci.bioinfo.smm.Subject;
import edu.utah.hci.bioinfo.smm.SubjectMatchMaker;
import util.gen.IO;
import util.gen.Misc;

public class AvatarDataWrangler {

	//user fields
	private File resourceDirectory = null;
	private File jobDir = null;
	private File dxCmdLineFile = null;

	//internal
	private File patientPhiFile = null;
	private File awsRepoListFile = null;
	private File linkageFile = null;
	private File qcMetricsFile = null;
	private File wesCramsFile = null;
	private File rnaCramsFile = null;
	private ClinicalMolLinkage linkage = null;
	private WesQCMetrics wes = null;
	private PersonIDConverter patientPhi = null;
	private HashMap<String, String> slIdWesCramName = new HashMap<String, String>();
	private HashMap<String, String> slIdRnaCramName = new HashMap<String, String>();
	private HashMap<String, PatientADW> patients = new HashMap<String, PatientADW>();
	private int numPatientsWithMissingData;
	private ArrayList<PatientADW> patientsWithJobsToRun = new ArrayList<PatientADW>();
	
	//subject match maker
	private File smmDirectory = null;
	private File smmRegistryDirectory = null;

	

	public AvatarDataWrangler (String[] args){
		long startTime = System.currentTimeMillis();
		try {
			processArgs(args);
			
			loadResourceFiles();

			loadM2GenFiles();
			
			loadAvailableCramFiles();

			buildPatients();
			
			makeAnalysisJobs();
			
			compareJobsToAlreadyProcessed();
			
			if (patientsWithJobsToRun.size()!=0) {
				
				//fetch the PHI
				if (loadOrRequestPhi()) {
					
					pullMolDataPatientIds();
					
					buildAnalysisJobDirs();
				}
			
				
				
			}
			else IO.pl("\nNo new jobs to run, exiting.");
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
	}



	private void buildAnalysisJobDirs() throws IOException {
		IO.pl("\nBuilding Job Directories...");
		
		ArrayList<String> dxCmds = new ArrayList<String>();
		
		//for each analysis job
		for (PatientADW p : patientsWithJobsToRun) {
			File jd = new File (jobDir, p.getSubjectMatchMaker().getCoreIdNewOrMatch()+"/Avatar/");
			jd.mkdirs();

			for (AvatarAnalysisJob aj: p.getAnalysisJobs()){
				String ajName = aj.getComparisonKey(p.getPatientId());
				File test = new File (jd, ajName);
				test.mkdir();
				//File testDir, ArrayList<String> dxCmds, HashMap<String, String> keyDiseaseType
				aj.makeAnalysisJob(ajName, test, dxCmds, linkage);
			}
			//write out the dxCmds
			if (IO.writeArrayList(dxCmds, dxCmdLineFile) == false) {
				dxCmdLineFile.delete();
				throw new IOException("\nFailed to write out the DNAnexus cmd line file. Deleting it and aborting.");
			}
			else dxCmdLineFile.setExecutable(true);
		}
		
	}



	private void pullMolDataPatientIds() throws IOException {
		IO.pl("\nMatching patient PHI for job directory generation...");
		
		//make IO for writing out queries
		File queryFile = new File (smmDirectory, "avatarPatientQueries_PHI.txt");
		File queryResDir = new File (smmDirectory, "SMMQueryResults_PHI");
		queryResDir.mkdirs();
		
		PrintWriter out = new PrintWriter (new FileWriter( queryFile));
		out.println("#lastName\tfirstName\tdobMonth(1-12)\tdobDay(1-31)\tdobYear(1900-2050)\tgender(M|F)\tmrn\tavatarId;thciId");
		
		//write out each for lookup
		HashMap<String, String[]> personIdDataLine = patientPhi.getPersonIDDataLine();
		for (PatientADW cp : patientsWithJobsToRun) {
			String[] atts = personIdDataLine.get(cp.getPatientId());
			out.println(patientPhi.fetchSubjectMatchMakerLine(atts));	
		}
		out.close();
		
		//run the SMM
		String[] args = {
				"-r", smmRegistryDirectory.getCanonicalPath(),
				"-q", queryFile.getCanonicalPath(),
				"-o", queryResDir.getCanonicalPath(),
				"-v", //turn off verbosity
				"-c", //making name matching case-insensitive for EDW
				"-a", //adding queries not found to registry
				"-u" //fill in info not found to registry entries from queries, e.g. missing mrns, otherIds
		};
		
		SubjectMatchMaker smm = new SubjectMatchMaker(args);
		//in same order as that which was printed out above from workingPatients
		Subject[] queries = smm.getQuerySubjects();
		for (int i=0; i< queries.length; i++) {
			patientsWithJobsToRun.get(i).setSubject(queries[i]);
		}
		
	}



	private void compareJobsToAlreadyProcessed() throws IOException {
		IO.pl("\nExcluding Analysis Jobs found in the AWS repo list....");
		
		//pull contents of the aws repo 
		AWSRepoList arl = new AWSRepoList(awsRepoListFile);
		HashSet<String> avatarJobDirNames = arl.parseSubDirNamesAfterParent("Avatar");
		IO.pl("\t"+ avatarJobDirNames.size()+ "\t# Analysis Jobs found on AWS");
		
		//Create a list of current jobs that should be run
		int numExcluded = 0;
		int numKept = 0;
		for (PatientADW p: patients.values()) {
			ArrayList<AvatarAnalysisJob> jobsToKeep = new ArrayList<AvatarAnalysisJob>();
			for (AvatarAnalysisJob j: p.getAnalysisJobs()) {
				String key = j.getComparisonKey(p.getPatientId());
				if (avatarJobDirNames.contains(key)) numExcluded++;
				else {
					jobsToKeep.add(j);
					numKept++;
				}
			}
			p.setAnalysisJobs(jobsToKeep);
			if (jobsToKeep.size()>0) patientsWithJobsToRun.add(p);
		}
		//assign filtered patients
		
		IO.pl("\t"+ numExcluded+ "\t# Analysis Jobs identified here that are already uploaded to AWS and are thus skipped.");
		IO.pl("\t"+ numKept+ "\t# New Analysis Jobs to run");
	}

	private boolean loadOrRequestPhi() throws IOException {
		IO.pl("\nLooking for and parsing PHI...");
		if (patientPhiFile == null) {
			StringBuilder sb = new StringBuilder(patientsWithJobsToRun.get(0).getPatientId());
			for (int i=1; i< patientsWithJobsToRun.size(); i++) {
				sb.append(", ");
				sb.append(patientsWithJobsToRun.get(i).getPatientId());
			}
			IO.pl("\tPlease pull PHI for the following ORIENAvatarPatientIDs from RISR's PersonID converter. Paste the following into \n"+
			"\thttps://hci-rdbs.hci.utah.edu/Reports/report/TCC/Get%20PersonID%20By%20MRN/Get%20PersonID%20By%20MRN , set the\n"+
			"\tInput Type to AvatarID, click View Report, and export/ download it as a csv file into the Resources folder naming it xxx_PatientPHI.csv\n"+
			"\tLastly, rerun this tool. AvatarIDs to enter: "+sb);
			return false;
		}
		
		else {
			patientPhi = new PersonIDConverter(patientPhiFile);
			if (patientPhi.isParsed()==false) throw new IOException("Failed to parse the patient PHI lookup file.");
			
			//check that all of the Patients have phi
			HashMap<String, String[]> avatarIDDataLine = patientPhi.getPersonIDDataLine();
			ArrayList<String> missingPHIForAvatarID = new ArrayList<String>();
			for (PatientADW patient: patientsWithJobsToRun) {
				String avatarId = patient.getPatientId();
				if (avatarIDDataLine.containsKey(avatarId) == false) missingPHIForAvatarID.add(avatarId);
			}
			if (missingPHIForAvatarID.size()!=0) {
				IO.el("Failed to find PHI for the following ORIENAvatarIds from "+ patientPhiFile +" correct and restart.\n" + missingPHIForAvatarID);
				return false;
			}
			else return true;
		}
	}

		
	private void loadResourceFiles() {
		
		IO.pl("Looking for resource files....");
		StringBuilder issues = new StringBuilder();
		

		linkageFile = fetchFromResource("ClinicalMolLinkage_V4.csv",issues);
		qcMetricsFile = fetchFromResource("WES_QC_Metrics.csv",issues);
		wesCramsFile = fetchFromResource("WesCramList.txt",issues);
		rnaCramsFile = fetchFromResource("RNACramList.txt",issues);
		awsRepoListFile = fetchFromResource("AWSRepoList.txt",issues);
		
		if (issues.length() !=0) Misc.printErrAndExit("One or more resource files are missing:\n"+issues);
				
		//assign the patient phi, this might be missing!
		patientPhiFile = fetchFromResource("PatientPHI.csv",issues);
	}
		
	
	private File fetchFromResource(String extension, StringBuilder issues) {
		File[] f = IO.extractFiles(resourceDirectory, extension);
		if (f.length==1) return f[0];
		issues.append("\tError: failed to find one xxx_"+extension+" file in "+resourceDirectory);
		return null;
	}

	private void loadAvailableCramFiles() throws IOException {
		IO.pl("\nLoading available cram files...");
		// FT-SA134847_st_g_markdup.cram 
		//    SL261633_st_g_markdup.cram 
		//    SL261681_st_t_markdup.cram 
		//      A59553_st_t_markdup.cram 
		//      A59554_st_g_markdup.cram
		String[] wesCrams = IO.loadFileIntoStringArray(wesCramsFile);
		for (String wc: wesCrams) {
			String[] f = Misc.UNDERSCORE.split(wc);
			if (f.length != 4) throw new IOException ("Failed to find 4 parts in this wes cram file entry "+wc);
			if (slIdWesCramName.containsKey(f[0])) throw new IOException ("Found duplicate SLID "+f[0]);
			slIdWesCramName.put(f[0], wc);
		}

		// FT-SA130920R.genome.cram 
		//     SL316725.genome.cram
		String[] rnaCrams = IO.loadFileIntoStringArray(rnaCramsFile);
		for (String wc: rnaCrams) {
			String[] f = Misc.PERIOD.split(wc);
			if (f.length != 3) throw new IOException ("Failed to find 3 parts in this RNA cram file entry "+wc);
			if (slIdRnaCramName.containsKey(f[0])) throw new IOException ("Found duplicate SLID "+f[0]);
			slIdRnaCramName.put(f[0], wc);
		}
	}

	private void makeAnalysisJobs() throws IOException {
		IO.pl("\nBuilding analysis jobs for...");
		int numAnalysisJobs = 0;
		
		for (PatientADW p: patients.values()) {
			buildAnalysisJobs(p);
			numAnalysisJobs+= p.getAnalysisJobs().size();
		}
		IO.pl("\t"+ numAnalysisJobs+"\t# Analysis jobs ready for processing");
	}

	private void buildAnalysisJobs(PatientADW p) throws IOException {

		IO.pl(p.getPatientId());

		//split normal samples by platform
		HashMap<String, ArrayList<NormalSampleADW>> platformNormalSamples = splitNormalSamplesByPlatform(p);

		//split tumor samples by trimmed generic specimineId
		HashMap<String, ArrayList<TumorSampleADW>> specimineIdTumorSamples = splitTumorSamplesBySpecimine(p.getTumorSamples());

		//for each tumor specimine
		for (String specimineId: specimineIdTumorSamples.keySet()) {

			/////////////
			//Define the Tumor exome and rna samples
			ArrayList<TumorSampleADW> tumorSamples = specimineIdTumorSamples.get(specimineId);

			//attempt to merge those with separate RNA and Exome ids
			tumorSamples = mergeSplitTumorExomeRNADatasets(tumorSamples);
			int numTumSamp = tumorSamples.size();

			//more than one?
			if (numTumSamp !=1) throw new IOException("\nFailed to merge tumor samples with the same trimmed specimineId for patient "+p.getPatientId());
			TumorSampleADW tumorSample = tumorSamples.get(0);

			////////////
			//Define the Normal? do so only for Tumor samples with a exome
			ArrayList<NormalSampleADW> normalSamplesToAdd = new ArrayList<NormalSampleADW>();
			boolean matchedPlatform = true;
			
			if (tumorSample.getTumorDnaName() != null) {
				ArrayList<NormalSampleADW> normalSamples = platformNormalSamples.get(tumorSample.getPlatformName());
				int numNorm = 0;
				if (normalSamples != null) numNorm = normalSamples.size();


				//more than one normal? 
				if (numNorm > 1) {
					IO.pl("\tWARNING: multiple normal files found in the same platform, these will all be added to the Fastq dir for merging: ");
					for (NormalSampleADW ns: normalSamples) {
						IO.pl("\t\t"+ns.getNormalWesCramFileNameToFetch());
						normalSamplesToAdd.add(ns);
					}
				}

				//one normal
				else if (numNorm == 1) normalSamplesToAdd.add(normalSamples.get(0));

				//for zero platform normals, will try to link in another platform one below
				else {
					//any normals?
					if (p.getNormalSamples().size() == 0) {
						IO.pl("\tWARNING: no normal found.");
					}
					else {
						IO.pl("\tWARNING: no normal found in the same platform, will add those from all other platforms for merging:");
						for (NormalSampleADW ns: p.getNormalSamples()) {
							normalSamplesToAdd.add(ns);
							IO.pl("\t\t"+ns.getPlatformName()+"\t"+ ns.getNormalWesCramFileNameToFetch());
							matchedPlatform = false;
						}
					}
				}
			}
			//add the job to the patient
			AvatarAnalysisJob adJob = new AvatarAnalysisJob(p, tumorSample, normalSamplesToAdd, matchedPlatform);
			p.getAnalysisJobs().add(adJob);

		}

		//just normal samples?
		if (p.getTumorSamples().size() == 0 && p.getNormalSamples().size()!=0) {
			IO.pl("\tWARNING: no tumor samples found, just normal");
			HashMap<String, ArrayList<NormalSampleADW>> platNorm = splitNormalSamplesByPlatform(p);
			boolean matchedPlatform = false;
			if (platNorm.size() == 1) matchedPlatform = true;
			AvatarAnalysisJob adJob = new AvatarAnalysisJob(p, null, p.getNormalSamples(), matchedPlatform);
			p.getAnalysisJobs().add(adJob);
		}
	}
	

	
	/**Attempts to merge tumor samples with the same trimmed specimine id.  These are all from the Heme group where the 
	 * tumor exome and tumor rna are on different clin link data lines */
	private ArrayList<TumorSampleADW> mergeSplitTumorExomeRNADatasets(ArrayList<TumorSampleADW> tumorSamples) {
		int num = tumorSamples.size();

		if (num == 2) {
			TumorSampleADW exomeOnly = null;
			TumorSampleADW rnaOnly = null;
			for (TumorSampleADW ts: tumorSamples) {
				if (ts.getTumorDnaName() != null && ts.getTumorRnaName() == null) exomeOnly = ts;
				else if (ts.getTumorDnaName() == null && ts.getTumorRnaName() != null) rnaOnly = ts;
			}
			//both found?
			if (exomeOnly != null && rnaOnly != null) {
				exomeOnly.setTumorRnaName(rnaOnly.getTumorRnaName());
				exomeOnly.setTumorRnaCramFileNameToFetch(rnaOnly.getTumorRnaCramFileNameToFetch());
				exomeOnly.getTumorLinkageDataLines().addAll(rnaOnly.getTumorLinkageDataLines());
				tumorSamples.clear();
				tumorSamples.add(exomeOnly);
			}
		}
		return tumorSamples;
	}

	private HashMap<String, ArrayList<NormalSampleADW>> splitNormalSamplesByPlatform(PatientADW p) {
		HashMap<String, ArrayList<NormalSampleADW>> platformNormalSamples = new HashMap<String, ArrayList<NormalSampleADW>>();
		for (NormalSampleADW ns: p.getNormalSamples()) {
			String platform = ns.getPlatformName();
			if (platform == null) platform= "NA";
			ArrayList<NormalSampleADW> al = platformNormalSamples.get(platform);
			if (al == null) {
				al = new ArrayList<NormalSampleADW>();
				platformNormalSamples.put(platform, al);
			}
			al.add(ns);
		}
		return platformNormalSamples;
	}

	private HashMap<String, ArrayList<TumorSampleADW>> splitTumorSamplesByPlatform(PatientADW p) {
		HashMap<String, ArrayList<TumorSampleADW>> platformTumorSamples = new HashMap<String, ArrayList<TumorSampleADW>>();
		for (TumorSampleADW ts: p.getTumorSamples()) {
			String platform = ts.getPlatformName();
			if (platform == null) platform= "NA";
			ArrayList<TumorSampleADW> al = platformTumorSamples.get(platform);
			if (al == null) {
				al = new ArrayList<TumorSampleADW>();
				platformTumorSamples.put(platform, al);
			}
			al.add(ts);
		}
		return platformTumorSamples;
	}
	
	private HashMap<String, ArrayList<TumorSampleADW>> splitTumorSamplesBySpecimine(ArrayList<TumorSampleADW> tss) {
		HashMap<String, ArrayList<TumorSampleADW>> specimineTumorSamples = new HashMap<String, ArrayList<TumorSampleADW>>();
		for (TumorSampleADW ts: tss) {
			String specimine = ts.getGenericSpecimineId();
			ArrayList<TumorSampleADW> al = specimineTumorSamples.get(specimine);
			if (al == null) {
				al = new ArrayList<TumorSampleADW>();
				specimineTumorSamples.put(specimine, al);
			}
			al.add(ts);
		}
		return specimineTumorSamples;
	}

	private void buildPatients() throws IOException {
		IO.pl("\nCollecting patient dataset info...");
		HashMap<String, Integer> map = linkage.getHeaderKeyIndex();

		Integer typeIndex = map.get("Tumor/Germline");
		Integer wesIndex = map.get("WES");
		Integer rnaIndex = map.get("RNASeq");
		Integer specimineIndex = map.get("ORIENSpecimenID");

		if (typeIndex == null || wesIndex == null || rnaIndex == null || specimineIndex == null) {
			Misc.printErrAndExit("\nERROR: failed to "
					+ "identify one or more of the following linkage file header columns, check that these are in your file: "
					+ "Tumor/Germline, WES, RNASeq, ORIENSpecimenID");
		}
		
		//fetch the platform and datalines
		HashMap<String, ArrayList<String[]>> idDataLines = linkage.getAvatarIdDataLines();
		HashMap<String, String> slidPlatform = wes.getSlidPlatform();

		//for each patient id
		for (String patientId: idDataLines.keySet()) {
			
			//fetch the patient
			PatientADW p = new PatientADW(patientId);
			boolean OK = true;
			
			//for each patient associated dataline
			for (String[] fields: idDataLines.get(patientId)) {
				
				//parse the specimine and trim off the trailing .xx    16-0061475.1b 17-0001137c
				String specimineId = fields[specimineIndex];
				if (specimineId.length()==0) throw new IOException("Failed to find a ORIENSpecimenID for one of the clinical linkage file data lines associated with patient "+patientId);
				int periodIndex = specimineId.indexOf(".");
				if (periodIndex!= -1) specimineId = specimineId.substring(0, periodIndex);
				
				String tumorGermline = fields[typeIndex];
				String wesId = fields[wesIndex];
				String rnaId = fields[rnaIndex];
				if (wesId.length()==0) wesId= null;
				if (rnaId.length()==0) rnaId= null;

				//pull platform, will be null for rna only samples
				String platform = null;
				if (wesId!=null) {
					platform = slidPlatform.get(wesId);
					if (platform == null) {
						Misc.printErrAndExit("\nERROR: failed to find a platform type for "+tumorGermline+" sample "+
								wesId+" for patient ID "+patientId+" Aborting.");
					}
				}
				
				
				//is it a tumor dataline? could be tumor exome, tumor rna, or both
				if (tumorGermline.equals("Tumor")) {
					TumorSampleADW ts = new TumorSampleADW(wesId, rnaId, platform, specimineId, fields);
					
					//any tumor wes?
					if (wesId != null) {
						String nameTumorWes = slIdWesCramName.get(wesId);
						if (nameTumorWes == null) {
							OK = false;
							IO.pl("\tFailed to find an available tumor WES cram file for "+ wesId+" in the M2Gen project, skipping patient "+patientId);
						}
						else ts.setTumorWesCramFileNameToFetch(nameTumorWes);
					}
					
					//any tumor rna
					if (rnaId != null) {
						String nameTumorRna = slIdRnaCramName.get(rnaId);
						if (nameTumorRna == null) {
							OK = false;
							IO.pl("\tFailed to find an available tumor RNA cram file for "+ rnaId+" in the M2Gen project, skipping patient "+patientId);
						}
						else ts.setTumorRnaCramFileNameToFetch(nameTumorRna);
					}
					p.getTumorSamples().add(ts);
				}
				
				//is it a germline normal dataline
				else if (tumorGermline.equals("Germline")) {
					String nameNormalWes = slIdWesCramName.get(wesId);
					if (nameNormalWes == null) {
						OK = false;
						IO.pl("\tFailed to find an available normal WES cram file for "+ wesId+" in the M2Gen project, skipping patient "+patientId);
					}
					else {
						NormalSampleADW ns = new NormalSampleADW(wesId, platform, fields);
						ns.setNormalWesCramFileNameToFetch(nameNormalWes);
						p.getNormalSamples().add(ns);
					}
					
				}
				else {
					Misc.printErrAndExit("\nERROR: didn't find 'Tumor' or 'Germline' in the type column for "+Misc.stringArrayToString(fields, " "));
				}
			}
			if (OK) patients.put(patientId, p);
			else numPatientsWithMissingData++;
		}
		IO.pl("\t"+ patients.size()+"\t# Patients ready for processing");
		IO.pl("\t"+ numPatientsWithMissingData+"\t# Patients missing cram files, see above.");
	}
	
	private void loadM2GenFiles() throws IOException {
		IO.pl("\nLoading M2GenFiles...");
		linkage = new ClinicalMolLinkage(linkageFile);
		if (linkage.isParsed() == false) throw new IOException("Problem parsing "+linkageFile);
		wes = new WesQCMetrics(qcMetricsFile);
		if (wes.isParsed() == false) throw new IOException("Problem parsing "+qcMetricsFile);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AvatarDataWrangler(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		Pattern pat = Pattern.compile("-[a-z]");
		File tmpDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'r': resourceDirectory = new File(args[++i]); break;
					case 'j': jobDir = new File(args[++i]); break;
					case 'p': printScript(); break;
					case 't': tmpDir = new File(args[++i]); break;
					case 'd': dxCmdLineFile = new File(args[++i]); break;
					case 's': smmRegistryDirectory = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//dx cmd file
		if (dxCmdLineFile == null || dxCmdLineFile.getParentFile().canWrite() == false) Misc.printErrAndExit("\nERROR: please provide a file path to save the DNAnexus download cmds. Aborting. "+dxCmdLineFile);
		
		
		//jobDir
		if (jobDir == null) Misc.printErrAndExit("\nERROR: cannot find your JobDir? Aborting. "+jobDir);
		jobDir.mkdirs();
		if (jobDir.exists() == false || jobDir.canWrite() == false) Misc.printErrAndExit("\nERROR: cannot write into your JobDir? Aborting. "+jobDir);
		
		//resourceDir
		if (resourceDirectory == null || resourceDirectory.exists()== false) Misc.printErrAndExit("\nERROR: cannot find your Resource directory? Aborting. "+resourceDirectory);
	
		//check smm registry directory
		if (smmRegistryDirectory == null || smmRegistryDirectory.exists() == false || smmRegistryDirectory.isDirectory()==false) Misc.printErrAndExit("Error: cannot find your Subject Match Maker registry directory, "+smmRegistryDirectory);
		
		//check tmpDir for the subject match maker
		if (tmpDir == null) Misc.printErrAndExit("Error: cannot find your temp subject ID directory "+tmpDir);
		tmpDir.mkdirs();
		if (tmpDir.exists() == false) Misc.printErrAndExit("Error: cannot make your temp subject ID directory "+tmpDir);
		smmDirectory = new File (tmpDir, System.currentTimeMillis()+"_DeleteMe");
		smmDirectory.mkdir();
	}


	private void printScript() {
		StringBuilder sb = new StringBuilder();
		sb.append("# Log in to the DNAnexus API\n");
		sb.append("t=$(cat ~/BioApps/DNAnexus/token) \n");
		sb.append("dx login --token $t \n");
		sb.append("\n");
		sb.append("# Download the linkage and metrics files \n");
		sb.append("mkdir -f ResourceFiles\n");
		sb.append("dx download HCI_Molecular_Data:/Avatar_MolecularData_hg38/Manifests_and_QC_Files/*_HCI_ClinicalMolLinkage_V4.csv -o ResourceFiles/\n");
		sb.append("dx download HCI_Molecular_Data:/Avatar_MolecularData_hg38/Manifests_and_QC_Files/*_HCI_WES_QC_Metrics.csv -o ResourceFiles/\n");
		sb.append("\n");
		sb.append("# Pull list of available WES and RNASeq crams \n");
		sb.append("d=$(date '+%d%b%Y')\n");
		sb.append("dx ls HCI_ORIEN_AVATAR_MOLECULAR_DATA:/Whole_Exome/alignment_crams/*cram | sort | uniq > ResourceFiles/$d'_WesCramList.txt'\n");
		sb.append("dx ls HCI_ORIEN_AVATAR_MOLECULAR_DATA:/RNAseq/alignment_crams/*cram | sort | uniq > ResourceFiles/$d'_RNACramList.txt'\n");
		sb.append("\n");
		sb.append("# Pull a list of all of the Avatar associated files in the AWS patient repo\n");
		sb.append("aws s3 ls s3://hcibioinfo-patient-molecular-repo/Patients/ --recursive | grep Avatar > ResourceFiles/$d'_AWSRepoList.txt'\n");
		sb.append("\n");
		sb.append("# Pull PHI for the indicated AvatarIds from running this app\n");
		sb.append("# Log in to https://hci-rdbs.hci.utah.edu/Reports/report/TCC/Get%20PersonID%20By%20MRN/Get%20PersonID%20By%20MRN\n");
		sb.append("# Change the Input Type to Avatar ID, copy paste the AvatarIDs, click View Report, \n");
		sb.append("# Export the report by selecting CSV, rename this xxx_PatientPHI.csv and place it in the Resource directory\n");
		Misc.printExit(sb.toString());
	}



	public static void printDocs(){

		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                           Avatar Data Wrangler : Sept 2022                       **\n" +
				"**************************************************************************************\n" +
				"Tool for assembling directories for TNRunner based on files provided by M2Gen via\n"+
				"download from DNAnexus, see\n"+
				"HCI_Molecular_Data:/Avatar_MolecularData_hg38/Manifests_and_QC_Files/ Handles patient\n"+
				"datasets from different exome capture platforms with multiple tumor samples.\n"+
				"Each TNRunner Job will be named as \n"+
				"MolDataPatientID-Platform-NormalExomeID-TumorExomeID-TumorRNASeqID, when missing they will be\n"+
				"labeled NA.\n"+
				"Only assembles dirs for processing where the data is complete.\n"+

				"\nRequired Parameters:\n"+
				"-r Path to a folder containing the following resource files.  Each must end with the\n"+
				"     given extensions:\n"+
				"        xxx_ClinicalMolLinkage_V4.csv - DNAnexus\n"+
				"        xxx_WES_QC_Metrics.csv - DNAnexus\n"+
				"        xxx_WesCramList.txt - DNAnexus\n"+
				"        xxx_RNACramList.txt - DNAnexus\n"+
				"        xxx_AWSRepoList.txt - AWS\n"+
				"        xxx_PatientPHI.csv - RISR's PatientIdCoverter (optional)\n"+
				"-p Print an example script to pull these resource files.\n"+
				"-j Job dir to build out the patient dir structure.\n"+
				"-t Directory to place temp files with PHI for Subject ID matching\n"+
				"-s Directory containing the SubjectMatchMaker 'currentRegistry_' file\n"+
				"-d Path to save a bash script for downloading the sequence read data.\n"+

                "\nExample: java -jar ~/USeqApps/AvatarDataWrangler ... \n\n"+


				"**************************************************************************************\n");
	}
}
