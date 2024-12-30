package edu.utah.seq.run.avproj.adw;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
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
	//private File cmdLineFile = null;
	private File asterDownloadDir = null;
	private String awsRepo = "s3://hcibioinfo-patient-molecular-repo/";

	//internal
	private File patientPhiFile = null;
	private File awsRepoListFile = null;
	private File linkageFile = null;
	private File qcMetricsFile = null;
	private File fastqFiles = null;
	private ClinicalMolLinkage linkage = null;
	private WesQCMetrics wes = null;
	private PersonIDConverter patientPhi = null;
	private HashMap<String, LinkedHashSet<String>> slIdFastqPaths = new HashMap<String, LinkedHashSet<String>>();
	private HashMap<String, LinkedHashSet<String>> slIdAwsPaths = new HashMap<String, LinkedHashSet<String>>();
	private HashMap<String, PatientADW> patients = new HashMap<String, PatientADW>();
	private int numPatientsWithMissingData;
	private ArrayList<PatientADW> patientsWithJobsToRun = new ArrayList<PatientADW>();
	AWSRepoList awsRepoList = null;;
	
	//Aster downloads, a moving target at present
	private String asterExecDir = "/uufs/chpc.utah.edu/common/HIPAA/u0028003/BioApps/Aster";
	private String asterProjectId = "project-F66v00Q0q4045Q4Y6PY2Xv7F"; 
	
	//subject match maker
	private File smmDirectory = null;
	private File smmRegistryDirectory = null;

	public AvatarDataWrangler (String[] args){
		long startTime = System.currentTimeMillis();
		try {
			processArgs(args);
			
			loadResourceFiles();

			loadM2GenFiles();
			
			loadAvailableFastqFiles();
			loadAvailableFastqFilesFromAws();

			buildPatients();
			
			makeAnalysisJobs();
			
			compareJobsToAlreadyProcessed();
			
			if (patientsWithJobsToRun.size()!=0) {
				//fetch the PHI
				if (loadOrRequestPhi()) {
					pullMolDataPatientIds();
					buildAnalysisJobDirsNew();
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

	private void buildAnalysisJobDirsNew() throws IOException {
		IO.pl("\nBuilding Job Directories...");
		
		ArrayList<String> linkCommandsToExecute = new ArrayList<String>();
		linkCommandsToExecute.add("set -e");
		
		LinkedHashSet<String> filesToDownload = new LinkedHashSet<String>();
		ArrayList<String> normalsToMerge = new ArrayList<String>();
		
		//for each analysis job
		for (PatientADW p : patientsWithJobsToRun) {
			File jd = new File (jobDir, p.getSubjectMatchMaker().getCoreIdNewOrMatch()+"/Avatar/");
			jd.mkdirs();
			for (AvatarAnalysisJob aj: p.getAnalysisJobs()){
				String ajName = aj.getComparisonKey(p.getPatientId());
				File test = new File (jd, ajName);
				test.mkdir();
				//File testDir, ArrayList<String> dxCmds, HashMap<String, String> keyDiseaseType
				aj.makeAnalysisJobAsterNew(asterDownloadDir, ajName, test, linkCommandsToExecute, linkage, filesToDownload);
				//normals to merge?
				if (aj.getNormalSamples().size()>1) normalsToMerge.add(aj.getNormalDNAFastqDir().toString());
			}
		}
		
		//write out ln commands
		File linkCommands = new File (resourceDirectory, "makeSoftLinks.sh");
		linkCommandsToExecute.add("echo COMPLETE");
		IO.writeArrayList(linkCommandsToExecute, linkCommands);
		linkCommands.setExecutable(true);
		
		String fullDownloadDirPath = asterDownloadDir.getCanonicalPath();
		
		//download commands, these are either from Aster or from AWS
		ArrayList<String> asterCmds = new ArrayList<String>();
		ArrayList<String> mvCmds = new ArrayList<String>();
		mvCmds.add("set -e");
		ArrayList<String> s3CopyPaths = new ArrayList<String>();
		int numAsterToDownload = 0;
		int numAwsToDownload = 0;
		
		for (String path: filesToDownload) {
			int lastSlashIndex = path.lastIndexOf('/');
			String fileName = path.substring(lastSlashIndex+1);
			String justPath = path.substring(0, lastSlashIndex);
			
			if (path.startsWith("/Avatar_Molecular")) {
				//add on a single line so can use parallel
				//aster scripts throw errors that are just bs so cannot test for completion
				asterCmds.add("echo STARTING_"+fileName+
						" && module load python3 && python3 "+asterExecDir+"/support_scripts/download_project.py"+
						" --no-dry-run --workers 1"+
						" --project-id "+ asterProjectId+
						" --exec "+asterExecDir+"/rwb.linux.x64"+
						" --include "+justPath+
						" --file-filter "+fileName+
						" --destination-path "+ fullDownloadDirPath+
						" &> /dev/null && echo OK_"+fileName+" || echo FAILED_"+fileName);
				mvCmds.add("mv "+fullDownloadDirPath+ path+ " "+fullDownloadDirPath+"/" );
				numAsterToDownload++;
			}
			else if (path.startsWith("Patients/")) {
				s3CopyPaths.add(awsRepo+path+" > "+asterDownloadDir+"/");
				numAwsToDownload++;
			}
			else Misc.printErrAndExit("\nERROR finding source for seq data download? "+path);
		}
		
		//write out Aster cmds
		if (numAsterToDownload !=0) {
			File asterFile = new File (resourceDirectory, "downloadAsterFiles.sh");
			IO.writeArrayList(asterCmds, asterFile);
			asterFile.setExecutable(true);
			
			File mvFile = new File (resourceDirectory, "moveAsterFiles.sh");
			mvCmds.add("echo COMPLETE");
			IO.writeArrayList(mvCmds, mvFile);
			mvFile.setExecutable(true);
		}

		
		//write out S3Copy cmds
		if (numAwsToDownload !=0) {
			File toCopyS3 = new File (resourceDirectory, "awsFilesToDownload.txt");
			IO.writeArrayList(s3CopyPaths, toCopyS3);
		}
		
		//any multiple normals?
		if (normalsToMerge.size()!=0) {
			File f = new File (resourceDirectory, "normalsToMerge.txt");
			IO.writeArrayList(normalsToMerge, f);
		}
		
		IO.pl("\tAster downloads: "+numAsterToDownload+"\tAws downloads: "+numAwsToDownload+"\tNormals to merge: "+normalsToMerge.size());
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
				"-c", //making name matching case-insensitive for EDW
				"-a", //adding queries not found to registry
				"-u" //fill in info not found to registry entries from queries, e.g. missing mrns, otherIds
				//"-v", //turn off verbosity
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
		HashSet<String> avatarJobDirNames = awsRepoList.parseSubDirNamesAfterParent("Avatar");
		IO.pl("\t"+ avatarJobDirNames.size()+ "\t# Analysis Jobs found on AWS");
		
		//Create a list of current jobs that should be run
		int numExcluded = 0;
		int numKept = 0;
		for (PatientADW p: patients.values()) {
			ArrayList<AvatarAnalysisJob> jobsToKeep = new ArrayList<AvatarAnalysisJob>();
			for (AvatarAnalysisJob j: p.getAnalysisJobs()) {
				String key = j.getComparisonKey(p.getPatientId());
				//doing an exact comparision
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

		
	private void loadResourceFiles() throws IOException {
		
		IO.pl("Looking for resource files....");
		StringBuilder issues = new StringBuilder();
		

		linkageFile = fetchFromResource("ClinicalMolLinkage_V4.csv",issues);
		qcMetricsFile = fetchFromResource("WES_QC_Metrics.csv",issues);
		fastqFiles = fetchFromResource("fastqFiles.txt",issues);
		awsRepoListFile = fetchFromResource("AWSRepoList.txt",issues);
		
		if (issues.length() !=0) Misc.printErrAndExit("One or more resource files are missing:\n"+issues);
				
		//assign the patient phi, this might be missing!
		patientPhiFile = fetchFromResource("PatientPHI.csv",issues);
		
		//load the aws objects
		awsRepoList = new AWSRepoList(awsRepoListFile);
	}
		
	
	private File fetchFromResource(String extension, StringBuilder issues) {
		File[] f = IO.extractFiles(resourceDirectory, extension);
		if (f.length==1) return f[0];
		issues.append("\tError: failed to find one xxx_"+extension+" file in "+resourceDirectory);
		return null;
	}
	
	private void loadAvailableFastqFiles() throws IOException {
		IO.pl("\nLoading available fastq files...");
		/* 
		/Avatar_MolecularData_hg38/2023_06_30/Whole_Exome/FASTq/FT-SA212052_R1.fastq.gz
		/Avatar_MolecularData_hg38/2023_06_30/Whole_Exome/FASTq/FT-SA212052_R2.fastq.gz
		/Avatar_MolecularData_hg38/2023_06_30/Whole_Exome/FASTq/SL600460_1.fastq.gz
		/Avatar_MolecularData_hg38/2023_06_30/Whole_Exome/FASTq/SL600460_2.fastq.gz
		/Avatar_MolecularData_hg38/2023_03_03/RNAseq/FASTq/FT-SA168005R_R1.fastq.gz
		/Avatar_MolecularData_hg38/2023_03_03/RNAseq/FASTq/FT-SA168005R_R2.fastq.gz
		/Avatar_MolecularData_hg38/2023_03_03/RNAseq/FASTq/SL526167_1.fastq.gz
		/Avatar_MolecularData_hg38/2023_03_03/RNAseq/FASTq/SL526167_2.fastq.gz
		 * /Avatar_MolecularData_hg38/2019_10_30/Whole_Exome/FASTq/18-0026130a_C046_0025_010931_PB_Whole_C1_K1ID2_A62010_R1.fastq.gz
		*/
		String[] wesFastq = IO.loadFileIntoStringArray(fastqFiles);
		for (String wc: wesFastq) {
			//  /Avatar_MolecularData_hg38/2023_06_30/Whole_Exome/FASTq/FT-SA212052_R1.fastq.gz
			//  /           1                   2          3        4        5
			//  /Avatar_MolecularData_hg38/2023_03_03/RNAseq/FASTq/SL526167_2.fastq.gz
			//  /          1                 2          3      4      5
			String[] l = Misc.FORWARD_SLASH.split(wc);

			
			String id = null;
			if (l[5].contains("_K1ID2_")) {
				String[] uSplit = Misc.UNDERSCORE.split(l[5]);
				//18-0026130a_C046_0025_010931_PB_Whole_C1_K1ID2_A62010_R1.fastq.gz
				//     0        1    2     3    4   5    6   7     8      9
				id = uSplit[8];
			}
			else {
				int indexLastUnderscore = l[5].lastIndexOf("_");
				id = l[5].substring(0, indexLastUnderscore);
			}
			LinkedHashSet<String> al = slIdFastqPaths.get(id);
			if (al == null) {
				al = new LinkedHashSet<String>();
					slIdFastqPaths.put(id, al);
			}
			al.add(wc);
		}
		//iterate through each keeping last pair, the most recent
		for (String id: slIdFastqPaths.keySet()) {
			LinkedHashSet<String> al = slIdFastqPaths.get(id);
			if (al.size()==2) {}
			else if (al.size()==4) {
				LinkedHashSet<String> justTwo = new LinkedHashSet<String>();
				Iterator<String> it = al.iterator();
				for (int i=0; i< 4; i++) {
					String f = it.next();
					if (i>1) justTwo.add(f);
				}
				slIdFastqPaths.put(id, justTwo);
			}
			else {
				Misc.printErrAndExit("ERROR with parsing two or four fastqs from the available Aster datasets for "+id);
			}
		}
	}
	
	private void loadAvailableFastqFilesFromAws() throws IOException {
		IO.pl("\nLoading available fastq files from the AWS HCI repo...");
		/*
		 Patients/znw5ZG9RJv/Avatar/A041582_FT-SA178581_FT-SA182621D_FT-SA182621R/Fastq/TumorDNA/FT-SA182621D_st_t_markdup.cram
  		 Patients/znw5ZG9RJv/Avatar/A041582_FT-SA178581_FT-SA182621D_FT-SA182621R/Fastq/TumorRNA/FT-SA182621R.genome.cram
  		 Patients/yj2jn3QQ/Avatar/A017808_SL531287_SL402396_SL405412/Fastq/NormalDNA/SL531287_2.fastq.gz
        Patients/yj2jn3QQ/Avatar/A017808_SL531287_SL402396_SL405412/Fastq/TumorDNA/SL402396_1.fastq.gz
        Patients/zY7KU9eD/Avatar/A030398_A79884_A79832_NA/Fastq/TumorDNA/18-0030048a_C046_0055_013056_BL_Whole_T1_K1ID2_A79832_R2.fastq.gz
           0         1       2              3               4      5                          6
		 */
		
		for (String wc: awsRepoList.getAwsObjectNames()) {
			if (wc.contains("/Fastq/") == false || wc.contains("/Avatar/") == false) continue;
			
			String[] l = Misc.FORWARD_SLASH.split(wc);
			String id = null;
			if (l[6].contains("_K1ID2_")) {
				String[] uSplit = Misc.UNDERSCORE.split(l[6]);
				//18-0026130a_C046_0025_010931_PB_Whole_C1_K1ID2_A62010_R1.fastq.gz
				//     0        1    2     3    4   5    6   7     8      9
				id = uSplit[8];
			}
			else {
				int index = l[6].indexOf("_");
				if (index == -1) index = l[6].indexOf(".");
				
				id = l[6].substring(0, index);
			}
			LinkedHashSet<String> al = slIdAwsPaths.get(id);
			if (al == null) {
				al = new LinkedHashSet<String>();
				slIdAwsPaths.put(id, al);
			}
			al.add(wc);
		}
		//for each id, keep the first cram, or two fastqs not both
		for (String id: slIdAwsPaths.keySet()) {
			LinkedHashSet<String> al = slIdAwsPaths.get(id);
			ArrayList<String> fastq = new ArrayList<String>();
			String cram = null;
			for (String path: al) {
				if (path.endsWith(".cram")) {
					cram = path;
					break;
				}
				else if (path.endsWith("q.gz")) fastq.add(path);
			}
			al.clear();
			if (cram!=null) al.add(cram);
			else {
				int numFastqs = fastq.size();
				if (numFastqs == 2) {
					al.addAll(fastq);
				}
				else if (numFastqs > 3) {
					//cases where the normal fastqs are available from multiple sample sets, take first set
					String firstPath = fastq.get(0);
					firstPath = firstPath.substring(0, firstPath.lastIndexOf('/')+1);
					
					String secondPath = null;
					for (int i=1; i< numFastqs; i++) {
						String test = fastq.get(i);
						if (test.startsWith(firstPath)) {
							secondPath = test;
							break;
						}
					}
					if (secondPath != null) {
						al.add(fastq.get(0));
						al.add(secondPath);
					}
					else {
						IO.pl("\tFAILED to find parse just two fastq for "+id);
						System.exit(0);
					}
				}
				else {
					IO.pl("\tFAILED to find one cram or two fastq for "+id);
					System.exit(0);
				}
			}
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
		
		//for each tumor specimen
		for (String specimineId: specimineIdTumorSamples.keySet()) {

			//Define the Tumor exome and rna samples
			ArrayList<TumorSampleADW> tumorSamples = specimineIdTumorSamples.get(specimineId);

			//attempt to merge those with separate RNA and Exome ids on different lines
			tumorSamples = mergeSplitTumorExomeRNADatasets(tumorSamples);

			//sometimes there is more than one due to multiple tumor exomes with multiple tumor rna
			for (TumorSampleADW tumorSample: tumorSamples) {

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
							IO.pl("\t\t"+ns.getNormalWesFastqPathsToFetch());
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
								IO.pl("\t\t"+ns.getPlatformName()+"\t"+ ns.getNormalWesFastqPathsToFetch());
								matchedPlatform = false;
							}
						}
					}
				}
				//add the job to the patient
				AvatarAnalysisJob adJob = new AvatarAnalysisJob(p, tumorSample, normalSamplesToAdd, matchedPlatform);
				p.getAnalysisJobs().add(adJob);
			}
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
				exomeOnly.setTumorRnaPathsToFetch(rnaOnly.getTumorRnaFastqPathsToFetch());
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
						IO.pl("\tFailed to find a platform type for "+tumorGermline+" sample "+
								wesId+", skipping patient "+patientId);

						OK = false;
					}
				}
				
				//is it a tumor dataline? could be tumor exome, tumor rna, or both
				if (tumorGermline.equals("Tumor")) {
					TumorSampleADW ts = new TumorSampleADW(wesId, rnaId, platform, specimineId, fields);
					
					//any tumor wes?
					if (wesId != null) {
						HashSet<String> nameTumorWes = slIdFastqPaths.get(wesId);
						if (nameTumorWes == null) nameTumorWes = slIdAwsPaths.get(wesId);
						if (nameTumorWes == null || checkFastqCram(nameTumorWes)==false) {
							OK = false;
							IO.pl("\tFailed to find tumor WES fastq or cram files '"+nameTumorWes+"' for "+ wesId+", skipping patient "+patientId);
						}
						else ts.setTumorWesFastqPathsToFetch(nameTumorWes);
					}
					
					//any tumor rna
					if (rnaId != null) {
						HashSet<String>  nameTumorRna = slIdFastqPaths.get(rnaId);
						if (nameTumorRna == null) nameTumorRna = slIdAwsPaths.get(rnaId);
						if (nameTumorRna == null || checkFastqCram(nameTumorRna)==false) {
							OK = false;
							IO.pl("\tFailed to find tumor RNA fastq or cram files '"+nameTumorRna+"' for "+ rnaId+", skipping patient "+patientId);
						}
						else {
							ts.setTumorRnaPathsToFetch(nameTumorRna);
						}
					}
					p.getTumorSamples().add(ts);
				}
				
				//is it a germline normal dataline
				else if (tumorGermline.equals("Germline")) {
					HashSet<String>  nameNormalWes = slIdFastqPaths.get(wesId);
					if (nameNormalWes == null) nameNormalWes = slIdAwsPaths.get(wesId);
					if (nameNormalWes == null || checkFastqCram(nameNormalWes)==false) {
						OK = false;
						IO.pl("\tFailed to find normal WES fastq or cram files '"+nameNormalWes+"' for "+ wesId+", skipping patient "+patientId);
					}
					else {
						NormalSampleADW ns = new NormalSampleADW(wesId, platform, fields);
						ns.setNormalWesFastqPathsToFetch(nameNormalWes);
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
		IO.pl("\t"+ numPatientsWithMissingData+"\t# Patients missing cram or fastq files, see above.");
	}
	
	/**checks if have just one cram or two fastq.gz files*/
	public static boolean checkFastqCram(HashSet<String> paths) {
		int numCram = 0;
		int numFastq = 0;
		for (String path: paths) {
			if (path.endsWith("q.gz")) numFastq++;
			else if (path.endsWith(".cram")) numCram++;
		}
		if (numCram == 1 && numFastq == 0) return true;
		if (numCram == 0 && numFastq == 2) return true;
		return false;
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
					case 'a': asterDownloadDir = new File(args[++i]); break;
					case 't': tmpDir = new File(args[++i]); break;
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
		//jobDir
		if (jobDir == null) Misc.printErrAndExit("\nERROR: cannot find your JobDir? Aborting. "+jobDir);
		jobDir.mkdirs();
		if (jobDir.exists() == false || jobDir.canWrite() == false) Misc.printErrAndExit("\nERROR: cannot write into your JobDir? Aborting. "+jobDir);
		
		//Aster download dir
		if (asterDownloadDir == null) Misc.printErrAndExit("\nERROR: cannot find your AsterDownloadDir? Aborting. "+asterDownloadDir);
		asterDownloadDir.mkdirs();
		if (asterDownloadDir.exists() == false || jobDir.canWrite() == false) Misc.printErrAndExit("\nERROR: cannot write into your AsterDownloadDir? Aborting. "+asterDownloadDir);

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

	public static void printDocs(){

		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                           Avatar Data Wrangler : Sept 2024                       **\n" +
				"**************************************************************************************\n" +
				"Tool for assembling directories for TNRunner based on files provided by Aster and AWS.\n"+
				"Handles patient datasets from different exome capture platforms with multiple tumor\n"+
				"and normal samples.\n"+
				"Each TNRunner Job will be named as AvatarID-NormalExomeID-TumorExomeID-TumorRNASeqID.\n"+
				"When missing they will be labeled NA. Uses the SubjectMatchMaker to fetch or make\n"+
				"molecular data patient ids shared with the other clinically processed datasets.\n"+
				"Four files are written to the -r dir: downloadAsterFiles.sh, awsFilesToDownload.txt,\n"+
				"moveAsterFiles.sh, makeSoftLinks.sh, and normalsToMerge.txt for additional processing.\n"+

				"\nRequired Parameters:\n"+
				"-r Path to a folder containing the following resource files.  Each must end with the\n"+
				"     given extensions:\n"+
				"        xxx_ClinicalMolLinkage_V4.csv - Aster\n"+
				"        xxx_WES_QC_Metrics.csv - Aster\n"+
				"        fastqFiles.txt - Aster\n"+
				"        xxx_AWSRepoList.txt - AWS\n"+
				"        (optional) xxx_PatientPHI.csv - RISR's PatientIdCoverter\n"+
				"   See ~/TNRunner/Workflows/Auto/Avatar/ for scripts to pull these.\n"+
				"-j Job dir to build out the patient dir structure.\n"+
				"-t Directory to place temp files with PHI for Subject ID matching\n"+
				"-s Directory containing the SubjectMatchMaker 'currentRegistry_' file\n"+
				"-a Path to where the primary seq datasets will be downloaded.\n"+

                "\nExample: java -jar ~/USeqApps/AvatarDataWrangler -r Resources/ -j AJobs -t SMM_PHI\n"+
                "   -s ~/PHI/SmmRegistry/ -d asterDownloadCmds.sh -a AsterDownloads\n\n"+


				"**************************************************************************************\n");
	}
}
