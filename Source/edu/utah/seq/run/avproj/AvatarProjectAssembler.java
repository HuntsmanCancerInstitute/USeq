package edu.utah.seq.run.avproj;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

public class AvatarProjectAssembler {

	//user fields
	private File genderFile = null;
	private File linkageFile = null;
	private File platformFile = null;
	private File jobDir = null;
	private File fastqDir = null;
	private boolean verbose = false;

	//internal
	private HashMap<String, String> patientGender = null;
	private HashMap<String, String> dataSetPlatform = null;
	private HashMap<String, File> nameFastqFile = new HashMap<String, File>();
	private HashMap<String, Patient> patients = new HashMap<String, Patient>();
	private int numJobDirs = 0;
	private StringBuilder normalsToDelete = new StringBuilder();

	

	public AvatarProjectAssembler (String[] args){
		long startTime = System.currentTimeMillis();
		try {
			processArgs(args);

			loadHashes();

			loadFastqCram();

			walkLinkageFile();
			
			makeJobs();
			
			if (normalsToDelete.length() !=0) IO.pl("\nExecute before running CNV, JointGenotyping, AgggregateQCStats, or any analysis that requires a capture design matched normal:\n"+ normalsToDelete);

		} catch (IOException e) {
			e.printStackTrace();
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
	}



	private void makeJobs() throws IOException {

		int numWithFastq = 0;
		for (Patient p: patients.values()) {
			if (p.getTumorSamples().size() !=0 || p.getNormalSamples().size() !=0) {
				numWithFastq++;
				buildJob(p);
			}
		}
		
		IO.pl("\nNumber patients "+ patients.size());
		IO.pl("Number patients with fastq/cram files "+ numWithFastq);
		IO.pl("Number patient job directories "+ numJobDirs);
		
	}



	private void buildJob(Patient p) throws IOException {
		
		IO.pl("Building folders for "+p.getPatientId());

		//split samples by platform
		HashMap<String, ArrayList<TumorSample>> platformTumorSamples = splitTumorSamplesByPlatform(p);
		HashMap<String, ArrayList<NormalSample>> platformNormalSamples = splitNormalSamplesByPlatform(p);
		HashSet<String> allPlatforms = new HashSet<String>();
		for (String plat: platformTumorSamples.keySet()) allPlatforms.add(plat);
		for (String plat: platformNormalSamples.keySet()) allPlatforms.add(plat);
		
		//for each platform
		for (String plat: allPlatforms) {
			ArrayList<NormalSample> normalSamples = platformNormalSamples.get(plat);
			ArrayList<TumorSample> tumorSamples = platformTumorSamples.get(plat);
			int numNorm = 0;
			if (normalSamples != null) numNorm = normalSamples.size();
			int numTum = 0;
			if (tumorSamples != null) numTum = tumorSamples.size();

			if (verbose) IO.pl("\tNumber normal samples "+numNorm);
			//more than one normal? throw error and skip
			if (numNorm > 1) {
				IO.el("\tERROR: skipping subject "+p.getPatientId()+ ", multiple normal files found in the same platform, cat these and modify the linkage file, then restart: ");
				for (NormalSample ns: normalSamples) {
					for (File f: ns.getNormalDnaFastqCram()) IO.el("\t"+f);
				}
				return;
			}

			NormalSample ns = null;
			if (numNorm !=0) {
				ns = normalSamples.get(0);
				for (NormalSample n: normalSamples) {
					if (n.getNormalDnaFastqCram().size() == 0) IO.el("\tWARNING: subject "+p.getPatientId()+" is missing the normal/ germline fastq/cram file(s) for "+n.getNormalDnaName()+", putting dataset name in job dir name but no files to link in! Find them and rerun." );
				}
			}
			if (verbose) IO.pl("\tNumber tumor samples "+numTum);
			
			//No tumor samples just yet?
			if (numTum == 0) {
				//create patient dir, patientID-Platform-NormalExomeID-TumorExomeID-TumorRNAID 
				File dir = new File(jobDir, p.getPatientId()+"_"+ns.getPlatformName()+"_"+ns.getNormalDnaName()+"_NA_NA");
				dir.mkdirs();
				addNormalSample(ns, dir);
				addGenderWithPlatform(p.getGender(),dir);
				numJobDirs++;
			}
			else {
				for (TumorSample ts: tumorSamples) {
					if (ts.getTumorDnaName() != null && ts.getTumorDnaFastqCram().size() ==0) {
						IO.el("\tWARNING: subject "+p.getPatientId()+" is missing Tumor DNA fastq/cram file(s) for "+ts.getTumorDnaName()+", putting dataset name in job dir name but no files to link in! Find them and rerun." );
					}
					if (ts.getTumorRnaName() != null && ts.getTumorRnaFastqCram().size() ==0) {
						IO.el("\tWARNING: subject "+p.getPatientId()+" is missing Tumor RNA fastq/cram file(s) for "+ts.getTumorRnaName()+", putting dataset name in job dir name but no files to link in! Find them and rerun." );
					}
					//create patient dir, patientID_Platform_NormalExomeID_TumorExomeID_TumorRNAID 
					String platform = "NA";
					if (ts.getPlatformName() != null) platform = ts.getPlatformName();
					String normDnaName = "NA";
					if (ns != null) {
						normDnaName = ns.getNormalDnaName();
					}
					String tumorDnaName = "NA";
					
					boolean normalWarning = false;
					if (ts.getTumorDnaName() != null) {
						tumorDnaName = ts.getTumorDnaName();
						//OK tumor exome present, is there a normal? if not try to find one from the other platforms
						if (ns == null) {
							for (String pl: allPlatforms) {
								ArrayList<NormalSample> nps = platformNormalSamples.get(pl);
								if (nps !=null && nps.size() == 1) {
									ns = nps.get(0);
									normalWarning = true;
									break;
								}
							}
						}
					}
					String tumorRnaName = "NA";
					if (ts.getTumorRnaName() != null) tumorRnaName = ts.getTumorRnaName();
						
					File dir = new File(jobDir, p.getPatientId()+"_"+platform+"_"+normDnaName+"_"+tumorDnaName+"_"+tumorRnaName);
					dir.mkdirs();
					addNormalSample(ns, dir);
					if (normalWarning) {
						IO.el("\tWARNING: added non platform matched Normal DNA sample for somatic variant calling.");
						normalsToDelete.append("rm -rf "+dir+"/Fastq/NormalDNA/\n");
						normalsToDelete.append("rm -rf "+dir+"/Alignment/*_NormalDNA/\n");
					}
					addTumorDnaSample(ts, dir);
					addTumorRnaSample(ts, dir);
					addGenderWithPlatform(p.getGender(),dir);
					numJobDirs++;
				}
				
			}
		}
		
	}
	
	private void addGenderWithPlatform(String gender, File dir) {
		File clinDir = new File (dir, "ClinicalReport");
		clinDir.mkdir();
		
		//From Tempus
		//TL-18-843E9B_XT.V1_2018-10-26_deid_Neeraj_Agarwal_M.json
		//	     0         1        2       3     4      5    6

		//For Avatar
		//AvatarID_Platform_Gender.json
		//ZBG7MKFOZ1_IDT.V1_M.json
		//File dir = new File(jobDir, p.getPatientId()+"_"+platform+"_"+normDnaName+"_"+tumorDnaName+"_"+tumorRnaName);
		//                                   0                 1
		String[] info = Misc.UNDERSCORE.split(dir.getName());
		//convert platform from NIM and IDT to be NIM.V1 ditto for IDT, there is an IDT.V2 coming!
		if (info[1].contains(".V") == false) info[1] = info[1]+".V1";
		
		String genderLetter = null;
		if (gender.toLowerCase().equals("female")) genderLetter = "F";
		else if (gender.toLowerCase().equals("male")) genderLetter = "M";

		File json = new File(clinDir, info[0]+"_"+ info[1]+"_"+ genderLetter+ ".json");
		IO.writeString("Gender "+gender+"\n", json);
		if (verbose) IO.pl("\tAdding gender info file to "+json);
	}



	private void addTumorRnaSample(TumorSample ts, File dir) throws IOException {
		if (ts.getTumorRnaFastqCram().size() ==0) return;
		//create fastq
		File fastqDir = new File (dir, "Fastq/TumorRNA");
		fastqDir.mkdirs();
		//link in files
		for (File f: ts.getTumorRnaFastqCram()) {
			File link = new File(fastqDir, f.getName());
			if (link.exists() == false) Files.createSymbolicLink(link.toPath(), f.toPath());
			if (verbose) IO.pl("\tLinking in Tumor Sample RNA "+f.getName()+" into "+fastqDir);
		}
		
	}

	private void addTumorDnaSample(TumorSample ts, File dir) throws IOException {
		if (ts.getTumorDnaFastqCram().size() ==0) return;
		//create fastq
		File fastqDir = new File (dir, "Fastq/TumorDNA");
		fastqDir.mkdirs();
		//link in files
		for (File f: ts.getTumorDnaFastqCram()) {
			File link = new File(fastqDir, f.getName());
			if (link.exists() == false) Files.createSymbolicLink(link.toPath(), f.toPath());
			if (verbose) IO.pl("\tLinking in Tumor Sample DNA "+f.getName()+" into "+fastqDir);
		}
		
	}

	private void addNormalSample(NormalSample ns, File dir) throws IOException {
		if (ns == null) return;
		//create fastq
		File fastqDir = new File (dir, "Fastq/NormalDNA");
		fastqDir.mkdirs();
		//link in files
		for (File f: ns.getNormalDnaFastqCram()) {
			File link = new File(fastqDir, f.getName());
			if (link.exists() == false) Files.createSymbolicLink(link.toPath(), f.toPath());
			if (verbose) IO.pl("\tLinking in NormalSample "+f.getName()+" into "+fastqDir);
		}
		
	}



	private HashMap<String, ArrayList<NormalSample>> splitNormalSamplesByPlatform(Patient p) {
		HashMap<String, ArrayList<NormalSample>> platformNormalSamples = new HashMap<String, ArrayList<NormalSample>>();
		for (NormalSample ns: p.getNormalSamples()) {
			String platform = ns.getPlatformName();
			if (platform == null) platform= "NA";
			ArrayList<NormalSample> al = platformNormalSamples.get(platform);
			if (al == null) {
				al = new ArrayList<NormalSample>();
				platformNormalSamples.put(platform, al);
			}
			al.add(ns);
		}
		return platformNormalSamples;
	}



	private HashMap<String, ArrayList<TumorSample>> splitTumorSamplesByPlatform(Patient p) {
		HashMap<String, ArrayList<TumorSample>> platformTumorSamples = new HashMap<String, ArrayList<TumorSample>>();
		for (TumorSample ts: p.getTumorSamples()) {
			String platform = ts.getPlatformName();
			if (platform == null) platform= "NA";
			ArrayList<TumorSample> al = platformTumorSamples.get(platform);
			if (al == null) {
				al = new ArrayList<TumorSample>();
				platformTumorSamples.put(platform, al);
			}
			al.add(ts);
		}
		return platformTumorSamples;
	}



	private void walkLinkageFile() throws IOException {
		BufferedReader in = IO.fetchBufferedReader(linkageFile);
		String headerLine = in.readLine();
		String[] fields = Misc.TAB.split(headerLine);
		int keyIndex = -1;
		int typeIndex = -1;
		int wesIndex = -1;
		int rnaIndex = -1;
		for (int i=0; i< fields.length; i++) {
			if (fields[i].equals("ORIENAvatarKey")) keyIndex = i;
			else if (fields[i].equals("Tumor/Germline")) typeIndex = i;
			else if (fields[i].equals("WES")) wesIndex = i;
			else if (fields[i].equals("RNASeq")) rnaIndex = i;
		}
		if (keyIndex == -1 || typeIndex == -1 || wesIndex == -1 || rnaIndex == -1) {
			in.close();
			Misc.printErrAndExit("\nERROR: failed to "
					+ "identify one or more of the following linkage file header columns, check that these are in your file: "
					+ "ORIENAvatarKey Tumor/Germline WES RNASeq");
		}

		//for each subsequent line
		String line = null;
		try {
		while ((line = in.readLine()) != null) {
			fields = Misc.TAB.split(line);
			String patientId = fields[keyIndex];
			String tumorGermline = fields[typeIndex];
			String wesId = fields[wesIndex];
			String rnaId = fields[rnaIndex];
			if (wesId.length()==0) wesId= null;
			if (rnaId.length()==0) rnaId= null;

			//fetch the patient, will create a new one if not present and fetch the gender.
			Patient p = fetchPatient(patientId);

			//pull platform, will be null for rna only samples
			String platform = null;
			if (wesId!=null) {
				platform = dataSetPlatform.get(wesId);
				if (platform == null) {
					in.close();
					Misc.printErrAndExit("\nERROR: failed to find a platform type for "+tumorGermline+" sample "+
						wesId+" for patient ID "+patientId+" Aborting.");
				}
			}
			
			if (tumorGermline.equals("Tumor")) {
				TumorSample ts = new TumorSample(wesId, rnaId, platform);
				if (ts.fetchFastqCrams(nameFastqFile, patientId) == true) p.getTumorSamples().add(ts);
			}
			else if (tumorGermline.equals("Germline")) {
				NormalSample ns = new NormalSample(wesId, platform);
				if (ns.fetchFastqCrams(nameFastqFile, patientId) == true) p.getNormalSamples().add(ns);
			}
			else {
				in.close();
				Misc.printErrAndExit("\nERROR: didn't find 'Tumor' or 'Germline' in the type column for "+line);
			}
		}
		} catch (ArrayIndexOutOfBoundsException e) {
			IO.el("\nERROR parsing line "+line);
			e.printStackTrace();
			System.exit(1);
		}
	}

	private Patient fetchPatient(String patientId) {
		Patient p = patients.get(patientId);
		if (p != null) return p;
		String gender = patientGender.get(patientId);
		if (gender == null) Misc.printErrAndExit("\nERROR: failed to find a gender for patient ID "+patientId+" Aborting.");
		p = new Patient(patientId, gender);
		patients.put(patientId, p);
		return p;
	}

	private void loadFastqCram() throws IOException {
		File[] allGz = IO.extractFiles(fastqDir, ".gz");
		for (File f: allGz) nameFastqFile.put(f.getName(), f.getCanonicalFile());
		File[] allCrams = IO.extractFiles(fastqDir, ".cram");
		for (File f: allCrams) nameFastqFile.put(f.getName(), f.getCanonicalFile());
	}

	private void loadHashes() {
		patientGender = IO.loadFileIntoHashMap(genderFile);
		dataSetPlatform = IO.loadFileIntoHashMap(platformFile);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AvatarProjectAssembler(args);
	}		

	/**This method will process each argument and assign new variables
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'g': genderFile = new File(args[++i]); break;
					case 'l': linkageFile = new File(args[++i]); break;
					case 'p': platformFile = new File(args[++i]); break;
					case 'j': jobDir = new File(args[++i]); break;
					case 'f': fastqDir = new File(args[++i]); break;
					case 'v': verbose = true; break;
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
		File[] files = new File[]{genderFile, linkageFile, platformFile, jobDir, fastqDir};
		for (File f: files) if (f== null) Misc.printErrAndExit("Error: missing one of the required five files (-g -l -p -j -f)");

		jobDir.mkdirs();
		if (jobDir.exists() == false || jobDir.canWrite() == false) Misc.printErrAndExit("\nERROR: cannot write into your JobDir? Aborting. "+jobDir);
	}


	public static void printDocs(){

		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                          Avatar Project Assembler : Nov 2022                    **\n" +
				"**************************************************************************************\n" +
				"Tool for assembling directories for TNRunner based on files provided by M2Gen.\n"+
				"Handles patient datasets from different exome capture platforms with multiple tumor\n"+
				"samples. Convert each M2Gen xxx.csv to tab delimited using Excel's save as option.\n"+
				"Each TNRunner Job will be named as \n"+
				"PatientID-Platform-NormalExomeID-TumorExomeID-TumorRNASeqID, when missing they will be\n"+
				"labeled NA. Fastq files are linked, not moved, so delete the Jobs dir to start over.\n"+


				"\nRequired Parameters:\n"+
				"-l Sample linkage file containing at least 4 tab delimited columns named\n"+
				"     'ORIENAvatarKey' 'Tumor/Germline' 'WES' 'RNASeq', see \n"+
				"     xxx_ClinicalMolLinkage_V4.txt \n" +
				"-g Patient gender file with just two tab delimited columns for the ORIENAvatarKey and\n"+
				"     Sex, see xxx__PatientMaster_V4.txt \n"+
				"-p Sample exome platform file with just two tab delimited columns, the sample IDs that\n"+
				"     match those in the 'WES' linkage file and the platform type NIM or IDT, see\n"+
				"     xxx_WES_Resource_table.txt \n"+
				"-f Directory containing all of the paired fastq and M2Gen cram files for the T/N exomes\n"+
				"     and T RNASeq, these should begin with the sample IDs that match those in the 'WES'\n"+
				"     and 'RNASeq' columns in the linkage file.\n"+
				"-j Job dir to build out the patient dir structure and linked fastqs. Existing files\n"+
				"     won't be overwritten.\n"+
				"-v Verbose debugging output.\n"+
				

                "\nExample: java -jar ~/USeqApps/AvatarProjectAssembler -l linkage.txt -g gender.txt -p\n"+
                "     platform.txt -f ~/UroProj/Fastq -j ~/UroProj/TNRunnerJobs/ -v\n\n"+


				"**************************************************************************************\n");
	}
}
