package edu.utah.seq.run;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.json.JSONObject;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class AvatarAssembler {
	
	private File gender = null;
	private File diagnosis = null;
	private File info = null;
	private File prep = null;
	private File sourceOrigin = null;
	private File path = null;
	private File results = null;
	private File jobDir = null;
	private String[] yearDirNames = new String[] {"2017", "2018", "2019", "2020", "2021", "2022"};
	private File[] yearDirs = null;
	private TreeMap<Integer, AvatarPatient> hciId2AvatarPatient = new TreeMap<Integer, AvatarPatient>();
	private HashMap<String, AvatarPatient> experimentId2AvatarPatient = new HashMap<String, AvatarPatient>(); 
	private HashMap<String, AvatarSample> sampleId2AvatarPatient = new HashMap<String, AvatarSample>();
	private int numTeNeTt = 0;
	private int numTeNe = 0;
	private int numTeTt = 0;
	private int numNeTt = 0;
	private int numNe = 0;
	private int numTe = 0;
	private int numTt = 0;
	private HashSet<String> links = new HashSet<String>();
	private ArrayList<AvatarPatient> toLink = new ArrayList<AvatarPatient>();
	private boolean linkOnlyTN = true;
	private String diagnosisFilter = null;
	private boolean printFastqErrors = true;
	
	private void incrementSampleCounters(AvatarPatient ap) {
		SamplesInfo is = ap.createSamplesInfo();
		if (is.tEPresent && is.nEPresent && is.tTPresent) {
			numTeNeTt++;
			toLink.add(ap);
		}
		else if (is.tEPresent && is.nEPresent) {
			toLink.add(ap);
			numTeNe++;
		}
		else if (is.tEPresent && is.tTPresent) {
			numTeTt++;
			if (linkOnlyTN == false) toLink.add(ap);
		}
		else if (is.nEPresent && is.tTPresent) {
			numNeTt++;
			if (linkOnlyTN == false) toLink.add(ap);
		}
		else if (is.nEPresent) {
			numNe++;
			if (linkOnlyTN == false) toLink.add(ap);
		}
		else if (is.tEPresent) {
			numTe++;
			if (linkOnlyTN == false) toLink.add(ap);
		}
		else if (is.tTPresent) {
			numTt++;
			if (linkOnlyTN == false) toLink.add(ap);
		}
		
	}
	
	public AvatarAssembler (String[] args){
		long startTime = System.currentTimeMillis();
		try {
			processArgs(args);
			
			loadHciIdHash();
			
			loadGender();
			
			loadDiagnosis();
			
			loadOrigin();
			
			loadPrep();
			
			filterPatients();

			outputStatPatients();
			
			for (AvatarPatient ap: toLink) ap.buildJobDirs();
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
	}
	
	private void filterPatients() {
		if (diagnosisFilter == null) return;
		String lcDF = diagnosisFilter.toLowerCase();
		TreeMap<Integer, AvatarPatient> keep = new TreeMap<Integer, AvatarPatient>();
		for (Integer i: hciId2AvatarPatient.keySet()) {
			AvatarPatient ap = hciId2AvatarPatient.get(i);
			for (AvatarSample as: ap.avatarSamples){
				String lcDig = as.diagnosis.toLowerCase();
				if (lcDig.contains(lcDF)){
					keep.put(i, ap);
					break;
				}
			}
		}
		hciId2AvatarPatient = keep;
	}

	private void outputStatPatients() throws FileNotFoundException, IOException {
		Gzipper out = new Gzipper(results);
		for (AvatarPatient ap: hciId2AvatarPatient.values()) {
			//write out dump
			out.println(ap);
			//increment sample info
			incrementSampleCounters(ap);
		}
		
		int total = numTeNeTt + numTeNe + numTeTt + numNeTt + numNe + numTe + numTt;
		out.println(hciId2AvatarPatient.size()+ "\tNum Patient Sets");
		out.println("\t"+ numTeNeTt + "\tTE NE TT ");
		out.println("\t"+ numTeNe  + "\tTE NE ");
		out.println("\t"+ numTeTt  + "\tTE TT ");
		out.println("\t"+ numNeTt + "\tNE TT ");
		out.println("\t"+numNe + "\tNE ");
		out.println("\t"+numTe + "\tTE ");
		out.println("\t"+numTt + "\tTT ");
		out.println("\t\t"+total + "\tTotal");
		
		out.closeNoException();
	}


	//these are all of the diagnosis from Foundation, Tempus, and Avatar so many won't have Avatar info
	private void loadDiagnosis() throws IOException {
		BufferedReader in = IO.fetchBufferedReader(diagnosis);
		String line = in.readLine();
		while ((line = in.readLine()) != null){
			String[] fields = Misc.TAB.split(line);
			if (fields.length != 2) {
				//often broken
				continue;
				//throw new IOException("Incorrect number "+fields.length+" of fields in "+line);
			}
			AvatarSample as = sampleId2AvatarPatient.get(fields[0]);
			if (as == null) {
				//IO.pl("WARNING failed to find a sample associated with this diagnosis "+line);
				continue;
				//throw new IOException("Failed to find a sample assoicated with "+line);
			}
			as.diagnosis = fields[1];
		}
		in.close();
	}
	
	//these are all of the origins from Foundation, Tempus, and Avatar so many won't have Avatar info
	private void loadOrigin() throws IOException {
		BufferedReader in = IO.fetchBufferedReader(sourceOrigin);
		String line = in.readLine();
		while ((line = in.readLine()) != null){
			String[] fields = Misc.TAB.split(line.trim());
			if (fields.length != 2) continue;
			AvatarSample as = sampleId2AvatarPatient.get(fields[0]);
			if (as == null) continue;
			as.origin = fields[1].trim();
		}
		in.close();
	}
	
	//these are all of the preps from Foundation, Tempus, and Avatar so many won't have Avatar info
	private void loadPrep() throws IOException {
		BufferedReader in = IO.fetchBufferedReader(prep);
		String line = in.readLine();
		while ((line = in.readLine()) != null){
			line = line.replaceAll("\"", "");
			String[] fields = Misc.TAB.split(line.trim());
			if (fields.length != 2) continue;
			AvatarSample as = sampleId2AvatarPatient.get(fields[0]);
			if (as == null) continue;
			as.prep = fields[1].trim();
		}
		in.close();
	}

	//these are all of the genders from Foundation, Tempus, and Avatar so many won't have Avatar info
	private void loadGender() throws IOException {
		BufferedReader in = IO.fetchBufferedReader(gender);
		String line = in.readLine();
		while ((line = in.readLine()) != null){
			String[] fields = Misc.TAB.split(line);
			if (fields.length != 2) {
				continue;
				//throw new IOException("Incorrect number "+fields.length+" of fields in "+line);
			}
			AvatarPatient ap = experimentId2AvatarPatient.get(fields[0]);
			if (ap == null) {
				continue;
				//throw new IOException("Failed to find a sample assoicated with "+line);
			}
			ap.gender = fields[1];
		}
		in.close();
	}

	private void loadHciIdHash() throws IOException {
		BufferedReader in = IO.fetchBufferedReader(info);
		String line = in.readLine();
		while ((line = in.readLine()) != null){
			String[] fields = Misc.TAB.split(line);
			if (fields.length != 8) throw new IOException("Incorrect number "+fields.length+" of fields in "+line);
			int hciId = Integer.parseInt(fields[5]);
			AvatarPatient ap = hciId2AvatarPatient.get(hciId);
			if (ap == null){
				ap = new AvatarPatient(fields);
				hciId2AvatarPatient.put(ap.hciId, ap);
				experimentId2AvatarPatient.put(ap.gExperimentId, ap);
			}
			ap.addSample(fields);
		}
		in.close();
	}

	private class AvatarPatient {
		int hciId = 0;
		String gExperimentId = null;
		String gAnalysisId = null;
		String gender = null;
		ArrayList<AvatarSample> avatarSamples = new ArrayList<AvatarSample>();
		SamplesInfo samplesInfo = null;
		HashMap<String, ArrayList<AvatarSample>> nameSamples = null;
		ArrayList<AvatarSample> normalExomes = new ArrayList<AvatarSample>();
		
		public AvatarPatient(String[] fields) {
			hciId = Integer.parseInt(fields[5]);
			gExperimentId = fields[1];
			gAnalysisId = fields[2];
		}
		
		public void buildSampleSets() {
			nameSamples = new HashMap<String, ArrayList<AvatarSample>>();
			for (AvatarSample as: avatarSamples) {
				String name = as.sampleSpecimen;
				//remove .xxx from each 14-0022632.15a to 14-0022632
				int index = name.lastIndexOf(".");
				if (index != -1) name = name.substring(0, index);
				
				
				ArrayList<AvatarSample> samples = nameSamples.get(name);
				if (samples == null) {
					samples = new ArrayList<AvatarSample>();
					nameSamples.put(name, samples);
				}
				samples.add(as);
				if (as.sampleType.equals("Exome") && as.sampleSource.equals("Normal")) normalExomes.add(as);
			}
		}
		
		public void buildJobDirs() {
			buildSampleSets();

			try {
				
			//for each specimen set
			for (String specimenName: nameSamples.keySet()) {
				ArrayList<AvatarSample> samples = nameSamples.get(specimenName);
				
				//find a tumor exome and normal exome
				AvatarSample te = null;
				int numTe = 0;
				AvatarSample ne = null;
				int numNe = 0;
				AvatarSample tt = null;
				int numTT = 0;
				for (AvatarSample as: samples) {
					if (as.sampleType.equals("Exome")) {
						if (as.sampleSource.equals("Tumor")) {
							te = as;
							numTe++;
						}
						else if (as.sampleSource.equals("Normal")) {
							ne = as;
							numNe++;
						}
					}
					else if (as.sampleType.equals("Transcriptome")) {
						tt = as;
						numTT++;
					}
				}
				
				//any double counts for this specimen?
				if (numTe > 1 || numNe > 1 || numTT > 1) {
					IO.pl("\nERROR: Multiple TE's ("+numTe+"), TN's("+numNe+"), or TT's("+numTT+") with the same specimen name("+specimenName+"), needs attention, skipping:\n"+this.toString());
					return;
				}
				
				//only build a set if a TE or TT is present, skip everything else
				if (numTe == 1 || numTT == 1) {
					
					//create dir to put linked files
					File patientDir = new File (jobDir, new Integer(hciId).toString()+"_"+specimenName);
					File patFastqDir = new File (patientDir, "Fastq");
					patFastqDir.mkdirs();
					
					//TE?
					if (numTe == 1) {
						//make the dir
						File fastqDir = new File (patFastqDir, "TumorDNA");
						fastqDir.mkdirs();
						//add the fastq links
						for (File oriF: te.fastq){
							File link = new File(fastqDir, "TumorDNA_"+oriF.getName());
							Files.createSymbolicLink(link.toPath(), oriF.toPath());
							links.add(link.toString());
						}
						
						//normal exome?
						//if null try to find one associated with the patient
						if (numNe == 0 && normalExomes.size() > 0) ne = normalExomes.get(0);
						if (ne != null) {
							//make the dir
							File fastqDirN = new File (patFastqDir, "NormalDNA");
							fastqDirN.mkdirs();
							//add the fastq links
							for (File oriF: ne.fastq){
								File link = new File(fastqDirN, "NormalDNA_"+oriF.getName());
								Files.createSymbolicLink(link.toPath(), oriF.toPath());
								links.add(link.toString());
							}
						}
					}
					
					//TT?
					if (numTT == 1) {
						//make the dir
						File fastqDir = new File (patFastqDir, "TumorRNA");
						fastqDir.mkdirs();
						//add the fastq links
						for (File oriF: tt.fastq){
							File link = new File(fastqDir, "TumorRNA_"+oriF.getName());
							Files.createSymbolicLink(link.toPath(), oriF.toPath());
							links.add(link.toString());
						}
					}
					
					//write out info
					File info = new File(patientDir, hciId+"_AvatarInfo.json.gz");
					Gzipper gz = new Gzipper(info);
					gz.println(this.toJson().toString(4));
					gz.close();
				}
			}
			


			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}

		/*
		public void writeSoftLinksOld() {
			File patFastqDir = null;
			File link = null;
			try {
				//check counts if weird then don't output, needs custom
				//if (samplesInfo.justSingle() == false) IO.pl("Needs custom:"+this.toString());

				//create dir to put linked files
				File patientDir = new File (jobDir, new Integer(hciId).toString());
				patFastqDir = new File (patientDir, "Fastq");
				patFastqDir.mkdirs();

				//for each sample
				for (AvatarSample as: avatarSamples){
					String type = null;
					if (as.sampleType == "Exome") type = "DNA";
					else if (as.sampleType == "Transcriptome") type = "RNA";
					else throw new IOException("\nFailed to identify the sampleType from -> "+as.sampleType);
					//Source? Tumor or Normal
					String dirName = as.sampleSource+ type;
					File fastqDir = new File (patFastqDir, dirName);
					fastqDir.mkdirs();
					//make soft links
					for (File oriF: as.fastq){
						link = new File(fastqDir, dirName+"_"+oriF.getName());
						//try to create it, catch errors here
						try {
							Files.createSymbolicLink(link.toPath(), oriF.toPath());
							links.add(link.toString());
						} catch (FileAlreadyExistsException ex) {}
					}
				}
				//write out info
				File info = new File(patientDir, hciId+"_AvatarInfo.json.gz");
				Gzipper gz = new Gzipper(info);
				gz.println(this.toJson().toString(4));
				gz.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}*/
		
		public void addSample(String[] fields) {
			//does it exist?
			AvatarSample as = sampleId2AvatarPatient.get(fields[0]);
			if (as == null){
				as = new AvatarSample(fields, this);
				sampleId2AvatarPatient.put(fields[0], as);
				avatarSamples.add(as);
			}
			else {
				//look for fastq in each year
				boolean notFound = true;
				for (File p: yearDirs) {
					//look for fastq.gz first
					String noS3Name = fields[7].replace(".S3.txt", "");
					File f = new File (p, noS3Name);
					if (f.exists()) {
						as.fastq.add(f);
						notFound = false;
						break;
					}
					//archived in Amazon S3? don't add if real file is present
					else {
						f = new File (p, noS3Name+".S3.txt");
						if (f.exists()) {
							as.fastq.add(f);
							notFound = false;
							break;
						}
					}
				}
				if (printFastqErrors && notFound) System.err.println("\nFailed to find the fastq file "+fields[7]+" and xxx.S3.txt associated with\n"+this);
			}
		}
		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append("HciId\t"+         hciId); sb.append("\n");
			sb.append("\tExperimentId\t"+ gExperimentId+"R"); sb.append("\n");
			sb.append("\tAnalysisId\tA"+   gAnalysisId); sb.append("\n");
			sb.append("\tGender\t"+        gender); sb.append("\n");
			for (AvatarSample as: avatarSamples) sb.append(as.toString());
			return sb.toString();
		}
		public JSONObject toJson(){
			JSONObject fo = new JSONObject();
			fo.put("HciId", new Integer(hciId).toString());
			fo.put("ExperimentId", gExperimentId+"R");
			fo.put("AnalysisId", "A"+gAnalysisId);
			fo.put("Gender", gender);
			
			ArrayList<JSONObject> al = new ArrayList<JSONObject>();
			for (AvatarSample as: avatarSamples) al.add(as.toJson());
			fo.put("Samples", al);
			return fo;
		}
		
		public SamplesInfo createSamplesInfo(){
			samplesInfo = new SamplesInfo(avatarSamples);
			return samplesInfo;
		}
	}
	
	private class SamplesInfo {
		boolean tEPresent = false;
		boolean nEPresent = false;
		boolean tTPresent = false;
		int numTE = 0;
		int numNE = 0;
		int numTT = 0;
		
		public SamplesInfo(ArrayList<AvatarSample> as){
			for (AvatarSample s: as){
				String type = s.sampleType; //Exome or Transcriptome
				String source = s.sampleSource; //Tumor or Normal
				if (type.equals("Exome")){
					if (source.equals("Tumor")) {
						tEPresent = true;
						numTE++;
					}
					else if (source.equals("Normal")) {
						nEPresent = true;
						numNE++;
					}
				}
				else if (type.equals("Transcriptome")){
					tTPresent = true;
					numTT++;
				}
			}
		}
		
		/*
		public boolean justSingle() {
			if (numTE > 1 || numNE > 1 || numTT >1) return false;
			return true;
		}*/
	}
	
	private class AvatarSample {
		AvatarPatient avatarPatient = null;
		String sampleId = null;
		String sampleName = null;
		String sampleSpecimen = null;
		String diagnosis = "";
		String origin = "";
		String prep = "";
		
		//Exome or Transcriptome
		String sampleType = null;
		//Tumor or Normal
		String sampleSource = null;
		
		TreeSet<File> fastq = new TreeSet<File>();

		public AvatarSample(String[] fields, AvatarPatient avatarPatient) {
			this.avatarPatient = avatarPatient;
			sampleId = fields[0];
			sampleSpecimen = fields[3];
			//int index = sampleSpecimen.lastIndexOf(".");
			//if (index != -1) sampleSpecimen = sampleSpecimen.substring(0, index);
			sampleName = fields[4];
			//set type
			if (fields[6].contains("WES")) {
				sampleType = "Exome";
				//set source
				if (fields[6].contains("Tumor")) sampleSource = "Tumor";
				else if (fields[6].contains("Germline")) sampleSource = "Normal";
			}
			else if (fields[6].contains("RNA")) {
				sampleType = "Transcriptome";
				sampleSource = "Tumor";
			}
			//look for fastq in each year
			boolean notFound = true;
			for (File p: yearDirs) {
				//look for fastq.gz first
				String noS3Name = fields[7].replace(".S3.txt", "");
				File f = new File (p, noS3Name);
				if (f.exists()) {
					fastq.add(f);
					notFound = false;
					break;
				}
				//archived in Amazon S3? don't add if real file is present
				else {
					f = new File (p, noS3Name+".S3.txt");
					if (f.exists()) {
						fastq.add(f);
						notFound = false;
						break;
					}
				}
			}
			if (printFastqErrors && notFound) {
				System.err.println("\nFailed to find the fastq file "+fields[7]+" associated with patient "+avatarPatient.hciId);
			}
		}

		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append("\tSampleName\t"+ sampleName); sb.append("\n");
			sb.append("\t\tId\t"+         sampleId); sb.append("\n");
			sb.append("\t\tSpecimen\t"+    sampleSpecimen); sb.append("\n");
			sb.append("\t\tDiagnosis\t"+   diagnosis); sb.append("\n");
			sb.append("\t\tType\t"+        sampleType); sb.append("\n");
			sb.append("\t\tSource\t"+sampleSource);  sb.append("\n");
			sb.append("\t\tOrigin\t"+   origin); sb.append("\n");
			sb.append("\t\tPrep\t"+   prep); sb.append("\n");
			
			for (File f: fastq) {
				sb.append("\t\t");
				sb.append(f.toString());
				sb.append("\n");
			}
			return sb.toString();
		}
		
		public JSONObject toJson(){
			JSONObject fo = new JSONObject();
			fo.put("Name", sampleName);
			fo.put("Id", sampleId);
			fo.put("Specimen", sampleSpecimen);
			fo.put("Diagnosis", diagnosis);
			fo.put("Type", sampleType);
			fo.put("Source", sampleSource);
			fo.put("Origin", origin);
			fo.put("Prep", prep);
			String[] fileNames = new String[fastq.size()];
			int counter = 0;
			for (File f: fastq) fileNames[counter++] = f.toString();
			fo.put("Fastq", fileNames);
			return fo;
		}
	}

	
	
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AvatarAssembler(args);
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
					case 'g': gender = new File(args[++i]); break;
					case 'd': diagnosis = new File(args[++i]); break;
					case 'i': info = new File(args[++i]); break;
					case 'p': path = new File(args[++i]); break;
					case 'r': results = new File(args[++i]); break;
					case 'o': sourceOrigin = new File(args[++i]); break;
					case 'e': prep = new File(args[++i]); break;
					case 'y': yearDirNames = Misc.splitString(args[++i], ",");
					case 'f': diagnosisFilter = args[++i]; break;
					case 'j': jobDir = new File(args[++i]); break;
					case 'l': linkOnlyTN = false; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		File[] files = new File[]{gender, diagnosis, info, path, results, jobDir, sourceOrigin, prep};
		for (File f: files) if (f== null) Misc.printErrAndExit("Error: missing one of the required eight files (-i -d -g -o -p -r -j -e)");
		
		int numYearDirs = yearDirNames.length;
		yearDirs = new File[numYearDirs];
		for (int i=0; i< numYearDirs; i++) yearDirs[i] = new File(path, yearDirNames[i]);
		
		jobDir.mkdirs();
	}


	public static void printDocs(){
		
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                              Avatar Assembler : July 2020                        **\n" +
				"**************************************************************************************\n" +
				"Tool for assembling fastq avatar datasets based on the results of several sql queries.\n"+
				"See https://ri-confluence.hci.utah.edu/x/KwBFAg  Login on hci-clingen1.hci.utah.edu\n"+

				"\nOptions:\n"+
				"-i Info\n" +
				"-d Diagnosis\n"+
				"-g Gender\n"+
				"-o Origin\n"+
				"-e Prep\n"+
				"-p Path to Exp dir w/o Year, e.g. /Repository/PersonData/\n"+
				"-y Year dirs to examine for fastq linking, defaults to 2017,2018,2019,2020,2021,2022\n"+
				"-j Job dir to place linked fastq\n"+
				"-f Only keep patients with a diagnosis containing this String, defaults to all.\n"+
				"-l Create Fastq links for all patient datasets, defaults to just those with both a\n"+
				"    Tumor and Normal exome.\n"+
				"-r Patient stats output file.\n\n"+
				
                "Example: java -jar -Xmx2G ~/USeqApps/AvatarAssembler -p /Repository/PersonData/\n"+
                "    -r avatarAssembler.log.gz -i sampleInfo.txt -d sampleDiagnosis.txt -g \n"+
                "    sampleGender.txt -s sampleOrigin.txt -e samplePrep.txt -f HEM \n"+
                "    > avatarAssemblerProblemSamples.txt \n\n"+


				"**************************************************************************************\n");
	}
}
