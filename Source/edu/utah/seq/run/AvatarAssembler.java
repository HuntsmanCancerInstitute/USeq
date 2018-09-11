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

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class AvatarAssembler {
	
	private File gender = null;
	private File diagnosis = null;
	private File info = null;
	private File path = null;
	private File results = null;
	private File linkDir = null;
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
		else if (is.tEPresent && is.tTPresent) numTeTt++;
		else if (is.nEPresent && is.tTPresent) numNeTt++;
		else if (is.nEPresent) numNe++;
		else if (is.tEPresent) numTe++;
		else if (is.tTPresent) numTt++;
	}
	
	public AvatarAssembler (String[] args){
		long startTime = System.currentTimeMillis();
		try {
			processArgs(args);
			
			loadHciIdHash();
			
			loadGender();
			
			loadDiagnosis();

			outputStatPatients();
			
			for (AvatarPatient ap: toLink) ap.writeSoftLinks();
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
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



	private void loadDiagnosis() throws IOException {
		BufferedReader in = IO.fetchBufferedReader(diagnosis);
		String line = in.readLine();
		while ((line = in.readLine()) != null){
			String[] fields = Misc.TAB.split(line);
			if (fields.length != 2) {
				continue;
				//throw new IOException("Incorrect number "+fields.length+" of fields in "+line);
			}
			AvatarSample as = sampleId2AvatarPatient.get(fields[0]);
			if (as == null) {
				continue;
				//throw new IOException("Failed to find a sample assoicated with "+line);
			}
			as.diagnosis = fields[1];
		}
		in.close();
	}

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
			if (fields.length != 7) throw new IOException("Incorrect number "+fields.length+" of fields in "+line);
			int hciId = Integer.parseInt(fields[4]);
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
		
		public AvatarPatient(String[] fields) {
			hciId = Integer.parseInt(fields[4]);
			gExperimentId = fields[1];
			gAnalysisId = fields[2];
		}
		
		public void writeSoftLinks() {
			File patFastqDir = null;
			File link = null;
			try {
				//check counts if weird then don't output, needs custom
				if (samplesInfo.justSingle() == false) IO.pl("Needs custom:"+this.toString());

				//create dir to put linked files
				patFastqDir = new File (linkDir, hciId+"/Fastq/");
				patFastqDir.mkdirs();
				//for each sample
				for (AvatarSample as: avatarSamples){
					//Source? Tumor or Normal
					String dirName = as.sampleSource+ as.sampleType;
					File fastqDir = new File (patFastqDir, dirName);
					fastqDir.mkdirs();
					//make soft links
					for (File oriF: as.fastq){
						link = new File(fastqDir, dirName+"_"+oriF.getName());
						if (links.contains(link.toString())){
							IO.pl("Already exists "+link.toString());
						}
						else {
							Files.createSymbolicLink(link.toPath(), oriF.toPath());
							links.add(link.toString());
						}
					}
				}
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
				//IO.pl("Duplicate/ broken links? "+this.toString()+"BadLink "+link+"\n");
				//IO.deleteDirectory(patFastqDir);
			}
		}
		
		public void addSample(String[] fields) {
			//does it exist?
			AvatarSample as = sampleId2AvatarPatient.get(fields[0]);
			if (as == null){
				as = new AvatarSample(fields, this);
				sampleId2AvatarPatient.put(fields[0], as);
				avatarSamples.add(as);
			}
			else as.fastq.add(new File(path, fields[6]));
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

		public boolean justSingle() {
			if (numTE > 1 || numNE > 1 || numTT >1) return false;
			return true;
		}
	}
	
	private class AvatarSample {
		AvatarPatient avatarPatient = null;
		String sampleId = null;
		String sampleName = null;
		String diagnosis = null;
		
		//Exome or Transcriptome
		String sampleType = null;
		//Tumor or Normal
		String sampleSource = null;
		
		TreeSet<File> fastq = new TreeSet<File>();

		public AvatarSample(String[] fields, AvatarPatient avatarPatient) {
			this.avatarPatient = avatarPatient;
			sampleId = fields[0];
			sampleName = fields[3];
			//set type
			if (fields[5].contains("WES")) {
				sampleType = "Exome";
				//set source
				if (fields[5].contains("Tumor")) sampleSource = "Tumor";
				else if (fields[5].contains("Germline")) sampleSource = "Normal";
			}
			else if (fields[5].contains("RNA")) {
				sampleType = "Transcriptome";
				sampleSource = "Tumor";
			}
			fastq.add(new File(path, fields[6]));
		}
		
		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append("\tSampleName\t"+ sampleName); sb.append("\n");
			sb.append("\t\tSampleId\t"+         sampleId); sb.append("\n");
			sb.append("\t\tDiagnosis\t"+   diagnosis); sb.append("\n");
			sb.append("\t\tSampleType\t"+        sampleType); sb.append("\n");
			sb.append("\t\tSampleSource\t"+sampleSource);  sb.append("\n");
			for (File f: fastq) {
				sb.append("\t\t");
				sb.append(f.toString());
				sb.append("\n");
			}
			return sb.toString();
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
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		File[] files = new File[]{gender, diagnosis, info, path, results};
		for (File f: files) if (f== null) Misc.printErrAndExit("Error: missing one of the required five files (-g -d -i -p -r).");
		linkDir = new File(results.getParentFile(), "PatientBatch"+Misc.getDateNoSpaces()).getCanonicalFile();
		IO.deleteDirectoryViaCmdLine(linkDir);
		linkDir.mkdirs();
	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                                TNRunner  September 2018                          **\n" +
				"**************************************************************************************\n" +
				"Tool for assembling fastq avatar datasets based on the results of three sql queries.\n"+
				"Beta!\n"+

				"\nOptions:\n"+
				"-i info.\n" +
				"-d diagnosis.\n"+
				"-g gender.\n"+
				"-p path to Exp dir, e.g. /Repository/PersonData/2018/\n"+
				"-r results file.\n"+
				


				"**************************************************************************************\n");
	}

}
