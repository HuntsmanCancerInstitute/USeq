package edu.utah.seq.qc;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.json.JSONObject;
import util.gen.IO;
import util.gen.Misc;

public class AggregateQCStats2 {

	//fields
	private File saveDirectory;
	private File jobDirectory;
	private String prependString = "";

	private String alignLogMatch = ".+AlignHg38.log";
	private String dupLogMatch = ".+Markdup.log";
	private String readCovJsonMatch = ".+UniObRC.json.gz";
	private String scJsonMatch = ".+SampleConcordance.json.gz";
	private String aiJsonMatch = ".+AvatarInfo.json.gz";
	private String normalDNAMatch = ".+NormalDNA.+";
	private String tumorDNAMatch = ".+TumorDNA.+";
	
	private Pattern alignLogPattern;
	private Pattern dupLogPattern;
	private Pattern readCovJsonPattern;
	private Pattern scJsonPattern;
	private Pattern aiJsonPattern;
	private Pattern normalDNAPattern;
	private Pattern tumorDNAPattern;

	private ArrayList<SampleQC2> samples = new ArrayList<SampleQC2>();


	//constructors
	public AggregateQCStats2(String[] args){

		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);
			makePatterns();

			System.out.println("Loading samples...");
			loadSamples();

			System.out.println("\nSaving aggregated data...");
			printStatsTxtSheet();

			System.out.println("\tRead coverage...");
			printReadCoverageTxtSheet(true);
			printReadCoverageTxtSheet(false);

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Sec\n");
			
		} catch (Exception e) {
			IO.el("Problem parsing files....");
			e.printStackTrace();
			System.exit(1);
		}
	}


	private void printStatsTxtSheet() {
		StringBuilder sb = new StringBuilder();

		//add header
		sb.append(SampleQC2.fetchTabbedHeader() ); sb.append("\n");

		//for each sample
		for (SampleQC2 s: samples){
			sb.append(s.fetchTabbedLine());
			//sb.append("\n");
		}

		//add descriptions 
		sb.append("\nDescriptions:\n");
		sb.append(SampleQC2.fetchDescriptions("", "\t", "\n"));

		File txt = new File(saveDirectory, prependString+ "qcStats.xls");
		IO.writeString(sb.toString(), txt);
	}

	private void printReadCoverageTxtSheet(boolean tumor) {
		//fetch samples with given type, must watch out for old samples with diff jsons
		ArrayList<SampleQC2> samplesAl = new ArrayList<SampleQC2>();
		for (SampleQC2 s: samples){
			if (tumor && s.getTumorDnaSample()!= null && s.getTumorDnaSample().getFractionTargetBpsWithIndexedCoverage()!= null) samplesAl.add(s);
			else if (tumor == false && s.getNormalDnaSample()!= null && s.getNormalDnaSample().getFractionTargetBpsWithIndexedCoverage() != null) samplesAl.add(s);
		}
		
		//build header
		StringBuilder sb = new StringBuilder();
		sb.append("Coverage");
		for (SampleQC2 s: samplesAl){
			sb.append("\t");
			sb.append(s.getSampleName());
		}
		sb.append("\n");

		//how many X axis numbers 0x, 1x, 2x, 3x, 4x,...
		int max25 = 0;
		//for each sample
		for (SampleQC2 s: samplesAl){
			int hit25 = 0;

			if (tumor) hit25 = s.getTumorDnaSample().whenHit25ReadCoverage();
			else hit25 = s.getNormalDnaSample().whenHit25ReadCoverage();			
			if (hit25> max25) max25= hit25;
		}
		if (max25 < 11) max25 = 11;

		//start at 1 to skip 0x coverage and graph better in excel
		for (int i=1; i< max25; i++){
			sb.append(i);
			for (SampleQC2 s: samplesAl){
				sb.append("\t");
				float[] f;
				if (tumor) f = s.getTumorDnaSample().getFractionTargetBpsWithIndexedCoverage();
				else f = s.getNormalDnaSample().getFractionTargetBpsWithIndexedCoverage();
				if (i>= f.length) sb.append(0f);
				else sb.append(f[i]);
			}
			sb.append("\n");
		}

		File txt;
		if (tumor) txt = new File(saveDirectory, "qcReadCoverageTumor.xls");
		else txt = new File(saveDirectory, "qcReadCoverageNormal.xls");
		IO.writeString(sb.toString(), txt);
	}
	
	private void makePatterns() {
		alignLogPattern = Pattern.compile(alignLogMatch);
		dupLogPattern = Pattern.compile(dupLogMatch);
		readCovJsonPattern = Pattern.compile(readCovJsonMatch);
		scJsonPattern = Pattern.compile(scJsonMatch);
		aiJsonPattern = Pattern.compile(aiJsonMatch);
		normalDNAPattern = Pattern.compile(normalDNAMatch);
		tumorDNAPattern = Pattern.compile(tumorDNAMatch);
	}

	private void loadSamples() throws Exception {

			//for each job dir
			File[] jobDirs = IO.extractOnlyDirectories(jobDirectory);
			for (File job: jobDirs){
				System.out.println("\t"+ job.getName());

				File normAlignLog = null;
				File normDupLog = null;
				File normReadCovJson = null;

				File tumAlignLog = null;
				File tumDupLog = null;
				File tumReadCovJson = null;

				File sampConcJson = null;
				File clinInfoJson = null;
				File clinInfoTxt = null;

				for (File f: IO.fetchAllFilesRecursively(job)) {
					String name = f.getName();
					//align log
					if (alignLogPattern.matcher(name).matches()) {
						int nto = fetchSource(f);
						if (nto == 0) normAlignLog = f;
						else tumAlignLog = f;
					}

					//dup log
					else if (dupLogPattern.matcher(name).matches()) {
						int nto = fetchSource(f);
						if (nto == 0) normDupLog = f;
						else tumDupLog = f;
					}

					//read coverage json
					else if (readCovJsonPattern.matcher(name).matches()) {
						int nto = fetchSource(f);
						if (nto == 0) normReadCovJson = f;
						else tumReadCovJson = f;
					}

					//sample concordance json
					else if (scJsonPattern.matcher(name).matches()) sampConcJson =f;

					//clin info json
					else if (aiJsonPattern.matcher(name).matches()) clinInfoJson =f;
				}
				
				//look for clinInfo?
				if (clinInfoJson == null) {
					File cr = new File(job, "ClinicalReport");
					if (cr.exists()) {
						File[] files = IO.extractFiles(cr);
						if (files.length == 1) clinInfoTxt = files[0];
					}
				}
				
				boolean filesFound = (normAlignLog!=null || normDupLog!=null || normReadCovJson!=null || tumAlignLog!=null || tumDupLog!=null || tumReadCovJson!=null);
				if (filesFound) samples.add(new SampleQC2(job.getName(), normAlignLog, normDupLog, normReadCovJson, tumAlignLog, tumDupLog, tumReadCovJson, sampConcJson, clinInfoJson, clinInfoTxt));
				else IO.pl("\t\tNo files, skipping!");
			}
		}

	private int fetchSource(File f) throws IOException {
		
		//attempt to pull from file name
		String name = f.getName();
		if (normalDNAPattern.matcher(name).matches()) return 0;
		if (tumorDNAPattern.matcher(name).matches()) return 1;
		
		//attempt to pull from the path, not canonical path because some of these use soft links
		name = f.getPath();
		if (normalDNAPattern.matcher(name).matches()) return 0;
		if (tumorDNAPattern.matcher(name).matches()) return 1;
		
		IO.el("\nERROR: failed to parse the tumor DNA or normal DNA source from "+f.getPath());
		System.exit(1);
		return -1;
	}

	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new AggregateQCStats2(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'j': jobDirectory = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'a': alignLogMatch = args[++i]; break;
					case 'd': dupLogMatch = args[++i]; break;
					case 'r': readCovJsonMatch = args[++i]; break;
					case 'c': scJsonMatch = args[++i]; break;
					case 'i': aiJsonMatch = args[++i]; break;
					case 'n': normalDNAMatch = args[++i]; break;
					case 't': tumorDNAMatch = args[++i]; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//save directory
		if (saveDirectory == null) Misc.printErrAndExit("\nPlease provide a directory to save the results.\n");
		saveDirectory.mkdirs();

	}	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Aggregate QC Stats2: March 2022                        **\n" +
				"**************************************************************************************\n" +
				"Parses and aggregates alignment quality statistics from log and json files produced by\n"+
				"the TNRunner2 DnaAlignQC and SampleConcordance workflows.\n"+

				"\nRequired Options:\n"+
				"-j Directory containing TNRunner2 jobs.\n" +
				"-s Directory for saving the AggQC results.\n"+

				"\nOptions:\n"+
				"-a Alignment log file match, defaults to '.+_AlignHg38.log'\n"+
				"-d Mark duplicates log file match, defaults to '.+_Markdup.log'\n"+
				"-r Read coverage json file match, defaults to '.+_UniObRC.json.gz'\n"+
				"-c Sample concordance json file match, defaults to '.+_SampleConcordance.json.gz'\n"+
				"-i Clinical info json file match, defaults to '.+_AvatarInfo.json.gz'\n"+ 
				"-n Normal DNA file match, defaults to '.+NormalDNA.+'\n"+
				"-t Tumor DNA file match, defaults to '.+TumorDNA.+'\n"+
				"\n"+

				"Example: java -Xmx1G -jar pathToUSeq/Apps/AggregateQCStats2 -j TNJobs/ -s QCStats/\n\n" +

				"**************************************************************************************\n");

	}	
}
