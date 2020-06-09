package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import util.gen.IO;
import util.gen.Misc;

public class TSample {
	
	private String id;
	private File rootDir;
	private TRunner tRunner;
	private boolean forceRestart;
	private boolean restartFailedJobs;
	private boolean softRestart;
	private ArrayList<String> info = new ArrayList<String>();
	
	//Input datasets
	private File[] inputDnaBam = null;
	private File inputRnaBam = null;
	private File inputXmlReport = null;
	
	//Alignment and passing region bed files
	private File[] dnaBamBai = null;
	private File dnaAlignDir = null;
	
	private File[] rnaBamBai = null;
	private File rnaAlignDir = null;
	
	//Somatic variant Vcf calls
	private File somatiCallDir = null;
	private File[] somaticVcfTbi = null;
	
	//Xml to VCF
	private File xmlIntegrationDir = null;
	private File[] mergedVariants = null;
	
	//Annotations
	private File annoDir = null;
	
	//Job status
	private boolean failed = false;
	private boolean running = false;
	private int numJobsLaunched = 0;
	

	public TSample(File rootDir, TRunner TRunner) throws IOException{
		this.rootDir = rootDir;
		this.tRunner = TRunner;
		this.forceRestart = TRunner.isForceRestart();
		this.restartFailedJobs = TRunner.isRestartFailed();
		this.softRestart = TRunner.isSoftRestart();
		
		//assign id
		id = rootDir.getName();
		info.add("\n"+id);
		
		//look for Bam and Xml and check they contain bam and bai, some may be empty
		checkForBamsXml();
		
		//launch available alignments
		checkAlignments();
		
		//launch t/n somatic capture?
		callSomaticVariants();
		
		//parse Xml and compare
		parseIntegrateFoundationVcf();
		
		//annotate merged somatic vcf?
		annotateSomaticVariants();
		
		//print messages?
		if (failed || TRunner.isVerbose()) {
			Misc.printArray(info);
			IO.pl("Failures found? "+failed);
			IO.pl("Running jobs? "+running);
		}
	}	
		
	private void parseIntegrateFoundationVcf() throws IOException {
		if (somaticVcfTbi == null || inputXmlReport == null) return;
		info.add("Checking Xml parsing and somatic variant integration...");		
		
		//make dir, ok if it already exists
		xmlIntegrationDir = new File (rootDir, "SomaticVariantCalls/"+id+"_Merged");
		xmlIntegrationDir.mkdirs();
		
		File[] toSoftLink = new File[]{inputXmlReport, somaticVcfTbi[0], somaticVcfTbi[1]};
		
		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(xmlIntegrationDir);
		if (nameFile.size() == 0) launch(xmlIntegrationDir, toSoftLink, tRunner.getVcfIntegrationDocs());

		//OK some files are present
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			//find the final vcf file
			File[] vcf = IO.extractFiles(new File(xmlIntegrationDir, "Vcfs"), "_final.vcf.gz");
			File[] index = IO.extractFiles(new File(xmlIntegrationDir, "Vcfs"), "_final.vcf.gz.tbi");
			if (vcf == null || vcf.length !=1 || index == null || index.length != 1) {
				info.add("\t\tFAILED The parsed and integrated variant Xml extraction was marked COMPLETE but failed to find the final vcf or index file in the Vcfs/ in "+xmlIntegrationDir);
				failed = true;
				new File(xmlIntegrationDir, "FAILED").createNewFile();
				return;
			}
			//remove softlinks
			for (File f: toSoftLink) new File(xmlIntegrationDir, f.getName()).delete();
			info.add("\t\tCOMPLETE "+xmlIntegrationDir);
			//set vars for annotator
			mergedVariants = new File[]{vcf[0], index[0]};
		}
		else checkJob(nameFile, xmlIntegrationDir, toSoftLink, tRunner.getVcfIntegrationDocs());
	}
	
	private void checkJob(HashMap<String, File> nameFile, File jobDir, File[] toSoftLink, File[] runDocs) throws IOException {
		//force a restart?
		if (nameFile.containsKey("FAILED")){
			if (forceRestart  || restartFailedJobs ){
				//cancel any slurm jobs and delete the directory
				TNSample.cancelDeleteJobDir(nameFile, jobDir, info, true);
				restart(jobDir, toSoftLink, runDocs);
			}
			else if (softRestart){
				//cancel any slurm jobs and delete the directory
				TNSample.cancelDeleteJobDir(nameFile, jobDir, info, false);
				restart(jobDir, toSoftLink, runDocs);
			}
			//FAILED but no restart
			else if (nameFile.containsKey("FAILED")){
				info.add("\t\tFAILED "+jobDir);
				failed = true;
			}
		}
		//QUEUED
		else if (nameFile.containsKey("QUEUED")){
			info.add("\t\tQUEUED "+jobDir);
			running = true;
		}
		//STARTED
		else if (nameFile.containsKey("STARTED")) {
			if (TNSample.checkQueue(nameFile, jobDir, info) == false) failed = true;
			else running = true;
		}
		//hmm no status files, probably something went wrong on the cluster? mark it as FAILED
		else {
			info.add("\t\tMarking as FAILED, no job status files in "+jobDir);
			new File(jobDir, "FAILED").createNewFile();
			failed = true;
		}
	}
	
	private void restart(File jobDir, File[] toSoftLink, File[] runDocs) throws IOException{
		//launch it
		launch(jobDir, toSoftLink, runDocs);
		new File(jobDir, "RESTARTED").createNewFile();
		info.add("\t\tRESTARTED");
	}
	
	private void launch(File jobDir, File[] toLink, File[] docs) throws IOException{
		info.add("\t\tLaunching "+jobDir);
		running = true;
		IO.createSymbolicLinks(toLink, jobDir);
		//replace any launch scripts with the current
		File shellScript = TNSample.copyInWorkflowDocs(docs, jobDir);
		//clear any progress files
		TNSample.removeProgressFiles(jobDir);
		//squeue the shell script
		new File(jobDir, "QUEUED").createNewFile();
		String alignDirPath = jobDir.getCanonicalPath();
		String[] cmd = null;
		if (tRunner.isNice()) cmd = new String[]{"sbatch", "--nice=10000", "-J", alignDirPath.replace(tRunner.getPathToTrim(), ""), "-D", alignDirPath, shellScript.getCanonicalPath()};
		else cmd = new String[]{"sbatch", "-J", alignDirPath.replace(tRunner.getPathToTrim(), ""), "-D", alignDirPath, shellScript.getCanonicalPath()};
		String[] output = IO.executeViaProcessBuilder(cmd, false);
		
		for (String o: output) info.add("\t\t"+o);
		numJobsLaunched++;
	}

	private void annotateSomaticVariants() throws IOException {
		if (mergedVariants == null) return;
		info.add("Checking merged somatic variant annotation...");		
		
		//make dir, ok if it already exists
		annoDir = new File (rootDir, "SomaticVariantCalls/"+id+"_Anno");
		annoDir.mkdirs();
		
		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(annoDir);
		if (nameFile.size() == 0) {
			//write out config file
			IO.writeString(tRunner.getSomaticAnnotatedVcfParser(), new File(annoDir, "annotatedVcfParser.config.txt"));
			launch(annoDir, mergedVariants, tRunner.getVarAnnoDocs());
		}
		
		//OK some files are present
		else if (nameFile.containsKey("COMPLETE")){
			//remove the links
			for (File f: mergedVariants) new File(annoDir, f.getName()).delete();
			new File(annoDir, "annotatedVcfParser.config.txt").delete();
			info.add("\t\tCOMPLETE "+annoDir);
		}
		else {
			checkJob(nameFile, annoDir, mergedVariants, tRunner.getVarAnnoDocs());
			IO.writeString(tRunner.getSomaticAnnotatedVcfParser(), new File(annoDir, "annotatedVcfParser.config.txt"));
		}
	}

	private void callSomaticVariants() throws IOException {
		if (dnaBamBai == null) return;
		info.add("Checking somatic variant calling...");		
		
		//make dir, ok if it already exists
		somatiCallDir = new File (rootDir,"SomaticVariantCalls/"+id+"_Illumina");
		somatiCallDir.mkdirs();
		
		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(somatiCallDir);
		if (nameFile.size() == 0) launch(somatiCallDir,dnaBamBai, tRunner.getSomaticVarCallDocs());
				
		//OK some files are present
		else if (nameFile.containsKey("COMPLETE")){
			//find the final vcf file
			File[] vcf = IO.extractFiles(new File(somatiCallDir, "Vcfs"), "_final.vcf.gz");
			File[] index = IO.extractFiles(new File(somatiCallDir, "Vcfs"), "_final.vcf.gz.tbi");
			if (vcf == null || vcf.length !=1 || index == null || index.length != 1) {
				info.add("\t\tFAILED The somatic variant calling was marked COMPLETE but failed to find the final vcf or index file in the Vcfs/ in "+somatiCallDir);
				new File(xmlIntegrationDir, "FAILED").createNewFile();
				failed = true;
				return;
			}
			somaticVcfTbi = new File[]{vcf[0], index[0]};
			//remove softlinks
			for (File f: dnaBamBai) new File(somatiCallDir, f.getName()).delete();
			info.add("\t\tCOMPLETE "+somatiCallDir);
		}
		else checkJob(nameFile, somatiCallDir, dnaBamBai, tRunner.getSomaticVarCallDocs());
	}

	private void checkAlignments() throws IOException {
		info.add("Checking alignments...");
		if (inputDnaBam != null) dnaAlignQC();
		if (inputRnaBam != null) rnaAlignQC();
	}

	/**For RNASeq Alignments and Picard CollectRNASeqMetrics*/
	private void rnaAlignQC() throws IOException {
		if (inputRnaBam == null) return;
		info.add("\tTumorRNA:");
		
		//make dir, ok if it already exists
		rnaAlignDir = new File (rootDir, "Alignment/"+id+"_TumorRNA");
		rnaAlignDir.mkdirs();
		
		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(rnaAlignDir);
		if (nameFile.size() == 0) {
			launch(rnaAlignDir, new File[]{inputRnaBam}, tRunner.getTranscriptomeAlignQCDocs());
		}
		//OK some files are present
		else if (nameFile.containsKey("COMPLETE")){
			//find the final bam and bed
			File bamDir = new File(rnaAlignDir, "Bam");
			rnaBamBai = checkBam(bamDir);
			
			if (rnaBamBai == null || rnaBamBai.length !=2) {
				info.add("\t\tFAILED The alignemnt was marked COMPLETE but failed to find the bam or bai files in "+bamDir);
				failed = true;
				new File(rnaAlignDir, "FAILED").createNewFile();
				return;
			}
			//remove the linked bam
			new File(rnaAlignDir, inputRnaBam.getName()).delete();
			info.add("\t\tCOMPLETE "+rnaAlignDir);
			return;
		}
		else checkJob(nameFile, rnaAlignDir, new File[]{inputRnaBam}, tRunner.getTranscriptomeAlignQCDocs());
	}

	/**For Captures, diff read coverage for passing bed generation
	 * @throws IOException */
	private void dnaAlignQC() throws IOException {
		if (inputDnaBam == null) return;
		info.add("\tTumorDNA:");
		
		//make dir, ok if it already exists
		dnaAlignDir = new File (rootDir, "Alignment/"+id+"_TumorDNA");
		dnaAlignDir.mkdirs();
		
		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(dnaAlignDir);
		if (nameFile.size() == 0) launch(dnaAlignDir, inputDnaBam, tRunner.getCaptureAlignQCDocs());
		
		//OK some files are present
		else if (nameFile.containsKey("COMPLETE")){
			//find the final bam and bed
			File bamDir = new File(dnaAlignDir, "Bam");
			dnaBamBai = checkBam(bamDir);
			//check for problems
			if (dnaBamBai == null) {
				info.add("\tFAILED The alignemnt was marked COMPLETE but failed to find the final bam and bai files in the Bam/ in "+dnaAlignDir);
				failed = true;
				new File(rnaAlignDir, "FAILED").createNewFile();
				return;
			}
			//remove the linked files
			for (File f: inputDnaBam) new File(dnaAlignDir, f.getName()).delete();
			info.add("\t\tCOMPLETE "+dnaAlignDir);
		}
		else checkJob(nameFile, dnaAlignDir, inputDnaBam, tRunner.getCaptureAlignQCDocs());
	}
		
	private File[] checkBam(File dir) throws IOException{
		File[] b = IO.extractFiles(dir, ".bam");
		if (b != null && b.length == 1) {
			File i = new File(b[0]+".bai");			
			if (i.exists() == false) {
				i = new File(Misc.removeExtension(b[0].getCanonicalPath())+".bai");			
			}
			if (i.exists()) {
				return new File[]{b[0], i};
			}
		}
		return null;
	}

	private void checkForBamsXml() throws IOException {
		info.add("Checking Input Bam and Xml file availability...");
		File bamDir = makeCheckFile(rootDir, "Bam");
		if (bamDir != null){
			String path = null;
			//DNA
			File[] bamBai = checkBam(new File(bamDir, "TumorDNA"));
			if (bamBai != null) {
				inputDnaBam = bamBai;
				path = inputDnaBam[0].getCanonicalPath();
			}
			info.add("\tTumorDNA\t"+path);
			path = null;
			
			//RNA
			File rDir = new File(bamDir, "TumorRNA");
			if (rDir.exists()) {
				File[] bams = IO.extractFiles(rDir, ".bam");
				if (bams.length ==1) {
					inputRnaBam = bams[0];
					path = inputRnaBam.getCanonicalPath();
				}
			}
			info.add("\tTumorRNA\t"+path);
		}
		else info.add("\tNo Bam dir?");

		File xmlDir = makeCheckFile(rootDir, "Xml");
		if (xmlDir != null){
			File[] x = IO.extractFiles(xmlDir, ".xml");
			if (x.length == 1) inputXmlReport = x[0];
			info.add("\tXmlReport\t"+inputXmlReport);
		}
		else info.add("\tNo Xml dir?");
	}

	public static File makeCheckFile(File parentDir, String fileName) throws IOException {
		File f = new File(parentDir, fileName);
		if (f.exists() == false) return null;
		return f;
	}

	public File getRootDir() {
		return rootDir;
	}

	public String getId() {
		return id;
	}

	public boolean isFailed() {
		return failed;
	}

	public boolean isRunning() {
		return running;
	}

	public int getNumJobsLaunched() {
		return numJobsLaunched;
	}
}
