package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.nio.file.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

/**Runs multiple workflows on paired T N DNA and transcriptiome datasets.*/
public class TNRunner2 {

	//user defined fields
	private File sampleDir = null;
	private File[] rootDirs = null;
	private boolean verbose = true;
	private File[] DNAAlignQCDocs = null;
	private File[] RNAAlignQCDocs = null;
	private File[] RNAFusionDocs = null;
	private File[] illuminaSomaticVarCallDocs = null;
	private File[] gatkSomaticVarCallDocs = null;
	private File[] sampleConcordanceDocs = null;
	private File[] varAnnoDocs = null;
	private File[] haplotypeCallDocs = null;
	private File[] gatkJointGenotypingDocs = null;
	private File[] illuminaJointGenotypingDocs = null;
	private File[] copyRatioDocs = null;
	private File[] clinicalVcfDocs = null;
	private File[] msiDocs = null;
	private File[] lohDocs = null;
	private File normalAlignmentDir = null;
	private File oncoKBConfig = null;
	private File[] nonMatchedNormal = null;
	private File copyRatioBkgDir = null;
	private File[] copyRatioBkgHdf5Files = null;
	private File[] copyRatioBkgIntervalListFiles = null;
	private File bpileupFileOrDir = null;
	private int minReadCoverageTumor = 12;
	private int minReadCoverageNormal = 12;
	private TNSample2[] tNSamples = null;
	private HashMap<String, TNSample2> idSample = new HashMap<String, TNSample2>();
	private String germlineAnnotatedVcfParser = null;
	private String somaticAnnotatedVcfParser = null;
	private String germlineVcfCallFreq = "-b Hg38/Germline/Avatar/Bed -v Hg38/Germline/Avatar/Vcf";
	private String somaticVcfCallFreq = "-b Hg38/Somatic/Avatar/Bed -v Hg38/Somatic/Avatar/Vcf";
	private ArrayList<String> info = new ArrayList<String>();
	private boolean gatkGroupProcessingComplete = false;
	private boolean gatkGroupProcessingFailed = false;
	private boolean illGroupProcessingComplete = false;
	private boolean illGroupProcessingFailed = false;
	private boolean restartFailed = true;
	private int maxNumJobs = 35;
	private boolean loop = false;
	private int numMinToSleep = 60;
	private boolean niceJobs = true;
	private String partition = "hci-rw";
	private boolean sbatch = true;
	private HashSet<String> panels2SkipForCopyRatio = null;

	private String pathToTrim = null;
	private String nice = null;
	private String accountOverride = null;

	public TNRunner2 (String[] args) {
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);

			//loop?			
			for (int i=0; i< 1000; i++){
				//check samples
				if (processIndividualSamples()){
					processGATKSampleGroup();
					processIlluminaSampleGroup();
					//check germline processing
					if (complete()) {
						IO.pl("\nALL COMPLETE!");
						loop = false;
					}
					else IO.pl("\nNOT COMPLETE!");
				}
				//samples still running
				else IO.pl("\nNOT COMPLETE!");

				//sleep or exit
				if (loop) sleep(i);
				else i=1000;
			}

			printSamplesWithFastqIssues();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			IO.pl("\nDone! "+Math.round(diffTime)+" Min\n");
		} catch (IOException e) {
			IO.el("\nERROR running TNRunner, aborting");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	private void sleep(int i) {
		try {
			IO.pl("\nWaiting... "+i);
			Thread.sleep(1000*60*numMinToSleep);
		} catch (InterruptedException e) {}
	}

	private void printSamplesWithFastqIssues() {
		ArrayList<String> issues = new ArrayList<String>();
		for (TNSample2 t: tNSamples){
			if (t!= null && t.isFastqIssue()) issues.add(t.getId());
		}
		if (issues.size() != 0 && verbose) IO.pl("\nThe following have been skipped due to fastq file count issues: "+issues);
	}

	private boolean complete() {
		if (gatkGroupProcessingComplete == false || gatkGroupProcessingFailed == true || illGroupProcessingComplete == false || illGroupProcessingFailed == true) return false;
		for (TNSample2 t: tNSamples){
			if (t == null) continue;
			if (t.isFailed() || t.isRunning()) return false;
		}
		return true;
	}

	private void processGATKSampleGroup() {
		try {
			//any joint genotyping workflow files? If none then they want to skip this, e.g. no matched normal
			if (gatkJointGenotypingDocs == null) {
				gatkGroupProcessingFailed = false;
				gatkGroupProcessingComplete = true;
				return;
			}
			IO.pl("\nChecking GATK joint genotype group processing...");


			//for each sample are the haplotype gvcf calls ready?
			ArrayList<File> toGeno = new ArrayList<File>();
			boolean allPresent = true;
			for (TNSample2 tns: tNSamples){
				//JointGenotyping already run?
				File gatkDir = new File(tns.getRootDir(), "GermlineVariantCalling/"+tns.getId()+"_GATK");
				if (gatkDir.exists()) {
					if (verbose) IO.pl("\tComplete");
					gatkGroupProcessingComplete = true;
					return;
				}

				//OK doesn't exist so wait until all gvcfs become available
				//does this sample have a normal sample?
				if (tns.getNormalDNADataset() != null) {
					File gvcfDir = new File(tns.getRootDir(), "GermlineVariantCalling/"+tns.getId()+"_GVCF/Vcfs");
					File[] res = IO.extractFiles(gvcfDir, ".g.vcf.gz");
					if (res == null || res.length !=1) {
						if (verbose) IO.pl("\tWaiting for germline GATK gvcf from "+tns.getId());
						allPresent = false;
					}
					//OK it's present, add it if it's not the mock normal
					else {
						if (res[0].getName().contains("NA12878") == false) {
							//add the vcf and matched index
							toGeno.add(res[0]);
							toGeno.add(new File(res[0].getCanonicalPath()+".tbi"));
						}
					}
				}
			}
			if (allPresent == false) return;

			//launch it
			checkJointGenotyping(toGeno);
			if (gatkGroupProcessingFailed || verbose) Misc.printArray(info);

			if (verbose) {
				IO.pl("\tFailed? "+gatkGroupProcessingFailed);
				IO.pl("\tComplete? "+gatkGroupProcessingComplete);
			}

		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem processing group analysis!");
		}

	}

	private void processIlluminaSampleGroup() throws IOException {
		//any joint genotyping workflow files? If none then they want to skip this
		if (illuminaJointGenotypingDocs == null) {
			illGroupProcessingFailed = false;
			illGroupProcessingComplete = true;
			return;
		}
		IO.pl("\nChecking Illumina joint genotype group processing...");

		//look for completed dirs
		boolean alreadyLaunched = lookAndDistributeIlluminaJointGenotypedSamples();
		if (alreadyLaunched) return;
		
		//for each sample are the normal cram files ready?
		ArrayList<TNSample2> done = new ArrayList<TNSample2>();
		ArrayList<TNSample2> toRun = new ArrayList<TNSample2>();
		ArrayList<TNSample2> waitingForAlignment = new ArrayList<TNSample2>();
		ArrayList<TNSample2> noNormal = new ArrayList<TNSample2>();

		for (TNSample2 tns: tNSamples){
			//Already run?
			File complete = new File(tns.getRootDir(), "GermlineVariantCalling/"+tns.getId()+"_Illumina/COMPLETE");
			boolean alreadyRun = complete.exists();
			if (alreadyRun) {
				done.add(tns);
				continue;
			}
			//Alignment complete?
			else if (tns.getNormalDNADataset() != null && tns.getNormalDNADataset().isComplete()) {
				//is it the NA12878 control?
				if (nonMatchedNormal != null) {
					String nonMatchName = nonMatchedNormal[0].getName();
					String normName = tns.getNormalDNADataset().getCramFile().getName();
					if (normName.equals(nonMatchName) == false) toRun.add(tns);
				}
				else toRun.add(tns);
			}
			//Fastq present?
			else if (tns.getNormalDNAFastq()!= null && tns.getNormalDNAFastq().isGoodToAlign()) waitingForAlignment.add(tns);
			//Fastq and alignment not present so skip
			else  {
				noNormal.add(tns);
				continue;
			}
		}
		
		//IO.pl("\nScaned samples done:"+done.size()+" torun:"+toRun.size()+" waitingForAlign:"+waitingForAlignment.size()+" noNorm:"+noNormal.size());

		//any normals awaiting?
		int numWaiting = waitingForAlignment.size();
		if (numWaiting != 0) {
			if (verbose) {
				for (TNSample2 ts2: waitingForAlignment){
					IO.pl("\tWaiting for normal alignment from "+ts2.getId());
				}
			}
		}
		//all ready so build and launch
		else if (toRun.size() !=0) {
			jointGenotypeIlluminaSamples(toRun);
		}
		//nothing to do 
		else {
			illGroupProcessingFailed = false;
			illGroupProcessingComplete = true;
		}	
	}


	private void jointGenotypeIlluminaSamples(ArrayList<TNSample2> toRun) throws IOException {
		//create a IlluminaJointGenotyping_XXX dir for every 10 or less normals

		//split into chunks with a max of 10 samples
		TNSample2[] tns = new TNSample2[toRun.size()];
		toRun.toArray(tns);
		Misc.randomize(tns, 1);
		Object[][] chunks = Misc.chunkWithMax(tns, 10);

		//for each chunk make a dir, link in the normal crams to call germline vars on, submit to the queue
		for (int i=0; i< chunks.length; i++) {
			File jjDir = new File(sampleDir, "IlluminaJointGenotyping_"+i);
			jjDir.mkdir();
			File toGenoDir = new File(jjDir, "ToGenotype");
			toGenoDir.mkdir();
			Object[] samples = chunks[i];
			ArrayList<File> filesToLink = new ArrayList<File>();
			for (int j=0; j< samples.length; j++) {
				TNSample2 tn = (TNSample2) samples[j];
				filesToLink.add(tn.getNormalDNADataset().getCramFile());
				filesToLink.add(tn.getNormalDNADataset().getCramIndexFile());
			}
			launchIlluminaJointGenotyping(filesToLink, jjDir);
		}
		illGroupProcessingComplete = false;
		illGroupProcessingFailed = false;
	}

	private boolean lookAndDistributeIlluminaJointGenotypedSamples() throws IOException {

		ArrayList<File> jjDirs = IO.extractOnlyDirectories(sampleDir, "IlluminaJointGenotyping");
		if (jjDirs.size() == 0) {
			return false;
		}

		illGroupProcessingComplete = true;
		illGroupProcessingFailed = false;
		boolean checkSamples = false;
		
		//check each Joint Geno dir
		for (File jobDir: jjDirs) {
			if (verbose) IO.pl("\tChecking "+jobDir);
			//any files?
			HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);

			//COMPLETE
			if (nameFile.containsKey("COMPLETE")){
				IO.pl("\tCOMPLETE "+jobDir.getCanonicalPath());
				//find the final vcfs
				File res = new File(jobDir, "Vcfs");
				File[] genoVcfs = IO.extractFiles(res, "JointGenotyped.vcf.gz");

				//any files? might have already been moved
				if (genoVcfs.length != 0) {
					//OK just completed, need to distribute files to individual Jobs
					File logs = new File(jobDir, "Logs");
					File runScripts = new File(jobDir, "RunScripts");
					File zip = new File(jobDir, "LogRunScripts.zip");
					if (zip.exists()==false && (logs.exists() == false || runScripts.exists() == false)) throw new IOException("\t"+jobDir.getName()+ 
							" was marked COMPLETE but failed to find the Logs or RunScripts directories or LogRunScripts.zip archive.");

					//for each separate genotyped vcf
					for (int i=0; i< genoVcfs.length; i++){

						//Extract the sample name, JZBCLZLXL6-IDT-SL409940-SL413507-SL416453_NormalDNA_Hg38_JointGenotyped.vcf.gz
						String fileName = genoVcfs[i].getName();			
						String sampleName = fileName.replace("_NormalDNA_Hg38_JointGenotyped.vcf.gz", "");
						TNSample2 tns = idSample.get(sampleName);
						if (tns == null) throw new IOException("\nERROR: Failed to find the sample dataset by parsing the JointGenotyping vcf names for "+fileName +" -> "+sampleName);

						//create a folder for the incoming genotyped vcf in the sample folder
						File resDir = new File(tns.getRootDir(), "GermlineVariantCalling/"+sampleName+"_Illumina");
						resDir.mkdirs();

						if (verbose) IO.pl("\tMoving genotyped vcfs to "+resDir);

						//copy and move over relevant files
						if (zip.exists()) {
							File newZip = new File(resDir,"LogRunScripts.zip");
							IO.copyViaFileChannel(zip, newZip);
						}
						else {
							File newLogDir = new File(resDir,"Logs");
							newLogDir.mkdir();
							File newRSDir = new File(resDir,"RunScripts");
							newRSDir.mkdir();
							IO.copyDirectoryRecursive(logs, newLogDir, null);
							IO.copyDirectoryRecursive(runScripts, newRSDir, null);
						}
						
						genoVcfs[i].renameTo(new File(resDir, genoVcfs[i].getName()));
						File tbi = new File (genoVcfs[i].getCanonicalPath()+".tbi");
						tbi.renameTo(new File(resDir, tbi.getName()));
						File complete = new File(resDir, "COMPLETE");
						complete.createNewFile();
						
					}
				}
				checkSamples = true;
			}

			//force a restart?
			else if (restartFailed && nameFile.containsKey("FAILED")){
				//cancel any slurm jobs and delete the directory
				TNSample.cancelDeleteJobDir(nameFile, jobDir, info, true);
				illGroupProcessingComplete = false;
				illGroupProcessingFailed = true;
			}
			//QUEUED
			else if (nameFile.containsKey("QUEUED")){
				info.add("\tQUEUED "+jobDir);
				illGroupProcessingComplete = false;
			}	
			//RUNME
			else if (nameFile.containsKey("RUNME")){
				info.add("\tRUNME "+jobDir);
				illGroupProcessingComplete = false;
			}
			//STARTED
			else if (nameFile.containsKey("STARTED")){
				if (TNSample.checkQueue(nameFile, jobDir, info) == false) illGroupProcessingFailed = true;
				illGroupProcessingComplete = false;
			}
			//FAILED but no forceRestart
			else if (nameFile.containsKey("FAILED")){
				if (verbose) IO.pl("\tFAILED "+jobDir);
				illGroupProcessingFailed = true;
				illGroupProcessingComplete = false;
			}
			//hmm no status files, probably something went wrong on the cluster? mark it as FAILED
			else {
				if (verbose) IO.pl("\tFAILED?, no job status files in "+jobDir);
				illGroupProcessingFailed = true;
				illGroupProcessingComplete = false;
			}
			//check samples to run annotator?
			//if (checkSamples) processIndividualSamples();
		}
		if (verbose) {
			IO.pl("\tAll complete? "+illGroupProcessingComplete);
			IO.pl("\tAny failed? "+illGroupProcessingFailed);
		}
		return true;
	}

	private void checkJointGenotyping(ArrayList<File> toGeno) throws IOException {

		//make dir, ok if it already exists
		File jobDir = new File (sampleDir, "GATKJointGenotyping");
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			launchGatkJointGenotyping(toGeno, jobDir);
			return;
		}

		//OK some files are present
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			IO.pl("\tCOMPLETE "+jobDir.getCanonicalPath());
			//find the final vcfs
			File res = new File(jobDir, "Vcfs");
			File[] genoVcfs = IO.extractFiles(res, "JointGenotyped.vcf.gz");
			//any files? might have already been moved
			if (genoVcfs.length == 0) {
				gatkGroupProcessingComplete = true;
				return;
			}

			//OK just completed, need to distribute files to individual Sample dirs
			File[] genoVcfIndexes = IO.extractFiles(res, "JointGenotyped.vcf.gz.tbi");
		
			
			File logs = new File(jobDir, "Logs");
			File runScripts = new File(jobDir, "RunScripts");
			File zip = new File(jobDir, "LogRunScripts.zip");
			if (zip.exists()==false && (logs.exists() == false || runScripts.exists() == false)) throw new IOException("\t"+jobDir.getName()+ 
					" was marked COMPLETE but failed to find the Logs or RunScripts directories or LogRunScripts.zip archive.");

			//for each separate genotyped vcf
			for (int i=0; i< genoVcfs.length; i++){
				//Extract the sample name 
				//03GCDX0K5N-IDT-SL416271-SL419057-SL420350_JointGenotyped.vcf.gz

				String fileName = genoVcfs[i].getName();			
				String sampleName = fileName.replace("_JointGenotyped.vcf.gz", "");
				if (sampleName.length() == fileName.length()) throw new IOException("\nERROR: Failed to parse the sample name from "+fileName);

				File sampleDir = null;
				String id = null;

				//find the sample
				TNSample2 toLaunch = null;
				for (TNSample2 tns: tNSamples){
					id = tns.getId();			
					if (id.equals(sampleName)){
						sampleDir = tns.getRootDir().getCanonicalFile();
						toLaunch = tns;
						break;
					}
				}
				if (sampleDir == null) throw new IOException("\nERROR: Failed to find the sampleDir by parsing the JointGenotyping vcf names. "+fileName);

				if (verbose) IO.pl("\tMoving genotyped vcfs to "+sampleDir);


				//create a folder for the incoming genotyped vcf in the sample folder
				File resDir = new File(sampleDir, "GermlineVariantCalling/"+sampleName+"_GATK");
				resDir.mkdirs();

				//copy and move over relevant files
				if (zip.exists()) {
					File newZip = new File(resDir,"LogRunScripts.zip");
					IO.copyViaFileChannel(zip, newZip);
				}
				else {
					File newLogDir = new File(resDir,"Logs");
					newLogDir.mkdir();
					File newRSDir = new File(resDir,"RunScripts");
					newRSDir.mkdir();
					IO.copyDirectoryRecursive(logs, newLogDir, null);
					IO.copyDirectoryRecursive(runScripts, newRSDir, null);
				}

				genoVcfs[i].renameTo(new File(resDir, genoVcfs[i].getName()));
				genoVcfIndexes[i].renameTo(new File(resDir, genoVcfIndexes[i].getName()));

				//launch germline annotator?
				if (getVarAnnoDocs()!= null) toLaunch.annotateGermlineVcf("GATK");
			}
			gatkGroupProcessingComplete = true;
		}

		//force a restart?
		else if (restartFailed && nameFile.containsKey("FAILED")){
			//cancel any slurm jobs and delete the directory
			TNSample.cancelDeleteJobDir(nameFile, jobDir, info, true);
			//launch it
			launchGatkJointGenotyping(toGeno, jobDir);
		}
		//QUEUED
		else if (nameFile.containsKey("QUEUED")){
			info.add("\tQUEUED "+jobDir);
		}	
		//RUNME
		else if (nameFile.containsKey("RUNME")){
			info.add("\tRUNME "+jobDir);
		}
		//STARTED
		else if (nameFile.containsKey("STARTED")){
			if (TNSample.checkQueue(nameFile, jobDir, info) == false) gatkGroupProcessingFailed = true;
		}
		//FAILED but no forceRestart
		else if (nameFile.containsKey("FAILED")){
			if (verbose) IO.pl("\tFAILED "+jobDir);
			gatkGroupProcessingFailed = true;
		}
		//hmm no status files, probably something went wrong on the cluster? mark it as FAILED
		else {
			if (verbose) IO.pl("\tMarking as FAILED, no job status files in "+jobDir);
			new File(jobDir, "FAILED").createNewFile();
			gatkGroupProcessingFailed = true;
		}
	}

	private void launchGatkJointGenotyping(ArrayList<File> toGeno, File jobDir) throws IOException {
		if (verbose) {
			if (sbatch) info.add("\tLAUNCHING "+jobDir);
			else info.add("\tSETTING UP "+jobDir);
		}
		
		//want to create/ replace the soft links
		createJointGenotypingLinks(toGeno, jobDir);

		//replace any launch scripts with the current
		File shellScript = TNSample.copyInWorkflowDocs(gatkJointGenotypingDocs, jobDir);

		//clear any progress files
		TNSample.removeProgressFiles(jobDir);

		//squeue the shell script
		String alignDirPath = jobDir.getCanonicalPath();
		String [] cmd = null;
		if (nice != null) cmd = new String[]{"sbatch", nice, "-J", alignDirPath.replace(pathToTrim, ""), "-D", alignDirPath, shellScript.getCanonicalPath()};
		else cmd = new String[]{"sbatch", "-J", alignDirPath.replace(pathToTrim, ""), "-D", alignDirPath, shellScript.getCanonicalPath()};
		if (sbatch) {
			String[] output = IO.executeViaProcessBuilder(cmd, false);
			if (verbose) for (String o: output) IO.pl("\t\t"+o);
			for (String o: output) {
				if (o.toLowerCase().contains("error")) {
					gatkGroupProcessingFailed = true;
					IO.pl("\tERROR launching "+Misc.stringArrayToString(cmd, " "));
					if (verbose == false) for (String x: output) IO.pl("\t\t"+x);
					new File(jobDir, "FAILED").createNewFile();
					return;
				}
			}
			new File(jobDir, "QUEUED").createNewFile();
		}
		else new File(jobDir, "RUNME").createNewFile();
	}

	private void launchIlluminaJointGenotyping(ArrayList<File> toGeno, File jobDir) throws IOException {
		if (verbose) {
			if (sbatch) info.add("\tLAUNCHING "+jobDir);
			else info.add("\tSETTING UP "+jobDir);
		}
		
		//want to create/ replace the soft links
		createJointGenotypingLinks(toGeno, jobDir);

		//replace any launch scripts with the current
		File shellScript = TNSample.copyInWorkflowDocs(illuminaJointGenotypingDocs, jobDir);

		//clear any progress files
		TNSample.removeProgressFiles(jobDir);

		//squeue the shell script
		String [] cmd = buildSBatchCommand (jobDir, shellScript);
		if (sbatch) {
			String[] output = IO.executeViaProcessBuilder(cmd, false);
			if (verbose) for (String o: output) IO.pl("\t\t"+o);
			for (String o: output) {
				if (o.toLowerCase().contains("error")) {
					illGroupProcessingFailed = true;
					IO.pl("\tERROR launching "+Misc.stringArrayToString(cmd, " "));
					if (verbose == false) for (String x: output) IO.pl("\t\t"+x);
					new File(jobDir, "FAILED").createNewFile();
					return;
				}
			}
			new File(jobDir, "QUEUED").createNewFile();
		}
		else new File(jobDir, "RUNME").createNewFile();

	}
	
	String[] buildSBatchCommand (File jobDir, File shellScript) throws IOException {
		ArrayList<String> cmds = new ArrayList<String>();
		String jobDirPath = jobDir.getCanonicalPath();
		//{"sbatch", nice, accountOverride "-J", alignDirPath.replace(pathToTrim, ""), "-D", alignDirPath, shellScript.getCanonicalPath()}
		cmds.add("sbatch");
		if (nice != null) cmds.add(nice);
		if (accountOverride != null) cmds.add(accountOverride);
		cmds.add("-J");
		cmds.add(jobDirPath.replace(pathToTrim, ""));
		cmds.add("-D");
		cmds.add(jobDirPath);
		cmds.add(shellScript.getCanonicalPath());
		String[] cmd = new String[cmds.size()];
		cmds.toArray(cmd);		
		return cmd;
	}

	private static int countNumberRunningJobs(String partition) throws IOException {
		String [] cmd = new String[]{"squeue", "-p", partition};
		String[] output = IO.executeViaProcessBuilder(cmd, false);
		return output.length -1;
	}


	private void createJointGenotypingLinks(ArrayList<File> toGeno, File jobDir) throws IOException {
		File linkDir = new File(jobDir.getCanonicalFile(), "ToGenotype");
		if (linkDir.exists()) IO.deleteDirectory(linkDir);
		linkDir.mkdir();

		for (File f: toGeno){
			File link = new File(linkDir, f.getName());
			Files.createSymbolicLink(link.toPath(), f.toPath());
		}
	}

	private boolean processIndividualSamples() {
		try {
			if (verbose) IO.pl("\nChecking individual samples...");
			else IO.pl("\nChecking individual samples (SampleID RunningJobs)");

			int numJobsLaunched = countNumberRunningJobs(partition);
			tNSamples = new TNSample2[rootDirs.length];
			for (int i=0; i< rootDirs.length; i++){
				if (sbatch == true && numJobsLaunched >= maxNumJobs) {
					IO.pl("\nMaximum number jobs launched, skipping remaining samples.");
					return false;
				}
				if (verbose == false) IO.p("\t"+rootDirs[i].getName());
				tNSamples[i] = new TNSample2(rootDirs[i].getCanonicalFile(), this);
				idSample.put(tNSamples[i].getId(), tNSamples[i]);
				if (verbose == false) IO.pl("\t"+tNSamples[i].isRunning());
				numJobsLaunched += tNSamples[i].getNumJobsLaunched();
			}
			return true;
		} catch (IOException e) {
			System.err.println("\n\nProblem processing individual sample:");
			e.printStackTrace();
			return false;
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new TNRunner2(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		try {
			IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
			Pattern pat = Pattern.compile("-[a-zA-Z]");
			File DNAWorkflowDir = null;
			File illSomVarCallWorkflowDir = null;
			File gatkSomVarCallWorkflowDir = null;
			File annoWorkflowDir = null;
			File sampleConWorkflowDir = null;
			File gatkJointGenoWorklfowDir = null;
			File illuminaJointGenoWorklfowDir = null;
			File haploWorklfowDir = null;
			File RNAWorkflowDir = null;
			File rnaFuseDir = null;
			File copyRatioDocsDir = null;
			File clinicalVcfDir = null;
			File msiWorkflowDir = null;
			File lohWorkflowDir = null;
			String panel2Skip = null;
			for (int i = 0; i<args.length; i++){
				Matcher mat = pat.matcher(args[i]);
				if (mat.matches()){
					char test = args[i].charAt(1);
					try{
						switch (test){
						case 'p': sampleDir = new File(args[++i]).getCanonicalFile(); break;
						case 'e': DNAWorkflowDir = new File(args[++i]); break;
						case 't': RNAWorkflowDir = new File(args[++i]); break;
						case 'f': rnaFuseDir = new File(args[++i]); break;
						case 'c': illSomVarCallWorkflowDir = new File(args[++i]); break;
						case 'C': gatkSomVarCallWorkflowDir = new File(args[++i]); break;
						case 'a': annoWorkflowDir = new File(args[++i]); break;
						case 'b': sampleConWorkflowDir = new File(args[++i]); break;
						case 'm': msiWorkflowDir = new File(args[++i]); break;
						case 'L': lohWorkflowDir = new File(args[++i]); break;
						case 'j': gatkJointGenoWorklfowDir = new File(args[++i]); break;
						case 'q': illuminaJointGenoWorklfowDir = new File(args[++i]); break;
						case 'h': haploWorklfowDir = new File(args[++i]); break;
						case 'y': copyRatioDocsDir = new File(args[++i]); break;
						case 'k': copyRatioBkgDir = new File(args[++i]); break;
						case 'o': oncoKBConfig = new File(args[++i]); break;
						case 'P': panel2Skip = args[++i]; break;
						case 'v': clinicalVcfDir = new File(args[++i]); break;
						case 'B': bpileupFileOrDir = new File(args[++i]); break;
						case 'w': normalAlignmentDir = new File(args[++i]); break;
						case 'g': germlineAnnotatedVcfParser = args[++i]; break;
						case 's': somaticAnnotatedVcfParser = args[++i]; break;
						case 'G': germlineVcfCallFreq = args[++i]; break;
						case 'S': somaticVcfCallFreq = args[++i]; break;
						case 'u': minReadCoverageTumor = Integer.parseInt(args[++i]); break;
						case 'i': minReadCoverageNormal = Integer.parseInt(args[++i]); break;
						case 'x': maxNumJobs = Integer.parseInt(args[++i]); break;
						case 'z': niceJobs = false; break;
						case 'J': sbatch = false; break;
						case 'l': loop = true; break;
						case 'd': restartFailed = false; break;
						case 'n': partition = args[++i]; break;
						case 'A': accountOverride = args[++i]; break;
						default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
						}
					}
					catch (Exception e){
						e.printStackTrace();
						Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
					}
				}
			}

			//root patient dirs? Looks for Fastq dirs then pulls their parent folder.
			if (sampleDir == null || sampleDir.exists() == false) Misc.printErrAndExit("Error: failed to find your starting data directory? "+sampleDir);
			
			ArrayList<File> fastqDirs = IO.fetchDirectoriesRecursively(sampleDir, "Fastq");
			rootDirs = IO.fetchParentDirectories(fastqDirs);
			rootDirs = removeThoseWithComplete(rootDirs);
			if (rootDirs == null || rootDirs.length == 0) Misc.printExit("WARNING: failed to find any job directories (those containing a 'Fastq' dir or those without a COMPLETE) to process? "+rootDirs);

			//remove JointGenotyping if present
			ArrayList<File> toKeep = new ArrayList<File>();
			for (File rd: rootDirs){
				String name = rd.getName();
				if (name.contains("JointGenotyping") == false) toKeep.add(rd);
			}
			rootDirs = new File[toKeep.size()];
			toKeep.toArray(rootDirs);

			//DNA align qc docs
			if (DNAWorkflowDir != null){
				if (DNAWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing DNA alignment workflow docs? "+DNAWorkflowDir);
				DNAAlignQCDocs = IO.extractFiles(DNAWorkflowDir);
			}

			//RNA align qc docs
			if (RNAWorkflowDir != null) {
				if (RNAWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing RNA alignment workflow docs? "+RNAWorkflowDir);
				RNAAlignQCDocs = IO.extractFiles(RNAWorkflowDir);
			}

			//RNA fusion docs
			if (rnaFuseDir != null) {
				if (rnaFuseDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing RNA fusion workflow docs? "+rnaFuseDir);
				RNAFusionDocs = IO.extractFiles(rnaFuseDir);
			}

			//variant calling docs using Manta and Strelka
			if (illSomVarCallWorkflowDir != null){
				if(illSomVarCallWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing Illumina somatic variant calling workflow docs? "+illSomVarCallWorkflowDir);
				illuminaSomaticVarCallDocs = IO.extractFiles(illSomVarCallWorkflowDir);
				if (bpileupFileOrDir == null || bpileupFileOrDir.exists() == false) Misc.printErrAndExit("Error: failed to find the bpileup file or dir for the Illumina somatic variant calling? "+bpileupFileOrDir);
				bpileupFileOrDir = bpileupFileOrDir.getCanonicalFile();
			}
			
			//variant calling docs using Mutect2 GATK
			if (gatkSomVarCallWorkflowDir != null){
				if(gatkSomVarCallWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing GATK somatic variant calling workflow docs? "+gatkSomVarCallWorkflowDir);
				gatkSomaticVarCallDocs = IO.extractFiles(gatkSomVarCallWorkflowDir);
			}

			//variant annotation
			if (annoWorkflowDir != null){
				if (annoWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing variant annotation workflow docs? "+annoWorkflowDir);
				varAnnoDocs = IO.extractFiles(annoWorkflowDir);
				if (oncoKBConfig == null || oncoKBConfig.exists()==false) Misc.printErrAndExit("Error: failed to find an OncoKB configuration file, see annotator.README.sh. "+ oncoKBConfig);
			}

			//sample concordance
			if (sampleConWorkflowDir != null){
				if (sampleConWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing sample concordance workflow docs? "+sampleConWorkflowDir);
				sampleConcordanceDocs = IO.extractFiles(sampleConWorkflowDir);
			}

			//haplotype calling
			if (haploWorklfowDir != null){
				if (haploWorklfowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing GATK haplotype calling workflow docs? "+haploWorklfowDir);
				haplotypeCallDocs = IO.extractFiles(haploWorklfowDir);
			}

			//gatk joint genotyping
			if (gatkJointGenoWorklfowDir != null){
				if (gatkJointGenoWorklfowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing GATK joint genotyping workflow docs? "+gatkJointGenoWorklfowDir);
				gatkJointGenotypingDocs = IO.extractFiles(gatkJointGenoWorklfowDir);
			}

			//strelka joint genotyping
			if (illuminaJointGenoWorklfowDir != null){
				if (illuminaJointGenoWorklfowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing strelka joint genotyping workflow docs? "+illuminaJointGenoWorklfowDir);
				illuminaJointGenotypingDocs = IO.extractFiles(illuminaJointGenoWorklfowDir);
			}

			//msi
			if (msiWorkflowDir != null){
				if (msiWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing msi workflow docs? "+msiWorkflowDir);
				msiDocs = IO.extractFiles(msiWorkflowDir);
			}
			
			//loh
			if (lohWorkflowDir != null){
				if (lohWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing LoH workflow docs? "+lohWorkflowDir);
				lohDocs = IO.extractFiles(lohWorkflowDir);
			}

			//clinical vcf merging workflow
			if (clinicalVcfDir != null){
				if (clinicalVcfDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing workflow docs for merging and comparing clinical test reports and vcfs? "+clinicalVcfDocs);
				clinicalVcfDocs = IO.extractFiles(clinicalVcfDir);
			}

			//copy ratio analysis?
			if (copyRatioDocsDir !=null){
				if (copyRatioDocsDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing copy ratio analysis workflow docs? "+copyRatioDocsDir);
				copyRatioDocs = IO.extractFiles(copyRatioDocsDir);
				//copy ratio analysis
				if (copyRatioBkgDir == null || copyRatioBkgDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing copy ratio background files? "+copyRatioBkgDir);
				extractCopyRtoFiles(copyRatioBkgDir);
				//any panels 2 skip?
				if (panel2Skip!= null) {
					panels2SkipForCopyRatio = new HashSet<String>();
					for (String x: Misc.COMMA.split(panel2Skip)) panels2SkipForCopyRatio.add(x);
				}
			}

			//non matched normal?
			if (normalAlignmentDir != null) {
				//need to generate cram files from these and rename the Pass to PassRC
				//find the final bam and bed
				File[] finalBam = IO.extractFiles(new File(normalAlignmentDir, "Alignment"), ".cram");
				File[] passingBed = IO.extractFiles(new File(normalAlignmentDir, "QC"), "PassRC.bed.gz");
				if (finalBam == null || finalBam.length !=1 )throw new IOException("\nFailed to find a single Alignment/xxx.cram in "+normalAlignmentDir);
				if (passingBed == null || passingBed.length !=1 )throw new IOException("\nFailed to find a single QC/xxxPassRC.bed.gz in "+normalAlignmentDir);
				nonMatchedNormal = new File[] {finalBam[0], passingBed[0]};
			}


			//AnnotatedVcfParser options for germline and somatic filtering
			if (germlineAnnotatedVcfParser == null) germlineAnnotatedVcfParser = "-d "+minReadCoverageNormal+" -m 0.075 -q 0.1 -p 0.01 -g D5S,D3S -n 4.4 -a HIGH -l -c "
					+ "Pathogenic,Likely_pathogenic,Conflicting_interpretations_of_pathogenicity,Drug_response -t 0.51 -e Benign,Likely_benign -o -b 0.1 -z 3 -u RYR1";
			if (somaticAnnotatedVcfParser == null) somaticAnnotatedVcfParser = "-d "+minReadCoverageTumor+" -f -m 0.025 -j ";

			//set path to trim
			pathToTrim = sampleDir.getCanonicalFile().getParentFile().getCanonicalPath()+"/";

			if (niceJobs) nice = "--nice=10000";
			
			if (accountOverride!=null) accountOverride = "--account="+accountOverride;

			if (verbose){
				IO.pl("Run parameters:");
				IO.pl("Patient sample directory\t"+rootDirs[0].getParent());
				IO.pl("DNA alignment workflow directory\t"+DNAWorkflowDir);
				IO.pl("RNA alignment workflow directory\t"+RNAWorkflowDir);
				IO.pl("Illumina somatic variant workflow directory\t"+illSomVarCallWorkflowDir);
				IO.pl("GATK somatic variant workflow directory\t"+ gatkSomVarCallWorkflowDir);
				if (normalAlignmentDir != null) IO.pl("Non matched normal alignment directory for somatic calling\t"+normalAlignmentDir);
				IO.pl("Variant annotation workflow directory\t"+annoWorkflowDir);
				if (annoWorkflowDir!=null) IO.pl("OncoKB configuration file\t"+ oncoKBConfig);
				IO.pl("Sample concordance workflow directory\t"+sampleConWorkflowDir);
				IO.pl("MSI workflow directory\t"+msiWorkflowDir);
				IO.pl("LoH workflow directory\t"+lohWorkflowDir);
				IO.pl("GATK haplotype calling workflow directory\t"+haploWorklfowDir);
				IO.pl("GATK Joint genotyping workflow directory\t"+gatkJointGenoWorklfowDir);
				IO.pl("Illumina Joint genotyping workflow directory\t"+illuminaJointGenoWorklfowDir);
				IO.pl("Copy ratio workflow directory\t"+copyRatioDocsDir);
				IO.pl("Copy ratio background directory\t"+copyRatioBkgDir);
				if (panels2SkipForCopyRatio != null) IO.pl("Copy ratio panels to skip\t"+panels2SkipForCopyRatio);

				IO.pl("Min tumor read coverage\t"+minReadCoverageTumor);
				IO.pl("Min normal read coverage\t"+minReadCoverageNormal);
				IO.pl("Germline Annotated Vcf Parser options\t"+germlineAnnotatedVcfParser);
				IO.pl("Somatic Annotated Vcf Parser options\t"+somaticAnnotatedVcfParser);

				IO.pl("Delete and restart failed jobs\t"+restartFailed);
				IO.pl("Just set up jobs, don't launch them\t"+ (sbatch==false));
				IO.pl("Verbose logging\t"+verbose);
				IO.pl("Max # jobs to launch\t"+maxNumJobs);
				IO.pl("Nice jobs\t"+niceJobs);
				IO.pl("Job Partition\t"+partition);
				IO.pl("Override cluster account\t"+accountOverride);
				IO.pl("Relaunch jobs until complete\t"+loop);
			}
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing arguments for TNRunner, aborting!");
		}
	}


	private void extractCopyRtoFiles(File dir) {
		ArrayList<File> hdf5 = new ArrayList<File>();
		ArrayList<File> il = new ArrayList<File>();
		File[] toCheck = IO.extractOnlyFiles(dir);
		//xxx.hdf5 and xxx.interval_list
		for (File f: toCheck) {
			String name = f.getName();
			if (name.endsWith(".hdf5")) hdf5.add(f);
			else if (name.endsWith(".interval_list")) il.add(f);
		}
		//minimum of one each!
		if (hdf5.size()< 1 || il.size()< 1) Misc.printErrAndExit("Error: failed to find at least one set of copy ratio bkgd files (xxx.hdf5 and xxx.interval_list in "+copyRatioBkgDir);
		copyRatioBkgHdf5Files = new File[hdf5.size()];
		hdf5.toArray(copyRatioBkgHdf5Files);
		copyRatioBkgIntervalListFiles = new File[il.size()];
		il.toArray(copyRatioBkgIntervalListFiles);
	}

	private File[] removeThoseWithComplete(File[] dirs) {
		IO.pl("\nChecking for job dirs for those marked COMPLETE...");
		ArrayList<File> toReturn = new ArrayList<File>();
		for (File d: dirs) {
			if (new File(d,"COMPLETE").exists()== false) toReturn.add(d);
			else IO.pl("\tSkipping "+d.getName());
		}
		File[] f = new File[toReturn.size()];
		toReturn.toArray(f);
		return f;
	}

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                                TNRunner2 : March 2025                            **\n" +
				"**************************************************************************************\n" +
				"TNRunner2 is designed to execute several containerized workflows on tumor normal\n"+
				"datasets via a slurm cluster.  Based on the availability of paired fastq datasets, \n"+
				"alignments are run, somatic and germline variants are called, and concordance measured\n"+
				"between sample crams. To execute TNRunner, create the following directory structure,\n"+
				"link or copy in the corresponding paired end Illumina gzipped fastq files. To run \n"+
				"GATK Joint Genotyping on sample sets where this has been run before, delete the\n "+
				"GermlineVariantCalling/GATK* and GATKJointGenotyping folders.\n"+

				"\nMyPatientSampleDatasets\n"+
				"   MyPatientA\n"+
				"      Fastq\n"+
				"         NormalDNA\n"+
				"         TumorDNA\n"+
				"         TumorRNA\n"+
				"   MyPatientB\n"+
				"      Fastq\n"+
				"         TumorDNA\n"+
				"         TumorRNA\n"+
				"   MyPatientC....\n"+

				"\nThe Fastq directory and sub directories must match this naming. Only include those\n"+
				"for which you have fastq.  Change the MyXXX to something relevant. TNRunner is\n"+
				"stateless so as more Fastq becomes available or issues are addressed, relaunch the\n"+
				"app. This won't effect running or queued slurm jobs. Relaunch periodically to assess\n"+
				"the current processing status and queue additional tasks or set option -l. Download \n"+
				"the latest workflows from \n"+
				"https://github.com/HuntsmanCancerInstitute/Workflows/tree/master/Hg38RunnerWorkflows \n"+
				"and the matching resource bundle from\n"+
				"https://hci-bio-app.hci.utah.edu/gnomex/gnomexFlex.jsp?analysisNumber=A5578 .\n"+
				"All workflow docs are optional although some require output from prior Analysis.\n"+

				"\nSetup Options:\n"+
				"-p Directory containing one or more patient data directories with a 'Fastq' folder.\n" +
				"      TNRunner looks for these recursively. Those patient dirs with a file named \n"+
				"      COMPLETE will be skipped.\n"+
				"-e Workflow docs for launching DNA alignments.\n"+
				"-t Workflow docs for launching RNA alignments.\n"+
				"-f Workflow docs for launching RNA fusion detection.\n"+
				"-c Workflow docs for launching Illumina somatic variant calling.\n"+
				"-C Workflow docs for launching GATK somatic variant calling.\n"+
				"-m Workflow docs for launching MSI status calling.\n"+
				"-a Workflow docs for launching variant annotation.\n"+
				"-o If -a provide an OncoKB configuration file, see the annotator.README.sh\n"+
				"-b Workflow docs for launching sample concordance.\n"+
				"-y Workflow docs for launching somatic copy analysis.\n"+
				"-q Workflow docs for launching Illumina germline joint genotyping.\n"+
				"-h Workflow docs for launching GATK haplotype calling.\n"+
				"-j Workflow docs for launching GATK germline joint genotyping.\n"+
				"-L Workflow docs for launching LoH analysis. Requires an annotated germline and\n"+
				"      copy ratio analysis.\n"+
				"-v Workflow docs for launching clinical test variant info. Add a ClinicalReport folder \n"+
				"      to each patient dir containing the json formatted clinical information. The file name\n"+
				"      should contain the unique panel name and end with the sex if running copy ratio analysis,\n"+
				"      e.g. TL-19-9DA8E2_XE.V2_M.json, \n"+
				"-k Directory containing the sex and platform matched xxx.hdf5 and xxx.interval_list files\n"+
				"      for somatic copy ratio analysis.\n"+
				"      files with unique panel names, e.g. XE.V2_M.hdf5 and XE.V2_F.hdf5\n"+
				"-P Comma delimited string of panels to skip for copy ratio analysis, no spaces.\n"+
				"-B BamPileup file or directory containing such for somatic variant VCFBkz scoring,\n"+
				"      These should end in xxx.bp.txt.gz with paired tbi indexes and contain unique\n"+
				"      panel names, e.g. XE.V2.bp.txt.gz\n"+
				"-g Germline AnnotatedVcfParser options, defaults to '-d 12 -m 0.075 -q 0.1 -p 0.01 -g\n"+
				"      D5S,D3S -n 4.4 -a HIGH -l -c Pathogenic,Likely_pathogenic,Conflicting_\n"+
				"      interpretations_of_pathogenicity,Drug_response -t 0.51 -e Benign,Likely_benign\n"+
				"      -o -b 0.1 -z 3 -u RYR1'\n"+
				"-s Somatic AnnotatedVcfParser options, defaults to '-f -d 12 -m 0.025 -j'\n"+
				"-G Germline VCFCallFrequency options, defaults to '-b Hg38/Germline/Avatar/Bed -v \n"+
				"      Hg38/Germline/Avatar/Vcf'\n"+
				"-S Somatic VCFCallFrequency options, defaults to '-b Hg38/Somatic/Avatar/Bed -v \n"+
				"      Hg38/Somatic/Avatar/Vcf'\n"+
				"-u Minimum read coverage for tumor sample vcf records, defaults to 12\n"+
				"-i Minimum read coverage for normal sample vcf records, defaults to 12\n"+
				"-w Non matched normal alignment directory (e.g. from NA12878) to use when no matched\n"+
				"      normal is available for somatic variant calling. Needs to contain \n"+
				"      Bam/xxx_final.bam and QC/xxx_Pass.bed.gz dirs and files with indexes.\n"+

				"\nJob Execution Options:\n"+
				"-d Do not restart FAILED jobs one time.\n"+
				"-J Just set up jobs, don't launch them via sbatch.\n"+
				"-x Maximum # jobs to run at any given time, defaults to 35.\n"+
				"-z Do not nice jobs (--nice=10000), run at maximum primority.\n"+
				"-l Check and launch jobs every hour until all are complete.\n"+
				"-n Cluster partitian, defaults to hci-rw\n"+
				"-A Cluster account override, e.g. -A hci-collab\n"+

				"\nExample: java -jar pathToUSeq/Apps/TNRunner -p PatientDirs -o ~/FoundationPatients/\n"+
				"     -e ~/Hg38/DNAAlignQC/ -c ~/Hg38/SomaticCaller/ -a ~/Hg38/Annotator/ -b \n"+
				"     ~/Hg38/BamConcordance/ -j ~/Hg38/JointGenotyping/ -t ~/Hg38/RNAAlignQC/\n"+
				"     -y ~/Hg38/CopyRatio/ -k /Hg38/CopyRatio/Bkg/ -s '-d 30 -r' -x 125 -l \n"+
				"     -v ~/Hg38/Tempus/TempusVcf -m ~/Hg38/Msi/ -f ~/Hg38/StarFusion/ -l -o oncoKB.txt\n"+


				"\n**************************************************************************************\n");
	}

	public boolean isVerbose() {
		return verbose;
	}
	public File[] getDNAAlignQCDocs() {
		return DNAAlignQCDocs;
	}
	public File[] getIlluminaSomaticVarCallDocs() {
		return illuminaSomaticVarCallDocs;
	}
	public File[] getGatkSomaticVarCallDocs() {
		return gatkSomaticVarCallDocs;
	}
	public File[] getVarAnnoDocs() {
		return varAnnoDocs;
	}
	public File[] getSampleConcordanceDocs() {
		return sampleConcordanceDocs;
	}
	public String getGermlineAnnotatedVcfParser() {
		return germlineAnnotatedVcfParser;
	}
	public String getSomaticAnnotatedVcfParser() {
		return somaticAnnotatedVcfParser;
	}
	public File[] getRNAAlignQCDocs() {
		return RNAAlignQCDocs;
	}
	public File[] getRNAFusionDocs() {
		return RNAFusionDocs;
	}
	public boolean isRestartFailed() {
		return restartFailed;
	}
	public File[] getCopyRatioDocs() {
		return copyRatioDocs;
	}
	public File[] getClinicalVcfDocs() {
		return clinicalVcfDocs;
	}

	public File[] getMsiDocs() {
		return msiDocs;
	}
	
	public File[] getLoHDocs() {
		return lohDocs;
	}

	public File[] getNonMatchedNormal() {
		return nonMatchedNormal;
	}

	public File[] getHaplotypeCallDocs() {
		return haplotypeCallDocs;
	}

	public File[] getGatkJointGenotypingDocs() {
		return gatkJointGenotypingDocs;
	}

	public File[] getStrelkaJointGenotypingDocs() {
		return illuminaJointGenotypingDocs;
	}

	public String getGermlineVcfCallFreq() {
		return germlineVcfCallFreq;
	}

	public String getSomaticVcfCallFreq() {
		return somaticVcfCallFreq;
	}

	public boolean isSbatch() {
		return sbatch;
	}

	public File getBpileupFileOrDir() {
		return bpileupFileOrDir;
	}

	public File[] getCopyRatioHdf5Files() {
		return copyRatioBkgHdf5Files;
	}
	public File[] getCopyRatioIntervalListFiles() {
		return copyRatioBkgIntervalListFiles;
	}

	public HashSet<String> getPanels2SkipForCopyRatio() {
		return panels2SkipForCopyRatio;
	}

	public String getAccountOverride() {
		return accountOverride;
	}

	public File getOncoKBConfig() {
		return oncoKBConfig;
	}
}
