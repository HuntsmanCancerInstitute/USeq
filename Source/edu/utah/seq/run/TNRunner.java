package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.nio.file.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

/**Runs multiple workflows on paired T N DNA and transcriptiome datasets.*/
public class TNRunner {

	//user defined fields
	private File sampleDir = null;
	private File[] rootDirs = null;
	private boolean verbose = true;
	private File[] DNAAlignQCDocs = null;
	private File[] RNAAlignQCDocs = null;
	private File[] RNAFusionDocs = null;
	private File[] somaticVarCallDocs = null;
	private File[] bamConcordanceDocs = null;
	private File[] varAnnoDocs = null;
	private File[] jointGenotypingDocs = null;
	private File[] copyRatioDocs = null;
	private File[] clinicalVcfDocs = null;
	private File[] msiDocs = null;
	private File normalAlignmentDir = null;
	private File[] nonMatchedNormal = null;
	private File maleBkg = null;
	private File femaleBkg = null;
	private boolean forceRestart = false;
	private boolean softRestart = false;
	private int minReadCoverageTumor = 20;
	private int minReadCoverageNormal = 15;
	private TNSample[] tNSamples = null;
	private String germlineAnnotatedVcfParser = null;
	private String somaticAnnotatedVcfParser = null;
	private ArrayList<String> info = new ArrayList<String>();
	private boolean groupProcessingComplete = false;
	private boolean groupProcessingFailed = false;
	private boolean restartFailed = false;
	private HashMap<String, File> otherBams = null;
	private int maxNumJobs = 25;
	private boolean loop = false;
	private int numMinToSleep = 60;
	private boolean niceJobs = true;
	private String partition = "hci-rw";

	private String pathToTrim = null;
	private String nice = "";

	public TNRunner (String[] args) {
		long startTime = System.currentTimeMillis();

		processArgs(args);

		//loop?
		int iterations = 1;
		int i=0;
		if (loop) iterations = 1000;
		
		for (; i< iterations; i++){
			if (processIndividualSamples()){
				processSampleGroup();
				if (complete()) {
					IO.pl("\nALL COMPLETE!");
					break;
				}
				else IO.pl("\nNOT COMPLETE!");
			}
			else IO.pl("\nNOT COMPLETE!");
			if (iterations !=1){
				try {
					IO.pl("\nWaiting... "+i);
					Thread.sleep(1000*60*numMinToSleep);
				} catch (InterruptedException e) {}
			}
		}
		
		printSamplesWithFastqIssues();
		

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
	}

	private void printSamplesWithFastqIssues() {
		ArrayList<String> issues = new ArrayList<String>();
		for (TNSample t: tNSamples){
			if (t!= null && t.isFastqIssue()) issues.add(t.getId());
		}
		if (issues.size() != 0 && verbose) IO.pl("\nThe following have been skipped due to fastq file count issues: "+issues);
	}

	private boolean complete() {
		if (groupProcessingComplete == false || groupProcessingFailed == true) return false;
		for (TNSample t: tNSamples){
			if (t.isFailed() || t.isRunning()) return false;
		}
		return true;
	}

	private void processSampleGroup() {
		try {
			//any joint genotyping workflow files? If none then they want to skip this, e.g. no matched normal
			if (jointGenotypingDocs == null) {
				groupProcessingFailed = false;
				groupProcessingComplete = true;
				return;
			}
			IO.pl("\nChecking group processing...");

			
			//for each sample with a fastq germline component are the alignments ready?
			ArrayList<File> toGeno = new ArrayList<File>();
			boolean allPresent = true;
			for (TNSample tns: tNSamples){
				//JointGenotyping already run?
				File gVcfDir = new File(tns.getRootDir(), "GermlineVariantCalling");
				if (gVcfDir.exists()) {
					if (verbose) IO.pl("\tComplete");
					groupProcessingComplete = true;
					return;
				}
				
				//OK doesn't exist so wait until all gvcfs become available
				//does this sample have a normal sample?
				if (tns.getNormalDNAFastq().isFastqDirExists()) {
					AlignmentDataset normalDNADataset = tns.getNormalDNADataset();
					if (normalDNADataset == null || normalDNADataset.isComplete() == false) {
						if (verbose) IO.pl("\tWaiting for germline gvcfs from "+tns.getId());
						allPresent = false;
					}
					//OK it's present, add it if it's not the mock normal
					else {
						if (normalDNADataset.getgVcfFile().getName().equals("NA12878_NormalDNA_Hg38_Haplo.g.vcf.gz") == false) {
							//add the vcf and matched index
							toGeno.add(normalDNADataset.getgVcfFile());
							toGeno.add(normalDNADataset.getgVcfIndexFile());
						}
					}
				}
			}
			if (allPresent == false) return;
			
			//launch it
			checkJointGenotyping(toGeno);
			if (groupProcessingFailed || verbose) Misc.printArray(info);

			if (verbose) {
				IO.pl("\tFailed? "+groupProcessingFailed);
				IO.pl("\tComplete? "+groupProcessingComplete);
			}
			
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem processing group analysis!");
		}

	}

	private void checkJointGenotyping(ArrayList<File> toGeno) throws IOException {

		//make dir, ok if it already exists
		File jobDir = new File (sampleDir, "JointGenotyping");
		jobDir.mkdirs();

		//any files?
		HashMap<String, File> nameFile = IO.fetchNamesAndFiles(jobDir);
		if (nameFile.size() == 0) {
			launchJointGenotyping(toGeno, jobDir);
			return;
		}

		//OK some files are present
		//COMPLETE
		else if (nameFile.containsKey("COMPLETE")){
			IO.pl("\tCOMPLETE "+jobDir.getCanonicalPath());
			//find the final vcfs
			File res = new File(jobDir, jobDir.getName()+"_Hg38_GenotypedVcfs");
			if (res.exists() == false) res = new File (jobDir, "Vcfs");
			File[] genoVcfs = IO.extractFiles(res, "Haplo_JointGenotyped.vcf.gz");
			
			//any files? might have already been moved
			if (genoVcfs.length == 0) {
				groupProcessingComplete = true;
				return;
			}
			
			//OK just completed, need to distribute files to individual Sample dirs
			File[] genoVcfIndexes = IO.extractFiles(res, "Haplo_JointGenotyped.vcf.gz.tbi");
			File logs = new File(jobDir, "Logs");
			File runScripts = new File(jobDir, "RunScripts");
			
			if (genoVcfs.length == 0 || genoVcfs.length != genoVcfIndexes.length) throw new IOException("\tThe JointGenotyping was marked COMPLETE but failed to find the final vcf files or their indexes in "+res);
			
			//for each separate genotyped vcf
			for (int i=0; i< genoVcfs.length; i++){

				//Extract the sample name 
				//HCI_P1_NormalDNA_Hg38_JointGenotyping_Hg38.vcf.gz
				//1218982_NormalDNA_Hg38_JointGenotyped.vcf.gz
				//1218982_Pancreas_NormalDNA_Hg38_JointGenotyped.vcf.gz
				
				//how about just remove _NormalDNA_Hg38_JointGenotyped.vcf.gz?
				
				
				String fileName = genoVcfs[i].getName();
				String sampleName = null;
				if (fileName.contains("_NormalDNA_Hg38_JointGenotyped.vcf.gz")) sampleName = fileName.replace("_NormalDNA_Hg38_JointGenotyped.vcf.gz", "");
				else sampleName = fileName.replace("_NormalDNA_Hg38_Haplo_JointGenotyped.vcf.gz", "");
				if (sampleName.length() == fileName.length()) throw new IOException("\nERROR: Failed to parse the sample name from "+fileName);
				
				File sampleDir = null;
				String id = null;
				
				//find the sample
				TNSample toLaunch = null;
				for (TNSample tns: tNSamples){
					id = tns.getId();			
					if (id.equals(sampleName)){
						sampleDir = tns.getRootDir().getCanonicalFile();
						toLaunch = tns;
						break;
					}
				}
				if (sampleDir == null) throw new IOException("\nERROR: Failed to find the sampleDir by parsing the JointGenotyping vcf names. "+fileName);
				
				if (verbose) IO.pl("\tMoving genotyped vcfs to "+sampleDir);
				
				
				//delete any existing germline folder
				IO.deleteDirectoryViaCmdLine(new File(sampleDir, "GermlineVariantCalling"));
				//create a folder for the incoming genotyped vcf in the sample folder
				File resDir = new File(sampleDir, "GermlineVariantCalling/"+fileName.replace("_Hg38_JointGenotyped.vcf.gz", ""));
				resDir.mkdirs();
				
				//copy and move over relevant files
				File newLogDir = new File(resDir,"Logs");
				newLogDir.mkdir();
				File newRSDir = new File(resDir,"RunScripts");
				newRSDir.mkdir();
				IO.copyDirectoryRecursive(logs, newLogDir, null);
				IO.copyDirectoryRecursive(runScripts, newRSDir, null);
				
				genoVcfs[i].renameTo(new File(resDir, genoVcfs[i].getName()));
				genoVcfIndexes[i].renameTo(new File(resDir, genoVcfIndexes[i].getName()));
				
				//launch germline annotator
				toLaunch.annotateGermlineVcf();
			}
			groupProcessingComplete = true;
		}

		//force a restart?
		else if (forceRestart || restartFailed && nameFile.containsKey("FAILED")){
			//cancel any slurm jobs and delete the directory
			TNSample.cancelDeleteJobDir(nameFile, jobDir, info, true);
			//launch it
			launchJointGenotyping(toGeno, jobDir);
		}
		//QUEUED
		else if (nameFile.containsKey("QUEUED")){
			info.add("\tQUEUED "+jobDir);
		}		
		//STARTED
		else if (nameFile.containsKey("STARTED")){
			if (TNSample.checkQueue(nameFile, jobDir, info) == false) groupProcessingFailed = true;
		}
		//FAILED but no forceRestart
		else if (nameFile.containsKey("FAILED")){
			if (verbose) IO.pl("\tFAILED "+jobDir);
			groupProcessingFailed = true;
		}
		//hmm no status files, probably something went wrong on the cluster? mark it as FAILED
		else {
			if (verbose) IO.pl("\tMarking as FAILED, no job status files in "+jobDir);
			new File(jobDir, "FAILED").createNewFile();
			groupProcessingFailed = true;
		}
	}

	private void launchJointGenotyping(ArrayList<File> toGeno, File jobDir) throws IOException {
		if (verbose) IO.pl("\tLAUNCHING "+jobDir);
		//want to create/ replace the soft links
		createJointGenotypingLinks(toGeno, jobDir);

		//replace any launch scripts with the current
		File shellScript = TNSample.copyInWorkflowDocs(jointGenotypingDocs, jobDir);

		//clear any progress files
		TNSample.removeProgressFiles(jobDir);

		//squeue the shell script
		String alignDirPath = jobDir.getCanonicalPath();
		String [] cmd = null;
		if (nice.length() !=0) cmd = new String[]{"sbatch", nice, "-J", alignDirPath.replace(pathToTrim, ""), "-D", alignDirPath, shellScript.getCanonicalPath()};
		else cmd = new String[]{"sbatch", "-J", alignDirPath.replace(pathToTrim, ""), "-D", alignDirPath, shellScript.getCanonicalPath()};
		String[] output = IO.executeViaProcessBuilder(cmd, false);
		if (verbose) for (String o: output) IO.pl("\t\t"+o);
		for (String o: output) {
			if (o.toLowerCase().contains("error")) {
				groupProcessingFailed = true;
				IO.pl("\tERROR launching "+Misc.stringArrayToString(cmd, " "));
				if (verbose == false) for (String x: output) IO.pl("\t\t"+x);
				new File(jobDir, "FAILED").createNewFile();
				return;
			}
		}
		new File(jobDir, "QUEUED").createNewFile();
		
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
			tNSamples = new TNSample[rootDirs.length];
			for (int i=0; i< rootDirs.length; i++){
				if (numJobsLaunched >= maxNumJobs) {
					IO.pl("\nMaximum number jobs launched, skipping remaining samples.");
					return false;
				}
				if (verbose == false) IO.p("\t"+rootDirs[i].getName());
				tNSamples[i] = new TNSample(rootDirs[i].getCanonicalFile(), this);
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
		new TNRunner(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		try {
			IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
			Pattern pat = Pattern.compile("-[a-z]");
			File DNAWorkflowDir = null;
			File somVarCallWorkflowDir = null;
			File annoWorkflowDir = null;
			File bamConWorkflowDir = null;
			File jointGenoWorklfowDir = null;
			File transWorkflowDir = null;
			File rnaFuseDir = null;
			File copyRatioDocsDir = null;
			File copyRatioBkgDir = null;
			File clinicalVcfDir = null;
			File msiWorkflowDir = null;
			File otherDir = null;
			for (int i = 0; i<args.length; i++){
				String lcArg = args[i].toLowerCase();
				Matcher mat = pat.matcher(lcArg);
				if (mat.matches()){
					char test = args[i].charAt(1);
					try{
						switch (test){
						case 'p': sampleDir = new File(args[++i]).getCanonicalFile(); break;
						case 'e': DNAWorkflowDir = new File(args[++i]); break;
						case 't': transWorkflowDir = new File(args[++i]); break;
						case 'f': rnaFuseDir = new File(args[++i]); break;
						case 'c': somVarCallWorkflowDir = new File(args[++i]); break;
						case 'a': annoWorkflowDir = new File(args[++i]); break;
						case 'b': bamConWorkflowDir = new File(args[++i]); break;
						case 'm': msiWorkflowDir = new File(args[++i]); break;
						case 'j': jointGenoWorklfowDir = new File(args[++i]); break;
						case 'y': copyRatioDocsDir = new File(args[++i]); break;
						case 'k': copyRatioBkgDir = new File(args[++i]); break;
						case 'v': clinicalVcfDir = new File(args[++i]); break;
						case 'o': otherDir = new File(args[++i]); break;
						case 'w': normalAlignmentDir = new File(args[++i]); break;
						case 'g': germlineAnnotatedVcfParser = args[++i]; break;
						case 's': somaticAnnotatedVcfParser = args[++i]; break;
						case 'u': minReadCoverageTumor = Integer.parseInt(args[++i]); break;
						case 'i': minReadCoverageNormal = Integer.parseInt(args[++i]); break;
						case 'x': maxNumJobs = Integer.parseInt(args[++i]); break;
						case 'z': niceJobs = true; break;
						case 'q': verbose = false; break;
						case 'h': forceRestart = true; break;
						case 'l': loop = true; break;
						case 'd': restartFailed = true; break;
						case 'r': softRestart = true; break;
						case 'n': partition = args[++i]; break;
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
			if (sampleDir == null || sampleDir.exists() == false) Misc.printErrAndExit("Error: failed to find your data directory? "+sampleDir);
			rootDirs = IO.extractOnlyDirectories(sampleDir);
			if (rootDirs == null || rootDirs.length == 0) Misc.printErrAndExit("Error: failed to find root directories to process? "+rootDirs);
			
			//remove JointGenotyping if present
			ArrayList<File> toKeep = new ArrayList<File>();
			for (File rd: rootDirs){
				String name = rd.getName();
				if (name.equals("JointGenotyping") == false) toKeep.add(rd);
			}
			rootDirs = new File[toKeep.size()];
			toKeep.toArray(rootDirs);
			
			//DNA align qc docs
			if (DNAWorkflowDir != null){
				if (DNAWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing DNA alignment workflow docs? "+DNAWorkflowDir);
				DNAAlignQCDocs = IO.extractFiles(DNAWorkflowDir);
			}

			//RNA align qc docs
			if (transWorkflowDir != null) {
				if (transWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing RNA alignment workflow docs? "+transWorkflowDir);
				RNAAlignQCDocs = IO.extractFiles(transWorkflowDir);
			}
			
			//RNA fusion docs
			if (rnaFuseDir != null) {
				if (rnaFuseDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing RNA fusion workflow docs? "+rnaFuseDir);
				RNAFusionDocs = IO.extractFiles(rnaFuseDir);
			}
			
			//variant calling docs
			if (somVarCallWorkflowDir != null){
				if(somVarCallWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing somatic variant calling workflow docs? "+somVarCallWorkflowDir);
				somaticVarCallDocs = IO.extractFiles(somVarCallWorkflowDir);
			}
			
			//variant annotation
			if (annoWorkflowDir != null){
				if (annoWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing variant annotation workflow docs? "+annoWorkflowDir);
				varAnnoDocs = IO.extractFiles(annoWorkflowDir);
			}
			
			//bam concordance
			if (bamConWorkflowDir != null){
				if (bamConWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing bam concordance workflow docs? "+bamConWorkflowDir);
				bamConcordanceDocs = IO.extractFiles(bamConWorkflowDir);
			}
			
			//joint genotyping
			if (jointGenoWorklfowDir != null){
				if (jointGenoWorklfowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing joint genotyping workflow docs? "+jointGenoWorklfowDir);
				jointGenotypingDocs = IO.extractFiles(jointGenoWorklfowDir);
			}
			
			//msi
			if (msiWorkflowDir != null){
				if (msiWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing msi workflow docs? "+msiWorkflowDir);
				msiDocs = IO.extractFiles(msiWorkflowDir);
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
				File[] f = IO.extractFiles(copyRatioBkgDir, "PoN.hdf5");
				if (f == null || f.length != 2) Misc.printErrAndExit("Error: failed to find two copy ratio background xxxPoN.hdf5 files in "+copyRatioBkgDir);
				String name = f[0].getName().toLowerCase();
				if (name.contains("female")) {
					femaleBkg = f[0];
					maleBkg = f[1];
				}
				else {
					femaleBkg = f[1];
					maleBkg = f[0];
				}
				if (maleBkg == null || femaleBkg == null) Misc.printErrAndExit("Error: failed to find male "+maleBkg+" and female "+femaleBkg+" background xxxFemalePoN.hdf5 or xxxMalePoN.hdf5 files in "+copyRatioBkgDir);
			}
			
			//non matched normal?
			if (normalAlignmentDir != null) {
				//find the final bam and bed
				File[] finalBam = IO.extractFiles(new File(normalAlignmentDir, "Bam"), "_final.bam");
				File[] passingBed = IO.extractFiles(new File(normalAlignmentDir, "QC"), "_Pass.bed.gz");
				if (finalBam == null || finalBam.length !=1 )throw new IOException("\nFailed to find a single Bam/xxx_final.bam in "+normalAlignmentDir);
				if (passingBed == null || passingBed.length !=1 )throw new IOException("\nFailed to find a single QC/xxx_Pass.bed.gz in "+normalAlignmentDir);
				nonMatchedNormal = new File[] {finalBam[0], passingBed[0]};
			}
			
			
			//AnnotatedVcfParser options for germline and somatic filtering
			if (germlineAnnotatedVcfParser == null) germlineAnnotatedVcfParser = "-d "+minReadCoverageNormal+" -m 0.1 -q 0.1 -p 0.01 -g D5S,D3S -n 4.4 -a HIGH -l -c "
					+ "Pathogenic,Likely_pathogenic,Conflicting_interpretations_of_pathogenicity,Drug_response -t 0.51 -e Benign,Likely_benign -o -b 0.1 -z 3 -u RYR1";
			if (somaticAnnotatedVcfParser == null) somaticAnnotatedVcfParser = "-d "+minReadCoverageTumor+" -f";

			//set path to trim
			pathToTrim = sampleDir.getCanonicalFile().getParentFile().getCanonicalPath()+"/";
			
			//any Foundation datasets?
			if (otherDir != null){
				File[] p = IO.extractOnlyDirectories(otherDir);
				if (p != null && p.length !=0){
					otherBams = new HashMap<String, File>();
					for (File i: p) otherBams.put(i.getName(), i);
				}
				
			}
			
			if (niceJobs) nice = "--nice=10000";
			
			if (verbose){
				IO.pl("Run parameters:");
				IO.pl("Patient sample directory\t"+rootDirs[0].getParent());
				IO.pl("DNA alignment workflow directory\t"+DNAWorkflowDir);
				IO.pl("RNA alignment workflow directory\t"+transWorkflowDir);
				IO.pl("Somatic variant workflow directory\t"+somVarCallWorkflowDir);
				if (normalAlignmentDir != null) IO.pl("Non matched normal alignment directory for somatic calling\t"+normalAlignmentDir);
				IO.pl("Variant annotation workflow directory\t"+annoWorkflowDir);
				IO.pl("Bam concordance workflow directory\t"+bamConWorkflowDir);
				IO.pl("MSI workflow directory\t"+msiWorkflowDir);
				IO.pl("Joint genotyping workflow directory\t"+jointGenoWorklfowDir);
				IO.pl("Copy ratio workflow directory\t"+copyRatioDocsDir);
				IO.pl("Copy ratio background directory\t"+copyRatioBkgDir);
				IO.pl("Other patient sample directory for concordance\t"+otherDir);
				IO.pl("Min tumor read coverage\t"+minReadCoverageTumor);
				IO.pl("Min normal read coverage\t"+minReadCoverageNormal);
				IO.pl("Germline Annotated Vcf Parser options\t"+germlineAnnotatedVcfParser);
				IO.pl("Somatic Annotated Vcf Parser options\t"+somaticAnnotatedVcfParser);
				IO.pl("Restart failed jobs\t"+softRestart);
				IO.pl("Delete and restart failed jobs\t"+restartFailed);
				IO.pl("Force restart\t"+forceRestart);
				IO.pl("Verbose logging\t"+verbose);
				IO.pl("Max # jobs to launch\t"+maxNumJobs);
				IO.pl("Nice jobs\t"+niceJobs);
				IO.pl("Job Partition\t"+partition);
				IO.pl("Relaunch jobs until complete\t"+loop);
			}
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing arguments for TNRunner, aborting!");
		}
	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                                  TNRunner : Nov 2021                             **\n" +
				"**************************************************************************************\n" +
				"TNRunner is designed to execute several containerized snakmake workflows on tumor\n"+
				"normal datasets via a slurm cluster.  Based on the availability of fastq, \n"+
				"alignments are run, somatic and germline variants are called, and concordance measured\n"+
				"between sample bams. To execute TNRunner, create the following directory structure and\n"+
				"link or copy in the corresponding paired end Illumina gzipped fastq files. To run \n"+
				"Joint Genotyping on sample sets where this has been run before, delete the\n "+
				"GermlineVariantCalling and JointGenotyping folders.\n"+
				
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
				"the current processing status and queue additional tasks or set option -l. Download the\n"+
				"latest workflows from \n"+
				"https://github.com/HuntsmanCancerInstitute/Workflows/tree/master/Hg38RunnerWorkflows \n"+
				"and the matching resource bundle from\n"+
				"https://hci-bio-app.hci.utah.edu/gnomex/gnomexFlex.jsp?analysisNumber=A5578 .\n"+
				"All workflow docs are optional although some require output from prior Analysis.\n"+
					
				"\nSetup Options:\n"+
				"-p Directory containing one or more patient data directories to process.\n" +
				"-e Workflow docs for launching DNA alignments.\n"+
				"-t Workflow docs for launching RNA alignments.\n"+
				"-f Workflow docs for launching RNA fusion detection.\n"+
				"-c Workflow docs for launching somatic variant calling.\n"+
				"-m Workflow docs for launching MSI status calling.\n"+
				"-a Workflow docs for launching variant annotation.\n"+
				"-b Workflow docs for launching bam concordance.\n"+
				"-j Workflow docs for launching joint genotyping.\n"+
				"-y Workflow docs for launching somatic copy analysis.\n"+
				"-v Workflow docs for launching clinical test variant info. Add a ClinicalReport folder to\n"+
				"      each patient dir containing the json formatted clinical information.\n"+
				"-o Other patient's directory, containing additional xxx_final.bam files to include in\n"+
				"   sample concordance. The patient directory naming must match.\n"+
				"-k Directory containing xxxMalePoN.hdf5 and xxxFemalePoN.hdf5 GATK copy ratio\n"+
				"      background files.\n"+
				"-g Germline AnnotatedVcfParser options, defaults to '-d 15 -m 0.1 -q 0.1 -p 0.01 -g\n"+
				"      D5S,D3S -n 4.4 -a HIGH -l -c Pathogenic,Likely_pathogenic,Conflicting_\n"+
				"      interpretations_of_pathogenicity,Drug_response -t 0.51 -e Benign,Likely_benign\n"+
				"      -o -b 0.1 -z 3 -u RYR1'\n"+
				"-s Somatic AnnotatedVcfParser options, defaults to '-d 20 -f'\n"+
				"-u Minimum read coverage for tumor DNA samples, defaults to 20\n"+
				"-i Minimum read coverage for normal DNA samples, defaults to 15\n"+
				"-w Non matched normal alignment directory (e.g. from NA12878) to use when no matched\n"+
				"      normal is available for somatic variant calling. Needs to contain Bam/xxx_final.bam\n"+
				"      and QC/xxx_Pass.bed.gz dirs and files with indexes.\n"+
				
				"\nJob Execution Options:\n"+
				"-r Attempt to restart FAILED jobs from last successfully completed rule.\n"+
				"-d Delete and restart FAILED jobs.\n"+
				"-h Force a restart of all running, queued, failed, and uncompleted jobs.\n"+
				"-q Quite output.\n"+
				"-x Maximum # jobs to run at any given time, defaults to 25.\n"+
				"-z Do not nice jobs (--nice=10000), run at maximum primority.\n"+
				"-l Check and launch jobs every hour until all are complete.\n"+
				"-n Cluster partitian, defaults to hci-rw\n"+

				"\nExample: java -jar pathToUSeq/Apps/TNRunner -p PatientDirs -o ~/FoundationPatients/\n"+
				"     -e ~/Hg38/DNAAlignQC/ -c ~/Hg38/SomaticCaller/ -a ~/Hg38/Annotator/ -b \n"+
				"     ~/Hg38/BamConcordance/ -j ~/Hg38/JointGenotyping/ -t ~/Hg38/RNAAlignQC/\n"+
				"     -y ~/Hg38/CopyRatio/ -k /Hg38/CopyRatio/Bkg/ -s '-d 30 -r' -x 10 -l \n"+
				"     -v ~/Hg38/Tempus/TempusVcf -m ~/Hg38/Msi/ -f ~/Hg38/StarFusion/ -l \n"+


				"\n**************************************************************************************\n");
	}

	public boolean isVerbose() {
		return verbose;
	}
	public File[] getDNAAlignQCDocs() {
		return DNAAlignQCDocs;
	}
	public boolean isForceRestart() {
		return forceRestart;
	}
	public int getMinReadCoverageTumor() {
		return minReadCoverageTumor;
	}
	public int getMinReadCoverageNormal() {
		return minReadCoverageNormal;
	}
	public String getPathToTrim() {
		return pathToTrim;
	}
	public File[] getSomaticVarCallDocs() {
		return somaticVarCallDocs;
	}
	public File[] getVarAnnoDocs() {
		return varAnnoDocs;
	}
	public File[] getBamConcordanceDocs() {
		return bamConcordanceDocs;
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
	public HashMap<String, File> getOtherBams() {
		return otherBams;
	}
	public File[] getCopyRatioDocs() {
		return copyRatioDocs;
	}
	public File getMaleBkg() {
		return maleBkg;
	}
	public File getFemaleBkg() {
		return femaleBkg;
	}
	public boolean isSoftRestart() {
		return softRestart;
	}

	public File[] getClinicalVcfDocs() {
		return clinicalVcfDocs;
	}

	public File[] getMsiDocs() {
		return msiDocs;
	}

	public File[] getNonMatchedNormal() {
		return nonMatchedNormal;
	}
}
