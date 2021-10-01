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
public class TNRunner2 {

	//user defined fields
	private File sampleDir = null;
	private File[] rootDirs = null;
	private boolean verbose = true;
	private File[] DNAAlignQCDocs = null;
	private File[] RNAAlignQCDocs = null;
	private File[] RNAFusionDocs = null;
	private File[] somaticVarCallDocs = null;
	private File[] sampleConcordanceDocs = null;
	private File[] varAnnoDocs = null;
	private File[] haplotypeCallDocs = null;
	private File[] gatkJointGenotypingDocs = null;
	private File[] strelkaJointGenotypingDocs = null;
	private File[] copyRatioDocs = null;
	private File[] clinicalVcfDocs = null;
	private File[] msiDocs = null;
	private File normalAlignmentDir = null;
	private File[] nonMatchedNormal = null;
	private File maleBkg = null;
	private File femaleBkg = null;
	private int minReadCoverageTumor = 12;
	private int minReadCoverageNormal = 12;
	private TNSample2[] tNSamples = null;
	private String germlineAnnotatedVcfParser = null;
	private String somaticAnnotatedVcfParser = null;
	private String germlineVcfCallFreq = "-b Hg38/Germline/Avatar/Bed -v Hg38/Germline/Avatar/Vcf";
	private String somaticVcfCallFreq = "-b Hg38/Somatic/Avatar/Bed -v Hg38/Somatic/Avatar/Vcf";
	private ArrayList<String> info = new ArrayList<String>();
	private boolean groupProcessingComplete = false;
	private boolean groupProcessingFailed = false;
	private boolean restartFailed = false;
	private int maxNumJobs = 35;
	private boolean loop = false;
	private int numMinToSleep = 60;
	private boolean niceJobs = true;
	private String partition = "hci-rw";

	private String pathToTrim = null;
	private String nice = "";

	public TNRunner2 (String[] args) {
		long startTime = System.currentTimeMillis();

		processArgs(args);

		//loop?
		int iterations = 1;
		int i=0;
		if (loop) iterations = 1000;
		
		for (; i< iterations; i++){
			if (processIndividualSamples()){
				processGATKSampleGroup();
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
		for (TNSample2 t: tNSamples){
			if (t!= null && t.isFastqIssue()) issues.add(t.getId());
		}
		if (issues.size() != 0 && verbose) IO.pl("\nThe following have been skipped due to fastq file count issues: "+issues);
	}

	private boolean complete() {
		if (groupProcessingComplete == false || groupProcessingFailed == true) return false;
		for (TNSample2 t: tNSamples){
			if (t.isFailed() || t.isRunning()) return false;
		}
		return true;
	}
	
	private void processGATKSampleGroup() {
		try {
			//any joint genotyping workflow files? If none then they want to skip this, e.g. no matched normal
			if (gatkJointGenotypingDocs == null) {
				groupProcessingFailed = false;
				groupProcessingComplete = true;
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
					groupProcessingComplete = true;
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
				groupProcessingComplete = true;
				return;
			}
			
			//OK just completed, need to distribute files to individual Sample dirs
			File[] genoVcfIndexes = IO.extractFiles(res, "JointGenotyped.vcf.gz.tbi");
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
				String sampleName = fileName.replace("_GVCF_Hg38_JointGenotyped.vcf.gz", "");
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
				
				
				//delete any existing germline folder
				IO.deleteDirectoryViaCmdLine(new File(sampleDir, "GermlineVariantCalling"));
				//create a folder for the incoming genotyped vcf in the sample folder
				File resDir = new File(sampleDir, "GermlineVariantCalling/"+sampleName+"_GATK");
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
				
				//launch germline annotator?
				if (this.getVarAnnoDocs()!= null) toLaunch.annotateGermlineVcf("GATK");
			}
			groupProcessingComplete = true;
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

	private void launchGatkJointGenotyping(ArrayList<File> toGeno, File jobDir) throws IOException {
		if (verbose) IO.pl("\tLAUNCHING "+jobDir);
		//want to create/ replace the soft links
		createJointGenotypingLinks(toGeno, jobDir);

		//replace any launch scripts with the current
		File shellScript = TNSample.copyInWorkflowDocs(gatkJointGenotypingDocs, jobDir);

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
			tNSamples = new TNSample2[rootDirs.length];
			for (int i=0; i< rootDirs.length; i++){
				if (numJobsLaunched >= maxNumJobs) {
					IO.pl("\nMaximum number jobs launched, skipping remaining samples.");
					return false;
				}
				if (verbose == false) IO.p("\t"+rootDirs[i].getName());
				tNSamples[i] = new TNSample2(rootDirs[i].getCanonicalFile(), this);
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
			File somVarCallWorkflowDir = null;
			File annoWorkflowDir = null;
			File sampleConWorkflowDir = null;
			File gatkJointGenoWorklfowDir = null;
			File strelkaJointGenoWorklfowDir = null;
			File haploWorklfowDir = null;
			File RNAWorkflowDir = null;
			File rnaFuseDir = null;
			File copyRatioDocsDir = null;
			File copyRatioBkgDir = null;
			File clinicalVcfDir = null;
			File msiWorkflowDir = null;
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
						case 'c': somVarCallWorkflowDir = new File(args[++i]); break;
						case 'a': annoWorkflowDir = new File(args[++i]); break;
						case 'b': sampleConWorkflowDir = new File(args[++i]); break;
						case 'm': msiWorkflowDir = new File(args[++i]); break;
						case 'j': gatkJointGenoWorklfowDir = new File(args[++i]); break;
//case 'q': strelkaJointGenoWorklfowDir = new File(args[++i]); break;
						case 'h': haploWorklfowDir = new File(args[++i]); break;
						case 'y': copyRatioDocsDir = new File(args[++i]); break;
						case 'k': copyRatioBkgDir = new File(args[++i]); break;
						case 'v': clinicalVcfDir = new File(args[++i]); break;
						case 'w': normalAlignmentDir = new File(args[++i]); break;
						case 'g': germlineAnnotatedVcfParser = args[++i]; break;
						case 's': somaticAnnotatedVcfParser = args[++i]; break;
						case 'G': germlineVcfCallFreq = args[++i]; break;
						case 'S': somaticVcfCallFreq = args[++i]; break;
						case 'u': minReadCoverageTumor = Integer.parseInt(args[++i]); break;
						case 'i': minReadCoverageNormal = Integer.parseInt(args[++i]); break;
						case 'x': maxNumJobs = Integer.parseInt(args[++i]); break;
						case 'z': niceJobs = true; break;
						case 'l': loop = true; break;
						case 'd': restartFailed = true; break;
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
			if (RNAWorkflowDir != null) {
				if (RNAWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing RNA alignment workflow docs? "+RNAWorkflowDir);
				RNAAlignQCDocs = IO.extractFiles(RNAWorkflowDir);
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
			if (strelkaJointGenoWorklfowDir != null){
				if (strelkaJointGenoWorklfowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing strelka joint genotyping workflow docs? "+strelkaJointGenoWorklfowDir);
				strelkaJointGenotypingDocs = IO.extractFiles(strelkaJointGenoWorklfowDir);
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
//need to generate cram files from these and rename the Pass to PassRC
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
			
			if (niceJobs) nice = "--nice=10000";
			
			if (verbose){
				IO.pl("Run parameters:");
				IO.pl("Patient sample directory\t"+rootDirs[0].getParent());
				IO.pl("DNA alignment workflow directory\t"+DNAWorkflowDir);
				IO.pl("RNA alignment workflow directory\t"+RNAWorkflowDir);
				IO.pl("Somatic variant workflow directory\t"+somVarCallWorkflowDir);
				if (normalAlignmentDir != null) IO.pl("Non matched normal alignment directory for somatic calling\t"+normalAlignmentDir);
				IO.pl("Variant annotation workflow directory\t"+annoWorkflowDir);
				IO.pl("Sample concordance workflow directory\t"+sampleConWorkflowDir);
				IO.pl("MSI workflow directory\t"+msiWorkflowDir);
				IO.pl("GATK haplotype calling workflow directory\t"+haploWorklfowDir);
				IO.pl("GATK Joint genotyping workflow directory\t"+gatkJointGenoWorklfowDir);
//IO.pl("Strelka Joint genotyping workflow directory\t"+strelkaJointGenoWorklfowDir);
				IO.pl("Copy ratio workflow directory\t"+copyRatioDocsDir);
				IO.pl("Copy ratio background directory\t"+copyRatioBkgDir);
	
				IO.pl("Min tumor read coverage\t"+minReadCoverageTumor);
				IO.pl("Min normal read coverage\t"+minReadCoverageNormal);
				IO.pl("Germline Annotated Vcf Parser options\t"+germlineAnnotatedVcfParser);
				IO.pl("Somatic Annotated Vcf Parser options\t"+somaticAnnotatedVcfParser);
	
				IO.pl("Delete and restart failed jobs\t"+restartFailed);
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
				"**                                TNRunner2 : Sept 2021                             **\n" +
				"**************************************************************************************\n" +
				"TNRunner is designed to execute several containerized snakmake workflows on tumor\n"+
				"normal datasets via a slurm cluster.  Based on the availability of fastq, \n"+
				"alignments are run, somatic and germline variants are called, and concordance measured\n"+
				"between sample bams. To execute TNRunner, create the following directory structure and\n"+
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
				"-b Workflow docs for launching sample concordance.\n"+
				"-y Workflow docs for launching somatic copy analysis.\n"+
				"-h Workflow docs for launching GATK haplotype calling.\n"+
				"-j Workflow docs for launching GATK joint genotyping.\n"+
//"-q Workflow docs for launching Strelka germline joint genotyping.\n"+
				"-v Workflow docs for launching clinical test variant info. Add a ClinicalReport folder to\n"+
				"      each patient dir containing the json formatted clinical information.\n"+
				"-k Directory containing xxxMalePoN.hdf5 and xxxFemalePoN.hdf5 GATK copy ratio\n"+
				"      background files.\n"+
				"-g Germline AnnotatedVcfParser options, defaults to '-d 12 -m 0.1 -q 0.1 -p 0.01 -g\n"+
				"      D5S,D3S -n 4.4 -a HIGH -l -c Pathogenic,Likely_pathogenic,Conflicting_\n"+
				"      interpretations_of_pathogenicity,Drug_response -t 0.51 -e Benign,Likely_benign\n"+
				"      -o -b 0.1 -z 3 -u RYR1'\n"+
				"-s Somatic AnnotatedVcfParser options, defaults to '-d 12 -f'\n"+
				"-G Germline VCFCallFrequency options, defaults to '-b Hg38/Germline/Avatar/Bed -v Hg38/Germline/Avatar/Vcf'\n"+
				"-S Somatic VCFCallFrequency options, defaults to '-b Hg38/Somatic/Avatar/Bed -v Hg38/Somatic/Avatar/Vcf'\n"+
				"-u Minimum read coverage for tumor sample vcf records, defaults to 12\n"+
				"-i Minimum read coverage for normal sample vcf records, defaults to 12\n"+
				"-w Non matched normal alignment directory (e.g. from NA12878) to use when no matched\n"+
				"      normal is available for somatic variant calling. Needs to contain Bam/xxx_final.bam\n"+
				"      and QC/xxx_Pass.bed.gz dirs and files with indexes.\n"+
				
				"\nJob Execution Options:\n"+
				"-d Delete and restart FAILED jobs.\n"+
				"-x Maximum # jobs to run at any given time, defaults to 35.\n"+
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
	public String getPathToTrim() {
		return pathToTrim;
	}
	public File[] getSomaticVarCallDocs() {
		return somaticVarCallDocs;
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
	public File getMaleBkg() {
		return maleBkg;
	}
	public File getFemaleBkg() {
		return femaleBkg;
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

	public File[] getHaplotypeCallDocs() {
		return haplotypeCallDocs;
	}

	public File[] getGatkJointGenotypingDocs() {
		return gatkJointGenotypingDocs;
	}

	public File[] getStrelkaJointGenotypingDocs() {
		return strelkaJointGenotypingDocs;
	}

	public String getGermlineVcfCallFreq() {
		return germlineVcfCallFreq;
	}

	public String getSomaticVcfCallFreq() {
		return somaticVcfCallFreq;
	}
}
