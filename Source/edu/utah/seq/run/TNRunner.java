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

/**Runs multiple workflows on paired T N exome and transcriptiome datasets.*/
public class TNRunner {

	//user defined fields
	private File sampleDir = null;
	private File[] rootDirs = null;
	private boolean verbose = true;
	private File[] exomeAlignQCDocs = null;
	private File[] transcriptomeAlignQCDocs = null;
	private File[] somaticVarCallDocs = null;
	private File[] bamConcordanceDocs = null;
	private File[] varAnnoDocs = null;
	private File[] jointGenotypingDocs = null;
	private File[] copyRatioDocs = null;
	private File maleBkg = null;
	private File femaleBkg = null;
	private boolean forceRestart = false;
	private int minReadCoverageTumor = 20;
	private int minReadCoverageNormal = 10;
	private TNSample[] tNSamples = null;
	private String germlineAnnotatedVcfParser = null;
	private String somaticAnnotatedVcfParser = null;
	private ArrayList<String> info = new ArrayList<String>();
	private boolean groupProcessingComplete = false;
	private boolean groupProcessingFailed = false;
	private boolean restartFailed = false;
	private HashMap<String, File> otherBams = null;
	private int maxNumJobsToSubmit = 40;

	private String pathToTrim = null;

	public TNRunner (String[] args) {
		long startTime = System.currentTimeMillis();

		processArgs(args);

		if (processIndividualSamples()){
			
			processSampleGroup();
			
			if (complete()) {
				//TODO: aggregateQCForGroup();
				IO.pl("\nALL COMPLETE!");
			}
			else IO.pl("\nNOT COMPLETE!");
		}

		

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
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
				
				//if it exists then wait until all gvcfs become available
				if (tns.getNormalExomeFastq().isFastqDirExists()) {
					File[] ab = tns.getNormalExomeBamBedGvcf();
					if (ab == null) {
						if (verbose) IO.pl("\tWaiting for germline gvcfs from "+tns.getId());
						allPresent = false;
					}
					else {
						//add the vcf and matched index
						toGeno.add(ab[2]);
						toGeno.add(ab[3]);
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
			System.err.println("\nProblem processing group analysis!");
			e.printStackTrace();
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
			IO.pl("\t\tCOMPLETE "+jobDir.getCanonicalPath());
			//find the final vcfs
			File res = new File(jobDir, jobDir.getName()+"_Hg38_GenotypedVcfs");
			File[] genoVcfs = IO.extractFiles(res, ".vcf.gz");
			//any files? might have already been moved
			if (genoVcfs.length == 0) {
				groupProcessingComplete = true;
				return;
			}
			
			//OK just completed, need to distribute files to individual Sample dirs
			File[] genoVcfIndexes = IO.extractFiles(res, ".vcf.gz.tbi");
			File logs = new File(jobDir, "Logs");
			File runScripts = new File(jobDir, "RunScripts");
			if (genoVcfs.length == 0 || genoVcfs.length != genoVcfIndexes.length) throw new IOException("\tThe JointGenotyping was marked COMPLETE but failed to find the final vcf files or their indexes in "+res);
			
			//for each separate genotyped vcf
			for (int i=0; i< genoVcfs.length; i++){
				
				//Extract the sample name 
				//HCI_P1_NormalExome_Hg38_JointGenotyping_Hg38.vcf.gz
				//1218982_NormalExome_Hg38_JointGenotyped.vcf.gz
				//1218982_Pancreas_NormalExome_Hg38_JointGenotyped.vcf.gz
				String fileName = genoVcfs[i].getName();
				File sampleDir = null;
				String id = null;
				
				//find the sample
				TNSample toLaunch = null;
				for (TNSample tns: tNSamples){
					id = tns.getId();
					if (fileName.startsWith(id)){
						sampleDir = tns.getRootDir().getCanonicalFile();
						toLaunch = tns;
						break;
					}
				}
				if (sampleDir == null) throw new IOException("\nERROR: Failed to find the sampleDir by parsing the JointGenotyping vcf names. "+fileName);
				if (verbose) IO.pl("\t\t\tMoving genotyped vcfs to "+sampleDir);
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
			TNSample.cancelDeleteJobDir(nameFile, jobDir, info);
			//launch it
			launchJointGenotyping(toGeno, jobDir);
		}
		//QUEUED
		else if (nameFile.containsKey("QUEUED")){
			info.add("\t\tQUEUED "+jobDir);
		}		
		//STARTED
		else if (nameFile.containsKey("STARTED")){
			if (TNSample.checkQueue(nameFile, jobDir, info) == false) groupProcessingFailed = true;
		}
		//FAILED but no forceRestart
		else if (nameFile.containsKey("FAILED")){
			if (verbose) IO.pl("\t\tFAILED "+jobDir);
			groupProcessingFailed = true;
		}
		//hmm no status files, probably something went wrong on the cluster? mark it as FAILED
		else {
			if (verbose) IO.pl("\t\tMarking as FAILED, no job status files in "+jobDir);
			new File(jobDir, "FAILED").createNewFile();
			groupProcessingFailed = true;
		}
	}

	private void launchJointGenotyping(ArrayList<File> toGeno, File jobDir) throws IOException {
		if (verbose) IO.pl("\t\tLaunching "+jobDir);
		//want to create/ replace the soft links
		createJointGenotypingLinks(toGeno, jobDir);

		//replace any launch scripts with the current
		File shellScript = TNSample.copyInWorkflowDocs(jointGenotypingDocs, jobDir);

		//clear any progress files
		TNSample.removeProgressFiles(jobDir);

		//squeue the shell script
		new File(jobDir, "QUEUED").createNewFile();
		String alignDirPath = jobDir.getCanonicalPath();
		String[] output = IO.executeViaProcessBuilder(new String[]{"sbatch", "-J", alignDirPath.replace(pathToTrim, ""), "-D", alignDirPath, shellScript.getCanonicalPath()}, false);
		if (verbose) for (String o: output) IO.pl("\t\t"+o);
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
			
			int numJobsLaunched = 0;
			tNSamples = new TNSample[rootDirs.length];
			for (int i=0; i< rootDirs.length; i++){
				if (verbose == false) IO.p("\t"+rootDirs[i].getName());
				tNSamples[i] = new TNSample(rootDirs[i].getCanonicalFile(), this);
				if (verbose == false) IO.pl("\t"+tNSamples[i].isRunning());
				numJobsLaunched += tNSamples[i].getNumJobsLaunched();
				if (numJobsLaunched > maxNumJobsToSubmit) {
					IO.pl("\nMaximum number jobs launched, skipping remaining samples.");
					return false;
				}
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
			File exomeWorkflowDir = null;
			File somVarCallWorkflowDir = null;
			File annoWorkflowDir = null;
			File bamConWorkflowDir = null;
			File jointGenoWorklfowDir = null;
			File transWorkflowDir = null;
			File copyRatioDocsDir = null;
			File copyRatioBkgDir = null;
			File otherDir = null;
			for (int i = 0; i<args.length; i++){
				String lcArg = args[i].toLowerCase();
				Matcher mat = pat.matcher(lcArg);
				if (mat.matches()){
					char test = args[i].charAt(1);
					try{
						switch (test){
						case 'p': sampleDir = new File(args[++i]).getCanonicalFile(); break;
						case 'e': exomeWorkflowDir = new File(args[++i]); break;
						case 't': transWorkflowDir = new File(args[++i]); break;
						case 'c': somVarCallWorkflowDir = new File(args[++i]); break;
						case 'a': annoWorkflowDir = new File(args[++i]); break;
						case 'b': bamConWorkflowDir = new File(args[++i]); break;
						case 'j': jointGenoWorklfowDir = new File(args[++i]); break;
						case 'y': copyRatioDocsDir = new File(args[++i]); break;
						case 'k': copyRatioBkgDir = new File(args[++i]); break;
						case 'o': otherDir = new File(args[++i]); break;
						case 'u': minReadCoverageTumor = Integer.parseInt(args[++i]); break;
						case 'n': minReadCoverageNormal = Integer.parseInt(args[++i]); break;
						case 'x': maxNumJobsToSubmit = Integer.parseInt(args[++i]); break;
						case 'q': verbose = false; break;
						case 'f': forceRestart = true; break;
						case 'r': restartFailed = true; break;
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
			
			//exome align qc docs
			if (exomeWorkflowDir == null || exomeWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing exome alignment workflow docs? "+exomeWorkflowDir);
			exomeAlignQCDocs = IO.extractFiles(exomeWorkflowDir);

			//transcriptome align qc docs
			if (transWorkflowDir == null || transWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing transcriptome alignment workflow docs? "+transWorkflowDir);
			transcriptomeAlignQCDocs = IO.extractFiles(transWorkflowDir);

			
			//variant calling docs
			if (somVarCallWorkflowDir == null || somVarCallWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing somatic variant calling workflow docs? "+somVarCallWorkflowDir);
			somaticVarCallDocs = IO.extractFiles(somVarCallWorkflowDir);

			//variant annotation
			if (annoWorkflowDir == null || annoWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing variant annotation workflow docs? "+annoWorkflowDir);
			varAnnoDocs = IO.extractFiles(annoWorkflowDir);

			//bam concordance
			if (bamConWorkflowDir == null || bamConWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing bam concordance workflow docs? "+bamConWorkflowDir);
			bamConcordanceDocs = IO.extractFiles(bamConWorkflowDir);

			//joint genotyping
			if (jointGenoWorklfowDir == null || jointGenoWorklfowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing joint genotyping workflow docs? "+jointGenoWorklfowDir);
			jointGenotypingDocs = IO.extractFiles(jointGenoWorklfowDir);

			//copy ratio analysis
			if (copyRatioDocsDir == null || copyRatioDocsDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing copy ratio analysis workflow docs? "+copyRatioDocsDir);
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
			
			//AnnotatedVcfParser options for germline and somatic filtering
			if (germlineAnnotatedVcfParser == null) germlineAnnotatedVcfParser = "-d "+minReadCoverageNormal+" -m 0.2 -x 1 -p 0.01 -g D5S,D3S -n 5 -a HIGH -c Pathogenic,Likely_pathogenic -o -e Benign,Likely_benign";
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
			
			if (verbose){
				IO.pl("Run parameters:");
				IO.pl("Patient sample directory\t"+rootDirs[0].getParent());
				IO.pl("Exome alignment workflow directory\t"+exomeWorkflowDir);
				IO.pl("Transcriptome alignment workflow directory\t"+transWorkflowDir);
				IO.pl("Somatic variant workflow directory\t"+somVarCallWorkflowDir);
				IO.pl("Variant annotation workflow directory\t"+annoWorkflowDir);
				IO.pl("Bam concordance workflow directory\t"+bamConWorkflowDir);
				IO.pl("Joint genotyping workflow directory\t"+jointGenoWorklfowDir);
				IO.pl("Copy ratio workflow directory\t"+copyRatioDocsDir);
				IO.pl("Copy ratio background directory\t"+copyRatioBkgDir);
				IO.pl("Other patient sample directory for concordance\t"+otherDir);
				IO.pl("Min tumor read coverage\t"+minReadCoverageTumor);
				IO.pl("Min normal read coverage\t"+minReadCoverageNormal);
				IO.pl("Germline Annotated Vcf Parser options\t"+germlineAnnotatedVcfParser);
				IO.pl("Somatic Annotated Vcf Parser options\t"+somaticAnnotatedVcfParser);
				IO.pl("Restart failed jobs\t"+restartFailed);
				IO.pl("Force restart\t"+forceRestart);
				IO.pl("Verbose logging\t"+verbose);
				IO.pl("Max # jobs to launch\t"+maxNumJobsToSubmit);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                               TNRunner  December 2018                            **\n" +
				"**************************************************************************************\n" +
				"TNRunner is designed to execute several dockerized snakmake workflows on human tumor\n"+
				"normal datasets via a slurm cluster.  Based on the availability of fastq, Hg38\n"+
				"alignments are run, somatic and germline variants are called, and concordance measured\n"+
				"between sample bams. To execute TNRunner, create the following directory structure and\n"+
				"link or copy in the corresponding paired end Illumina gzipped fastq files.\n"+
				"MyPatientSampleDatasets\n"+
				"   MyPatientA\n"+
				"      Fastq\n"+
				"         NormalExome\n"+
				"         TumorExome\n"+
				"         TumorTranscriptome\n"+
				"   MyPatientB\n"+
				"      Fastq\n"+
				"         TumorExome\n"+
				"         TumorTranscriptome\n"+
				"   MyPatientC....\n"+
				"The Fastq directory and sub directories must match this naming. Only include those\n"+
				"for which you have fastq.  Change the MyXXX to something relevant. TNRunner is\n"+
				"stateless so as more Fastq becomes available or issues are addressed, relaunch the\n"+
				"app. This won't effect running or queued slurm jobs. Relaunch periodically to assess\n"+
				"the current processing status and queue additional tasks. Download the latest\n"+
				"workflows from https://github.com/HuntsmanCancerInstitute/Workflows/tree/master/Hg38 \n"+
				"and the matching resource bundle from\n"+
				"https://hci-bio-app.hci.utah.edu/gnomex/gnomexFlex.jsp?analysisNumber=A5578\n"+
					
				"\nOptions:\n"+
				"-p Directory containing one or more patient data directories to process.\n" +
				"-o Other patient's directory, containing additional xxx_final.bam files to include in\n"+
				"   sample concordance. The patient directory naming must match.\n"+
				"-k Directory containing xxxMalePoN.hdf5 and xxxFemalePoN.hdf5 GATK copy ratio\n"+
				"      background files.\n"+
				"-y Workflow docs for launching copy ratio analysis.\n"+
				"-e Workflow docs for launching exome alignments.\n"+
				"-t Workflow docs for launching transcriptome alignments.\n"+
				"-c Workflow docs for launching somatic variant calling.\n"+
				"-a Workflow docs for launching variant annotation.\n"+
				"-b Workflow docs for launching bam concordance.\n"+
				"-j Workflow docs for launching joint genotyping.\n"+
				"-g Germline AnnotatedVcfParser options, defaults to '-d 10 -m 0.2 -x 1 -p 0.01 -g \n"+
				"      D5S,D3S -n 5 -a HIGH -c Pathogenic,Likely_pathogenic -o -e Benign,Likely_benign'\n"+
				"-s Somatic AnnotatedVcfParser options, defaults to '-d 20 -f'\n"+
				"-r Restart FAILED jobs.\n"+
				"-f Force a restart of all uncompleted jobs.\n"+
				"-q Quite output.\n"+
				"-x Maximum # jobs to launch, defaults to 40.\n"+

				"\nExample: java -jar pathToUSeq/Apps/TNRunner -p AvatarPatients -o ~/FoundationPatients/\n"+
				"     -e ~/Hg38/ExomeAlignQC/ -c ~/Hg38/SomExoCaller/ -a ~/Hg38/Annotator/ -b \n"+
				"     ~/Hg38/BamConcordance/ -j ~/Hg38/JointGenotyping/ -t ~/Hg38/TranscriptomeAlignQC/\n"+
				"     -y /Hg38/CopyRatio/ -k /Hg38/CopyRatio/Bkg/ -s '-d 30 -r'  \n\n"+


				"**************************************************************************************\n");
	}

	public boolean isVerbose() {
		return verbose;
	}
	public File[] getExomeAlignQCDocs() {
		return exomeAlignQCDocs;
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
	public File[] getTranscriptomeAlignQCDocs() {
		return transcriptomeAlignQCDocs;
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
}
