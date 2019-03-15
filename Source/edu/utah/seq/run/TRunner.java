package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

/**Runs multiple workflows on Tumore capture, tumor transcriptiome, and xml reports, e.g. Foundation*/
public class TRunner {

	//user defined fields
	private File sampleDir = null;
	private File[] rootDirs = null;
	private boolean verbose = true;
	private File[] captureAlignQCDocs = null;
	private File[] transcriptomeAlignQCDocs = null;
	private File[] somaticVarCallDocs = null;
	private File[] vcfIntegrationDocs = null;
	private File[] varAnnoDocs = null;
	private boolean forceRestart = false;
	private boolean softRestart = false;
	private TSample[] tSamples = null;
	private String somaticAnnotatedVcfParser = "-d 50 -f";
	private boolean restartFailed = false;
	private String pathToTrim = null;
	private int maxNumJobsToSubmit = 0;

	public TRunner (String[] args) {
		long startTime = System.currentTimeMillis();

		processArgs(args);

		if (processIndividualSamples()){
			if (complete()) IO.pl("\nALL COMPLETE!");
			else IO.pl("\nNOT COMPLETE!");
		}
		else IO.pl("\nNOT COMPLETE!");

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
	}

	private boolean complete() {
		for (TSample t: tSamples){
			if (t.isFailed() || t.isRunning()) return false;
		}
		return true;
	}

	private boolean processIndividualSamples() {
		try {
			if (verbose) IO.pl("\nChecking individual samples...");
			else IO.pl("\nChecking individual samples (SampleID RunningJobs)");
			
			int numJobsLaunched = 0;
			tSamples = new TSample[rootDirs.length];
			for (int i=0; i< rootDirs.length; i++){
				if (verbose == false) IO.p("\t"+rootDirs[i].getName());
				tSamples[i] = new TSample(rootDirs[i].getCanonicalFile(), this);
				if (verbose == false) IO.pl("\t"+tSamples[i].isRunning());
				numJobsLaunched += tSamples[i].getNumJobsLaunched();
				//hit max?
				if (maxNumJobsToSubmit!=0 && numJobsLaunched > maxNumJobsToSubmit) {
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
		new TRunner(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		try {
			IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
			Pattern pat = Pattern.compile("-[a-z]");
			File captureWorkflowDir = null;
			File somVarCallWorkflowDir = null;
			File annoWorkflowDir = null;
			File vcfIntegrationDir = null;
			File transWorkflowDir = null;
			for (int i = 0; i<args.length; i++){
				String lcArg = args[i].toLowerCase();
				Matcher mat = pat.matcher(lcArg);
				if (mat.matches()){
					char test = args[i].charAt(1);
					try{
						switch (test){
						case 'p': sampleDir = new File(args[++i]).getCanonicalFile(); break;
						case 'e': captureWorkflowDir = new File(args[++i]); break;
						case 't': transWorkflowDir = new File(args[++i]); break;
						case 'c': somVarCallWorkflowDir = new File(args[++i]); break;
						case 'a': annoWorkflowDir = new File(args[++i]); break;
						case 'v': vcfIntegrationDir = new File(args[++i]); break;
						case 's': somaticAnnotatedVcfParser = args[++i]; break;
						case 'x': maxNumJobsToSubmit = Integer.parseInt(args[++i]); break;
						case 'q': verbose = false; break;
						case 'f': forceRestart = true; break;
						case 'd': restartFailed = true; break;
						case 'r': softRestart = true; break;
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
			
			//capture align qc docs
			if (captureWorkflowDir == null || captureWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing capture alignment workflow docs? "+captureWorkflowDir);
			captureAlignQCDocs = IO.extractFiles(captureWorkflowDir);

			//transcriptome align qc docs
			if (transWorkflowDir == null || transWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing transcriptome alignment workflow docs? "+transWorkflowDir);
			transcriptomeAlignQCDocs = IO.extractFiles(transWorkflowDir);

			//variant calling docs
			if (somVarCallWorkflowDir == null || somVarCallWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing somatic variant calling workflow docs? "+somVarCallWorkflowDir);
			somaticVarCallDocs = IO.extractFiles(somVarCallWorkflowDir);

			//variant annotation
			if (annoWorkflowDir == null || annoWorkflowDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing variant annotation workflow docs? "+annoWorkflowDir);
			varAnnoDocs = IO.extractFiles(annoWorkflowDir);

			//vcf integration
			if (vcfIntegrationDir == null || vcfIntegrationDir.exists() == false) Misc.printErrAndExit("Error: failed to find a directory containing vcf integrationDocs workflow docs? "+vcfIntegrationDir);
			vcfIntegrationDocs = IO.extractFiles(vcfIntegrationDir);

			//set path to trim
			pathToTrim = sampleDir.getCanonicalFile().getParentFile().getCanonicalPath()+"/";
			
			if (verbose){
				IO.pl("Run parameters:");
				IO.pl("Patient sample directory\t"+rootDirs[0].getParent());
				IO.pl("Capture alignment workflow docs\t"+captureWorkflowDir);
				IO.pl("Transcriptome alignment workflow docs\t"+transWorkflowDir);
				IO.pl("Somatic variant workflow docs\t"+somVarCallWorkflowDir);
				IO.pl("Somatic Annotated Vcf Parser options\t"+somaticAnnotatedVcfParser);
				IO.pl("Vcf integration workflow directory\t"+vcfIntegrationDir);
				IO.pl("Variant annotation workflow docs\t"+annoWorkflowDir);
				IO.pl("Restart failed jobs\t"+softRestart);
				IO.pl("Delete and restart failed jobs\t"+restartFailed);
				IO.pl("Force restart\t"+forceRestart);
				IO.pl("Verbose logging\t"+verbose);
				if (maxNumJobsToSubmit!=0) IO.pl("Max # jobs to launch\t"+maxNumJobsToSubmit);
				else IO.pl("Max # jobs to launch\tno limit");
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                                 TRunner December 2018                            **\n" +
				"**************************************************************************************\n" +
				"TRunner is designed to execute several dockerized snakmake workflows on human tumor\n"+
				"only datasets via a slurm cluster.  Based on the availability of the bams and xml,\n"+
				"fastq is extracted, Hg38 alignments run, somatic variants called, and xml variant info\n"+
				"parsed and compared. To execute TRunner, create the following directory structure and\n"+
				"link or copy in the corresponding raw bam, bai, and xml files.\n"+
				"MyPatientSampleDatasets\n"+
				"   MyPatient1\n"+
				"      Bam\n"+
				"         TumorDNA\n"+
				"         TumorRNA\n"+
				"      Xml\n"+
				"   MyPatient2\n"+
				"      Bam\n"+
				"         TumorDNA\n"+
				"      Xml\n"+
				"   MyPatient3....\n"+
				
				"\nThe Bam directory and sub directories must match this naming. Change MyXXX to \n"+
				"something relevant. As more files becomes available or issues are addressed, relaunch\n"+
				"the app. This won't effect running slurm jobs. Relaunch periodically to assess\n"+
				"the current processing status and queue additional tasks. Download the latest\n"+
				"workflows from https://github.com/HuntsmanCancerInstitute/Workflows/tree/master/\n"+
				"Hg38RunnerWorkflows/Foundation and the matching resource bundle from\n"+
				"https://hci-bio-app.hci.utah.edu/gnomex/gnomexFlex.jsp?analysisNumber=A5578\n"+
					
				"\nOptions:\n"+
				"-p Directory containing one or more patient data directories to process.\n" +
				"-e Workflow docs for launching DNA capture alignments.\n"+
				"-t Workflow docs for launching transcriptome alignments.\n"+
				"-c Workflow docs for launching somatic variant calling.\n"+
				"-a Workflow docs for launching variant annotation.\n"+
				"-v Workflow docs for launching xml variant integration.\n"+
				"-s Somatic AnnotatedVcfParser options, defaults to '-d 50 -f'\n"+
				"-r Attempt to restart FAILED jobs from last successfully completed rule.\n"+
				"-d Delete and restart FAILED jobs.\n"+
				"-f Force a restart of all running and uncompleted jobs.\n"+
				"-q Quite output.\n"+
				"-x Maximum # jobs to launch, defaults to 0, no limit\n"+

				"\nExample: java -jar pathToUSeq/Apps/TRunner -p FoundationPatients -e \n"+
				"     ~/Hg38/CaptureAlignQC/WorkflowDocs/ -c ~/Hg38/SomExoCaller/WorkflowDocs/ -a \n"+
				"     ~/Hg38/Annotator/WorkflowDocs/ -b ~/Hg38/BamConcordance/WorkflowDocs/ -j\n"+
				"     ~/Hg38/JointGenotyping/WorkflowDocs/ -t ~/Hg38/TranscriptomeAlignQC/WorkflowDocs/ \n"+
				"     -s '-d 100 -f' -x 20 \n\n"+


				"**************************************************************************************\n");
	}

	public boolean isVerbose() {
		return verbose;
	}
	public File[] getCaptureAlignQCDocs() {
		return captureAlignQCDocs;
	}
	public boolean isForceRestart() {
		return forceRestart;
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
	public String getSomaticAnnotatedVcfParser() {
		return somaticAnnotatedVcfParser;
	}
	public File[] getTranscriptomeAlignQCDocs() {
		return transcriptomeAlignQCDocs;
	}
	public boolean isRestartFailed() {
		return restartFailed;
	}
	public File[] getVcfIntegrationDocs() {
		return vcfIntegrationDocs;
	}

	public boolean isSoftRestart() {
		return softRestart;
	}
}
