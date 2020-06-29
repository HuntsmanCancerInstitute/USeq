package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;

/**Tool to assist in identifying new Avatar datasets for processing.*/
public class AvatarComparator {
	
	private File jobDir;
	private HashMap<String, File> patientDirs = null;
	private TreeSet<File> toProcess = new TreeSet<File>();
	private TreeSet<File> needsAttention = new TreeSet<File>();
	private TreeSet<File> complete = new TreeSet<File>();
	private ArrayList<File> s3FilesToRestore = new ArrayList<File>();
	private boolean verbose = false;
	
	public AvatarComparator (String[] args){
		long startTime = System.currentTimeMillis();
		try {
			processArgs(args);
			
			//fetch dirs
			patientDirs = IO.extractDirectories(jobDir);
			patientDirs.remove("AggregateQCStats");
			patientDirs.remove("JointGenotyping");
			
			walkPatientDirs();
			
			printResults();
			
		} catch (IOException e) {
			e.printStackTrace();
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" Sec\n");
	}
	
    private void printResults() throws IOException {
		IO.pl("\nReady for additional analysis: "+toProcess.size());
		for (File f: toProcess) IO.pl("ln -s "+f.getCanonicalPath()+" .");
		
		IO.pl("\nNeed to fetch these from Amazon S3.  As root, rename them, and then run GSync:");
		for (File f: s3FilesToRestore) IO.pl("mv "+f.getCanonicalPath()+" "+f.getCanonicalPath()+".restore");
		
		IO.pl("\nNeeds attention: "+needsAttention.size());
		for (File f: needsAttention) IO.pl("ls -1 "+f.getCanonicalPath()+"/*/*");
		
		if (verbose){
			IO.pl("\nCurrent: "+complete.size());
			for (File f: complete) IO.pl("\t"+f);
		}
		else IO.pl("\nNumber COMPLETE "+complete.size());
		
	}

	private void walkPatientDirs() throws IOException {
		for (String name: patientDirs.keySet()){
			if (verbose) IO.pl("Checking "+name);
			File patientDir = patientDirs.get(name).getCanonicalFile();
			
			//see if this has been marked COMPLETE
			File c = new File(patientDir, "COMPLETE");
			if (c.exists()){
				if (verbose) IO.pl(name+ "\tMarked COMPLETE");
				complete.add(patientDir);
				continue;
			}
			
			//pull Fastq dir
			File fastqDir = new File(patientDir, "Fastq");
			if (fastqDir.exists() == false) {
				IO.pl(name+"\tERROR: Failed to find "+name+"/Fastq dir");
				needsAttention.add(patientDir);
				continue;
			}
			
			//pull dirs in Fastq: NormalDNA  TumorDNA  TumorRNA
			File[] fastqDirs = IO.extractOnlyDirectories(fastqDir);
			if (fastqDirs == null || fastqDirs.length == 0) {
				IO.pl(name+"\tERROR: Failed to find one or more FastqDirs in "+name+"/Fastq/");
				needsAttention.add(patientDir);
				continue;
			}
			
			//check each Fastq dir for just two fastq files, could be xxx.S3.txt placeholders too, fetch all the fastq .S3.txt files too
			HashMap<String, ArrayList<File>> nameS3Fastq = new HashMap<String, ArrayList<File>>();
			for (File fqd: fastqDirs){
				File[] fastq = IO.extractFiles(fqd);
				
				if (fastq == null || fastq.length < 2){
					IO.pl(name+"\tERROR: Failed to find at least 2 fastq in "+fqd);
					needsAttention.add(patientDir);
					continue;
				}
				//collapse the names
				HashSet<String> collapsedNames = new HashSet<String>();
				for (File f: fastq) {
					if (f.getName().endsWith(".S3.txt")) {
						collapsedNames.add(f.getName().replace(".S3.txt", ""));
						ArrayList<File> toFetch = nameS3Fastq.get(fqd.getName());
						if (toFetch == null) {
							toFetch = new ArrayList<File>();
							nameS3Fastq.put(fqd.getName(), toFetch);
						}
						toFetch.add(f);
					}
					else collapsedNames.add(f.getName());
				}
				if (collapsedNames.size() != 2){
					IO.pl(name+"\tERROR: Failed to find just 2 fastq in "+fqd);
					needsAttention.add(patientDir);
					continue;
				}
			}
			
			//pull alignment dir and dirs
			File alignDir = new File(patientDir, "Alignment");
			if (alignDir.exists() == false) {
				if (verbose) IO.pl(name+"\tNo "+name+"/Alignment thus new dataset for processing");
				toProcess.add(patientDir);
				//pull all fastq and look for those that need to be restored
				for (String fastqDirName: nameS3Fastq.keySet()) {
					//only add to restore list if the .gz isn't present
					for (File fq: nameS3Fastq.get(fastqDirName)) {
						String conPath = fq.getCanonicalPath().replace(".S3.txt", "");
						if (new File (conPath).exists() == false) s3FilesToRestore.add(fq);
					}
				}
				continue;
			}
			
			//pull AlignmentDirs: 1218022_NormalDNA  1218022_TumorDNA  1218022_TumorRNA
			File[] alignDirs = IO.extractOnlyDirectories(alignDir);
			if (alignDirs == null || alignDirs.length == 0) {
				IO.pl(name+"\tERROR: Failed to find one or more Alignment dirs in "+name+"/Alignment/, delete the Alignment dir?");
				needsAttention.add(patientDir);
				continue;
			}
			
			if (fastqDirs.length < alignDirs.length){
				IO.pl(name+"\tERROR: # Fastq dirs < # Alignment dirs?");
				needsAttention.add(patientDir);
				continue;
			}
			
			//compare fastq dirs against alignment dirs, fetch all the bam and cram .S3.txt files too, with new Fastq, need to redo the SampleConcordance
			if (fastqDirs.length > alignDirs.length){
				if (verbose) IO.pl(name+"\tNew "+name+"/Fastq dir thus new dataset for processing");
				toProcess.add(patientDir);
				//for each alignment dir, remove any fastq S3, not needed
				for (File ad: alignDirs) {

					int index = ad.getName().lastIndexOf("_");
					if (index != -1) nameS3Fastq.remove(ad.getName().substring(index+1));
					File bam = new File(ad, "Bam");
					if (bam.exists() == false) {
						IO.pl(name+"\tERROR: Failed to find the Bam/ folder in "+ad+" alignment dir?!");
						needsAttention.add(patientDir);
						continue;
					}
					File[] s3 = IO.extractFiles(bam, ".S3.txt");
					for (File f: s3) {
						String conPath = f.getCanonicalPath().replace(".S3.txt", "");
						if (new File (conPath).exists() == false) s3FilesToRestore.add(f);
					}
				}
				//pull all the remaining S3 fastq that remain
				for (String fastqDirName: nameS3Fastq.keySet()) {
					for (File fq: nameS3Fastq.get(fastqDirName)) s3FilesToRestore.add(fq);
				}
				continue;
			}
			

			
			//tumor and normal present? If so then these folders need to be present
			String alignDirNames = Misc.stringArrayToString(IO.fetchFileNames(alignDirs), ",");
			if (alignDirNames.contains("_NormalDNA") && alignDirNames.contains("_TumorDNA")){

				File cr = new File(patientDir, "CopyAnalysis"); 
				File gc = new File(patientDir, "GermlineVariantCalling");
				File sc = new File(patientDir, "SampleConcordance");
				File sv = new File(patientDir, "SomaticVariantCalls");
				File[] toCheck = new File[]{cr, gc, sc, sv};
				for (File dir: toCheck){
					if (dir == null || dir.exists() == false) {
						toProcess.add(patientDir);
						continue;
					}
					else {
						File[] f = IO.extractFiles(dir);
						if (f == null || f.length == 0){
							toProcess.add(patientDir);
							continue;
						}
					}
				}
			}
			
			//nothing to be done
			if(verbose) IO.pl(name+"\tNo additional processing needed");
			complete.add(patientDir);
		}
		
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AvatarComparator(args);
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
					case 'j': jobDir = new File(args[++i]); break;
					case 'v': verbose = true; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					e.printStackTrace();
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                           Avatar Comparator : June 2020                          **\n" +
				"**************************************************************************************\n" +
				"Tool for identifying AVATAR datasets that are ready for analysis or need attention.\n"+

				"\nOptions:\n"+
				"-j Patient job directory\n" +
				"-v Verbose output\n"+
				
                "\nExample: java -jar -Xmx2G ~/USeqApps/AvatarComparator -j AJobs/ \n\n"+


				"**************************************************************************************\n");
	}

}
