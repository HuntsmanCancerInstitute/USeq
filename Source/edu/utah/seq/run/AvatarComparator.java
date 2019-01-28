package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
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
		
		IO.pl("\nNeeds attention: "+needsAttention.size());
		for (File f: needsAttention) IO.pl("ls -1 "+f.getCanonicalPath()+"/*/*");
		
		if (verbose){
			IO.pl("\nCurrent: "+complete.size());
			for (File f: complete) IO.pl("\t"+f);
		}
		
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
			
			//pull dirs in Fastq: NormalExome  TumorExome  TumorTranscriptome
			File[] fastqDirs = IO.extractOnlyDirectories(fastqDir);
			if (fastqDirs == null || fastqDirs.length == 0) {
				IO.pl(name+"\tERROR: Failed to find one or more FastqDirs in "+name+"/Fastq/");
				needsAttention.add(patientDir);
				continue;
			}
			
			//check each Fastq dir for just two fastq files
			for (File fqd: fastqDirs){
				File[] fastq = IO.extractFiles(fqd, ".gz");
				if (fastq == null || fastq.length !=2){
					IO.pl(name+"\tERROR: Failed to find just two fastq in "+fqd);
					needsAttention.add(patientDir);
					continue;
				}
			}
			
			//pull alignment dir and dirs
			File alignDir = new File(patientDir, "Alignment");
			if (alignDir.exists() == false) {
				if (verbose) IO.pl(name+"\tNo "+name+"/Alignment thus new dataset for processing");
				toProcess.add(patientDir);
				continue;
			}
			
			//pull AlignmentDirs: 1218022_NormalExome  1218022_TumorExome  1218022_TumorTranscriptome
			File[] alignDirs = IO.extractOnlyDirectories(alignDir);
			if (alignDirs == null || alignDirs.length == 0) {
				IO.pl(name+"\tERROR: Failed to find one or more Alignment dirs in "+name+"/Alignment/");
				needsAttention.add(patientDir);
				continue;
			}
			
			//compare fastq dirs against alignment dirs
			if (fastqDirs.length > alignDirs.length){
				if (verbose) IO.pl(name+"\tNew "+name+"/Fastq dir thus new dataset for processing");
				toProcess.add(patientDir);
				continue;
			}
			
			if (fastqDirs.length < alignDirs.length){
				IO.pl(name+"\tERROR: # Fastq dirs < # Alignment dirs?");
				needsAttention.add(patientDir);
				continue;
			}
			
			//tumor and normal present?
			String alignDirNames = Misc.stringArrayToString(IO.fetchFileNames(alignDirs), ",");
			if (alignDirNames.contains("_NormalExome") && alignDirNames.contains("_TumorExome")){

//need to update this when you make the rename
File cr = new File(patientDir, "CopyRatio"); 

				File gc = new File(patientDir, "GermlineVariantCalling");
				File sc = new File(patientDir, "SampleConcordance");
				File sv = new File(patientDir, "SomaticVariantCalls");
				File[] toCheck = new File[]{cr, gc, sc, sv};
				boolean fail = false;
				for (File dir: toCheck){
					if (dir.exists() == false) {
						IO.pl(name+"\tERROR: Failed to find "+dir);
						fail = true;
					}
					else {
						File[] f = IO.extractFiles(dir);
						if (f == null || f.length == 0){
							IO.pl(name+"\tERROR: Failed to find any files in "+dir);
							fail = true;
						}
					}
				}
				if (fail){ 
					needsAttention.add(patientDir);
					continue;
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
				"**                           Avatar Comparator : Jan 2019                           **\n" +
				"**************************************************************************************\n" +
				"Tool for identifying AVATAR datasets that are ready for analysis or need attention.\n"+

				"\nOptions:\n"+
				"-j Patient job directory\n" +
				"-v Verbose output\n"+
				
                "\nExample: java -jar -Xmx2G ~/USeqApps/AvatarComparator -j AJobs/ \n\n"+


				"**************************************************************************************\n");
	}

}
