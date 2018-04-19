package edu.utah.seq.vcf;

import java.io.*;
import java.util.regex.*;

import edu.utah.seq.parsers.MergeSams;
import util.bio.annotation.Bed;
import util.gen.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * @author david.nix@hci.utah.edu 
 * Takes a bed file of desired regions, splits it, writes out each, then calls GATK on each, finally merges the results.
 **/
public class GatkRunner {

	//user defined fields
	private File bedFile;
	private File saveDirectory;
	private String gatkArgs = null;
	private int numberThreads = 0;
	boolean useLowerCaseL = false;
	boolean bamOut = false;
	
	//internal fields
	GatkRunnerChunk[] runners;
	
	//constructors
	public GatkRunner(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running GatkRunner!");
		}
	}

	public void doWork() throws Exception{
		
		System.out.println("Parsing and spliting bed file...");
		
		//parse regions and randomize
		Bed[] regions = Bed.parseFile(bedFile, 0, 0);
		Misc.randomize(regions, 0);
		int numRegionsPerChunk = 1+ (int)Math.round( (double)regions.length/ (double)numberThreads );
		
		//split regions
		ArrayList<Bed[]> splitBed = Bed.splitByNumber(regions, numRegionsPerChunk);
		
		//create runners and start
		System.out.println("Launching...\n");
		runners = new GatkRunnerChunk[splitBed.size()];
		ExecutorService executor = Executors.newFixedThreadPool(runners.length);
		for (int i=0; i< runners.length; i++) {
			runners[i] = new GatkRunnerChunk(splitBed.get(i), this, ""+i);
			executor.execute(runners[i]);
		}
		executor.shutdown();
		//spins here until the executer is terminated, e.g. all threads complete
        while (!executor.isTerminated()) {
        }
		
		//check runners and delete temp files
        for (GatkRunnerChunk c: runners) {
			if (c.isFailed()) Misc.printErrAndExit("\nERROR: Failed runner, aborting! \n"+c.getCmd());
			c.getTempBed().deleteOnExit();
			c.getTempVcf().deleteOnExit();
			c.getTempLog().deleteOnExit();
			new File (c.getTempVcf()+".idx").deleteOnExit();
		}
		
		//merge vcf
		System.out.println("\nMerging vcfs...");
		mergeVcfs();
		
		//merge bams?
		if (bamOut) {
			System.out.println("\nMerging bams...");
			mergeBams();
		}
	}

	private void mergeVcfs() throws Exception{
		ArrayList<VCFRecord> recordsAL = new ArrayList<VCFRecord>();
		VCFParser vp = null;
		for (GatkRunnerChunk c: runners) {
			vp = new VCFParser(c.getTempVcf(), true, false, false);
			VCFRecord[] records = vp.getVcfRecords();
			if (records == null) continue;
			for (VCFRecord r: records) recordsAL.add(r);
		}
		VCFRecord[] all = new VCFRecord[recordsAL.size()];
		recordsAL.toArray(all);
		Arrays.sort(all);
		PrintWriter out = new PrintWriter( new FileWriter( new File( saveDirectory, "gatk.raw.vcf")));
		
		//write out header
		ArrayList<String> header = vp.getVcfComments().getOriginalComments();
		for (String h : header) out.println(h);
		
		//write out unique records
		HashSet<String> unique = new HashSet<String>();
		for (VCFRecord r : all) {
			String ori = r.getOriginalRecord();
			if (unique.contains(ori) == false) {
				out.println(ori);
				unique.add(ori);
			}
		}
		out.close();
		System.out.println("\t"+unique.size()+"\tvariants found");
	}
	
	private void mergeBams() throws Exception{
		File[] bams = new File[runners.length];
		for (int i=0; i< runners.length; i++) bams[i] = runners[i].getBamOut();
		File mergedBam = new File (saveDirectory, "mutect.realigned.bam");
		MergeSams ms = new MergeSams(bams, mergedBam, false, true);
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new GatkRunner(args);
	}		
	
	/**This method will process each argument and assign new varibles*/
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
					case 'r': bedFile = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'c': {
						StringInt si = parseCmd(args, ++i); 
						gatkArgs = si.cmd;
						i = si.i-1;
						break;
					}
					case 'l': useLowerCaseL= true; break;
					case 'b': bamOut= true; break;
					case 't': numberThreads = Integer.parseInt(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//threads to use
		double gigaBytesAvailable = ((double)Runtime.getRuntime().maxMemory())/ 1073741824.0;
		int numPossCores = (int)Math.round(gigaBytesAvailable/5.0);
		if (numPossCores < 1) numPossCores = 1;
		int numPossThreads = Runtime.getRuntime().availableProcessors();
		if (numberThreads == 0){
			if (numPossCores <= numPossThreads) numberThreads = numPossCores;
			else numberThreads = numPossThreads;
			System.out.println("Core usage:\n\tGB available to Java:\t"+ Num.formatNumber(gigaBytesAvailable, 1));
			System.out.println("\tAvailable cores:\t"+numPossThreads);
			System.out.println("\tNumber cores to use @ 5GB/core:\t"+numberThreads+"\n");
		}

		//check save dir
		if (saveDirectory == null) Misc.printErrAndExit("\nError: cannot find your save directory!\n"+saveDirectory);
		saveDirectory.mkdirs();
		if (saveDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: your save directory does not appear to be a directory?\n");

		//check bed
		if (bedFile == null || bedFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find your bed file of regions to interrogate?\n"+bedFile);
		
		//check args
		if (gatkArgs == null) Misc.printErrAndExit("\nError: please provide a gatk launch cmd similar to the example.\n\n");
		if (gatkArgs.contains("~") || gatkArgs.contains("./")) Misc.printErrAndExit("\nError: provide full paths in the GATK command.\n"+gatkArgs);
		if (gatkArgs.contains(" -o ") || gatkArgs.contains(" -L ")) Misc.printErrAndExit("\nError: please don't provide a -o or -L argument to the cmd.\n"+gatkArgs);	
		
		
	}	
	
	private StringInt parseCmd(String[] args, int i) {
		ArrayList<String> al = new ArrayList<String>();
		
		for (; i< args.length; i++){
			if (args[i].equals("-r") || args[i].equals("-s") || args[i].equals("-l") || args[i].equals("-b") || args[i].equals("-t")) break;
			al.add(args[i]);
		}
		String cmd = Misc.stringArrayListToString(al, " ");
		return new StringInt(i,cmd);
	}
	
	private class StringInt{
		int i;
		String cmd;
		public StringInt(int i, String cmd){
			this.i = i;
			this.cmd = cmd;
		}
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Gatk Runner: March 2018                            **\n" +
				"**************************************************************************************\n" +
				"Takes a bed file of target regions, splits it by the number of threads, writes out\n"+
				"each, executes the GATK Gatktype caller, and merges the results. Set the -Xmx to the\n"+
				"maximum available on the machine to enable correct cpu thread usage.\n"+

				"\nOptions:\n"+
				"-r A regions bed file (chr, start, stop,...) to intersect, see\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1 , gz/zip OK.\n"+
				"-s Path to a directory for saving the results.\n"+
				"-t Number concurrent thread override. Sets itself based on the memory and cpus \n"+
				"     available to the JVM.\n"+
				"-c GATK command to execute, see the example below, modify to match your enviroment.\n"+
				"     Most resources require full paths. Don't set -o or -L\n"+
				"-l Use lowercased l for Lofreq compatability.\n"+
				"-b Add a -bamout argument and merge bam chunks.\n"+
				
				"\nExample: java -Xmx24G -jar pathToUSeq/Apps/GatkRunner -b -r /SS/targets.bed -s\n" +
				"     /SS/HC/ -c 'java -Xmx4G -jar /SS/GenomeAnalysisTK.jar -T MuTect2 \n"+
				"    -R /SS/human_g1k_v37.fasta --dbsnp /SS/dbsnp_138.b37.vcf \n"+
				"    --cosmic /SS/v76_GRCh37_CosmicCodingMuts.vcf.gz -I:tumor /SS/sar.bam -I:normal \n"+
				"    /SS/normal.bam'\n\n" +

				"**************************************************************************************\n");

	}


	public File getSaveDirectory() {
		return saveDirectory;
	}

	public String getGatkArgs() {
		return gatkArgs;
	}

	
}
