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
	boolean useUpperCaseO = false;
	boolean bamOut = false;
	
	//internal fields
	GatkRunnerChunk[] runners;
	private int maxBpPerRegion = 1000;
	
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
		
		//parse regions 
		Bed[] regions = Bed.parseFile(bedFile, 0, 0);
		regions = Bed.splitBigRegions(regions, maxBpPerRegion);
	
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
        while (!executor.isTerminated()) {}
		
		//check runners and delete temp files
        File[] bams = new File[runners.length];
        File[] vcfs = new File[runners.length];
        for (int i=0; i< runners.length; i++) {
			if (runners[i].isFailed()) Misc.printErrAndExit("\nERROR: Failed runner, aborting! \n"+runners[i].getCmd());
			runners[i].getTempBed().deleteOnExit();
			runners[i].getTempVcf().deleteOnExit();
			runners[i].getTempLog().deleteOnExit();
			new File (runners[i].getTempVcf()+".idx").deleteOnExit();
			vcfs[i] = runners[i].getTempVcf();
			bams[i] = runners[i].getBamOut();
		}
		
		//merge vcf
		System.out.println("\nMerging vcfs...");
		mergeVcfs(vcfs, new File( saveDirectory, "gatk.raw.vcf.gz") );
		
		//merge bams?
		if (bamOut) {
			System.out.println("\nMerging bams...");
			mergeBams(bams, new File (saveDirectory, "realigned.bam"));
		}
	}
	
	private void mergeVcfs(File[] vcfFiles, File gzippedVcfOut) throws Exception {
		Gzipper out = new Gzipper(gzippedVcfOut);
		long totalRecords = 0;
		for (int i=0; i< vcfFiles.length; i++) {
			File vcfFile = vcfFiles[i];
			if (vcfFile.exists() == false) throw new IOException("ERROR: failed to find the tmp vcf for "+vcfFile);
			BufferedReader in = new BufferedReader( new FileReader(vcfFile));
			String line = null;
			while ((line = in.readLine())!= null) {
				if (line.startsWith("#")) {
					if (i==0) out.println(line);
				}
				else {
					totalRecords++;
					out.println(line);
				}
			}
			in.close();
		}
		out.close();
		System.out.println("\t"+totalRecords+ " processed vcf records.");
	}
	
	private void mergeBams(File[] bams, File mergedBam) throws Exception{
		new MergeSams(bams, mergedBam, false, true);
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
		int numberProcPerChunk = 4;
		double numberGBPerChunk = 5.0;
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
					case 'u': useUpperCaseO= true; break;
					case 'b': bamOut= true; break;
					case 'p': numberProcPerChunk = Integer.parseInt(args[++i]); break;
					case 'g': numberGBPerChunk = Double.parseDouble(args[++i]); break;
					
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//threads to use
		double gbAvail = ((double)Runtime.getRuntime().maxMemory())/ 1073741824.0;
		float procAvail = Runtime.getRuntime().availableProcessors();
		
		int numChunkGB = (int)Math.round((gbAvail-1.0)/numberGBPerChunk);
		int numChunkProc = Math.round(procAvail/(float)numberProcPerChunk);
		
		if (numChunkGB <= numChunkProc) numberThreads = numChunkGB;
		else numberThreads = numChunkProc;
		
		
		System.out.println("Core usage:\n\tGB available to Java:\t"+ Num.formatNumber(gbAvail, 1));
		System.out.println("\tAvailable cores:\t"+procAvail);
		System.out.println("\tNumber jobs to launch with a minimum "+numberGBPerChunk+"GB/Job and "+numberProcPerChunk+"Proc/Job: "+numberThreads+"\n");
		

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
				"**                               Gatk Runner: Jan 2022                              **\n" +
				"**************************************************************************************\n" +
				"Takes a bed file of target regions, splits it by the number of threads, writes out\n"+
				"each, executes the GATK commands, and merges the results. Set the -Xmx to the\n"+
				"maximum available on the machine to enable correct thread usage.\n"+

				"\nOptions:\n"+
				"-r A regions bed file (chr, start, stop,...) to intersect, see\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1 , gz/zip OK.\n"+
				"-s Path to a directory for saving the results.\n"+
				"-p Minimum number processors per job, defaults to 4\n"+
				"-g Minimum GB RAM per job, defaults to 5\n"+
				"-c Command to execute, see the example below, modify to match your enviroment.\n"+
				"     Most resources require full paths. Don't set -o, -O, -l, or -L\n"+
				"-l Use lowercased -l for Lofreq compatability, defaults to -L\n"+
				"-u Use uppercase -O for GATK 4 compatability, defaults to -o\n"+
				"-b Add a -bamout argument and merge bam chunks.\n"+
				
				"\nExample: java -Xmx128G -jar pathToUSeq/Apps/GatkRunner -b -r /SS/targets.bed -s\n" +
				"     /SS/HC/ -c 'java -Xmx4G -jar /SS/GenomeAnalysisTK.jar -T MuTect2 \n"+
				"    -R /SS/human_g1k_v37.fasta --dbsnp /SS/dbsnp_138.b37.vcf \n"+
				"    --cosmic /SS/v76_GRCh37_CosmicCodingMuts.vcf.gz -I:tumor /SS/sar.bam -I:normal \n"+
				"    /SS/normal.bam' -u \n\n" +

				"**************************************************************************************\n");

	}


	public File getSaveDirectory() {
		return saveDirectory;
	}

	public String getGatkArgs() {
		return gatkArgs;
	}

	public boolean isUseUpperCaseO() {
		return useUpperCaseO;
	}

	
}
