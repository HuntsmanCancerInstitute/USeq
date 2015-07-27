package edu.utah.seq.vcf;

import java.io.*;
import java.util.regex.*;
import util.bio.annotation.Bed;
import util.gen.*;
import java.util.*;

/**
 * @author david.nix@hci.utah.edu 
 * Takes a bed file of desired regions, splits it, writes out each, then cals the haplotype caller on each, finally merges the results.
 **/
public class HaploRunner {

	//user defined fields
	private File bedFile;
	private File saveDirectory;
	private String haploArgs = null;
	private int numberConcurrentThreads = 4;
	
	//internal fields
	HaploRunnerChunk[] runners;
	
	//constructors
	public HaploRunner(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			for (HaploRunnerChunk c: runners) c.interrupt();
			Misc.printErrAndExit("\nProblem running HaploRunner!");
		}
	}

	public void doWork() throws Exception{
		
		System.out.println("Parsing and spliting bed file...");
		
		//parse regions and randomize
		Bed[] regions = Bed.parseFile(bedFile, 0, 0);
		Misc.randomize(regions, 0);
		int numRegionsPerChunk = 1+ (int)Math.round( (double)regions.length/ (double)numberConcurrentThreads );
		
		//split regions
		ArrayList<Bed[]> splitBed = Bed.splitByNumber(regions, numRegionsPerChunk);
		
		//create runners and start
		System.out.println("Launching...");
		runners = new HaploRunnerChunk[splitBed.size()];
		for (int i=0; i< runners.length; i++) runners[i] = new HaploRunnerChunk(splitBed.get(i), this, ""+i);
		
		//check threads, will throw exception if something explodes
		checkThreads();
		
		//merge vcf
		System.out.println("\nMerging vcfs...");
		mergeVcfs();
		
		//delete temp files
		for (HaploRunnerChunk c: runners) {
			c.getTempBed().deleteOnExit();
			c.getTempVcf().deleteOnExit();
		}
	}

	private void mergeVcfs() throws Exception{
		ArrayList<VCFRecord> recordsAL = new ArrayList<VCFRecord>();
		VCFParser vp = null;
		for (HaploRunnerChunk c: runners) {
			vp = new VCFParser(c.getTempVcf(), true, false, false);
			VCFRecord[] records = vp.getVcfRecords();
			for (VCFRecord r: records) recordsAL.add(r);
		}
		VCFRecord[] all = new VCFRecord[recordsAL.size()];
		recordsAL.toArray(all);
		Arrays.sort(all);
		PrintWriter out = new PrintWriter( new FileWriter( new File( saveDirectory, "haplo.raw.vcf")));
		
		//write out header
		ArrayList<String> header = vp.getVcfComments().getOriginalComments();
		for (String h : header) out.println(h);
		
		//write out records
		for (VCFRecord r : all) out.println(r.getOriginalRecord());
		out.close();
		
		System.out.println("\t"+all.length+"\tvariants found");
	}

	private void checkThreads() throws Exception {
		int priorNumComplete = 0;
		while (true){ 
				int numComplete = 0;
				for (HaploRunnerChunk c: runners) {
					if (c.isFailed()) throw new Exception ("Failed runner "+c.getName());
					if (c.isComplete()) numComplete++;
				}
				if (priorNumComplete != numComplete){
					System.out.println("\t"+(runners.length- numComplete)+" remaining to complete");
					priorNumComplete = numComplete;
				}
				if (numComplete == runners.length) break;
				//sleep, 15 sec
				Thread.sleep(15000);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new HaploRunner(args);
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
					case 'c': haploArgs = args[++i]; break;
					case 't': numberConcurrentThreads = Integer.parseInt(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//check save dir
		if (saveDirectory == null) Misc.printErrAndExit("\nError: cannot find your save directory!\n"+saveDirectory);
		saveDirectory.mkdirs();
		if (saveDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: your save directory does not appear to be a directory?\n");

		//check bed
		if (bedFile == null || bedFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find your bed file of regions to interrogate?\n"+bedFile);
		
		//check args
		if (haploArgs == null) Misc.printErrAndExit("\nError: please provide a gatk haplotype caller launch cmd similar to the following where you "
				+ "replace the $xxx with the correct path to these resources on your system:\n'java -Xmx4G -jar $GenomeAnalysisTK.jar -T "
				+ "HaplotypeCaller -stand_call_conf 5 -stand_emit_conf 5 --min_mapping_quality_score 13 -R $fasta --dbsnp $dbsnp -I $bam'\n");
		if (haploArgs.contains("~") || haploArgs.contains("./")) Misc.printErrAndExit("\nError: full paths in the GATK command.\n"+haploArgs);
		if (haploArgs.contains("-o") || haploArgs.contains("-L")) Misc.printErrAndExit("\nError: please don't provide a -o or -L argument to the cmd.\n"+haploArgs);	
	
		//determine number of threads
		double gigaBytesAvailable = ((double)Runtime.getRuntime().maxMemory())/ 1073741824.0;
		
	
	}	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Haplo Runner: July 2015                             **\n" +
				"**************************************************************************************\n" +
				"Takes a bed file of target regions, splits it by the number of threads, writes out\n"+
				"each, executes the GATK Haplotype caller, and merges the vcf results. \n"+

				"\nRequired Options:\n"+
				"-r A regions bed file (chr, start, stop,...) to intersect, see\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1 , gz/zip OK.\n"+
				"-s Path to a directory for saving the results.\n"+
				"-t Number concurrent threads to run, will need > 4G RAM each.\n"+
				"-c GATK command to execute, see the example below, modify to match your enviroment.\n"+
				"     Most resources require full paths.\n"+

				"\nExample: java -Xmx24G -jar pathToUSeq/Apps/HaploRunner -r /SS/targets.bed -t 4 -s\n" +
				"     /SS/HC/ -c 'java -Xmx4G -jar /SS/GenomeAnalysisTK.jar -T HaplotypeCaller \n"+
				"    --downsample_to_coverage 500 -stand_call_conf 5 -stand_emit_conf 5 --min_mapping_quality_score 13 -R \n"+
				"    /SS/human_g1k_v37.fasta --dbsnp /SS/dbsnp_138.b37.vcf -I /SS/HC/sarc.bam'\n\n" +

				"**************************************************************************************\n");

	}


	public File getSaveDirectory() {
		return saveDirectory;
	}

	public String getHaploArgs() {
		return haploArgs;
	}

	
}
