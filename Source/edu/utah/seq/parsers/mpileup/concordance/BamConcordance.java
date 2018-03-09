package edu.utah.seq.parsers.mpileup.concordance;

import java.io.*;
import java.util.regex.*;

import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
import util.bio.annotation.Bed;
import util.bio.parsers.UCSCGeneLine;
import util.gen.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class BamConcordance {

	//user defined fields
	private File bedFile;
	private File commonSnvBed;
	private File samtools;
	private File fasta;
	private File[] bamFiles;
	private int minSnvDP = 25;
	private double minAFForHom = 0.95;
	private double minAFForMatch = 0.90;
	private double minAFForHis = 0.05;
	private int minBaseQuality = 20;
	private int minMappingQuality = 20;

	//internal fields
	private File tempDirectory;
	private int numberThreads = 0;
	private double maxIndel = 0.01;
	private Histogram[] afHist;
	private Histogram[] chrXAfHist;
	private Similarity[] similarities;
	private String[] sampleNames = null;
	private String bamNames = null;
	
	//internal fields
	ConcordanceChunk[] runners;
	

	
	//constructors
	public BamConcordance(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			printThresholds();
			
			doWork();

			printStats();
			
			printGenderRatios();

			printHist();
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running Bam Concordance Calculator!");
		}
	}

	public void doWork() throws Exception{
		System.out.println("\nParsing and spliting bed files...");
		
		//parse regions and sort
		Bed[] regions = Bed.parseFile(bedFile, 0, 0);
		Arrays.sort(regions);
		int numRegionsPerChunk = 1+ (int)Math.round( (double)regions.length/ (double)numberThreads );
		
		//split regions
		ArrayList<Bed[]> splitBed = Bed.splitByNumber(regions, numRegionsPerChunk);
		
		//create runners and start
		System.out.println("\nLaunching "+splitBed.size() +" jobs...");
		runners = new ConcordanceChunk[splitBed.size()];
		ExecutorService executor = Executors.newFixedThreadPool(runners.length);
		
		for (int i=0; i< runners.length; i++) {
			runners[i] = new ConcordanceChunk(splitBed.get(i), this, ""+i);
			executor.execute(runners[i]);
		}
		executor.shutdown();
		//spins here until the executer is terminated, e.g. all threads complete
        while (!executor.isTerminated()) {
        }
		
		//check runners and pull bed gzippers
        ArrayList<File> toMerge = new ArrayList<File>();
        for (ConcordanceChunk c: runners) {
			if (c.isFailed()) Misc.printErrAndExit("\nERROR: Failed runner, aborting! \n"+c.getCmd());
			c.getTempBed().delete();
			toMerge.add(c.getMisMatchBed().getGzipFile());
		}
        
        //write out mismatch bed file
        String name= Misc.removeExtension(bamFiles[0].getParentFile().getName());
        File misMatchBed = new File(bamFiles[0].getParentFile(), name+"_MisMatch.bed.gz");
        IO.concatinateFiles(toMerge, misMatchBed);
        
        aggregateStats(runners);
        
		//delete temp dir and any remaining files
		IO.deleteDirectory(tempDirectory);
	}

	private void aggregateStats(ConcordanceChunk[] runners) throws Exception {
        //collect stats
    	afHist = runners[0].getAfHist();
    	chrXAfHist = runners[0].getChrXAfHist();
    	similarities = runners[0].getSimilarities();
    	    	
    	//for each runner
    	for (int i=1; i< runners.length; i++){
    		
    		//histograms
    		Histogram[] a = runners[i].getAfHist();
    		Histogram[] b = runners[i].getChrXAfHist();
    		for (int x=0; x< afHist.length; x++) {
    			afHist[x].addCounts(a[x]);
    			chrXAfHist[x].addCounts(b[x]);
    		}
    	
    		//similarities
    		Similarity[] s = runners[i].getSimilarities();
    		for (int x=0; x< s.length; x++) similarities[x].add(s[x]);
    		
    		runners[i].getMisMatchBed().getGzipFile().delete();
    	}
	}

	public void printThresholds(){
		IO.p("Settings:");
		IO.p(minSnvDP+"\tMinSnvDP");
		IO.p(minAFForHom+"\tMinAFForHom");
		IO.p(minAFForMatch+"\tminAFForMatch");
		IO.p(minAFForHis+"\tMinAFForHis");
		IO.p(minBaseQuality+"\tMinBaseQuality");
		IO.p(minMappingQuality+"\tMinMappingQuality");
		IO.p(maxIndel+"\tMaxIndel");
		if (commonSnvBed!= null) IO.p(commonSnvBed.getName()+ "\tCommon SNV exclusion file");
		else IO.p("All SNVs counted, common and uncommon.");
	}

	public void printStats(){
		//calculate max match and sort
		for (int i=0; i< similarities.length; i++) similarities[i].calculateMaxMatch();
		Arrays.sort(similarities);
		System.out.println("\nStats:");
		for (int i=0; i< similarities.length; i++) System.out.println(similarities[i].toString(sampleNames));
	}

	public void printGenderRatios(){
		IO.p("Het/Hom histogram AF count ratios for AllChrs, ChrX, log2(All/X)");
		for (int i=0; i< afHist.length; i++){
			double allChr = ratioCenterVsLast(afHist[i]);
			double chrX = ratioCenterVsLast(chrXAfHist[i]);
			double lgrto = Num.log2(allChr/chrX);
			IO.p(sampleNames[i] +"\t"+Num.formatNumber(allChr, 3)+"\t"+Num.formatNumber(chrX, 3)+"\t"+Num.formatNumber(lgrto, 3));
		}
	}
	
	/**Sums the num of middle and end bins, returns the ratio. Hard coded!*/
	public double ratioCenterVsLast(Histogram his){
		long[] counts = his.getBinCounts();
		
		//0.895 - 1.0
		double end = 0;
		end+= counts[22];
		end+= counts[23];
		end+= counts[24];
		if (end ==0) end++;
		
		//0.396 - 0.588
		double middle = 0;
		middle+= counts[9];
		middle+= counts[10];
		middle+= counts[11];
		middle+= counts[12];
		
		return middle/end;
	}
	
	public void printHist(){
		IO.p("\nHistograms of AFs observed:");
		for (int i=0; i< afHist.length; i++){
			System.out.println("Sample\t"+sampleNames[i]);
			System.out.println("All NonRef AFs >= "+minAFForHis+":");
			afHist[i].printScaledHistogram();
			System.out.println("ChrX NonRef AFs >= "+minAFForHis+" - gender check:");
			chrXAfHist[i].printScaledHistogram();
			System.out.println();
		}
		IO.p("\nHistograms of the AFs for mis matched homozygous variants - contamination check:");
		for (int i=0; i< similarities.length; i++) {
			similarities[i].toStringMismatch(sampleNames);
			System.out.println();
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BamConcordance(args);
	}		
	
	/**This method will process each argument and assign new varibles
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
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
					case 'c': commonSnvBed = new File(args[++i]); break;
					case 's': samtools = new File(args[++i]); break;
					case 'f': fasta = new File(args[++i]); break;
					case 'b': bamFiles = IO.extractFiles(new File(args[++i]), ".bam"); break;
					case 'd': minSnvDP = Integer.parseInt(args[++i]); break; 
					case 'a': minAFForHom = Double.parseDouble(args[++i]); break;
					case 'm': minAFForMatch = Double.parseDouble(args[++i]); break;
					case 'q': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'u': minMappingQuality = Integer.parseInt(args[++i]); break;
					case 't': numberThreads = Integer.parseInt(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//check bed
		if (bedFile == null || bedFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find your bed file of regions to interrogate? "+bedFile);
		//check samtools
		if (samtools == null || samtools.canExecute() == false) Misc.printErrAndExit("\nError: cannot find your samtools executable? "+samtools);
		//check bed
		if (fasta == null || fasta.canRead() == false) Misc.printErrAndExit("\nError: cannot find your indexed fasta reference file "+fasta);
		//check bams
		if (bamFiles == null || bamFiles.length == 0) Misc.printErrAndExit("\nError: cannot find your bam files to parse? ");
		//threads to use
		double gigaBytesAvailable = ((double)Runtime.getRuntime().maxMemory())/ 1073741824.0;
		int numPossCores = (int)Math.round(gigaBytesAvailable/5.0);
		if (numPossCores < 1) numPossCores = 1;
		int numPossThreads = Runtime.getRuntime().availableProcessors();
		if (numberThreads == 0){
			if (numPossCores <= numPossThreads) numberThreads = numPossCores;
			else numberThreads = numPossThreads;
			System.out.println("Core usage:\n\tTotal GB available to Java:\t"+ Num.formatNumber(gigaBytesAvailable, 1));
			System.out.println("\tTotal available cores:\t"+numPossThreads);
			System.out.println("\tNumber cores to use @ 5GB/core:\t"+numberThreads+"\n");
		}
				
		//make temp dir
		File workingDir = new File (System.getProperty("user.dir"));
		tempDirectory = new File (workingDir, "TempDirBCC_"+Misc.getRandomString(6));
		tempDirectory.mkdirs();
		
		//sample names
		sampleNames = new String[bamFiles.length];
		StringBuilder sb = new StringBuilder();
		for (int i=0; i< bamFiles.length; i++) {
			sampleNames[i] = bamFiles[i].getName().replace(".bam", "");
			sb.append(bamFiles[i].getCanonicalPath()); sb.append(" ");
		}
		bamNames = sb.toString();
	
	}	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Bam Concordance: Feb  2018                         **\n" +
				"**************************************************************************************\n" +
				"BC calculates sample level concordance based on uncommon homozygous SNVs found in bam\n"+
				"files. Samples from the same person will show high similarity (>0.9). Run BC on\n"+
				"related sample bams (e.g tumor & normal exomes, tumor RNASeq) plus an unrelated bam\n"+
				"for comparison. Mismatches passing filters are written to file. BC also generates a\n"+
				"variety of AF histograms for checking gender and sample contamination. Although\n"+
				"threaded, BC runs slowly with more that a few bams. Use the USeq ClusterMultiSampleVCF\n"+
				"app to check large batches of vcfs to identify the likely mismatched sample pairs.\n\n"+
				
				"WARNING! Mpileup does not check that the chr order is the same across samples and the\n"+
				"fasta reference, be certain they are or mpileup will produce incorrect counts. Use\n"+
				"Picard's ReorderSam app if in doubt.\n\n"+
				
				"Note re FFPE derived RNASeq data: A fair bit of systematic error is found in these\n"+
				"datasets.  As such, the RNA-> DNA contrasts are low. Yet the DNA->RNA are > 0.9\n\n"+

				"Options:\n"+
				"-r Path to a bed file of regions to interrogate.\n"+
				"-s Path to the samtools executable.\n"+
				"-f Path to an indexed reference fasta file.\n"+
				"-b Path to a directory containing indexed bam files.\n"+
				"-c Path to a tabix indexed bed file of common dbSNPs. Download 00-common_all.vcf.gz \n"+
				"       from ftp://ftp.ncbi.nih.gov/snp/organisms/, run VCF2Bed, bgzip and tabix it\n"+
				"       with https://github.com/samtools/htslib, defaults to no exclusion from calcs. \n"+
				"-d Minimum read depth, defaults to 25.\n"+
				"-a Minimum allele frequency to count as a homozygous variant, defaults to 0.95\n"+
				"-m Minimum allele frequency to count a homozygous match, defaults to 0.9\n"+
				"-q Minimum base quality, defaults to 20.\n"+
				"-u Minimum mapping quality, defaults to 20.\n"+
				"-t Number of threads to use.  If not set, determines this based on the number of\n"+
				"      threads and memory available to the JVM so set the -Xmx value to the max.\n\n"+

				"Example: java -Xmx120G -jar pathTo/USeq/Apps/BamConcordance -r ~/exomeTargets.bed\n"+
				"      -s ~/Samtools1.3.1/bin/samtools -b ~/Patient7Bams -d 10 -a 0.9 -m 0.8 -f\n"+
				"      ~/B37/human_g1k_v37.fasta -c ~/B37/b38ComSnps.bed.gz\n\n" +

				"**************************************************************************************\n");

	}

	public File getSamtools() {
		return samtools;
	}

	public int getMinSnvDP() {
		return minSnvDP;
	}

	public double getMinAFForHom() {
		return minAFForHom;
	}

	public double getMinAFForMatch() {
		return minAFForMatch;
	}

	public double getMinAFForHis() {
		return minAFForHis;
	}

	public int getMinBaseQuality() {
		return minBaseQuality;
	}

	public File getTempDirectory() {
		return tempDirectory;
	}

	public double getMaxIndel() {
		return maxIndel;
	}

	public File getFasta() {
		return fasta;
	}

	public String getBamNames() {
		return bamNames;
	}

	public File[] getBamFiles() {
		return bamFiles;
	}

	public int getMinMappingQuality() {
		return minMappingQuality;
	}

	public File getCommonSnvBed() {
		return commonSnvBed;
	}

	
}
