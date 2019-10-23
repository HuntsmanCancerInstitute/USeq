package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.parsers.mpileup.MpileupTabixLoader;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**
 * Depreciated, moved and updated this app using internal BamPileup.  Mpileup was causing issues.
 * 
 * 
 * 
 * 
 * 
 * 
 * Given a vcf file compares each record against a set of normal samples from a tabix indexed mpileup file and marks those with evidence of the variant in the normals.
 * tabix-0.2.6/bgzip ~/repro3Ct8Norm.mpileup
 * tabix-0.2.6/tabix -s 1 -b 2 -e 2 ~/repro3Ct8Norm.mpileup.gz
 * @author Nix*/
public class VCFBkzDepreciated {
	
	//user defined fields
	private File forExtraction = null;
	private File[] vcfFiles;
	private File mpileup;
	private int[] sampleIndexesToUse = null;
	private File saveDir;
	private int minBaseQuality = 20;
	private int minReadCoverage = 20;
	private int minNumSamples = 3;
	private boolean verbose = false;
	private boolean removeNonZScoredRecords = false;
	private double minimumZScore = 0;
	private double maxSampleAF = 0.1;
	private boolean replaceQualScore = false;
	private int numberThreads = 0;
	private String afInfoName = "T_AF";
	private String dpInfoName = "T_DP";
	private int bpPad = 0;
	private boolean excludeBKAFs = false;
	
	//internal
	private static final int numVcfToProc = 100;
	private MpileupTabixLoader[] loaders;
	public static final double zscoreForInfinity = 1000;
	public static final double  zscoreLessThanZero = 0.001;

	//working
	private File vcfFile;
	private BufferedReader vcfIn;
	private int numRecords = 0;
	private int numNotScored = 0;
	private int numFailingZscore = 0;
	private int numWithBKAF = 0;
	private int numSaved = 0;
	private ArrayList<String> tooFewSamples = new ArrayList<String>();
	private ArrayList<String> vcfHeader = new ArrayList<String>();
	private ArrayList<String> vcfRecords = new ArrayList<String>();
	
	//constructor
	public VCFBkzDepreciated(String[] args){
		long startTime = System.currentTimeMillis();
		try {	
			processArgs(args);
			
			//make Loaders
			loaders = new MpileupTabixLoader[numberThreads];
			for (int i=0; i< loaders.length; i++) loaders[i] = new MpileupTabixLoader(mpileup, this);

			//for each vcf file
			IO.pl("\nParsing vcf files...");
			IO.pl("Name\tRecords\tSaved\tNotScored\tFailingBKZScore\tFailingBKAFFilter");
			for (int i=0; i< vcfFiles.length; i++){
				//set values and clear past data
				vcfFile = vcfFiles[i];
				numRecords = 0;
				numNotScored = 0;
				numFailingZscore = 0;
				numWithBKAF = 0;
				numSaved = 0;
				tooFewSamples.clear();
				vcfRecords.clear();
				
				IO.p(vcfFile.getName());
				String name = Misc.removeExtension(vcfFile.getName());
				createReaderSaveHeader();

				
				//load mpileup data
				ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
				for (MpileupTabixLoader l: loaders) executor.execute(l);
				executor.shutdown();
				while (!executor.isTerminated()) {}

				//check loaders 
				for (MpileupTabixLoader l: loaders) {
					if (l.isFailed()) throw new IOException("ERROR: File Loader issue! \n");
				}
				vcfIn.close();
				
				VcfSorter[] finalVcfs = fetchVcfSorterArray();
				Arrays.sort(finalVcfs);
				
				writeOutModifiedVcf (finalVcfs, new File(saveDir, name+".bkz.vcf.gz"));
				
				IO.pl("\t"+numRecords+"\t"+numSaved+"\t"+numNotScored+"\t"+numFailingZscore+"\t"+numWithBKAF);
			}
			
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("Problem encountered, aborting!");
		} finally {
			//shut down loaders
			for (MpileupTabixLoader m : loaders) m.getTabixReader().close();
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			IO.pl("\nDone! "+Math.round(diffTime)+" min\n");
		}
	}

	private void writeOutModifiedVcf(VcfSorter[] finalVcfs, File file) throws IOException {
		Gzipper out = new Gzipper(file);
		//add header
		for (String s: vcfHeader) out.println(s);
		for (VcfSorter v : finalVcfs) out.println(v.fields, "\t");
		out.close();
	}
	
	private VcfSorter[] fetchVcfSorterArray() throws Exception {
		int num= vcfRecords.size();

		VcfSorter[] v = new VcfSorter[num];
		for (int x = 0; x < num; x++){
			String vcf = vcfRecords.get(x);
			v[x] = new VcfSorter(Misc.TAB.split(vcf));
		}
		return v;
	}
	
	private class VcfSorter implements Comparable<VcfSorter>{
		private String chr;
		private int pos;
		//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLES
		String[] fields;
		
		public VcfSorter(String[] fields){
			this.fields = fields;
			chr = fields[0];
			pos = Integer.parseInt(fields[1]);
		}
		
		public int compareTo(VcfSorter other) {
			//sort by chromosome
			int x = chr.compareTo(other.chr);
			if (x !=0) return x;
			//sort by position
			if (pos < other.pos) return -1;
			if (pos > other.pos) return 1;
			return 0;
		}
	}

	public synchronized void save(ArrayList<String> vcf, ArrayList<int[][]> baseCounts) {
		int num = vcf.size();
		numSaved += num;
		
		//save vcf records
		vcfRecords.addAll(vcf);
	}

	/**Loads the AL with vcf lines up to the numVcfToProc, if none could be read, returns false, otherwise true.*/
	public synchronized boolean loadVcfRecords(ArrayList<String> chunk) throws IOException{
		String record;
		int count = 0;
		while ((record = vcfIn.readLine()) != null){
			record = record.trim();
			if (record.length() == 0) continue;
			chunk.add(record);
			if (count++ > numVcfToProc) break;
		}
		
		if (chunk.size() == 0) return false;
		numRecords += chunk.size();
		return true;
	}
	
	public synchronized void update(int numNotScored, int numFailingZscore, ArrayList<String> tooFewSamples, int numWithBKAF) {
		this.numNotScored += numNotScored;
		this.numFailingZscore+= numFailingZscore;
		this.tooFewSamples.addAll(tooFewSamples);
		this.numWithBKAF += numWithBKAF;
	}
	
	private void createReaderSaveHeader() throws Exception {
		vcfIn = IO.fetchBufferedReader(vcfFile);
		vcfHeader.clear();
		String record;
		boolean addInfo = true;
		boolean endFound = false;
		while ((record = vcfIn.readLine()) != null){
			record = record.trim();
			if (record.length() == 0) continue;
			if (record.startsWith("#")){
				if (addInfo && record.startsWith("##INFO=")) {
					vcfHeader.add("##FILTER=<ID=BKAF,Description=\"One or more background sample AFs are >= variant AF.\">");
					vcfHeader.add("##INFO=<ID=BKZ,Number=1,Type=Float,Description=\"Smallest AF z-score calculated from background AFs over effected bases. "
							+ "Values < ~4 are suspicous, non reference observations are likely present in the background samples.\">");
					vcfHeader.add("##INFO=<ID=BKAF,Number=1,Type=String,Description=\"Sorted list (largest to smallest) of background non-reference AFs used to calculate the BKZ.\">");					
					addInfo = false;
				}
				vcfHeader.add(record);
				if (record.startsWith("#CHROM")){
					endFound = true;
					break;
				}
			}
			else {
				endFound = true;
				break;
			}
		}
		if (endFound == false) Misc.printErrAndExit("\nERROR: failed to find the #CHROM line, aborting\n");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFBkzDepreciated(args);
	}		
	


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 's': saveDir = new File(args[++i]); break;
					case 'm': mpileup = new File(args[++i]); break;
					case 'z': minimumZScore = Double.parseDouble(args[++i]); break;
					case 'q': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'c': minReadCoverage = Integer.parseInt(args[++i]); break;
					case 'a': minNumSamples = Integer.parseInt(args[++i]); break;
					case 'b': bpPad = Integer.parseInt(args[++i]); break;
					case 'e': removeNonZScoredRecords = true; break;
					case 'd': verbose = true; break;
					case 'l': excludeBKAFs = true; break;
					case 'i': sampleIndexesToUse = Num.parseInts(args[++i], Misc.COMMA); break;
					case 'u': replaceQualScore = true; break;
					case 't': numberThreads = Integer.parseInt(args[++i]); break;
					case 'f': afInfoName = args[++i]; break;
					case 'p': dpInfoName = args[++i]; break;
					case 'x': printPoNExample(); System.exit(0);
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.fetchFilesRecursively(forExtraction, ".vcf");
		tot[1] = IO.fetchFilesRecursively(forExtraction,".vcf.gz");
		tot[2] = IO.fetchFilesRecursively(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");

		//save dir?
		if (saveDir == null) saveDir = vcfFiles[0].getParentFile();
		else saveDir.mkdirs();
		if (saveDir.exists()==false || saveDir.canWrite()==false) Misc.printExit("\nError: failed to find or create a writeable save directory? Aborting.\n");
		
		//tabix indexed mpileup?
		if (mpileup == null || mpileup.canRead() == false) Misc.printExit("\nError: please provide a path to a bgzipped tabix indexed mpileup file.\n");
		File index = new File (mpileup.toString()+".tbi");
		if (index.exists() == false) Misc.printExit("\nError: cannot find the '"+index.getName()+"' index file corresponding to this indexed mpileup file "+mpileup.getName());	
		
		//threads
		int numProc = Runtime.getRuntime().availableProcessors() - 1;
		if (numberThreads == 0 || numberThreads > numProc) numberThreads = numProc;
		
		printSettings();
		
	}	
	
	private static void printPoNExample() {
		IO.pl("#!/bin/bash");
		IO.pl("bed=/scratch/mammoth/serial/Bed/HSV1_GBM_IDT_Probes_Hg38Pad150bps_91K.bed");
		IO.pl("ref=/uufs/chpc.utah.edu/TNRunner/Indexes/Hg38IndexForBwa-0.7.17/hs38DH.fa");
		IO.pl("");
		IO.pl("echo \"#SampleOrder-index; 15352X1-0 15352X2-1 15352X3-2 15352X4-3 15352X6-4 15352X8-5 15352X9-6 15352X10-7 \\");
		IO.pl("15352X11-8 15352X12-9 15352X14-10 15352X16-11 15352X17-12 15352X18-13 15352X19-14 15352X20-15 15352X22-16 15352X24-17 \\");
		IO.pl("> 15352R.mpileup");
		IO.pl("# Be sure to use samtools 1.8 with these settings, any changes should be evaluated carefully.");
		IO.pl("~/BioApps/Samtools/1.8/bin/samtools mpileup \\");
		IO.pl("--count-orphans \\");
		IO.pl("--max-depth 100000 \\");
		IO.pl("--max-idepth 100000 \\");
		IO.pl("--gap-frac 0.001 \\");
		IO.pl("--per-sample-mF \\");
		IO.pl("--ignore-RG \\");
		IO.pl("--min-MQ 13 \\");
		IO.pl("--no-BAQ \\");
		IO.pl("--fasta-ref $ref \\");
		IO.pl("--ff UNMAP,SECONDARY,QCFAIL \\");
		IO.pl("-l $bed \\");
		IO.pl("-o 15352R.mpileup \\");
		IO.pl("15352X1_Hg38_final.bam 15352X2_Hg38_final.bam 15352X3_Hg38_final.bam 15352X4_Hg38_final.bam 15352X6_Hg38_final.bam \\");
		IO.pl("15352X8_Hg38_final.bam 15352X9_Hg38_final.bam 15352X10_Hg38_final.bam 15352X11_Hg38_final.bam 15352X12_Hg38_final.bam \\");
		IO.pl("15352X14_Hg38_final.bam 15352X16_Hg38_final.bam 15352X17_Hg38_final.bam 15352X18_Hg38_final.bam 15352X19_Hg38_final.bam \\");
		IO.pl("15352X20_Hg38_final.bam 15352X22_Hg38_final.bam 15352X24_Hg38_final.bam >> 15352R.mpileup");
		IO.pl("");
		IO.pl("~/BioApps/HTSlib/1.3/bin/bgzip --threads 20 15352R.mpileup");
		IO.pl("~/BioApps/HTSlib/1.3/bin/tabix -s 1 -b 2 -e 2 15352R.mpileup.gz");
	}

	public void printSettings(){
		IO.pl("Settings:");
		
		IO.pl(" -v Vcf file(s) "+ forExtraction);
		IO.pl(" -m Mpileup bkg "+ IO.getCanonicalPath(mpileup));
		IO.pl(" -s Save dir    "+ IO.getCanonicalPath(saveDir));
		if (sampleIndexesToUse != null) IO.pl(" -x Sample idxs "+ Misc.intArrayToString(sampleIndexesToUse, " "));
		IO.pl(" -f "+ afInfoName+ "\tTumor AF INFO name");  
		IO.pl(" -p "+ dpInfoName+ "\tTumor DP INFO name");
		IO.pl(" -c "+ minReadCoverage+"\tMin mpileup sample read coverage");
		IO.pl(" -q "+ minBaseQuality+"\tMin mpileup sample base quality");
		IO.pl(" -f "+ maxSampleAF+"\tMax mpileup sample AF");
		IO.pl(" -a "+ minNumSamples+"\tMin # samples for z-score calc");
		IO.pl(" -z "+ minimumZScore+"\tMin vcf AF z-score to save");
		IO.pl(" -b "+ bpPad+"\tBP padding +/- to vcf record for mpileup scanning");
		IO.pl(" -t "+ numberThreads+"\tCPUs");
		IO.pl(" -e "+ removeNonZScoredRecords+ "\tExclude vcf records that could not be z-scored");
		IO.pl(" -l "+ excludeBKAFs+ "\tExclude vcf records with a bkg AF >= record AF");
		IO.pl(" -u "+ replaceQualScore+"\tReplace QUAL score with z-score and set non scored records to 0");
		IO.pl(" -d "+ verbose+"\tVerbose output");
	}
	
	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                                VCF Bkz : Oct 2019                                **\n" +
				"**************************************************************************************\n" +
				"VCFBkz calculates non-reference allele frequencies (AF) from a background multi-sample \n"+
				"mpileup file over each vcf record. It then calculates a z-score for the vcf AF and \n"+
				"appends it to the INFO field. If multiple bps are affected (e.g. INDELs) or bp padding\n"+
				"provided, the lowest bp z-score is appended. Z-scores < ~3 are indicative of non\n"+
				"reference bps in the background samples. A flag is appended to the FILTER field if a\n"+
				"background AF was found >= vcf AF. Note, VCFBkz requires tumor AF and\n"+
				"DP tags in the INFO field of each record to use in the z-score calculation, see -f\n"+
				"and -p. \n"+

				"\nRequired:\n"+
				"-v Path to a xxx.vcf(.gz/.zip OK) file or directory containing such.\n" +
				"-m Path to a bgzip compressed and tabix indexed multi-sample mpileup file. See -x\n"+
						
				"\nOptional:\n" +
				"-s Path to a directory to save the modified vcf file(s), defaults to vcf parent dir\n"+
				"-f Tumor AF INFO name, defaults to T_AF\n"+
				"-p Tumor DP INFO name, defaults to T_DP\n"+
				"-z Minimum vcf z-score, defaults to 0, no filtering. Unscored vars are kept.\n"+
				"-q Minimum mpileup sample bp quality, defaults to 20\n"+
				"-c Minimum mpileup sample read coverge, defaults to 20\n"+
				"-f Maximum mpileup sample AF to include, defaults to 0.1, set to exclude germline.\n"+
				"-a Minimum # mpileup samples for z-score calculation, defaults to 3\n" +
				"-b Pad the size of each vcf record, defaults to 0 bp\n"+
				"-i Sample indexes to use, comma delimited, no spaces, zero based, defaults to all\n"+
				"-e Exclude vcf records that could not be z-scored\n"+
				"-l Exclude vcf records with a bkg AF >= the record AF.\n"+
				"-u Replace QUAL value with z-score. Un scored vars will be assigned 0\n"+
				"-d Print verbose debugging output\n" +
				"-t Number of threads to use, defaults to all\n"+
				"-x Print an example mpileup PoN bash script\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/VCFBkz -v SomaticVcfs/ -z 3\n"+
				"-m bkg.mpileup.gz -s BkgFiltVcfs/ -q 13 -u -b 2\n\n"+

		        "**************************************************************************************\n");
	}

	//getters and setters

	public int getMinBaseQuality() {
		return minBaseQuality;
	}

	public int getMinReadCoverage() {
		return minReadCoverage;
	}

	public int getMinNumSamples() {
		return minNumSamples;
	}

	public boolean isRemoveNonZScoredRecords() {
		return removeNonZScoredRecords;
	}

	public double getMinimumZScore() {
		return minimumZScore;
	}

	public double getMaxSampleAF() {
		return maxSampleAF;
	}

	public boolean isReplaceQualScore() {
		return replaceQualScore;
	}

	public boolean isVerbose() {
		return verbose;
	}

	public String getAFInfoName() {
		return afInfoName;
	}

	public String getDPInfoName() {
		return dpInfoName;
	}

	public int getBpPad() {
		return bpPad;
	}

	public int[] getSampleIndexesToUse() {
		return sampleIndexesToUse;
	}

	public boolean isExcludeBKAFs() {
		return excludeBKAFs;
	}
}
