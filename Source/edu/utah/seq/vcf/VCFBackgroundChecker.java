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

/**Given a vcf file compares each record against a set of normal samples from a tabix indexed mpileup file and marks those with evidence of the variant in the normals.
 * tabix-0.2.6/bgzip ~/repro3Ct8Norm.mpileup
 * tabix-0.2.6/tabix -s 1 -b 2 -e 2 ~/repro3Ct8Norm.mpileup.gz
 * @author Nix*/
public class VCFBackgroundChecker {
	
	//user defined fields
	private File[] vcfFiles;
	private File mpileup;
	private File saveDir;
	private int minBaseQuality = 20;
	private int minReadCoverage = 20;
	private int minNumSamples = 3;
	private boolean verbose = false;
	private boolean removeNonZScoredRecords = false;
	private double minimumZScore = 0;
	private double maxSampleAF = 0.3;
	private boolean replaceQualScore = false;
	private int numberThreads = 0;
	private String afInfoName = "T_AF";
	private String dpInfoName = "T_DP";
	private int bpPad = 0;
	
	//internal
	private static final int numVcfToProc = 100;
	private MpileupTabixLoader[] loaders;
	public static final double zscoreForInfinity = 1000;
	public static final double  zscoreLessThanZero = 0.001;

	//working
	private File vcfFile;
	private BufferedReader vcfIn;
	private PrintWriter rOut;
	private int numRecords = 0;
	private int numNotScored = 0;
	private int numFailingZscore = 0;
	private int numSaved = 0;
	private ArrayList<String> tooFewSamples = new ArrayList<String>();
	private ArrayList<String> vcfHeader = new ArrayList<String>();
	private ArrayList<String> vcfRecords = new ArrayList<String>();
	
	//constructor
	public VCFBackgroundChecker(String[] args){
		long startTime = System.currentTimeMillis();
		try {	
			processArgs(args);
			
			//make Loaders
			loaders = new MpileupTabixLoader[numberThreads];
			for (int i=0; i< loaders.length; i++) loaders[i] = new MpileupTabixLoader(mpileup, this);
			
			//make Alun's BM caller
//BinomialMixture bm = new BinomialMixture(fullPathToR, saveDir, deleteTempFiles);

			//for each vcf file
			IO.pl("\nParsing vcf files...");
			for (int i=0; i< vcfFiles.length; i++){
				//set values and clear past data
				vcfFile = vcfFiles[i];
				numRecords = 0;
				numNotScored = 0;
				numFailingZscore = 0;
				numSaved = 0;
				tooFewSamples.clear();
				vcfRecords.clear();
				
				IO.pl(vcfFile.getName());
				String name = Misc.removeExtension(vcfFile.getName());
				createReaderSaveHeader();
				File tempRData = new File(saveDir, name+"_RObs.txt");
//if (deleteTempFiles) tempRData.deleteOnExit();
				tempRData.deleteOnExit();
				rOut = new PrintWriter( new FileWriter(tempRData));
				
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
				rOut.close();
				
				//estimate pvalue of tumor obs vs panel of normals data
				//disable for now, just sort
				
//double[] pvals = bm.estimatePValues(tempRData);
//Misc.printArray(pvals);
//if (pvals == null) throw new Exception("Error estimating PON pvalues!");
//VcfSorter[] finalVcfs = addPVals(pvals);
				
				VcfSorter[] finalVcfs = fetchVcfSorterArray();
				Arrays.sort(finalVcfs);
				
				writeOutModifiedVcf (finalVcfs, new File(saveDir, name+".vcf.gz"));
				
				//print summary stats
				/*if (tooFewSamples.size() !=0){
					System.err.println("WARNING, the following did not have any mpileup lines (e.g. outside the compiled regions?) or enough background "
							+ "samples passing thresholds to calculate a z-score. ");
					System.err.println(Misc.stringArrayListToString(tooFewSamples, "\n"));
				}*/
				IO.pl("\t#Rec="+numRecords+" #Saved="+numSaved+" #NotScored="+numNotScored+" #FailingBKZ="+numFailingZscore+"\n");
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

	/*
	private VcfSorter[] addPVals(double[] pvals) throws Exception {
		if (pvals.length != this.vcfRecords.size()) throw new Exception ("The number of pvalues and number of modified vcf records do not match?");

		VcfSorter[] v = new VcfSorter[pvals.length];
		for (int x = 0; x < pvals.length; x++){
			String vcf = vcfRecords.get(x);
			String pval = Num.formatNumberJustMax(pvals[x], 1);
			
			//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLES
			String[] fields = Misc.TAB.split(vcf);
			fields[7] = "BKP="+ pval+ ";" + fields[7];
			v[x] = new VcfSorter(fields);
		}
		return v;
	}*/
	
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
		
		//write out data for R
//disable!
		/*
		for (int i =0; i< num; i++){
			rOut.println("#"+vcf.get(i));
			int[][] bc = baseCounts.get(i);
			String nonRef = Misc.intArrayToString(bc[0], "\t");
			String total = Misc.intArrayToString(bc[1], "\t");
			rOut.println(nonRef);
			rOut.println(total);
		}*/
		
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
	
	public synchronized void update(int numNotScored, int numFailingZscore, ArrayList<String> tooFewSamples) {
		this.numNotScored += numNotScored;
		this.numFailingZscore+= numFailingZscore;
		this.tooFewSamples.addAll(tooFewSamples);
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
					vcfHeader.add("##FILTER=<ID=BKAF,Description=\"One or more background control sample AFs are >= 0.9 x variant AF.\">");
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
		new VCFBackgroundChecker(args);
	}		
	


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
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
					case 'u': replaceQualScore = true; break;
					case 't': numberThreads = Integer.parseInt(args[++i]); break;
					case 'f': afInfoName = args[++i]; break;
					case 'p': dpInfoName = args[++i]; break;
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
		if (saveDir == null) Misc.printExit("\nError: please provide a directory in which to save the marked vcf files\n");
		saveDir.mkdirs();
		if (saveDir.exists()==false || saveDir.canWrite()==false) Misc.printExit("\nError: failed to find or create a writeable save directory? Aborting.\n");
		
		//tabix indexed mpileup?
		if (mpileup == null || mpileup.canRead() == false) Misc.printExit("\nError: please provide a path to a bgzipped tabix indexed mpileup file.\n");
		File index = new File (mpileup.toString()+".tbi");
		if (index.exists() == false) Misc.printExit("\nError: cannot find the '"+index.getName()+"' index file corresponding to this indexed mpileup file "+mpileup.getName());	
		
		//threads to use
		int numAvail = Runtime.getRuntime().availableProcessors();
		if (numberThreads < 1 || numberThreads > numAvail) numberThreads =  numAvail - 1;
		
		printSettings();
		
	}	
	

	
	public void printSettings(){
		IO.pl("Settings:\nBackground\t"+mpileup);
		IO.pl("Save dir\t"+saveDir);
		IO.pl(afInfoName+ "\tTumor AF INFO name");
		IO.pl(dpInfoName+ "\tTumor DP INFO name");
		IO.pl(minReadCoverage+"\tMin mpileup sample read coverage");
		IO.pl(minBaseQuality+"\tMin mpileup sample base quality");
		IO.pl(maxSampleAF+"\tMax mpileup sample AF");
		IO.pl(minNumSamples+"\tMin # samples for z-score calc");
		IO.pl(minimumZScore+"\tMin vcf AF z-score to save");
		IO.pl(bpPad+"\tBP padding +/- to vcf record for mpileup scanning");
		IO.pl(numberThreads+"\tCPUs");
		IO.pl(removeNonZScoredRecords+ "\tExclude vcf records that could not be z-scored");
		IO.pl(verbose+"\tVerbose");
		IO.pl(replaceQualScore+"\tReplace QUAL score with z-score and set non scored records to 0");
	}

	
	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                         VCF Background Checker : May 2019                        **\n" +
				"**************************************************************************************\n" +
				"VBC calculates non-reference allele frequencies (AF) from a background multi-sample \n"+
				"mpileup file over each vcf record. It then calculates a z-score for the vcf AF and \n"+
				"appends it to the INFO field. If multiple bps are affected (e.g. INDELs) or bp padding\n"+
				"provided, the lowest bp z-score is appended. Z-scores < ~4 are indicative of non\n"+
				"reference bps in the background samples. A flag is appended the FILTER field if a\n"+
				"background AF was found within 10% of the vcf AF. Note, VBC requires AF and DP tags\n"+
				"in the INFO field of each record to use in the z-score calculation, see -f and -p. \n"+

				"\nRequired:\n"+
				"-v Path to a xxx.vcf(.gz/.zip OK) file or directory containing such.\n" +
				"-m Path to a bgzip compressed and tabix indexed multi-sample mpileup file. e.g.:\n"+
				"      1) Mpileup: 'echo \"#SampleOrder: \"$(ls *bam) > bkg.mpileup; samtools mpileup\n"+
				"             -B -q 13 -d 1000000 -f $fastaIndex -l $bedFile *.bam >> bkg.mpileup'\n"+
				"      2) (Optional) MpileupRandomizer: java -jar -Xmx10G ~/USeqApps/MpileupRandomizer\n"+
				"             -r 20 -s 3 -m bkg.mpileup\n"+
				"      3) Bgzip: 'tabix-0.2.6/bgzip bkg.mpileup_DP20MS3.txt'\n"+
                "         Tabix: 'tabix-0.2.6/tabix -s 1 -b 2 -e 2 bkg.mpileup_DP20MS3.txt.gz'\n"+
				"-s Path to directory in which to save the modified vcf file(s)\n"+
						
				"\nOptional:\n" +
				"-f Tumor AF INFO name, defaults to T_AF\n"+
				"-p Tumor DP INFO name, defaults to T_DP\n"+
				"-z Minimum vcf z-score, defaults to 0, no filtering. Unscored vars are kept.\n"+
				"-q Minimum mpileup sample bp quality, defaults to 20\n"+
				"-c Minimum mpileup sample read coverge, defaults to 20\n"+
				"-f Maximum mpileup sample AF, defaults to 0.3\n"+
				"-a Minimum # mpileup samples for z-score calculation, defaults to 3\n" +
				"-b Pad the size of each vcf record, defaults to 0 bp\n"+
				"-e Exclude vcf records that could not be z-scored\n"+
				"-u Replace QUAL value with z-score. Un scored vars will be assigned 0\n"+
				"-d Print verbose debugging output\n" +
				"-t Number of threads to use, defaults to all\n"+
				

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/VCFBackgroundChecker -v SomaticVcfs/ -z 3\n"+
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
}
