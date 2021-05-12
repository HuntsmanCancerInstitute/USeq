package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.parsers.jpileup.BamPileupTabixLoader;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Given a vcf file compares each record against a set of normal samples from a tabix indexed BamPileup 
 * USeq app file and marks those with evidence of the variant in the normals.
 * @author Nix*/
public class VCFBkz {
	
	//user defined fields
	private File forExtraction = null;
	private File[] vcfFiles;
	private File bpileup;
	private HashSet<Integer> sampleIndexesToExclude = null;
	private File saveDir;
	private int minReadCoverage = 100;
	private int minNumSamples = 3;
	private boolean removeNonZScoredRecords = false;
	private double minimumZScore = 0;
	private double maxSampleAF = 0.1;
	private boolean replaceQualScore = false;
	private int numberThreads = 0;
	private String afInfoName = "T_AF";
	private boolean excludeBKAFs = false;
	private boolean calcAllBkgAfs = false;
	private int bpIndelPad = 0;
	private boolean includeIndelCountsForSnvs = false;
	
	//internal
	private static final int numVcfToProc = 100;
	private BamPileupTabixLoader[] loaders;
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
	public VCFBkz(String[] args){
		long startTime = System.currentTimeMillis();
		try {	
			processArgs(args);
			
			//make Loaders
			loaders = new BamPileupTabixLoader[numberThreads];
			for (int i=0; i< loaders.length; i++) loaders[i] = new BamPileupTabixLoader(bpileup, this);

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

				
				//load bpileup data
				ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
				for (BamPileupTabixLoader l: loaders) executor.execute(l);
				executor.shutdown();
				while (!executor.isTerminated()) {}

				//check loaders 
				for (BamPileupTabixLoader l: loaders) {
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
		} finally {
			//shut down loaders
			for (BamPileupTabixLoader m : loaders) m.getTabixReader().close();
			
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

	public synchronized void save(ArrayList<String> vcf) {
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
			if (record.length() == 0 || record.startsWith("#")) continue;
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
					vcfHeader.add("##FILTER=<ID=BKAF,Description=\"Two or more background sample AFs are >= variant AF.\">");
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
		new VCFBkz(args);
	}		
	


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		String sampleIndexes = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 's': saveDir = new File(args[++i]); break;
					case 'b': bpileup = new File(args[++i]); break;
					case 'z': minimumZScore = Double.parseDouble(args[++i]); break;
					case 'f': maxSampleAF = Double.parseDouble(args[++i]); break;
					case 'c': minReadCoverage = Integer.parseInt(args[++i]); break;
					case 'n': minNumSamples = Integer.parseInt(args[++i]); break;
					case 'x': bpIndelPad = Integer.parseInt(args[++i]); break;
					case 'e': removeNonZScoredRecords = true; break;
					case 'l': excludeBKAFs = true; break;
					case 'i': sampleIndexes = args[++i]; break;
					case 'u': replaceQualScore = true; break;
					case 'y': includeIndelCountsForSnvs = true; break;
					case 'a': calcAllBkgAfs = true; break;
					case 'p': numberThreads = Integer.parseInt(args[++i]); break;
					case 't': afInfoName = args[++i]; break;
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
		tot[0] = IO.extractFiles(forExtraction, ".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");

		//save dir?
		if (saveDir == null) saveDir = vcfFiles[0].getParentFile();
		else saveDir.mkdirs();
		if (saveDir.exists()==false || saveDir.canWrite()==false) Misc.printExit("\nError: failed to find or create a writeable save directory? Aborting.\n");
		
		//tabix indexed bpileup?
		if (bpileup == null || bpileup.canRead() == false) Misc.printExit("\nError: please provide a path to a bgzipped tabix indexed bpileup file.\n");
		File index = new File (bpileup.toString()+".tbi");
		if (index.exists() == false) Misc.printExit("\nError: cannot find the '"+index.getName()+"' index file corresponding to this indexed bpileup file "+bpileup.getName());	
		
		//threads
		int numProc = Runtime.getRuntime().availableProcessors() - 1;
		if (numberThreads == 0 || numberThreads > numProc) numberThreads = numProc;
		
		//sample indexes?
		if (sampleIndexes != null) {
			int[] idx = Num.stringArrayToInts(sampleIndexes, ",");
			sampleIndexesToExclude = new HashSet<Integer>();
			for (Integer i: idx) sampleIndexesToExclude.add(i);
		}
		
		printSettings();
		
	}	
	
	
	public void printSettings(){
		IO.pl("Settings:");
		
		IO.pl(" -v Vcf file(s) "+ forExtraction);
		IO.pl(" -b Bpileup bkg "+ IO.getCanonicalPath(bpileup));
		IO.pl(" -s Save dir    "+ IO.getCanonicalPath(saveDir));
		if (sampleIndexesToExclude != null) IO.pl(" -i Sample idxs exclude "+ sampleIndexesToExclude);
		IO.pl(" -t "+ afInfoName+ "\tTumor AF INFO name");  
		IO.pl(" -c "+ minReadCoverage+"\tMin bpileup sample read coverage");
		IO.pl(" -f "+ maxSampleAF+"\tMax bpileup sample AF");
		IO.pl(" -n "+ minNumSamples+"\tMin # samples for z-score calc");
		IO.pl(" -z "+ minimumZScore+"\tMin vcf AF z-score to save");
		IO.pl(" -p "+ numberThreads+"\tCPUs");
		IO.pl(" -a "+ calcAllBkgAfs+ "\tCalculate z-scores using the sum of all classes that pass maxSampleAF");
		IO.pl(" -e "+ removeNonZScoredRecords+ "\tExclude vcf records that could not be z-scored");
		IO.pl(" -l "+ excludeBKAFs+ "\tExclude vcf records with 2 or more bkg AF >= record AF");
		IO.pl(" -u "+ replaceQualScore+"\tReplace QUAL score with z-score and set non scored records to 0");
		IO.pl(" -x "+ bpIndelPad+"\tBP padding for scanning INDEL backgrounds");
		IO.pl(" -y "+ includeIndelCountsForSnvs+"\tInclude INDEL counts when scoring SNV bkgs");
	}

	
	//FP flags with tumor sample:
	//High N frequency - zscore, 3 vs 5 mer - snvs, what depth? 1000
	//Both Ins Del  - for ins #ins/#indels, for del #del/#indels
	//Proximal calls 
	
	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                                VCF Bkz : May 2021                                **\n" +
				"**************************************************************************************\n" +
				"VCFBkz uses a panel of normals to calculate non-germline allele frequencies (AF) from \n"+
				"each sample that intersect a vcf record. It then calculates a z-score based on the vcf\n"+
				"AF and appends it to the INFO field. If multiple bps are affected (e.g. deletion), the \n"+
				"lowest z-score is appended. Z-scores < ~3 are indicative of false positives. A flag is \n"+
				"appended to the FILTER field if 2 or more background AF were found >= vcf AF*.9 . Note,\n"+
				"VCFBkz requires a tumor AF tag in the INFO field of each record to use in the z-score\n"+
				"calculation. Be sure to vt normalize, decompose, decompose_blocksub your vcf files, see\n"+
				"https://genome.sph.umich.edu/wiki/Vt#Normalization\n"+

				"\nRequired:\n"+
				"-v Path to a vt normalized xxx.vcf(.gz/.zip OK) file or directory containing such.\n" +
				"-b Path to a bgzip compressed and tabix indexed multi normal sample bpileup file. See\n"+
				"      the USeq BamPileup app.\n"+
						
				"\nOptional:\n" +
				"-s Path to a directory to save the modified vcf file(s), defaults to vcf parent dir.\n"+
				"-t Tumor AF INFO name, defaults to T_AF.\n"+
				"-a Calculate bkg AFs using a sum of all of the base and indel observations < max \n"+
				"      sample AF, defaults to just calculating bkg AFs using the vcf alt allele.\n"+
				"-z Minimum vcf z-score, defaults to 0, no filtering. Unscored vars are kept.\n"+
				"-c Minimum bpileup sample read coverge, defaults to 100\n"+
				"-f Maximum bpileup sample AF to include, defaults to 0.1, set to exclude germline.\n"+
				"-n Minimum # bpileup samples for z-score calculation, defaults to 3\n" +
				"-i Sample indexes to exclude, comma delimited, no spaces, zero based, defaults to all,\n"+
				"      see the bpileup file header for bam sample order.\n"+
				"-e Exclude vcf records that could not be z-scored, recommended, indicative of high bkg.\n"+
				"-l Exclude vcf records with two or more bkg AFs >= the record AF*0.9, recommended.\n"+
				"-u Replace QUAL value with z-score. Un scored vars will be assigned 0\n"+
				"-x Bp padding for scanning INDEL bkg, defaults to 0\n"+
				"-y Include INDEL counts when scoring SNVs to down weight SNV calls that overlap.\n"+
				"-p Number of processors to use, defaults to all\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/VCFBkz -v SomaticVcfs/ -z 3 -u -l -e -n 8\n"+
				"      -b bkgPoN.bpileup.gz -s BkgFiltVcfs/ -x 3 -y\n\n"+

		        "**************************************************************************************\n");
	}

	//getters and setters
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
	public String getAFInfoName() {
		return afInfoName;
	}
	public HashSet<Integer> getSampleIndexesToExclude() {
		return sampleIndexesToExclude;
	}
	public boolean isExcludeBKAFs() {
		return excludeBKAFs;
	}

	public boolean isCalcAllBkgAfs() {
		return calcAllBkgAfs;
	}

	public int getBpIndelPad() {
		return bpIndelPad;
	}

	public boolean isIncludeIndelCountsForSnvs() {
		return includeIndelCountsForSnvs;
	}
}
