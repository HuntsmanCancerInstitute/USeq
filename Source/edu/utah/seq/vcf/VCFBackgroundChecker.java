package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.parsers.mpileup.MpileupLine;
import edu.utah.seq.parsers.mpileup.MpileupSample;
import edu.utah.seq.parsers.mpileup.MpileupTabixLoader;
import edu.utah.seq.query.QueryIndexFileLoader;
import htsjdk.tribble.readers.TabixReader;
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
	private HashSet<Integer> sampleIndexesToExamine = null;
	private int bpBuffer = 0;
	private int minBaseQuality = 20;
	private int minReadCoverage = 20;
	private int minNumSamples = 3;
	private boolean verbose = false;
	private boolean removeNonZScoredRecords = false;
	private double minimumZScore = 0;
	private double maxSampleAF = 0.3;
	private boolean replaceQualScore = false;
	private int numberThreads = 0;
	
	//internal
	private static final int numVcfToProc = 100;
	private MpileupTabixLoader[] loaders;

	//working
	private File vcfFile;
	private BufferedReader vcfIn;
	private Gzipper vcfOut;
	private int numRecords = 0;
	private int numNotScored = 0;
	private int numFailingZscore = 0;
	private int numSaved = 0;
	private ArrayList<String> tooFewSamples = new ArrayList<String>();
	
	//constructor
	public VCFBackgroundChecker(String[] args){
		try {	
			processArgs(args);
			
			//make Loaders
			loaders = new MpileupTabixLoader[numberThreads];
			for (int i=0; i< loaders.length; i++) loaders[i] = new MpileupTabixLoader(mpileup, this);
			

			//for each vcf file
			System.out.println("\nParsing vcf files...");
			for (int i=0; i< vcfFiles.length; i++){
				vcfFile = vcfFiles[i];
				numRecords = 0;
				numNotScored = 0;
				numFailingZscore = 0;
				numSaved = 0;
				tooFewSamples.clear();
				System.out.println(vcfFile.getName());
				String name = Misc.removeExtension(vcfFile.getName());
				vcfOut = new Gzipper(new File(saveDir, name+"_BKZed.vcf.gz"));
				createReaderSaveHeader();
				
				ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
				for (MpileupTabixLoader l: loaders) executor.execute(l);
				executor.shutdown();

				//spins here until the executer is terminated, e.g. all threads complete
				while (!executor.isTerminated()) {}

				//check loaders 
				for (MpileupTabixLoader l: loaders) {
					if (l.isFailed()) throw new IOException("ERROR: File Loader issue! \n");
				}
				
				vcfOut.close();
				vcfIn.close();
				
				//print summary stats
				if (tooFewSamples.size() !=0){
					System.err.println("WARNING, the following did not have any mpileup lines (e.g. outside the compiled regions?) or enough background "
							+ "samples passing thresholds to calculate a z-score. ");
					System.err.println(Misc.stringArrayListToString(tooFewSamples, "\n"));
				}
				System.out.println("\t#Rec="+numRecords+" #Saved="+numSaved+" #NotScored="+numNotScored+" #FailingBKZ="+numFailingZscore+"\n");
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			//shut down loaders
			for (MpileupTabixLoader m : loaders) m.getTabixReader().close();
		}
	}
	
	/**Just prints out records.*/
	public synchronized void saveModifiedVcf(ArrayList<String> vcf) throws IOException{
		for (String s: vcf) vcfOut.println(s);
		numSaved+= vcf.size();
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
		String record;
		boolean addInfo = true;
		boolean endFound = false;
		while ((record = vcfIn.readLine()) != null){
			record = record.trim();
			if (record.length() == 0) continue;
			if (record.startsWith("#")){
				if (addInfo && record.startsWith("##INFO=")) {
					vcfOut.println("##INFO=<ID=BKZ,Number=1,Type=Float,Description=\"Smallest AF z-score calculated from background AFs over effected bases. "
							+ "Values < ~4 are suspicous, non reference observations are likely present in the background samples.\">");
					vcfOut.println("##INFO=<ID=BKAF,Number=1,Type=String,Description=\"Non-reference AFs from the background samples used to calculate the BKZ.\">");
					addInfo = false;
				}
				vcfOut.println(record);
				if (record.startsWith("#CHROM")){
					endFound = true;
					break;
				}
			}
			else {
				endFound = true;
				break;
				//numRecords++;
				//score(record);
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
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		int[] sampleIndexes = null;
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
					case 'i': sampleIndexes = Num.parseInts(args[++i], Misc.COMMA); break;
					case 'z': minimumZScore = Double.parseDouble(args[++i]); break;
					case 'b': bpBuffer = Integer.parseInt(args[++i]); break;
					case 'q': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'c': minReadCoverage = Integer.parseInt(args[++i]); break;
					case 'a': minNumSamples = Integer.parseInt(args[++i]); break;
					case 'e': removeNonZScoredRecords = true; break;
					case 'd': verbose = true; break;
					case 'r': replaceQualScore = true; break;
					case 't': numberThreads = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (sampleIndexes!= null){
			sampleIndexesToExamine = new HashSet<Integer>(sampleIndexes.length);
			for (Integer i: sampleIndexes) sampleIndexesToExamine.add(i);
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
		System.out.println("Settings:\nBackground\t"+mpileup);
		System.out.println("Save dir\t"+saveDir);
		System.out.println(bpBuffer+"\tBP buffer");
		System.out.println(minReadCoverage+"\tMin mpileup sample read coverage");
		System.out.println(minBaseQuality+"\tMin mpileup sample base quality");
		System.out.println(maxSampleAF+"\tMax mpileup sample AF");
		System.out.println(minNumSamples+"\tMin # samples for z-score calc");
		System.out.println(minimumZScore+"\tMin vcf AF z-score to save");
		System.out.println(numberThreads+"\tCPUs");
		System.out.println(removeNonZScoredRecords+ "\tExclude vcf records that could not be z-scored");
		System.out.println(verbose+"\tVerbose");
		System.out.println(replaceQualScore+"\tReplace QUAL score with z-score and set non scored records to 0");
		
		if (sampleIndexesToExamine!=null) System.out.println(Misc.hashSetToString(sampleIndexesToExamine, ",")+"\tMpileup samples to examine");
		else System.out.println("All\tMpileup samples to examine");
	}

	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         VCF Background Checker : March 2016                      **\n" +
				"**************************************************************************************\n" +
				"VBC calculates non-reference allele frequencies (AF) from a background multi-sample \n"+
				"mpileup file over each vcf record. It then calculates a z-score for the vcf AF and \n"+
				"appends it to the INFO field. If multiple bps are affected (e.g. INDELs) or bp padding\n"+
				"provided, the lowest bp z-score is appended. Z-scores < ~4 are indicative of non\n"+
				"reference bps in the background samples. Note, VBC requires an AF tag in the INFO\n"+
				"field of each record. \n"+

				"\nRequired:\n"+
				"-v Path to a xxx.vcf(.gz/.zip OK) file or directory containing such.\n" +
				"-m Path to a bgzip compressed and tabix indexed multi-sample mpileup file. e.g.:\n"+
				"      1) Mpileup: 'echo \"#SampleOrder: \"$(ls *bam) > bkg.mpileup; samtools mpileup\n"+
				"             -B -q 13 -d 1000000 -f $fastaIndex -l $bedFile *.bam >> bkg.mpileup'\n"+
				"      2) Bgzip: 'tabix-0.2.6/bgzip bkg.mpileup'\n"+
                "         Tabix: 'tabix-0.2.6/tabix -s 1 -b 2 -e 2 bkg.mpileup.gz'\n"+
				"-s Path to directory in which to save the modified vcf file(s)\n"+
						
				"\nOptional:\n" +
				"-i Comma delimited list (zero is 1st sample, no spaces) of mpileup sample indexes to\n"+
				"     examine, defaults to all.\n"+
				"-b BP padding, defaults to 0\n"+
				"-z Minimum vcf z-score, defaults to 0, no filtering\n"+
				"-q Minimum mpileup sample bp quality, defaults to 20\n"+
				"-c Minimum mpileup sample read coverge, defaults to 20\n"+
				"-f Maximum mpileup sample AF, defaults to 0.3\n"+
				"-a Minimum # mpileup samples for z-score calculation, defaults to 3\n" +
				"-e Exclude vcf records that could not be z-scored\n"+
				"-r Replace QUAL value with z-score. Un scored vars will be assigned 0\n"+
				"-d Print verbose debugging output\n" +
				"-t Number of threads to use, defaults to all\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/VCFBackgroundChecker -v SomaticVcfs/ -z 3\n"+
				"-m bkg.mpileup.gz -s BkgFiltVcfs/ -b 1 -q 13 -e -r \n\n"+

		        "**************************************************************************************\n");
	}

	//getters and setters
	public int getBpBuffer() {
		return bpBuffer;
	}

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

	public HashSet<Integer> getSampleIndexesToExamine() {
		return sampleIndexesToExamine;
	}



}
