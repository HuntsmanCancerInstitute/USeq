package edu.utah.seq.vcf.sim;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.parsers.SamAlignmentExtractor;
import edu.utah.seq.vcf.xml.foundation.SimpleVcf;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Generates mixes of alignments for tumor normal simulations with BamBlaster
 * @author Nix
 * */
public class ExactBamMixer {

	//fields
	private File injectedBamFile;
	private File unModifiedMatchingBamFile;
	private File unModifiedNoMatchBamFile;
	private File saveDirectory;
	private File vcfFile;
	private String readGroupName = "EBM";
	private double[] fractions = {0.025, 0.05, 0.1, 0.2};
	private int minNumAltReads = 2;

	//internal fields
	private SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	private BufferedReader vcfIn;
	private int numberThreads = 0;
	private static final int numVcfToProc = 100;
	private BamMixerLoader[] loaders = null;
	private String vcfHeader = null;
	private File unModMatchingNameSortedBam = null;
	private File[] modNameSortedBams = null;


	public ExactBamMixer (String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);

		doWork();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}

	public void doWork(){

		try {
			createVcfReader();

			parseAlignments();

			concatinateTxtResults();
			
			concatinateVcfs();

			concatinateNameSortSams();
			
			nameSortUnModified();
			
			walkWriteFinalBams();

			vcfIn.close();
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void walkWriteFinalBams() throws IOException {

		WalkSortSam[] wss = new WalkSortSam[fractions.length];

		//for each fraction, make a walker
		for (int i=0; i< fractions.length; i++){
			File finalBam = new File (saveDirectory, fractions[i]+".bam");
			wss[i] = new WalkSortSam(modNameSortedBams[i], unModMatchingNameSortedBam, finalBam, unModifiedNoMatchBamFile, i, readGroupName+"_"+fractions[i]);
		}

		//run threads
		int num = fractions.length;
		if (numberThreads < num) num = numberThreads;
		IO.pl("Walking unmodified adding alignments for each mix with "+num+" threads...");
		ExecutorService executor = Executors.newFixedThreadPool(num);
		for (WalkSortSam l: wss) executor.execute(l);
		executor.shutdown();
		while (!executor.isTerminated()) {}
		for (WalkSortSam l: wss) {
			if (l.isFailed()) throw new IOException("ERROR: Walk Sort Sam issue! \n");
		}
		
	}

	private void nameSortUnModified() {
		IO.pl("Name sorting unmodified matching bam...");
		unModMatchingNameSortedBam = new File (saveDirectory, Misc.removeExtension(unModifiedMatchingBamFile.getName())+".nameSorted.unmod.bam");
		unModMatchingNameSortedBam.deleteOnExit();
		new PicardSortSam (unModifiedMatchingBamFile, unModMatchingNameSortedBam, false);
	}

	private void concatinateTxtResults() throws IOException {
		IO.pl("Saving VCF AF stats...");
		ArrayList<File> toMerge = new ArrayList<File>();
		//for each loader pull appropriate gzipped file
		for (int j=0; j< loaders.length; j++) toMerge.add(loaders[j].getTargetVarResults().getGzipFile());
		File f = new File (saveDirectory, Misc.removeExtension(vcfFile.getName())+".txt.gz");
		IO.concatinateFiles(toMerge, f);
	}
	private void concatinateVcfs() throws IOException {
		IO.pl("Saving VCFs...");
		for (int i=0; i< fractions.length; i++){
			ArrayList<File> toMerge = new ArrayList<File>();
			//for each loader pull appropriate gzipped file
			for (int j=0; j< loaders.length; j++) toMerge.add(loaders[j].getVcfWriters()[i].getGzipFile());
			File unsortedVcf = new File (saveDirectory, "unsorted_"+fractions[i]+".vcf.gz");
			unsortedVcf.deleteOnExit();
			IO.concatinateFiles(toMerge, unsortedVcf);
			
			//sort it
			File vcf = new File (saveDirectory, fractions[i]+".vcf.gz");
			SimpleVcf.sortVcf(unsortedVcf, vcf);
		}
	}

	private void parseAlignments() throws IOException {
		IO.pl("Loading and mixing alignments with "+numberThreads+" threads...");

		//make Loaders
		loaders = new BamMixerLoader[numberThreads];
		for (int i=0; i< loaders.length; i++) loaders[i] = new BamMixerLoader(injectedBamFile, unModifiedMatchingBamFile, this, i);

		//parse alignments
		ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
		for (BamMixerLoader l: loaders) executor.execute(l);
		executor.shutdown();
		while (!executor.isTerminated()) {}
		for (BamMixerLoader l: loaders) {
			if (l.isFailed()) throw new IOException("ERROR: File Loader issue! \n");
		}
	}

	private void concatinateNameSortSams() throws IOException {

		CatSortSam[] css = new CatSortSam[fractions.length];
		modNameSortedBams = new File[fractions.length];

		//for each fraction
		for (int i=0; i< fractions.length; i++){
			ArrayList<File> toMerge = new ArrayList<File>();
			//for each loader pull appropriate gzipped file
			for (int j=0; j< loaders.length; j++) toMerge.add(loaders[j].getSamWriters()[i].getGzipFile());
			modNameSortedBams[i] = new File (saveDirectory, fractions[i]+".nameSorted.mod.bam");
			modNameSortedBams[i].deleteOnExit();
			css[i] = new CatSortSam(toMerge, modNameSortedBams[i], i);
		}

		//run threads
		int num = fractions.length;
		if (numberThreads < num) num = numberThreads;
		IO.pl("Combining and sorting sams with "+num+" threads...");
		ExecutorService executor = Executors.newFixedThreadPool(num);
		for (CatSortSam l: css) executor.execute(l);
		executor.shutdown();
		while (!executor.isTerminated()) {}
		for (CatSortSam l: css) {
			if (l.isFailed()) throw new IOException("ERROR: Cat Sort Sam issue! \n");
		}
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
		return true;
	}

	private void createVcfReader() throws Exception {
		vcfIn = IO.fetchBufferedReader(vcfFile);
		ArrayList<String> vcfHeaderAl = new ArrayList<String>();
		String record;
		boolean endFound = false;
		//advance to just before the first vcf record
		while ((record = vcfIn.readLine()) != null){
			record = record.trim();
			if (record.length() == 0) continue;
			if (record.startsWith("#")){
				vcfHeaderAl.add(record);
				if (record.startsWith("#CHROM")){
					endFound = true;
					break;
				}
			}
		}
		if (endFound == false) Misc.printErrAndExit("\nERROR: failed to find the #CHROM line in the vcf file, aborting\n");
		vcfHeader = Misc.stringArrayListToString(vcfHeaderAl, "\n");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ExactBamMixer(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		String frac = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'i': injectedBamFile = new File(args[++i]); break;
					case 'u': unModifiedMatchingBamFile = new File(args[++i]); break;
					case 'f': unModifiedNoMatchBamFile = new File(args[++i]); break;
					case 'r': saveDirectory = new File(args[++i]); break;
					case 'v': vcfFile = new File(args[++i]); break;
					case 't': numberThreads = Integer.parseInt(args[++i]); break;
					case 'm': frac = args[++i]; break;
					case 'n': readGroupName = args[++i]; break;
					case 'a': minNumAltReads = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//diff fractions than defaults?
		if (frac != null) fractions = Num.parseDoubles(Misc.COMMA.split(frac));

		//check bams
		if (injectedBamFile == null || injectedBamFile.canRead() == false) Misc.printErrAndExit("Error: please proved a sorted bam alignment file containing the merged paired and single end alignments from aligning the fastq.gz files from your BamBlaster run.\n");
		if (unModifiedMatchingBamFile == null || unModifiedMatchingBamFile.canRead() == false) Misc.printErrAndExit("Error: please proved a path to the xxx_unmodified.bam from your BamBlaster run.\n");
		if (unModifiedNoMatchBamFile == null || unModifiedNoMatchBamFile.canRead() == false) Misc.printErrAndExit("Error: please proved a path to the xxx_filtered.bam from your BamBlaster run.\n");
		SamAlignmentExtractor.lookForBaiIndexes(new File[]{injectedBamFile,unModifiedMatchingBamFile,unModifiedMatchingBamFile,unModifiedNoMatchBamFile}, true);

		//check save directory
		if (saveDirectory == null) Misc.printErrAndExit("Error: please provide a path to a directory for saving the results.");
		saveDirectory.mkdirs();

		//check vcf
		if (vcfFile == null || vcfFile.exists() == false) Misc.printErrAndExit("Error: please provide a path to the vcf file containing variants to use in extracting alignments to mix.");

		//threads to use
		int numAvail = Runtime.getRuntime().availableProcessors();
		if (numberThreads < 1 || numberThreads > numAvail) numberThreads =  numAvail - 1;
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Exact Bam Mixer : May 2020                          **\n" +
				"**************************************************************************************\n" +
				"Combines bam alignment files in different fractions to simulate multiple variant\n"+
				"frequencies. Run BamBlaster first. Threaded, so provide almost all the memory available\n"+
				"to java. The ExactBamMixer attempts to create bam files containing variants will very\n"+
				"similar AFs.  The BamMixer produces more of a spread of AFs. Ignore Picard and log4j\n"+
				"NOTE and ERROR output. These are just warnings from Picard that can't be silenced.\n\n"+

				"Required:\n"+
				"-r Path to a directory to save the results\n" +
				"-u Path to the xxx_unmodified.bam from your BamBlaster run\n"+
				"-f Path to the xxx_filtered.bam from your BamBlaster run\n"+
				"-i Path to your realigned bam containing injected variants, merge the single and \n"+
				"     paired end alignment files with MergeSams USeq app.\n"+
				"-v Path to the vcf file containing variants used to modify the BamBlaster alignments\n"+

				"\nOptional:\n"+
				"-m Fractions to mix in the variant alignments, comma delimited, no spaces, defaults to\n"+
				"     0.025,0.05,0.1,0.2\n"+
				"-t Number of threads to use, defaults to all\n"+
				"-a Minimum number alt read pairs to include an injected variant in a particular mixed\n"+
				"     bam, defaults to 2\n"+
				"-n Name to prepend onto the read group, defaults to EBM\n"+

				"\nExample: java -Xmx100G -jar pathTo/USeq/Apps/ExactBamMixer -r ~/TumorSim/ -v snv.vcf\n"+
				"    -u bb_unmodified.bam -f bb_filtered.bam -i bb_mergedReAlign.bam -n Snv\n\n" +

				"**************************************************************************************\n");
	}

	public SamReaderFactory getReaderFactory() {
		return readerFactory;
	}

	public double[] getFractions() {
		return fractions;
	}

	public File getSaveDirectory() {
		return saveDirectory;
	}

	public File getUnModifiedNoMatchBamFile() {
		return unModifiedNoMatchBamFile;
	}

	public int getMinNumAltReads() {
		return minNumAltReads;
	}

	public String getVcfHeader() {
		return vcfHeader;
	}
}
