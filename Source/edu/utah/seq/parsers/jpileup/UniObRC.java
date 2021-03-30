package edu.utah.seq.parsers.jpileup;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.*;
import htsjdk.samtools.*;
import htsjdk.samtools.cram.ref.ReferenceSource;
import util.apps.MergeRegions;
import util.bio.annotation.Bed;
import util.gen.*;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream; //needed by cram

/** Stripped down version of BamPileup to calculate unique observation read coverage
 * Set MinBaseQual to 1, MinMapQual to 0, includeOverlaps to true, to replicate IGV.
 * @author Nix
 * */
public class UniObRC {

	//user defined fields
	private File bedFile;
	private File bamFile;
	private File results;
	private File fastaFile;
	private File tabix;
	private File bgzip;
	private int minMappingQuality = 13;
	private int minBaseQuality = 13;
	private int minimumReadDepth = 10;
	private SamReaderFactory samFactory;
	private int numCpu = 25;
	private int maxBpOfRegion = 1000;
	private int maximumCoverageCalculated = 1001;
	private File tempDir = null;
	private ArrayList<File> tempBedFiles = new ArrayList<File>();
	private Histogram histogram = null;
	private ArrayList<Double> fractionTargetBpsAL = null;
	private int coverageAt95 = 0;
	private int coverageAt90 = -1;

	//constructor for stand alone use
	public UniObRC(String[] args){
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Num.formatNumber(diffTime, 2)+" minutes");
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR running UniObRC! Correct and restart.");
		}
	}
	
	public void doWork() throws Exception {
		//create reader
		samFactory = SamReaderFactory.makeDefault().referenceSource(new ReferenceSource(fastaFile)).validationStringency(ValidationStringency.SILENT);
		
		//load and chunk the bed file of regions to scan
		int numRegionsPerWorker = 0;
		Bed[] regions = Bed.parseFile(bedFile, 0, 0);
		Arrays.sort(regions);
		int numOriRegions = regions.length;
		
		//split big regions
		regions = Bed.splitBigRegions(regions, maxBpOfRegion);
		IO.pl("Splitting "+numOriRegions+" regions into "+regions.length+" with a max size of "+maxBpOfRegion+" bps...");
		
		if (regions.length <= numCpu) {
			numCpu = regions.length;
			numRegionsPerWorker = 1;
		}
		else numRegionsPerWorker = (int)Math.round( (double)regions.length/(double)numCpu );
		Object[][] chunks = Misc.chunk(regions, numRegionsPerWorker);
		
		//make and execute loaders
		IO.p("Launching "+chunks.length+" loaders...");
		ExecutorService executor = Executors.newFixedThreadPool(chunks.length);
		RCBamLoader[] loaders = new RCBamLoader[chunks.length];
		for (int i=0; i< chunks.length; i++) {
			Bed[] b = new Bed[chunks[i].length];
			for (int x = 0; x< b.length; x++) b[x] = (Bed)chunks[i][x];
			loaders[i] = new RCBamLoader(this, i, b);
			executor.execute(loaders[i]);
		}
		executor.shutdown();
		while (!executor.isTerminated()) {}
		IO.pl();
		
		//check loaders fetch files
		histogram = new Histogram(0, maximumCoverageCalculated, maximumCoverageCalculated);
		for (RCBamLoader l: loaders) {
			if (l.isFailed()) throw new IOException("ERROR: File Loader issue! \n");
			tempBedFiles.add(l.getCoverageFile());
			histogram.addCounts(l.getHistogram());
		}
		
		printCoverageStats(histogram);
		procTempFiles();
		printJson();

	}
	
	private void printJson() throws Exception {

		String jPath = Misc.removeExtension(results.getName()) + ".json.gz";
		Gzipper gz = new Gzipper(new File (results.getParentFile(), jPath));
		gz.println("{");
		gz.printJson("meanCoverage", Num.formatNumber(histogram.getStandardDeviation().getMean(), 1), true);
		gz.printJson("coverageAt0.95OfTargetBps", coverageAt95, true);
		gz.printJson("coverageAt0.90OfTargetBps", coverageAt90, true);
		gz.printJson("minimumMappingQuality", minMappingQuality, true);
		gz.printJson("minimumBaseQuality", minBaseQuality, true);
		gz.printJson("minimumPassingBedCoverageThreshold", minimumReadDepth, true);
		gz.printJson("passingBedFilePath", results.getCanonicalPath(), true);
		gz.printJson("targetRegionFilePath", bedFile.getCanonicalPath(), true);
		gz.printJson("alignmentFilePath", bamFile.getCanonicalPath(), true);
		gz.printJson("fractionTargetBpsWithIndexedCoverage", Num.arrayListToDoubles(fractionTargetBpsAL), true);
		gz.println("}");
		gz.close();
	}

	private void procTempFiles() throws IOException {
		//strip off any gz
		String name = results.getCanonicalPath();
		if (name.endsWith(".gz")) name = name.substring(0, name.length()-3);
		
		//merge the temp files
		IO.pl("\nMerging uniOb read coverge beds...");
		File mergedBed = new File(name);
		File[] beds = new File[tempBedFiles.size()];
		tempBedFiles.toArray(beds);
		new MergeRegions (beds, mergedBed, false, false);
		
		//add header
		addHeader(mergedBed);
			
		//tabix?
		if (bgzip != null) BamPileupMerger.compressAndIndex(bgzip, tabix, mergedBed, new String[] {"-f", "-s", "1", "-b", "2", "-e", "3"}, false);
		
	}
	
	public void addHeader(File source) throws IOException{
		BufferedReader in = IO.fetchBufferedReader(source);
		File temp = new File(tempDir, Misc.getRandomString(10)+"_merged.bed");
		PrintWriter out = new PrintWriter (new FileWriter (temp));
		
		//add header
		out.println("# BamCram "+ IO.getCanonicalPath(bamFile));
		out.println("# Bed "+IO.getCanonicalPath(bedFile));
		out.println("# MinUniObRC\t"+minimumReadDepth);
		out.println("# MinMapQual\t"+minMappingQuality);
		out.println("# MinBaseQual\t"+minBaseQuality);
		out.println("# MeanCoverage\t" +Num.formatNumber(histogram.getStandardDeviation().getMean(), 1));
		out.println("# CoverageAt0.95OfTargetBps\t"+ coverageAt95);
		out.println("# CoverageAt0.90OfTargetBps\t"+ coverageAt90);

		
		String line;
		while ((line = in.readLine()) != null) out.println(line);
		in.close();
		out.close();
		
		//rename
		temp.renameTo(source);
				
	}

	
	public void printCoverageStats(Histogram histogram){

			long[] counts = histogram.getBinCounts();
			double total = histogram.getTotalBinCounts();
			fractionTargetBpsAL = new ArrayList<Double>();

			IO.pl("\nUniObAlignmentCoverage\tFractionTargetBPs");
			double numCounts = 0;
			String zero = "0";
			for (int i=0; i< counts.length; i++){
				double cumFrac = (total-numCounts)/ total;
				//runout up to just two decimals for calculating 0.90 and 0.95
				String formNum = Num.formatNumber(cumFrac, 2);
				if (formNum.equals(zero)) {
					fractionTargetBpsAL.add(0.0);
					break;
				}
				IO.pl(i+"\t"+formNum);
				fractionTargetBpsAL.add(Double.parseDouble(formNum));
				numCounts += counts[i];
				if (numCounts == total) break;
			}
			
			//calc 0.95 and 0.9
			for (int i=fractionTargetBpsAL.size()-1; i >=0; i--){
				double cov = fractionTargetBpsAL.get(i).doubleValue();
				if (cov >= 0.9 && coverageAt90 == -1) coverageAt90 = i;
				if (cov >= 0.95) {
					coverageAt95 = i;
					break;
				}
			}
			if (coverageAt90 == -1) coverageAt90 = 0;
			
			//print summary stats
			IO.pl("\nTotal interrogated bases\t"+ (int)total);
			IO.pl("Mean Coverage\t"+Num.formatNumber(histogram.getStandardDeviation().getMean(), 1));
			IO.pl("CoverageAt0.95OfTargetBps\t" + coverageAt95);
			IO.pl("CoverageAt0.90OfTargetBps\t" + coverageAt90);
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new UniObRC(args);
	}		
	
	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args) throws Exception{
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File tabixBinDirectory = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bamFile = new File(args[++i]).getCanonicalFile(); break;
					case 'r': bedFile = new File(args[++i]); break;
					case 's': results = new File(args[++i]).getCanonicalFile(); break;
					case 'f': fastaFile = new File(args[++i]); break;
					case 'q': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'm': minMappingQuality = Integer.parseInt(args[++i]); break;
					case 'x': maxBpOfRegion = Integer.parseInt(args[++i]); break;
					case 'd': minimumReadDepth = Integer.parseInt(args[++i]); break;
					case 'p': numCpu = Integer.parseInt(args[++i]); break;
					case 'c': maximumCoverageCalculated  = (Integer.parseInt(args[++i]) + 1); break;
					case 't': tabixBinDirectory = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull bam and cram files
		if (bamFile == null || bamFile.canRead() == false) Misc.printExit("\nError: cannot find your xxx.bam or xxx.cram file!\n");		
		
		//Create fasta fetcher
		if (bamFile.getName().endsWith(".cram") && (fastaFile == null || fastaFile.canRead() == false))  Misc.printErrAndExit("\nError: please provide an reference genome fasta file and it's index for cram conversion");	
		if (bedFile == null ||  bedFile.canRead() == false) Misc.printErrAndExit("\nError: please provide a file of regions in bed format.");
		if (results == null ) Misc.printErrAndExit("\nError: please provide a results file that ends with xxx.bed.gz");
		
		//pull tabix and bgzip
		if (tabixBinDirectory != null) {
			bgzip = new File (tabixBinDirectory, "bgzip");
			tabix = new File (tabixBinDirectory, "tabix");
			if (bgzip.canExecute() == false || tabix.canExecute() == false) Misc.printExit("\nCannot find or execute bgzip or tabix executables from "+tabixBinDirectory);
		}
				
		//number of workers
		int numProc = Runtime.getRuntime().availableProcessors();
		if (numCpu == 0 || numCpu > numProc) numCpu = numProc;
		
		tempDir = results.getParentFile();
		tempDir.mkdirs();
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Unique Observation Read Coverage:  March 2021             **\n" +
				"**************************************************************************************\n" +
				"UniObRC calculates read coverage statistics and generates a bed file of regions\n"+
				"that pass thresholds for read depth. Overlapping pairs, secondary, supplementary,\n"+
				"duplicate, and failing vendorQC alignments are excluded. INDEL counts are included.\n"+
				"Summary statistics and app setting are saved in json format.\n"+

				
				"\nRequired Options:\n"+
				"-b Path to a coordinate sorted bam/cram file with index.\n"+
				"-r Bed file of non overlapping regions to interrogate. Run the USeq MergeRegions app\n"+
				"      if unsure. xxx.bed.gz/.zip OK\n"+
				"-s File path to write the passing region bed file.\n"+

				"\nDefault Options:\n"+
				"-q Minimum base quality, defaults to 13\n"+
				"-m Minimum alignment mapping quality, defaults to 13\n"+
				"-d Minimum read depth for passing coverge bed, defaults to 10\n"+
				"-i Include counts from overlapping paired reads, defaults to excluding them.\n"+
				"-p Number processors to use, defaults to 25, reduce if out of memory errors occur.\n"+
				"-x Max length of region chunk, defaults to 1000, set smaller if out of memory errors\n"+
				"      occur.\n"+
				"-c Maximum read coverage stats calculated, defaults to 1000.\n"+
				"-f For xxx.cram alignments, provide a path to the reference fasta with index.\n"+
				"-t Path to a directory containing the bgzip and tabix executables to compress and\n"+
				"      index the passing bed file, see htslib.org \n"+

				"\nExample: java -Xmx100G -jar pathTo/USeq/Apps/UniObRC -b my.cram -r target.bed\n"+
				"-f Ref/human_g1k_v37_decoy.fasta -s passing.bed.gz -t ~/BioApps/HTSlib/bin/ -d 20\n\n" +

				"**************************************************************************************\n");
	}

	public File getBedFile() {
		return bedFile;
	}
	public File getBamFile() {
		return bamFile;
	}
	public File getTempDir() {
		return tempDir;
	}
	public File getFastaFile() {
		return fastaFile;
	}
	public int getMinMappingQuality() {
		return minMappingQuality;
	}
	public int getMinBaseQuality() {
		return minBaseQuality;
	}
	public SamReaderFactory getSamFactory() {
		return samFactory;
	}
	public double getMaximumCoverageCalculated() {
		return maximumCoverageCalculated;
	}
	public int getMinimumReadDepth() {
		return minimumReadDepth;
	}
}
