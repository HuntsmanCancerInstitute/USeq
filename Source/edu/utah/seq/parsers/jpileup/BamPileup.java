package edu.utah.seq.parsers.jpileup;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.*;
import htsjdk.samtools.*;
import htsjdk.samtools.cram.ref.ReferenceSource;
import util.bio.annotation.Bed;
import util.gen.*;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream; //needed by cram

/** Set MinBaseQual to 1, MinMapQual to 0, includeOverlaps to true, to replicate IGV.
 * /scratch/mammoth/serial/u0028003/Underhill/SpikeBuild/SpikedAlignments/FirstSet/15352X1_9/15352X1_BB_Indels/BamMixer/0.01.bam
#chr1	2479985	2479988	DelSmall	0	.
#chr17	31214456	31214475	DelLarge	0	.
#chr17	31218969	31218970	InsSmall-39T,35I	0	.
#chr7	140753338	140753339	InsLarge-5301G,1N,50I	0	.
#chr11	102952198	102952201	MixINDEL	0	.
#chr13	48380133	48380135	Ns	0	.
#chr13	48379517	48379522	Mess	0	.

#Chr	1BasePos	Ref	A,C,G,T,N,Del,Ins,FailBQ
chr1	2479986	G	1,0,6914,3,2,6,0,0
chr17	31214457	A	4191,1,0,1,0,0,0,0
chr17	31214458	C	3,4179,0,0,0,32,0,0
chr17	31214459	A	4180,0,0,0,4,32,0,0
chr17	31214460	T	0,2,0,4173,2,32,0,0
chr17	31214461	T	0,0,0,4162,0,32,0,0
chr17	31214462	T	1,1,2,4111,1,32,0,0
chr17	31214463	A	4085,1,0,0,0,32,0,0
chr17	31214464	A	4111,0,0,1,1,32,0,0
chr17	31214465	A	4127,0,0,0,0,32,0,0
chr17	31214466	G	1,0,4138,1,1,32,0,0
chr17	31214467	A	4163,0,0,0,0,38,0,0
chr17	31214468	A	4219,0,0,0,1,32,0,0
chr17	31214469	A	4258,0,0,0,1,32,0,0
chr17	31214470	A	4283,0,0,0,0,32,0,0
chr17	31214471	A	4299,0,0,0,1,32,0,0
chr17	31214472	G	0,0,4257,0,5,32,0,0
chr17	31214473	T	1,0,0,4242,3,32,0,0
chr17	31214474	A	4260,1,0,0,0,2,0,0
chr17	31214475	A	4262,0,0,0,0,0,0,0
chr17	31218970	T	0,0,0,39,0,0,35,0
chr7	140753339	G	0,0,5201,0,1,0,50,0
chr11	102952199	G	3,1,6224,1,0,0,0,0
chr11	102952200	A	6136,0,0,0,1,69,18,0
chr11	102952201	A	6137,0,1,1,2,1,0,0
chr13	48380134	G	25,0,1813,6,71,0,0,0
chr13	48380135	T	2,0,1,1914,18,0,0,0
chr13	48379518	T	4,3,0,1564,53,0,1,0
chr13	48379519	C	156,1417,0,7,137,4,76,0
chr13	48379520	A	2206,19,0,1,4,141,970,0 - diff T from IGV
chr13	48379521	A	2373,8,0,0,0,47,0,0
chr13	48379522	A	2480,0,0,1,1,13,1,0


 * @author Nix
 * */
public class BamPileup {

	//user defined fields
	private File bedFile;
	private File bamFiles[];
	private File tempDir;
	private File results;
	private File fastaFile;
	private int minMappingQuality = 13;
	private int minBaseQuality = 10;
	private SamReaderFactory samFactory;
	private boolean includeOverlaps = false;
	private int numCpu = 0;

	//constructor
	public BamPileup(String[] args){
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			//create reader
			samFactory = SamReaderFactory.makeDefault().referenceSource(new ReferenceSource(fastaFile)).validationStringency(ValidationStringency.SILENT);
			
			//load and chunk the bed file of regions to scan
			int numRegionsPerWorker = 0;
			Bed[] regions = Bed.parseFile(bedFile, 0, 0);
			Arrays.sort(regions);
			if (regions.length <= numCpu) {
				numCpu = regions.length;
				numRegionsPerWorker = 1;
			}
			else numRegionsPerWorker = (int)Math.round( (double)regions.length/(double)numCpu );
			Object[][] chunks = Misc.chunk(regions, numRegionsPerWorker);
			
			//make and execute loaders
			IO.p("Launching "+chunks.length+" loaders");
			ExecutorService executor = Executors.newFixedThreadPool(chunks.length);
			BamPileupLoader[] loaders = new BamPileupLoader[chunks.length];
			for (int i=0; i< chunks.length; i++) {
				Bed[] b = new Bed[chunks[i].length];
				for (int x = 0; x< b.length; x++) b[x] = (Bed)chunks[i][x];
				loaders[i] = new BamPileupLoader(this, i, b);
				executor.execute(loaders[i]);
			}
			executor.shutdown();
			while (!executor.isTerminated()) {}
			IO.pl();
			
			//check loaders fetch files
			ArrayList<File> tempFiles = new ArrayList<File>();
			for (BamPileupLoader l: loaders) {
				if (l.isFailed()) throw new IOException("ERROR: File Loader issue! \n");
				tempFiles.add(l.getResultsFile());
			}
			
			//concatinate gzipped files
			IO.concatinateFiles(tempFiles, results);

			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Num.formatNumber(diffTime, 2)+" minutes");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BamPileup(args);
	}		
	
	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args) throws Exception{
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': forExtraction = new File(args[++i]); break;
					case 'r': bedFile = new File(args[++i]); break;
					case 'f': fastaFile = new File(args[++i]); break;
					case 's': results = new File(args[++i]); break;
					case 'q': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'm': minMappingQuality = Integer.parseInt(args[++i]); break;
					case 'p': numCpu = Integer.parseInt(args[++i]); break;
					case 'i': includeOverlaps = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//Create fasta fetcher
		if (fastaFile == null || fastaFile.canRead() == false)  Misc.printErrAndExit("\nError: please provide an reference genome fasta file and it's index.");		
		if (bedFile == null ||  bedFile.canRead() == false) Misc.printErrAndExit("\nError: please provide a file of regions in bed format.");
		if (results == null ) Misc.printErrAndExit("\nError: please provide a results file that ends with xxx.gz");
		
		//pull bam and cram files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to directory containing xxx.bam and xxx.cram files with their indexes.\n");
		File[][] tot = new File[2][];
		tot[0] = IO.extractFiles(forExtraction, ".bam");
		tot[1] = IO.extractFiles(forExtraction,".cram");
		bamFiles = IO.collapseFileArray(tot);
		if (bamFiles == null || bamFiles.length ==0 || bamFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.bam or xxx.cram files!\n");		
		
		//number of workers
		int numProc = Runtime.getRuntime().availableProcessors();
		if (numCpu == 0 || numCpu > numProc) numCpu = numProc;
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Bam Pileup:  March 2020                           **\n" +
				"**************************************************************************************\n" +
				"BP extracts pileup information for each bam file over a list of regions. This includes\n"+
				"the # A,C,G,T,N,Del,Ins,FailingBQ bps for each bam. Provide the max memory available to\n"+
				"JVM, multithreaded. BGZip and Tabix index the results to make it searchable, e.g.\n"+
				"      gunzip 15352R.bp.txt.gz\n"+
				"      ~/BioApps/HTSlib/1.3/bin/bgzip --threads 20 15352R.bp.txt\n"+
				"      ~/BioApps/HTSlib/1.3/bin/tabix -s 1 -b 2 -e 2 15352R.bp.txt.gz\n"+
				"The header of the results file contains the order of the bams. To approximate IGV\n"+
				"pileup info, set -q 0 -m 0 -i\n"+
				
				"\nRequired Options:\n"+
				"-b Path to a directory of coordinate sorted bams or crams with indexes.\n"+
				"-r Bed file of regions to extract pileup information. MUST BE NON OVERLAPPING. Run the\n"+
				"      USeq MergeRegions app if unsure. xxx.bed.gz/.zip OK\n"+
				"-f Path to the reference fasta with and xxx.fai index.\n"+
				"-s Path to a gzip file to save the pileup information, must end in xxx.gz\n"+

				"\nDefault Options:\n"+
				"-q Minimum base quality, defaults to 10\n"+
				"-m Minimum alignment mapping quality, defaults to 13\n"+
				"-i Include counts from overlapping paired reads, defaults to excluding them.\n"+
				"-p Number processors to use, defaults to all.\n"+

				"\nExample: java -Xmx100G -jar pathTo/USeq/Apps/BamPileup -b CramFiles/ -r target.bed\n"+
				"-f Ref/human_g1k_v37_decoy.fasta -s 15352RX.bp.txt.gz \n\n" +

				"**************************************************************************************\n");
	}

	public File getBedFile() {
		return bedFile;
	}
	public File[] getBamFiles() {
		return bamFiles;
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
	public boolean isIncludeOverlaps() {
		return includeOverlaps;
	}

}
