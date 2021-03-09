package edu.utah.seq.vcf.splice;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.*;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
import edu.utah.seq.mes.*;
import edu.utah.seq.parsers.BamLoader;
import edu.utah.seq.vcf.VCFLookUp;
import edu.utah.seq.vcf.VCFParser;
import edu.utah.seq.vcf.VCFRecord;
import edu.utah.seq.vcf.xml.foundation.SimpleVcf;
import util.bio.annotation.ExonIntron;
import util.bio.annotation.ExportIntergenicRegions;
import util.bio.parsers.UCSCGeneLine;
import util.bio.parsers.UCSCGeneModelTableReader;
import util.bio.seq.Seq;
import util.gen.*;

/**Scans known splices
 * @author Nix
 * */
public class KnownSpliceJunctionScanner {

	//user fields
	private File spliceModelDirectory;
	private IndexedFastaSequenceFile fasta; 
	private File transcriptSeqFile;
	private File tempDirectory;
	private int numberThreads = 0;
	private File bedResultsFile;


	//internal fields
	private HashMap<String,UCSCGeneLine[]> chromTranscripts;
	private KnownSpliceLoader[] loaders;
	private ArrayList<File> gzippedToConcat = new ArrayList<File>();
	private Histogram spliceHistogram5 = new Histogram(-14, 14, 100);
	private Histogram spliceHistogram3 = new Histogram(-14, 14, 100);
	private long numberUniqueSpliceJunctionsScored = 0;

	//constructor
	public KnownSpliceJunctionScanner(String[] args){
		try {
			//start clock
			long startTime = System.currentTimeMillis();

			//process args
			processArgs(args);

			//load transcripts
			System.out.println("Loading transcripts...");
			loadTranscripts();
			annotateKnownSplices();
			
			System.out.println("Concatinating beds...");
			IO.concatinateFiles(gzippedToConcat, bedResultsFile);
			IO.deleteDirectory(tempDirectory);
			
			//print stats
			System.out.println("\nScore stats:\n"+numberUniqueSpliceJunctionsScored+"\t# Unique splice junctions scored");
			System.out.println(Num.formatNumber(spliceHistogram5.getStandardDeviation().getMean(), 3)+"\tMean 5' splice junctions");
			System.out.println(Num.formatNumber(spliceHistogram3.getStandardDeviation().getMean(), 3)+"\tMean 3' splice junctions");
			System.out.println(Num.formatNumber(spliceHistogram5.getStandardDeviation().getStandardDeviation(), 3)+"\tStndDev 5' splice junctions");
			System.out.println(Num.formatNumber(spliceHistogram3.getStandardDeviation().getStandardDeviation(), 3)+"\tStndDev 3' splice junctions");
			System.out.println("\nHistogram 5' splices:");
			spliceHistogram5.printScaledHistogram();
			System.out.println("\nHistogram 3' splices:");
			spliceHistogram3.printScaledHistogram();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
			System.out.println("\nDone "+Math.round(diffTime)+" min!\n");
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public void annotateKnownSplices() {

		//create some loaders
		loaders = new KnownSpliceLoader[chromTranscripts.size()];
		int i = 0; 
		for (String chromName: chromTranscripts.keySet()) {
			loaders[i++] = new KnownSpliceLoader(this, chromName, new File(tempDirectory, i+"_.bed.gz"));
		}

		try {
			//launch em!
			ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
			for (KnownSpliceLoader l: loaders) executor.execute(l);
			executor.shutdown();
			while (!executor.isTerminated()) {}  //wait here until complete

			//check loaders 
			for (KnownSpliceLoader l: loaders) {
				if (l.isFailed()) throw new IOException("ERROR: Failed to annotate splices with one of the loaders?");
				int num = l.getNumUniqueJunctions();
				if (num !=0){
					gzippedToConcat.add(l.getBedOut().getGzipFile());
					spliceHistogram3.addCounts(l.getSpliceHistogram3());
					spliceHistogram5.addCounts(l.getSpliceHistogram5());
					numberUniqueSpliceJunctionsScored += num;
				}
			}
		}catch (Exception e) {
			e.printStackTrace();
			IO.deleteDirectory(tempDirectory);
			System.exit(1);
		}
	}



	/*TODO: replace with IntervalST*/
	public void loadTranscripts(){
		//main transcripts
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(transcriptSeqFile, 0);
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your transcript's coordinates are reversed. Check that each start is less than the stop.\n");
		chromTranscripts = reader.getChromSpecificGeneLines();
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new KnownSpliceJunctionScanner(args);
	}		


	/**This method will process each argument and assign new variables
	 * @throws FileNotFoundException */
	public void processArgs(String[] args) throws FileNotFoundException{
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File indexedFasta = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': indexedFasta = new File(args[++i]); break;
					case 'm': spliceModelDirectory = new File(args[++i]); break;
					case 'u': transcriptSeqFile = new File(args[++i]); break;
					case 'r': bedResultsFile = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//checkfiles
		//Create fasta fetcher
		if (indexedFasta == null || indexedFasta.exists() == false)  Misc.printErrAndExit("\nError: cannot find indexed fasta file? -> "+ indexedFasta);
		fasta = new IndexedFastaSequenceFile(indexedFasta);
		if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n"+ indexedFasta);

		if (spliceModelDirectory == null || spliceModelDirectory.isDirectory() == false ) {
			Misc.printErrAndExit("\nError: please provide a path to a directory containing the splice models.\n");
		}
		if (transcriptSeqFile == null || transcriptSeqFile.canRead()== false){
			Misc.printErrAndExit("\nPlease enter a transcript table file in ucsc refflat format.\n");
		}

		//temp directory
		if (bedResultsFile == null) Misc.printErrAndExit("\nPlease provide the path and name of a gzipped vcf to use in writing the results.\n");
		if (bedResultsFile.getName().endsWith(".gz") == false) Misc.printErrAndExit("\nSorry your output vcf file must end in xxx.gz\n");
		File parentDir = bedResultsFile.getParentFile();
		if (parentDir == null) parentDir = new File (System.getProperty("user.dir"));
		tempDirectory = new File (parentDir, "TempKnownSpliceJunctionScanner_DeleteMe_"+Misc.getRandomString(6));
		if (tempDirectory.mkdirs() == false) Misc.printErrAndExit("\nFailed to make the temp directory? "+tempDirectory);;

		//threads to use
		double totalGbAvailable = (double)(Runtime.getRuntime().maxMemory()/1000000000.0);
		int numPossCores = (int)Math.round(totalGbAvailable/7.0);
		if (numPossCores < 1) numPossCores = 1;
		int numPossThreads = Runtime.getRuntime().availableProcessors();

		if (numPossCores <= numPossThreads) numberThreads = numPossCores;
		else numberThreads = numPossThreads;


		System.out.println("Core usage:\n\tTotal GB available to Java:\t"+ Num.formatNumber(totalGbAvailable, 1));
		System.out.println("\tTotal available cores:\t"+numPossThreads);
		System.out.println("\tNumber cores to use @ 7GB/core:\t"+numberThreads+"\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                      Known Splice Junction Scanner : Sept 2017                   **\n" +
				"**************************************************************************************\n" +
				"Scores know splice junctions using the MaxEntScan algorithms. See Yeo and\n"+
				"Burge 2004, http://www.ncbi.nlm.nih.gov/pubmed/15285897 for details. \n\n" +

				"Required Options:\n"+
				"-r Name of a gzipped bed file to use in saving the results, will over write.\n"+
				"-f Path to the reference fasta with associated xxx.fai index\n"+
				"-u UCSC RefFlat or RefSeq transcript (not merged genes) file, full path. See RefSeq \n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (uniqueName1 name2(optional) chrom\n" +
				"       strand txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts\n" +
				"       (commaDelimited)exonEnds). Example: ENSG00000183888 C1orf64 chr1 + 16203317\n" +
				"       16207889 16203385 16205428 2 16203317,16205000 16203467,16207889 .\n"+
				"-m Full path directory name containing the me2x3acc1-9, splice5sequences and me2x5\n"+
				"       splice model files. See USeqDocumentation/splicemodels/ or \n"+
				"       http://genes.mit.edu/burgelab/maxent/download/ \n"+
				"\n"+

				"Example: java -Xmx10G -jar ~/USeq/Apps/KnownSpliceJunctionScanner -f ~/Hg19/hg19.fasta\n"+
				"       -r ~/exm2.bed.gz -m ~/USeq/Documentation/splicemodels -u ~/hg19EnsTrans.ucsc.zip\n"+
	

				"\n**************************************************************************************\n");
	}




	public IndexedFastaSequenceFile getFasta() {
		return fasta;
	}


	public HashMap<String, UCSCGeneLine[]>  getChromTranscripts() {
		return chromTranscripts;
	}


	public File getSpliceModelDirectory() {
		return spliceModelDirectory;
	}

}
