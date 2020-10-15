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
import edu.utah.seq.vcf.xml.SimpleVcf;
import util.bio.annotation.ExonIntron;
import util.bio.annotation.ExportIntergenicRegions;
import util.bio.parsers.UCSCGeneLine;
import util.bio.parsers.UCSCGeneModelTableReader;
import util.bio.seq.Seq;
import util.gen.*;

/**Simplified and threaded version of the VCFSpliceAnnotator
 * @author Nix
 * */
public class VCFSpliceScanner {

	//user fields
	private File vcfFile;
	private File spliceModelDirectory;
	private IndexedFastaSequenceFile fasta; 
	private File transcriptSeqFile;
	private File annotatedVcfFile = null;
	private File tempDirectory;

	private double minNewAltScore = 3;
	private double minNewScoreDelta = 1;
	private double maxDamagedAltScore = 3;
	private double minDamagedScoreDelta = 1;
	
	private boolean removeInfoDropNonAffected = false;
	private boolean scoreNovelIntronJunctions = false;
	private boolean scoreNovelExonJunctions = true;
	private boolean scoreNovelSpliceJunctionsInSplice = true;
	private short vcfExportCategory = 2;
	private int numberThreads = 0;
	private int chunkSize = 50;
	
	//internal fields
	private HashMap<String,UCSCGeneLine[]> chromGenes;
	private String workingChromosomeName = "";
	private String workingSequence = null;
	private ArrayList<SimpleVcf> workingVcfRecords = new ArrayList<SimpleVcf>();
	private int numWorkingRecords;
	private int workingVcfRecordIndex;
	private boolean allVcfsProcessed;
	private SimpleVcf simpleVcf;
	private HashSet<String> processedChromosomes = new HashSet<String>();
	private SpliceAnnotationLoader[] loaders;
	
	private UCSCGeneLine[] workingTranscripts;
	private Gzipper vcfOut;
	
	private int numTranscripts = 0; 
	private HashSet<String> intersectingTranscriptNames = new HashSet<String>();
	private int numVariantsScanned = 0;
	private int numVariantsIntersectingTranscripts = 0;
	private int numVariantsIntersectingExons = 0;
	private int numVariantsIntersectingIntrons = 0;
	private int numVariantsIntersectingSJs = 0;
	private int numExonSJsGained = 0;
	private int numIntronSJsGained = 0;
	private int numSpliceJunctionSJsGained = 0;
	private int numSJsLost = 0;
	private int numVariantsEffectingSJs = 0;
	private ArrayList<File> gzippedToConcat = new ArrayList<File>();

	//constructor
	public VCFSpliceScanner(String[] args){
		try {
			//start clock
			long startTime = System.currentTimeMillis();

			//process args
			processArgs(args);
			printThresholds();

			//load transcripts
			System.out.println("Loading transcripts...");
			loadTranscripts();

			annotateVCFWithSplices();
			concatTempVcf();
			printSummary();
			IO.deleteDirectory(tempDirectory);
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
			System.out.println("\nDone "+Math.round(diffTime)+" min!\n");
		} catch (Exception e){
			e.printStackTrace();
		}
	}


	private void concatTempVcf() throws IOException {
		File toSortVcf = new File(tempDirectory, "toSort.vcf.gz");
		IO.concatinateFiles(gzippedToConcat, toSortVcf);
		System.out.println("Sorting vcf...");
		SimpleVcf.sortVcf(toSortVcf, annotatedVcfFile);
	}

	private void printThresholds(){
		StringBuilder sb = new StringBuilder();
		sb.append("Threholds:\n");
		sb.append(minNewAltScore +"\tMinimum new splice junction threshold.\n");
		sb.append(minNewScoreDelta +"\tMinimum new splice junction score difference, new score - refseq score.\n");
		sb.append(maxDamagedAltScore +"\tMaximum damaged splice junction score.\n");
		sb.append(minDamagedScoreDelta +"\tMinimum  damaged score difference, refseq score - new score.\n");
		sb.append(scoreNovelIntronJunctions +"\tLook for novel intron splice junctions, outside of known splice junctions.\n");
		sb.append(scoreNovelExonJunctions +"\tLook for novel splice junctions in exons.\n");
		sb.append(scoreNovelSpliceJunctionsInSplice +"\tLook for novel splice junctions inside known splice junctions.\n");
		sb.append(this.removeInfoDropNonAffected +"\tPrint only variants that effect splice junctions and remove unnecessary info for downstream annotators.\n");
		sb.append(vcfExportCategory +"\tExport catagory for adding info to vcf records.\n");
		System.out.println(sb);
	}
	
	private void printSummary(){
		int total = numVariantsIntersectingIntrons+ numVariantsIntersectingExons + numVariantsIntersectingSJs;
		StringBuilder sb = new StringBuilder();
		sb.append("\nSummary stats:\n");
		sb.append(numTranscripts+"\tTranscripts\n"); 
		sb.append(intersectingTranscriptNames.size()+"\tTranscripts intersecting variants\n"); 
		sb.append(numVariantsScanned+"\tVariants (including alternates)\n");
		sb.append(numVariantsIntersectingTranscripts+"\tVariants intersecting transcripts, no repeat counting.\n");
		sb.append(numVariantsEffectingSJs+"\tVariants effecting splice junctions.\n");
		sb.append("\nCounts including repeats with overlapping transcripts:\n");
		sb.append(total + "\tVariants intersecting annotations\n");
		if (scoreNovelExonJunctions) sb.append(numVariantsIntersectingExons+"\tVariants intersecting exons (non splice junction)\n"); 
		if (scoreNovelIntronJunctions) sb.append(numVariantsIntersectingIntrons+"\tVariants intersecting introns (non splice junction)\n"); 
		sb.append(numVariantsIntersectingSJs+"\tVariants intersecting splice junctions\n"); 
		if (scoreNovelExonJunctions) sb.append(numExonSJsGained+"\tExonic splice junctions gained\n"); 
		if (scoreNovelIntronJunctions) sb.append(numIntronSJsGained+"\tIntronic splice junctions gained\n"); 
		if (scoreNovelSpliceJunctionsInSplice) sb.append(numSpliceJunctionSJsGained+"\tKnown splice junctions with novel splice junction gained\n"); 
		sb.append(numSJsLost+"\tKnown splice junctions damaged\n");
		System.out.println(sb);
	}
	
	/*
	 * 	private boolean scoreNovelExonJunctions = true;
	private boolean scoreNovelSpliceJunctionsInSplice = true;
	 */
	
	
	private boolean loadWorkingVcfRecords(BufferedReader in) throws Exception {
		//clear old
		workingVcfRecords.clear();
		workingVcfRecordIndex = 0;
		 
		//add initial vcf
		if (simpleVcf == null) {
			simpleVcf = new SimpleVcf(in.readLine(), 0);
		}
		
		//add new chrom
		workingChromosomeName = simpleVcf.getChr();
		workingVcfRecords.add(simpleVcf);
		
		//already processed?
		if (processedChromosomes.contains(workingChromosomeName)){
			IO.deleteDirectory(tempDirectory);
			throw new Exception("Error: "+workingChromosomeName+" chromosome has already been processed, is your vcf file sorted by chromosome and position? Aborting.");
		}
		processedChromosomes.add(workingChromosomeName);
		
		//load em until chrom changes
		String vcfLine;
		allVcfsProcessed = true;
		while ((vcfLine = in.readLine()) != null){
			simpleVcf = new SimpleVcf(vcfLine, 0);
			if (workingChromosomeName.equals(simpleVcf.getChr())) workingVcfRecords.add(simpleVcf);
			else {
				allVcfsProcessed = false;
				break;
			}
		}
		numWorkingRecords = workingVcfRecords.size();
		if (numWorkingRecords == 0) return false;
		return true;
	}

	public void annotateVCFWithSplices() {
		BufferedReader in = null;
		
		int numVcfRecords = (int) IO.countNonBlankOrCommentLines(vcfFile);
		
		//create some loaders
		int numLoadersToCreate = (int)Math.round((double)numVcfRecords/ (double)chunkSize) +1;
		if (numLoadersToCreate > numberThreads) numLoadersToCreate = numberThreads;
		loaders = new SpliceAnnotationLoader[numLoadersToCreate];
		for (int i=0; i< loaders.length; i++) loaders[i] = new SpliceAnnotationLoader(this, new File(tempDirectory, i+"_.vcf.gz"));
		
		System.out.println(numLoadersToCreate+" Annotators created...");
		
		try {
			in  = IO.fetchBufferedReader(vcfFile);
			
			//add ##INFO line and find "#CHROM" line 
			loadAndModifyHeader(in);

			
			//for each chromosome of vcf records
			
			while (allVcfsProcessed == false && loadWorkingVcfRecords(in)){
				
				//load chrom specific data
				loadChromosomeData(workingChromosomeName);
				
				//check if any transcripts were found
				if (workingTranscripts == null) {
					for (SimpleVcf r: workingVcfRecords) vcfOut.println(r.getOriginalRecord());
					continue;
				}
				
				//launch em!
				ExecutorService executor = Executors.newFixedThreadPool(numLoadersToCreate);
				for (SpliceAnnotationLoader l: loaders) executor.execute(l);
				executor.shutdown();
				while (!executor.isTerminated()) {}  //wait here until complete

				//check loaders 
				for (SpliceAnnotationLoader l: loaders) {
					if (l.isFailed()) throw new IOException("ERROR: Failed to annotate splices with one of the loaders?");
				}
				System.out.println();
			}
			
			//clean up
			in.close();
			vcfOut.close();
			
			//shut down the writers, fetch their files, and pull stats data
			for (SpliceAnnotationLoader l: loaders) {
				l.getVcfOut().close();
				gzippedToConcat.add(l.getVcfOut().getGzipFile());
				appendFields(l);
			}
			
		}catch (Exception e) {
			e.printStackTrace();
			IO.deleteDirectory(tempDirectory);
			System.exit(1);
		}
	}
	
	private void appendFields(SpliceAnnotationLoader l) {
		numVariantsScanned +=  l.getNumVariantsScanned();
		numVariantsIntersectingTranscripts +=  l.getNumVariantsIntersectingTranscripts();
		numVariantsIntersectingIntrons +=  l.getNumVariantsIntersectingIntrons();
		numVariantsIntersectingSJs +=  l.getNumVariantsIntersectingSJs();
		numVariantsIntersectingExons +=  l.getNumVariantsIntersectingExons();
		numSJsLost +=  l.getNumSJsLost();
		numExonSJsGained +=  l.getNumExonSJsGained();
		numSpliceJunctionSJsGained +=  l.getNumSpliceJunctionSJsGained();
		numIntronSJsGained +=  l.getNumIntronSJsGained();
		intersectingTranscriptNames.addAll(l.getIntersectingTranscriptNames());
		numVariantsEffectingSJs += l.getNumVariantsEffectingSJs();
	}


	/*The threads use this to pull vcf records to process.*/
	public synchronized boolean loadRecords(ArrayList<SimpleVcf> al){	
		int count = 0;
		for (; workingVcfRecordIndex< numWorkingRecords; workingVcfRecordIndex++){
			al.add(workingVcfRecords.get(workingVcfRecordIndex));
			if (++count >= chunkSize) {
				workingVcfRecordIndex++;
				break;
			}
		}
		
		if (al.size() != 0) {
			System.out.print(".");
			return true;
		}
		return false;
	}
	
	/*TODO: replace with IntervalST*/
	public void loadTranscripts(){
		//main transcripts
		UCSCGeneModelTableReader reader = new UCSCGeneModelTableReader(transcriptSeqFile, 0);
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your transcript's coordinates are reversed. Check that each start is less than the stop.\n");
		chromGenes = reader.getChromSpecificGeneLines();
		numTranscripts = reader.getGeneLines().length;
	}

	private void loadChromosomeData(String chrom) {
		try {
		//set new name
		workingChromosomeName = chrom;
		
		//find and load sequence
		ReferenceSequence p = fasta.getSequence(workingChromosomeName);
		if (p == null ) throw new IOException ("\n\nFailed to find or load a fasta sequence for '"+workingChromosomeName+"', aborting.\n");
		workingSequence = new String(p.getBases());
		workingSequence = workingSequence.toUpperCase();
		
		//fetch transcripts
		workingTranscripts = chromGenes.get(workingChromosomeName);
		if (workingTranscripts == null) {
			System.out.println("\tWARNING: no transcripts found for "+workingChromosomeName+", skipping all associated vcf records.");
		}
		else System.out.print("\tAnnotating: "+workingChromosomeName+"\tLen: "+workingSequence.length()+"\tTrans: "+workingTranscripts.length);

		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nDo your chromosome names in the gene table and vcf file match? Watch for 'chr'. "+workingChromosomeName);
		}
	}

	public void loadAndModifyHeader(BufferedReader in) throws IOException{
		File header = new File(tempDirectory, "headerPlusUnAnnotated.vcf.gz").getCanonicalFile();
		header.createNewFile();
		gzippedToConcat.add(header);
		
		vcfOut = new Gzipper(header);
		boolean addedInfo = false;
		boolean foundChrom = false;
		String line;
		while ((line=in.readLine()) != null){
			//comments
			if (line.startsWith("#")){
				//add info lines?
				if (line.startsWith("##INFO=")){
					if (addedInfo == false) { 
					addedInfo = true;
					vcfOut.println("##INFO=<ID=VCFSS,Number=.,Type=String,Description=\"USeq VCFSpliceScanner output. "
							+ "One or more splice junction (SJ) effects delimited by a &. Each unit consists of one or "
							+ "more comma delimited transcript names, a colon, and the affect info. The affect info is comma "
							+ "delimited containing Type, Positon, RefSeq, AltSeq, RefScore, AltScore, and ScoreDelta. "
							+ "SJ Types are a three char string with Gain or Damage; 3 or 5 prime; Intron, Exon or Splice "
							+ "(e.g. D3S, G5E, G3S). The SJ Position is the interbase coordinate position of the damaged or "
							+ "gained SJ for SNVs.  For INDELS, it is the variant position. The scores are the MaxEntScan "
							+ "values for the reference and alternate sequences. See Yeo and Burge 2004, "
							+ "http://www.ncbi.nlm.nih.gov/pubmed/15285897\">");
					}
					//print the info?
					if (removeInfoDropNonAffected == false) vcfOut.println(line);
				}
				//Format lines
				else if (line.startsWith("##FORMAT") || line.startsWith("##FILTER")){
					if(removeInfoDropNonAffected == false) vcfOut.println(line);
				}
				
				//Chrome line
				else if (line.startsWith("#CHROM")){
					foundChrom = true;
					if (removeInfoDropNonAffected) vcfOut.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
					else vcfOut.println(line);
					break;
				}

				//must be something else so just print it.
				else vcfOut.println(line);
			}
		}
		vcfOut.flush();
		if (foundChrom == false) throw new IOException("\tError: Failed to find the #CHROM header line? Aborting.\n");
		if (addedInfo == false) throw new IOException("\tError: Failed to find any ##INFO header lines? Aborting.\n");
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFSpliceScanner(args);
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
					case 'v': vcfFile = new File(args[++i]); break;
					case 'f': indexedFasta = new File(args[++i]); break;
					case 'm': spliceModelDirectory = new File(args[++i]); break;
					case 'u': transcriptSeqFile = new File(args[++i]); break;
					case 'r': annotatedVcfFile = new File(args[++i]); break;
					case 'i': scoreNovelIntronJunctions = true; break;
					case 'x': vcfExportCategory = Short.parseShort(args[++i]); break;
					case 'a': minNewAltScore = Double.parseDouble(args[++i]); break;
					case 'b': minNewScoreDelta = Double.parseDouble(args[++i]); break;
					case 'c': maxDamagedAltScore = Double.parseDouble(args[++i]); break;
					case 'd': minDamagedScoreDelta = Double.parseDouble(args[++i]); break;
					case 's': removeInfoDropNonAffected = true; break;
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

		if (vcfFile == null || vcfFile.canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz) file!\n");
		
		//temp directory
		if (annotatedVcfFile == null) Misc.printErrAndExit("\nPlease provide the path and name of a gzipped vcf to use in writing the results.\n");
		if (annotatedVcfFile.getName().endsWith(".gz") == false) Misc.printErrAndExit("\nSorry your output vcf file must end in xxx.gz\n");
		File parentDir = annotatedVcfFile.getParentFile();
		if (parentDir == null) parentDir = new File (System.getProperty("user.dir"));
		tempDirectory = new File (parentDir, "TempVCFSpliceScannerDir_DeleteMe_"+Misc.getRandomString(6));
		tempDirectory.mkdir();
		
		//threads to use
		double totalGbAvailable = (double)(Runtime.getRuntime().maxMemory()/1000000000.0);
		int numPossCores = (int)Math.round(totalGbAvailable/10.0);
		if (numPossCores < 1) numPossCores = 1;
		int numPossThreads = Runtime.getRuntime().availableProcessors();
		
		if (numPossCores <= numPossThreads) numberThreads = numPossCores;
		else numberThreads = numPossThreads;
		
		
		System.out.println("Core usage:\n\tTotal GB available to Java:\t"+ Num.formatNumber(totalGbAvailable, 1));
		System.out.println("\tTotal available cores:\t"+numPossThreads);
		System.out.println("\tNumber cores to use @ 10GB/core:\t"+numberThreads+"\n");
		
		
		//flip booleans?
		setExportBooleans();
	}	



	private void setExportBooleans() {
		if (vcfExportCategory == 0){
			scoreNovelExonJunctions = true;
			scoreNovelSpliceJunctionsInSplice = true;
			scoreNovelIntronJunctions = true;
		}
		else if (vcfExportCategory == 1){
			scoreNovelExonJunctions = false;
			scoreNovelSpliceJunctionsInSplice = false;
			scoreNovelIntronJunctions = false;
		}
		else if (vcfExportCategory == 2){
			scoreNovelExonJunctions = true;
			scoreNovelSpliceJunctionsInSplice = true;
			scoreNovelIntronJunctions = false;
		}
		else Misc.printErrAndExit("\nLo siento, yo no comprendo su Export Category?\n");
		
	}
	
	private static void scoreKnownSplices(HashMap<String,UCSCGeneLine[]> chromGenes, String workingSequence, String workingChromosomeName, File spliceModelDirectory) {

		int lenMinOne = workingSequence.length() - 1;
		HashSet<String> scored5Plus = new HashSet<String>();
		HashSet<String> scored3Plus = new HashSet<String>();
		HashSet<String> scored5Minus = new HashSet<String>();
		HashSet<String> scored3Minus = new HashSet<String>();
		Histogram spliceHistogram5 = new Histogram(-14, 14, 2500);
		Histogram spliceHistogram3 = new Histogram(-14, 14, 2500);

		MaxEntScanScore5 score5 = new MaxEntScanScore5(spliceModelDirectory);
		MaxEntScanScore3 score3 = new MaxEntScanScore3(spliceModelDirectory);

		//for each transcript
		for (UCSCGeneLine g: chromGenes.get(workingChromosomeName)){
			ExonIntron[] introns = g.getIntrons();
			if (introns == null) continue;
			boolean plusStrand = g.getStrand().equals("+");
			//mask both 5' and 3' junctions
			int startJunc = 0;
			int endJunc = 0;
			for (int i=0; i< introns.length; i++){
				//plus strand 
				if (plusStrand){
					//5'
					startJunc = introns[i].getStart()-3;     
					endJunc = introns[i].getStart()+6;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					String coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored5Plus.contains(coor) == false){
						scored5Plus.add(coor);

						String seq = workingSequence.substring(startJunc, endJunc);
						double score = score5.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE) {
							spliceHistogram5.count(score);
							System.out.println(workingChromosomeName+"\t"+startJunc+"\t"+endJunc+"\t5_"+seq+"\t"+score+"\t+");
						}


					}
					//3'
					startJunc = introns[i].getEnd()-20;
					endJunc = introns[i].getEnd()+3;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					coor = startJunc+"-"+endJunc;
					if (scored3Plus.contains(coor) == false){
						scored3Plus.add(coor);

						String seq = workingSequence.substring(startJunc, endJunc);
						double score = score3.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE) spliceHistogram3.count(score);
						//System.out.println("+ 3\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);
						System.out.println(workingChromosomeName+"\t"+startJunc+"\t"+endJunc+"\t3_"+seq+"\t"+score+"\t+");

					}
				}
				//minus strand
				else {
					//3'
					startJunc = introns[i].getStart()-3;
					endJunc = introns[i].getStart()+20;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					String coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored3Minus.contains(coor) == false){
						scored3Minus.add(coor);
						String seq = workingSequence.substring(startJunc, endJunc);
						seq = Seq.reverseComplementDNA(seq);
						double score = score3.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE) spliceHistogram3.count(score);
						//System.out.println("- 3\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);
						System.out.println(workingChromosomeName+"\t"+startJunc+"\t"+endJunc+"\t3_"+seq+"\t"+score+"\t-");

					}

					//5'
					startJunc = introns[i].getEnd()-6;
					endJunc = introns[i].getEnd()+3;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored5Minus.contains(coor) == false){
						scored5Minus.add(coor);
						String seq = workingSequence.substring(startJunc, endJunc);
						seq = Seq.reverseComplementDNA(seq);
						double score = score5.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE) {
							spliceHistogram5.count(score);
							//System.out.println("- 5\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);
							System.out.println(workingChromosomeName+"\t"+startJunc+"\t"+endJunc+"\t5_"+seq+"\t"+score+"\t-");
						}
					}
				}
			}
		}
	}



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            VCF Splice Scanner : Sept 2020                        **\n" +
				"**************************************************************************************\n" +
				"Scores variants for changes in splicing using the MaxEntScan algorithms. See Yeo and\n"+
				"Burge 2004, http://www.ncbi.nlm.nih.gov/pubmed/15285897 for details. Known splice\n"+
				"acceptors and donors are scored for loss of a junction.  Exonic, intronic, and splice\n"+
				"bases are scanned for novel junctions in a window around each variant. See the vcf\n"+
				"INFO header for a description of the output. Use this information to identify \n"+
				"snv and indel variants that may effect splicing.\n\n" +

				"Required Options:\n"+
				"-v Path to a vcf file to annotate (xxx.vcf(.gz/.zip OK)).\n"+
				"-r Name of a gzipped vcf file to use in saving the results, will over write.\n"+
				"-f Path to the reference fasta with xxx.fai index\n"+
				"-u UCSC RefFlat or RefSeq transcript (not merged genes) file, full path. See RefSeq \n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (uniqueName1 name2(optional) chrom\n" +
				"       strand txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts\n" +
				"       (commaDelimited)exonEnds). Example: ENSG00000183888 C1orf64 chr1 + 16203317\n" +
				"       16207889 16203385 16205428 2 16203317,16205000 16203467,16207889 .\n"+
				"-m Full path directory name containing the me2x3acc1-9, splice5sequences and me2x5\n"+
				"       splice model files. See USeqDocumentation/splicemodels/ or \n"+
				"       http://genes.mit.edu/burgelab/maxent/download/ \n"+
				"\n"+
				
				"Optional options:\n"+
				"-x Export category for adding info to the vcf file, defaults to 2:\n"+
				"       0 All types (gain or damaged in exon, intron, and splice)\n"+
				"       1 Just report damaged splices\n"+
				"       2 Report damaged splices and novel splices in exons and splice junctions.\n"+
				"-a Minimum new splice junction score in Alt, max score in Ref, defaults to 3.\n"+
				"-b Minimum new splice junction score difference, new - refseq, defaults to 1.\n"+
				"-c Maximum damaged splice junction score, defaults to 3.\n"+
				"-d Minimum damaged score difference, refseq - new, defaults to 1.\n"+
				"-s Format vcf with minimal output for downstream annotators, e.g. snpEff.\n"+

				"\n"+
				
				"Example: java -Xmx20G -jar ~/USeq/Apps/VCFSpliceAnnotator -f ~/Hg19/Fa/ -v ~/exm2.vcf\n"+
				"       -m ~/USeq/Documentation/splicemodels -i -u ~/Hg19/hg19EnsTrans.ucsc.zip -r\n"+
				"       ~/ExmSJAnno/exm2VSSAnno.vcf.gz -x 0\n"+

				"**************************************************************************************\n");
	}
	
	//getters and setters
	public boolean isScoreNovelIntronJunctions() {
		return scoreNovelIntronJunctions;
	}

	public boolean isScoreNovelExonJunctions() {
		return scoreNovelExonJunctions;
	}

	public boolean isScoreNovelSpliceJunctionsInSplice() {
		return scoreNovelSpliceJunctionsInSplice;
	}

	public String getWorkingSequence() {
		return workingSequence;
	}
	
	public File getSpliceModelDirectory() {
		return spliceModelDirectory;
	}

	public UCSCGeneLine[] getWorkingTranscripts() {
		return workingTranscripts;
	}

	public short getVcfExportCategory() {
		return vcfExportCategory;
	}


	public double getMinNewAltScore() {
		return minNewAltScore;
	}
	public double getMinNewScoreDelta() {
		return minNewScoreDelta;
	}
	public double getMaxDamagedAltScore() {
		return maxDamagedAltScore;
	}
	public double getMinDamagedScoreDelta() {
		return minDamagedScoreDelta;
	}


	public boolean isRemoveInfoDropNonAffected() {
		return removeInfoDropNonAffected;
	}
	
}
