package edu.utah.seq.barcodes;

import htsjdk.samtools.SAMRecord;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.data.sam.SamAlignmentFlags;
import util.gen.Gzipper;
import util.gen.Histogram;
import util.gen.Misc;
import util.gen.Num;

public class BarcodeChromMerger implements Runnable {

	//constants for this worker
	private Consensus consensus;
	private BarcodeLoader barcodeLoader;
	private BarcodeClusterEngine clusterEngine;
	private ConsensusEngine consensusEngine;
	private int threadId;
	private Gzipper passingSamWriter;
	private Gzipper failingSamWriter;
	private Gzipper fastqReadOneWriter;
	private Gzipper fastqReadTwoWriter;
	private Gzipper fastqUnpairedWriter;
	private long numBadBarcodes = 0;
	private long numFastqPairsOut = 0;
	private long numFastqUnpairedOut = 0;
	private long numPassingSamOut = 0;
	private long numFailingSamOut = 0;
	private long totalPassingBps = 0;
	private Histogram famSizeHist = new Histogram(0, 100, 100);
	private static final Pattern alignmentScore = Pattern.compile("\\s+AS:i:(\\d+)\\s*");
	
	//reset with each chunk
	private String chromosome;
	private int start;
	private int stop;
	private ArrayList<SAMRecord>[] orderedRecords;

	//working objects
	private ArrayList<SAMRecord> workingRecords = null;
	private ArrayList<SAMRecord> workingMates = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> firstForward = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> firstReverse = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> secondForward = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> secondReverse = new ArrayList<SAMRecord>();
	private boolean workFailed;
	private File passingSamFile;
	private File failingSamFile;
	private File fastqFileReadOne;
	private File fastqFileReadTwo;
	private File fastqFileUnpaired;
	
	public BarcodeChromMerger(Consensus consensus, int threadId) throws Exception {
		this.consensus = consensus;
		barcodeLoader = consensus.getBarcodeLoader();
		this.threadId = threadId;
		
		//make engines
		clusterEngine = new BarcodeClusterEngine(consensus.getMinBarcodeBaseQuality(), 
				consensus.getMinNumPassingBarcodeBases(), consensus.getMinFractionBarcodeIdentity());
		consensusEngine = new ConsensusEngine(consensus.getMinReadBaseQuality(), consensus.getMinReadBaseFractionIdentity());
		
		makeReadersAndWriters();
	}
	
	public void setChunkFields(BarcodeChunk chunk){
		//fetch fields from chunk
		chromosome = chunk.getChromosome();
		start = chunk.getStart();
		stop = chunk.getStop();
		orderedRecords = chunk.getOrderedRecords();
	}
	
	public void run() {
		try {
			while (true){
				long startTime = System.currentTimeMillis();
				//get and load a BarcodeChunk
				BarcodeChunk chunk = barcodeLoader.loadNextChunk();
				if (chunk == null) {
					incrementCounters();
					closeIO();
					return;
				}
				setChunkFields(chunk);

				processRecords();
				
				//print run info
				System.out.println("\t"+toString()+" "+memoryUsed(startTime));
			}
		} catch (Throwable ex){
			workFailed = true;
			System.err.println("\nError encountered when processing chromosome "+toString()+"\n");
			ex.printStackTrace();
		} 
	}
	
	/**Returns the amount of memory being used.*/
	public static String memoryUsed(long startTime){
		System.gc();
		Runtime rt = Runtime.getRuntime();
		double usedGB = (double)(rt.totalMemory() - rt.freeMemory()) / 1024.0 / 1024.0/ 1024.0;
		long endTime = System.currentTimeMillis();
		return Num.formatNumberOneFraction(usedGB)+"GB\t"+(endTime-startTime)/1000+"S";
	}
	
	private void incrementCounters() throws Exception {
		barcodeLoader.incrementNumFastqPairsOut(numFastqPairsOut);
		barcodeLoader.incrementNumFastqUnpairedOut(numFastqUnpairedOut);
		barcodeLoader.incrementNumPassingSamOut(numPassingSamOut);
		barcodeLoader.incrementNumFailingSamOut(numFailingSamOut);
		barcodeLoader.incrementNumBadBarcodes(numBadBarcodes);
		barcodeLoader.incrementFamSizeHist(famSizeHist);
		barcodeLoader.incrementNumModifiedBps(consensusEngine.getBpsModified());
		barcodeLoader.incrementTotalPassingBps(totalPassingBps);
	}
	
	public void closeIO() throws IOException{
		passingSamWriter.close();
		failingSamWriter.close();
		fastqReadOneWriter.close();
		fastqReadTwoWriter.close();
		fastqUnpairedWriter.close();
	}
	
	private void makeReadersAndWriters() throws Exception {
		//bam writers
		passingSamFile = new File (consensus.getSaveDirectory(), threadId+"_passing.sam.gz");
		passingSamWriter = new Gzipper(passingSamFile);
		failingSamFile = new File (consensus.getSaveDirectory(), threadId+"_failing.sam.gz");
		failingSamWriter = new Gzipper(failingSamFile);

		//fastq
		fastqFileReadOne = new File (consensus.getSaveDirectory(), threadId+"_merged_1.fastq.gz");
		fastqFileReadTwo = new File (consensus.getSaveDirectory(), threadId+"_merged_2.fastq.gz");
		fastqFileUnpaired = new File (consensus.getSaveDirectory(), threadId+"_merged_unpaired.fastq.gz");
		fastqReadOneWriter = new Gzipper(fastqFileReadOne);
		fastqReadTwoWriter = new Gzipper(fastqFileReadTwo);
		fastqUnpairedWriter = new Gzipper(fastqFileUnpaired);
		
		passingSamFile.deleteOnExit();
		failingSamFile.deleteOnExit();
		fastqFileReadOne.deleteOnExit();
		fastqFileReadTwo.deleteOnExit();
		fastqFileUnpaired.deleteOnExit();
	}

	private void processRecords() throws Exception {

		//walk through each base position
		for (int i=0; i< orderedRecords.length; i++) {
			//any records?
			workingRecords = orderedRecords[i];
			if (workingRecords == null) continue;
			
			/*
			for (SAMRecord sam: workingRecords){
				if (sam.getReadName().equals("NS500690:33:H5W2KBGXX:4:12606:23383:4972:BMF:TGGATAGACTGCGATC@CC645F3@@F@F3C/")) {
					System.out.println("Found it in working\n"+sam.getSAMString());
				}
			}*/
			
			//process the records found at this particular unclipped start position
			processRecordsAtOnePosition();

		}
	}

	private void processRecordsAtOnePosition() throws Exception {

		//just one record? or a group that might need to be split for clustering
		if (workingRecords.size() == 1) {
			SAMRecord sam = workingRecords.get(0);
//if (sam.getReadName().equals("NS500690:33:H5W2KBGXX:1:21203:7637:18198:BMF:AGACTTGTGCTTCTGT@CCCCFFFFFEFFFFF")) System.out.println("proc solo");
			processSolo(sam);
		}
		else processGroup();
	}

	private void processGroup() throws Exception {
		//clear and split records into the four types of alignments
		clearAndSplit();
		
		//for each group, are there enough to cluster collapse or are they solo?
		//first forward
		int numRec = firstForward.size();
		if (numRec > 0){
			if (numRec == 1) processSolo (firstForward.get(0));
			else clusterCollapseFirstOfPair( firstForward );
		}
		
		//first reverse
		numRec = firstReverse.size();
		if (numRec > 0){
			if (numRec == 1) processSolo (firstReverse.get(0));
			else clusterCollapseFirstOfPair( firstReverse );
		}
		
		//second forward
		numRec = secondForward.size();
		if (numRec > 0){
			if (numRec == 1) processSolo (secondForward.get(0));
			else clusterCollapseSecondOfPair( secondForward );
		}
		
		//second reverse
		numRec = secondReverse.size();
		if (numRec > 0){
			if (numRec == 1) processSolo (secondReverse.get(0));
			else clusterCollapseSecondOfPair( secondReverse );
		}
	}

	
	
	private void clusterCollapseSecondOfPair(ArrayList<SAMRecord> records) throws Exception {
		//for (SAMRecord sam: records){
			//if (sam.getReadName().equals("NS500690:33:H5W2KBGXX:4:12602:2299:19865:BMF:TAGAGGGGAGAGTGCC@CCCCFFFFFCFFFFF")) {
				//System.out.println("Found it in cluster collapse SecondOfPair\n"+sam.getSAMString());
			//}
		//}
		
		//the only time we want to collapse secondOfPairs are if they have no aligned mate in the file
		//if they do then skip them since they'll get processed with the first of pair
		parseNoMateRecordsIntoWorkingMates(records);
		
		//any remaining?
		int num = workingMates.size();
		if (num != 0){
			if (num == 1) {
				writeOutSingle(workingMates.get(0));
				famSizeHist.count(1);
			}
			else {
				//fetch SAMRecord[cluster number][SAMRecords in a given cluster]
				SAMRecord[][] clusters = clusterEngine.cluster(records);
				
				//any with bad barcodes?
				processBadBarcodeRecords(clusterEngine.getSkippedRecords());
				
				//collapse clusters
				for (SAMRecord[] c2Col : clusters){
					int familySize = c2Col.length;
					//just one, thus no cluster
					if (familySize == 1) writeOutSingle(c2Col[0]);
					
					else {
						//call consensus on these secondOfPairs with no mapped mates in file
						String[] fastq = consensusEngine.callConsensus(c2Col);
						
						//fetch unmapped mates, these would be first of pair
						SAMRecord[] unmapped = collectMates(c2Col);
						// no mates? write out as single end, otherwise cluster and write as paired
						if (unmapped == null) writeOutSingleEndFastq(fastq);
						else {
							if (unmapped.length != familySize) Misc.printErrAndExit("Error: mate numbers don't match? "+familySize+" vs "+unmapped.length);
							writeOutPairedEndFastq(consensusEngine.callConsensus(unmapped), fastq);
						}
					}
					//increment family size in histogram
					famSizeHist.count(familySize);
				}
			}
		}	
	}


	private void clusterCollapseFirstOfPair(ArrayList<SAMRecord> records) throws Exception {
		
		//for (SAMRecord sam: records){
			//if (sam.getReadName().equals("NS500690:33:H5W2KBGXX:4:12606:23383:4972:BMF:TGGATAGACTGCGATC@CC645F3@@F@F3C/")) {
				//System.out.println("Found it in cluster collapse\n"+sam.getSAMString());
			//}
		//}
		
		//this is a first of pair so cluster
		//fetch SAMRecord[cluster number][SAMRecords in a given cluster]
		SAMRecord[][] clusters = clusterEngine.cluster(records);		
		
		//any bad barcodes?
		processBadBarcodeRecords(clusterEngine.getSkippedRecords());

		//for each cluster, call consensus?
		for (SAMRecord[] c2Col : clusters){
			
			//just one record, thus no collapsing
			if (c2Col.length == 1) processSolo (c2Col[0]);
	
			else {
				//call consensus on these firstOfPairs
				String[] fastq = consensusEngine.callConsensus(c2Col);
//if (fastq[0].contains("GGCGGGAGTCGGTGTA")) Misc.printArray(fastq);
				//attempt to collect mates of these first of pairs 
				SAMRecord[] mates = collectMates(c2Col);
				if (mates == null) writeOutSingleEndFastq(fastq);
				else {
					if (mates.length != c2Col.length) Misc.printErrAndExit("Error: mate numbers don't match for first of pair? "+c2Col.length+" vs "+mates.length);
					writeOutPairedEndFastq(fastq, consensusEngine.callConsensus(mates));
				}
				//increment family size in histogram
				famSizeHist.count(c2Col.length);
			}
		}
	}

	private SAMRecord[] collectMates(SAMRecord[] c2Col) throws Exception {
		workingMates.clear();
		for (SAMRecord sam: c2Col){
			//any mapped mate?
			String mate = restoreSam(sam);
			if (mate != null) {
				insertMateInfo(sam, mate);
				workingMates.add(sam);
			}
		}
		if (workingMates.size() == 0) return null;
		
		SAMRecord[] mates = new SAMRecord[workingMates.size()];
		workingMates.toArray(mates);
		return mates;
	}
	
	private void parseNoMateRecordsIntoWorkingMates(ArrayList<SAMRecord> secondOfPairs) throws Exception {
		workingMates.clear();
		for (SAMRecord sam: secondOfPairs){
			String mate = fetchMate(sam);
			//no mate?
			if (mate == null) workingMates.add(sam);
			else {
				//OK, it has a mate, but is it mapped?
				SamAlignmentFlags saf = fetchBitFlagObj(mate);
				if (saf.isUnmapped()) workingMates.add(sam);
			}
		}
	}
	
	private SamAlignmentFlags fetchBitFlagObj(String trunkSamRecord) throws Exception{
		String[] fields = Misc.WHITESPACE.split(trunkSamRecord);
		Short flag = Short.parseShort(fields[0]);
		return new SamAlignmentFlags(flag.shortValue());
	}

	/**Writes to fail the records and their mates if present.*/
	private void processBadBarcodeRecords(ArrayList<SAMRecord> skippedRecords) throws Exception {
		int numBad = skippedRecords.size();
		if (numBad !=0) {
			//write out skipped records and possibly their mates as bam
			for (SAMRecord sam: skippedRecords){
				String mateString = restoreSam(sam);
				writeOutFailingSam(sam);
				if (mateString != null){
					insertMateInfo(sam, mateString);
					writeOutFailingSam(sam);
				}
			}
			numBadBarcodes+= numBad;
		}
	}
	

	private void clearAndSplit() {
		firstForward.clear();
		firstReverse.clear();
		secondForward.clear();
		secondReverse.clear();
		
		//split into 4 catagories
		for (SAMRecord sam : workingRecords){
			if (sam.getFirstOfPairFlag()){
				if (sam.getReadNegativeStrandFlag()) firstReverse.add(sam);
				else firstForward.add(sam);
			}
			else {
				if (sam.getReadNegativeStrandFlag()) secondReverse.add(sam);
				else secondForward.add(sam);
			}
		}
	}

	/**This writes out the firstOfPair and it's mate, if it is a second of pair, it will do likewise provided the mate is missing or is unmapped.*/
	private void processSolo( SAMRecord sam) throws Exception {
		//is it firstOfPair
		if (sam.getFirstOfPairFlag()) {
			writeOutSingle(sam);
			//increment family size histogram
			famSizeHist.count(1);
		}
		//nope secondOfPair
		else {
			String mateString = fetchMate(sam);
			//is it null or is the mate unmapped? if not then write this out, otherwise do nothing since it will be processed when hitting the first of pair
			if (mateString == null || fetchBitFlagObj(mateString).isUnmapped()) {
				writeOutSingle(sam);
				//increment family size histogram
				famSizeHist.count(1);
			}
		}
	}
	
	/**This writes out the sam and its mate, if present, to passing.*/
	private void writeOutSingle(SAMRecord sam) throws Exception{
		String mate = restoreSam(sam);
		writeOutPassingSam(sam);
		
		//does it have a mate? if so write it out too.
		if (mate != null) {
			insertMateInfo(sam, mate);
			writeOutPassingSam(sam);
		}
	}
	
	/**This uses the sam as a shell and inserts info from the mate into it making it the mate SAMRecord.*/
	private void insertMateInfo(SAMRecord sam, String mateString) throws Exception {
		//<FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL>  [<TAG>:<VTYPE>:<VALUE> [...]]
		//  0       1      2     3       4       5     6       7      8     9       10 +
		String[] fields = Misc.WHITESPACE.split(mateString);
		//bitflag
		sam.setFlags(Integer.parseInt(fields[0]));
		//chrom name
		sam.setReferenceName(fields[1]);
		//position
		sam.setAlignmentStart(Integer.parseInt(fields[2]));
		//mapping qual
		sam.setMappingQuality(Integer.parseInt(fields[3]));
		//cigar
		sam.setCigarString(fields[4]);
		//mate reference
		sam.setMateReferenceName(fields[5]);
		//set mate position
		sam.setMateAlignmentStart(Integer.parseInt(fields[6]));
		//insert size
		sam.setInferredInsertSize(Integer.parseInt(fields[7]));
		//sequence
		sam.setReadString(fields[8]);
		//qualities
		sam.setBaseQualityString(fields[9]);
		//AS attribute
		sam.setAttribute("AS", fetchAS(mateString));
	}
	
	private Integer fetchAS(String samString) throws Exception{
		Matcher mat = alignmentScore.matcher(samString);
		if (mat.find())return Integer.parseInt(mat.group(1));
		else throw new Exception ("\nFailed to find the AS:i:xxx from "+samString);
	}

	/**This removes any MT attribute and converts the position from the unclipped value from MatchMates back to the original position
	 * @return the truncated mate if present otherwise null.*/
	private String restoreSam(SAMRecord sam){
		//reset position
		int numClippedBases = sam.getAlignmentStart() - sam.getUnclippedStart();
		if (numClippedBases != 0) sam.setAlignmentStart(sam.getAlignmentStart() + numClippedBases);
				
		Object mate = sam.getAttribute("MT");
		//remove MT 
		if (mate !=null) {
			sam.setAttribute("MT", null);
			return (String)mate;
		}
		return null;
	}
	
	/**Fetches the truncated mate if present, otherwise null, makes no changes to the SAMRecord.*/
	public String fetchMate(SAMRecord sam){
		Object mate = sam.getAttribute("MT");
		if (mate != null) return (String) mate;
		return null;
	}
	
	public String toString(){
		return threadId+" "+chromosome+":"+start+"-"+stop;
	}

	public boolean isWorkFailed() {
		return workFailed;
	}
	
	public void writeOutPassingSam(SAMRecord sam) throws IOException{
		passingSamWriter.print(sam.getSAMString());
		numPassingSamOut++;
		totalPassingBps += sam.getReadBases().length;
	}
	public void writeOutFailingSam(SAMRecord sam) throws IOException{
		failingSamWriter.print(sam.getSAMString());
		numFailingSamOut++;
	}
	public void writeOutPairedEndFastq(String[] first, String[] second) throws IOException{
		fastqReadOneWriter.println(first);
		fastqReadTwoWriter.println(second);
		numFastqPairsOut++;
		totalPassingBps += first.length;
		totalPassingBps += second.length;
	}
	public void writeOutSingleEndFastq(String[] single) throws IOException{
		fastqUnpairedWriter.println(single);
		numFastqUnpairedOut++;
		totalPassingBps += single.length;
	}

	public File getPassingSamFile() {
		return passingSamFile;
	}

	public File getFailingSamFile() {
		return failingSamFile;
	}

	public File getFastqFileReadOne() {
		return fastqFileReadOne;
	}

	public File getFastqFileReadTwo() {
		return fastqFileReadTwo;
	}

	public File getFastqFileUnpaired() {
		return fastqFileUnpaired;
	}

	public ConsensusEngine getConsensusEngine() {
		return consensusEngine;
	}

}

