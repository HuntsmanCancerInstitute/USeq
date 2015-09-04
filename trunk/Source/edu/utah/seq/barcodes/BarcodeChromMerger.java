package edu.utah.seq.barcodes;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import util.gen.Gzipper;
import util.gen.Histogram;
import util.gen.Misc;

public class BarcodeChromMerger implements Runnable {

	private String chromosome;
	private int start;
	private int stop;
	private Consensus consensus;
	private BarcodeClusterEngine clusterEngine;
	private ConsensusEngine consensusEngine;
	private SAMFileWriter passingBamWriter;
	private SAMFileWriter failingBamWriter;
	private Gzipper fastqReadOneWriter;
	private Gzipper fastqReadTwoWriter;
	private Gzipper fastqUnpairedWriter;
	private SamReader bamReaderForLoading;
	private ArrayList<SAMRecord>[] orderedRecords;
	private HashMap<String, SAMRecord[]> nameMatePair;
	private long numRecords = 0;
	private long numFastqPairsOut = 0;
	private long numFastqUnpairedOut = 0;
	private long numPassingBamOut = 0;
	private long numFailingBamOut = 0;
	private long numPropPairMatesNotFound = 0;
	private int numBasesSubsampled = 0;
	private int maxAlignmentsToCluster;
	private Random random = new Random(0);
	
	private File passingBamFile;
	private File failingBamFile;
	private File fastqFileReadOne;
	private File fastqFileReadTwo;
	private File fastqFileUnpaired;
	private Histogram famSizeHist = new Histogram(0, 100, 100);
	private HashMap<String, HashMap<String, SAMRecord>[]> chromDiscordAlign;

	//working objects
	private ArrayList<SAMRecord> workingRecords = null;
	private ArrayList<SAMRecord> workingMates = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> firstForward = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> firstReverse = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> secondForward = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> secondReverse = new ArrayList<SAMRecord>();
	private boolean workFailed = false;
	
	/*Must be very careful to not reference something external that could cause this object to avoid garbage collection after it has finished running.  
	 * Huge pain with memory leaks!*/
	public BarcodeChromMerger(String chromosome, int start, int stop, Consensus consensus) throws Exception {
		this.chromosome = new String (chromosome);
		this.start = start;
		this.stop = stop;
		this.consensus = consensus;
		this.chromDiscordAlign = consensus.getChromDiscordAlign();
		maxAlignmentsToCluster = consensus.getMaxAlignmentsToCluster();
		
		//make engines
		clusterEngine = new BarcodeClusterEngine(consensus.getMinBarcodeBaseQuality(), 
				consensus.getMinNumPassingBarcodeBases(), consensus.getMinFractionBarcodeIdentity());
		consensusEngine = new ConsensusEngine(consensus.getMinReadBaseQuality(), consensus.getMinReadBaseFractionIdentity());
		
		makeReadersAndWriters();	
	}
	
	// TODO work on chunking so that we resolve the edges.  load more that bounds yes process starting at bounds?
	
	public void run() {
		try {
			//load records, if false then none found and work is complete
			if (loadRecords() == false) closeAndDeleteIO();
			else {
				processRecords();
				closeIO();
				incrementCounters();
			}
			System.out.println(chromosome+":"+start+"-"+stop+" complete "+ numRecords);
		} catch (Exception e){
			System.err.println("\nError encountered when processing chromosome "+chromosome+":"+start+"-"+stop+"\n");
			e.printStackTrace();
			workFailed = true;
			
		} 
	}
	

	private void incrementCounters() throws Exception {
		consensus.incrementNumRecords(numRecords);
		consensus.incrementNumFastqPairsOut(numFastqPairsOut);
		consensus.incrementNumFastqUnpairedOut(numFastqUnpairedOut);
		consensus.incrementNumPassingBamOut(numPassingBamOut);
		consensus.incrementNumFailingBamOut(numFailingBamOut);
		consensus.incrementNumBasesSubsampled(numBasesSubsampled);
		consensus.incrementNumPropPairMatesNotFound(numPropPairMatesNotFound);
		consensus.incrementFamSizeHist(famSizeHist);
	}

	private void processRecords() throws IOException {
//int counter = 0;

		for (ArrayList<SAMRecord> al : orderedRecords){
			//any records?
			if (al == null) continue;
			workingRecords = al;
			
			//check to see if over the limit and subsample
			checkSizeWorkingRecordsAndSubsample();
			
			//process the records found at this particular unclipped start position
			processRecordsAtOnePosition();

			//if (++counter > 2000) break;
//System.out.println("Rec: "+workingRecords.size()+" FPO: "+numFastqPairsOut+" FUO: "+numFastqUnpairedOut+" BPO: "+numPassingBamOut+" BFO: "+numFailingBamOut);
			//if (position == 9413061){
				//for (SAMRecord s: workingRecords) System.out.println(s.getSAMString().trim());
			//}
//if (position == 9533516) break;

		}
	}

	private void processRecordsAtOnePosition() throws IOException {
		//just one record? or a group that might need to be split for clustering
		if (workingRecords.size() == 1) {
			SAMRecord sam = workingRecords.get(0);
			processSolo(sam, sam.getFirstOfPairFlag());
		}
		else processGroup();
	}

	private void processGroup() throws IOException {
		//clear and split records into the four types of alignments
		clearAndSplit();
		
		//for each group, are there enough to cluster collapse or are they solo?
		//first forward
		int numRec = firstForward.size();
		if (numRec > 0){
			if (numRec == 1) processSolo (firstForward.get(0), true);
			else clusterCollapseFirstOfPair( firstForward );
		}
		
		//first reverse
		numRec = firstReverse.size();
		if (numRec > 0){
			if (numRec == 1) processSolo (firstReverse.get(0), true);
			else clusterCollapseFirstOfPair( firstReverse );
		}
		
		//second forward
		numRec = secondForward.size();
		if (numRec > 0){
			if (numRec == 1) processSolo (secondForward.get(0), false);
			else clusterCollapseSecondOfPair( secondForward );
		}
		
		//second reverse
		numRec = secondReverse.size();
		if (numRec > 0){
			if (numRec == 1) processSolo (secondReverse.get(0), false);
			else clusterCollapseSecondOfPair( secondReverse );
		}
	}

	private void clusterCollapseSecondOfPair(ArrayList<SAMRecord> records) throws IOException {
		//the only time we want to collapse these are if they have no mates in the file
		//if they do then skip them since they'll get processed with the first of pair
		//soo
		parseNoMateRecordsIntoWorkingMates(records);
		//any remaining?
		int num = workingMates.size();
		if (num != 0){
			if (num == 1) writeOutPassingSam(workingMates.get(0));
			else {
				//fetch SAMRecord[cluster number][SAMRecords in a given cluster]
				SAMRecord[][] clusters = clusterEngine.cluster(records);
				
				//any with bad barcodes?
				processBadBarcodeRecords(clusterEngine.getSkippedRecords(), false);
				
				//collapse clusters
				for (SAMRecord[] c2Col : clusters){
					//just one, thus no cluster? or some to collapse
					if (c2Col.length == 1) writeOutPassingSam(c2Col[0]);
					else {
						//call consensus on these secondOfPairs with no mates in file
						String[] fastq = consensusEngine.callConsensus(c2Col);
						// write out as single end 
						writeOutSingleEndFastq(fastq);
					}
					//increment family size in histogram
					famSizeHist.count(c2Col.length);
				}
			}
		}	
	}

	private void clusterCollapseFirstOfPair(ArrayList<SAMRecord> records) throws IOException {

		//fetch SAMRecord[cluster number][SAMRecords in a given cluster]
		SAMRecord[][] clusters = clusterEngine.cluster(records);		
		
		//any with bad barcodes?
		processBadBarcodeRecords(clusterEngine.getSkippedRecords(), true);

		//collapse clusters
		for (SAMRecord[] c2Col : clusters){
			//just one, thus no cluster? or some to collapse
			if (c2Col.length == 1) {
				processSolo (c2Col[0], true);
			}
			else {
				//call consensus on these firstOfPairs
				String[] fastq = consensusEngine.callConsensus(c2Col);
				//attempt to collect mates of these first of pairs then write out as single end or as paired
				SAMRecord[] mates = collectMates(c2Col);
				if (mates == null) writeOutSingleEndFastq(fastq);
				else writeOutPairedEndFastq(fastq, consensusEngine.callConsensus(mates));
				//increment family size in histogram
				famSizeHist.count(c2Col.length);
			}
		}
	}

	private SAMRecord[] collectMates(SAMRecord[] c2Col) {
		workingMates.clear();
		for (SAMRecord sam: c2Col){
			//any mapped mate?
			SAMRecord mate = queryMate(sam);
			if (mate != null) workingMates.add(mate);
		}
		if (workingMates.size() == 0) return null;
		
		SAMRecord[] mates = new SAMRecord[workingMates.size()];
		workingMates.toArray(mates);
		return mates;
	}
	
	private void parseNoMateRecordsIntoWorkingMates(ArrayList<SAMRecord> secondOfPairs) {
		workingMates.clear();
		for (SAMRecord sam: secondOfPairs){
			SAMRecord mate = queryMate(sam);
			if (mate == null) workingMates.add(sam);
		}
	}

	private void processBadBarcodeRecords(ArrayList<SAMRecord> skippedRecords,  boolean findMatesAndPrint) {
		//write out skipped records and possibly their mates as bam
		for (SAMRecord sam: skippedRecords){
			writeOutFailingSam(sam);
			if (findMatesAndPrint){
				SAMRecord mate = queryMate(sam);
				if (mate != null) writeOutFailingSam(mate);
			}
		}
	}
	
	private SAMRecord queryMate(SAMRecord sam){
		SAMRecord mate = null;
		
		//is it a proper pair? if so then query local hash
		if (sam.getProperPairFlag()){
//StringBuilder sb = new StringBuilder();
//sb.append("\nSAM: "+sam.getSAMString());
			SAMRecord[] pair = nameMatePair.get(sam.getReadName());
			if (pair != null){
				if (sam.getFirstOfPairFlag()) mate = pair[1];
				else mate = pair[0];
				//sb.append("Pair found\n");
			}
			else {
				//sb.append("Pair not found\n");
			}
			
			if (mate == null) {
				numPropPairMatesNotFound++;
				//sb.append("Mate null\n");
				//Misc.printErrAndExit(sb.toString());
			}
			//else sb.append("Mate found\n"+mate.getSAMString());
			
			//if (sam.getReadName().equals("NS500690:33:H5W2KBGXX:4:12608:14678:16772:BMF:ATAACGAGGGCCTCAG@CCCCFFFFFFFEFFF")) Misc.printErrAndExit(sb.toString());
			
			return mate;
		}
		
		//not a proper pair so attempt to fetch from main thread
		else {
			//is mate unmapped?
			if (sam.getMateUnmappedFlag()) return null;
			
			//attempt to fetch from the discordant alignments
			return fetchDiscordantMate(sam);
		}
	}

	private SAMRecord fetchDiscordantMate(SAMRecord sam) {
		//fetch HashMap<String, SAMRecord>[] for the mate chromosome
		HashMap<String, SAMRecord>[] positionHash = chromDiscordAlign.get(sam.getMateReferenceName());
		//it might not exist due to having no discordant alignments
		if (positionHash == null) return null;
		
		//fetch HashMap<String, SAMRecord> for a particular base position/ index
		HashMap<String, SAMRecord> recordsHash = positionHash[sam.getMateAlignmentStart()];
		//might not have any records at that position
		if (recordsHash == null) return null;
		
		//make name for the mate and see if in that hash
		String nameOrder;
		if (sam.getFirstOfPairFlag()) nameOrder = sam.getReadName()+ Consensus.secondOfPair;
		else nameOrder = sam.getReadName()+ Consensus.firstOfPair;
		
		//return it, could be null
		return recordsHash.get(nameOrder);
		
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

	private void processSolo( SAMRecord sam, boolean firstOfPair ) {
		//is it firstOfPair
		if (firstOfPair) {
			//any mapped mate? if not just write it out
			if (sam.getMateUnmappedFlag()) writeOutPassingSam(sam);
			else {
				//attempt to fetch mate and write out
				SAMRecord mate = queryMate(sam);
				if (mate != null) writeOutPassingSam(mate);
			}
		}
		//nope second
		else {
			//is there a mate? if not then write this out, 
			//otherwise do nothing since it will be processed when hitting the first of pair
			SAMRecord mate = queryMate(sam);
			if (mate == null) writeOutPassingSam(sam);
		}
		famSizeHist.count(1);
	}

	private boolean loadRecords() throws Exception {
		SAMRecordIterator iterator = bamReaderForLoading.queryOverlapping(chromosome, start, stop);
		
		//any records to load?
		if (iterator.hasNext() == false) return false;
		
		//make array and hash to hold them
		orderedRecords = new ArrayList[(stop - start)];
		nameMatePair = new HashMap<String, SAMRecord[]>();

		//load all the records for this chunk

		while (iterator.hasNext()){
			SAMRecord sam = iterator.next();
			numRecords++;
			
			//is it a failing alignment?
			if (sam.getReadUnmappedFlag() || sam.isSecondaryOrSupplementary()) writeOutFailingSam(sam);
				
			//nope, OK
			else{
				//check start
				int si = sam.getUnclippedStart()- start;
				if (si < 0 || si >= orderedRecords.length) continue;
				
				//get position and any existing ArrayList containing records
				ArrayList<SAMRecord> positionRecords = orderedRecords[si];
				if (positionRecords == null) orderedRecords[si] = new ArrayList<SAMRecord>();

				//add record
				orderedRecords[si].add(sam);
				
				//add record to mate hash
				SAMRecord[] pair = nameMatePair.get(sam.getReadName());
				if (pair == null){
					pair = new SAMRecord[2];
					nameMatePair.put(sam.getReadName(), pair);
				}
				if (sam.getFirstOfPairFlag()) pair[0] = sam;
				else pair[1] = sam;
			}
		}
		iterator.close();
		return true;
	}

	private void makeReadersAndWriters() throws Exception {
		String name = Misc.removeExtension( consensus.getBamFile().getName() );

		//bam readers
		bamReaderForLoading = consensus.getReaderfactory().open(consensus.getBamFile());

		//bam writers
		passingBamFile = new File (consensus.getSaveDirectory(), "passingNotMerged_"+chromosome+"_"+name+".bam");
		passingBamWriter = consensus.getWriterFactory().makeBAMWriter(bamReaderForLoading.getFileHeader(), false, passingBamFile);
		failingBamFile = new File (consensus.getSaveDirectory(), "failingNotMerged_"+chromosome+"_"+name+".bam");
		failingBamWriter = consensus.getWriterFactory().makeBAMWriter(bamReaderForLoading.getFileHeader(), false, failingBamFile);

		//fastq
		fastqFileReadOne = new File (consensus.getSaveDirectory(), "merged_"+ chromosome+ "_"+name+"_1.fastq.gz");
		fastqFileReadTwo = new File (consensus.getSaveDirectory(), "merged_"+ chromosome+ "_"+name+"_2.fastq.gz");
		fastqFileUnpaired = new File (consensus.getSaveDirectory(), "merged_"+ chromosome+ "_"+name+"_unpaired.fastq.gz");
		fastqReadOneWriter = new Gzipper(fastqFileReadOne);
		fastqReadTwoWriter = new Gzipper(fastqFileReadTwo);
		fastqUnpairedWriter = new Gzipper(fastqFileUnpaired);
	}
	
	private void closeIO() throws IOException{
		bamReaderForLoading.close();
		passingBamWriter.close();
		failingBamWriter.close();
		fastqReadOneWriter.close();
		fastqReadTwoWriter.close();
		fastqUnpairedWriter.close();
	}
	
	private void closeAndDeleteIO() throws IOException {
		closeIO();
		passingBamFile.deleteOnExit();
		failingBamFile.deleteOnExit();
		fastqFileReadOne.deleteOnExit();
		fastqFileReadTwo.deleteOnExit();
		fastqFileUnpaired.deleteOnExit();
	}

	private void checkSizeWorkingRecordsAndSubsample() {
		if (workingRecords.size() > maxAlignmentsToCluster){
			System.out.println("\nSubsampling chr"+chromosome+":"+workingRecords.get(0).getAlignmentStart()+"  # rec: "+workingRecords.size());
			while (workingRecords.size() > maxAlignmentsToCluster){
				int indexToDelete = random.nextInt(workingRecords.size());
				writeOutFailingSam( workingRecords.remove(indexToDelete) );
			}
			numBasesSubsampled++;
		}
	}
	
	private void writeOutFailingSam(SAMRecord sam){
		failingBamWriter.addAlignment(sam);
		numFailingBamOut++;
	}
	
	private void writeOutPassingSam(SAMRecord sam) {
		passingBamWriter.addAlignment(sam);
		numPassingBamOut++;
	}
	
	private void writeOutPairedEndFastq(String[] fastqFirst, String[] fastqSecond) throws IOException {
		fastqReadOneWriter.println(fastqFirst);
		fastqReadTwoWriter.println(fastqSecond);
		numFastqPairsOut++;
	}

	private void writeOutSingleEndFastq(String[] fastq) throws IOException {
		fastqUnpairedWriter.println(fastq);
		numFastqUnpairedOut++;
	}
	
	public String toString(){
		return chromosome+":"+start+"-"+stop;
	}
	
	public String getChromosome() {
		return chromosome;
	}
	public boolean getWorkFailed() {
		return workFailed;
	}

	public long getNumRecords() {
		return numRecords;
	}

	public long getNumFastqPairsOut() {
		return numFastqPairsOut;
	}

	public long getNumFastqUnpairedOut() {
		return numFastqUnpairedOut;
	}

	public long getNumPassingBamOut() {
		return numPassingBamOut;
	}

	public long getNumFailingBamOut() {
		return numFailingBamOut;
	}

	public Histogram getFamSizeHist() {
		return famSizeHist;
	}

	public long getNumPropPairMatesNotFound() {
		return numPropPairMatesNotFound;
	}

	public int getNumBasesSubsampled() {
		return numBasesSubsampled;
	}

	/*Don't need this?!*/
	private void cleanup() {
		
		
		for (ArrayList<SAMRecord> al : orderedRecords){
			if (al != null) for (SAMRecord s : al) s = null;
		}
		orderedRecords = null;
		for (String s: nameMatePair.keySet()){
			SAMRecord[] sams = nameMatePair.get(s);
			for (SAMRecord toDel : sams) toDel = null;
			s = null;
		}
		
		consensus = null;
		chromDiscordAlign = null;
		chromosome = null;
		clusterEngine = null;
		consensusEngine = null;
		orderedRecords = null;
		nameMatePair = null;
		passingBamWriter = null;
		failingBamWriter = null;
		fastqReadOneWriter = null;
		fastqReadTwoWriter = null;
		fastqUnpairedWriter = null;
		bamReaderForLoading = null;

		random = null;
		passingBamFile = null;
		failingBamFile = null;
		fastqFileReadOne = null;
		fastqFileReadTwo = null;
		fastqFileUnpaired = null;
		famSizeHist = null;
		workingRecords = null;
		workingMates = null;
		firstForward = null;
		firstReverse = null;
		secondForward = null;
		secondReverse = null;
		System.gc();
	}
}

