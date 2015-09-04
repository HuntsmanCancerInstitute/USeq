package edu.utah.seq.barcodes;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.Random;
import java.util.TreeMap;

import util.gen.Gzipper;
import util.gen.Histogram;
import util.gen.Misc;

public class OldBarcodeChromMerger implements Runnable {

	private Consensus consensus;
	private String chromosome;
	private BarcodeClusterEngine clusterEngine;
	private ConsensusEngine consensusEngine;
	private SAMFileWriter passingBamWriter;
	private SAMFileWriter failingBamWriter;
	private Gzipper fastqReadOneWriter;
	private Gzipper fastqReadTwoWriter;
	private Gzipper fastqUnpairedWriter;
	private SAMRecordIterator iterator;
	private SamReader bamReaderForLoading;
	private SamReader bamReaderForMates;
	private TreeMap<Integer, ArrayList<SAMRecord>> orderedRecords = new TreeMap<Integer, ArrayList<SAMRecord>>();
	private LinkedHashMap<String, SAMRecord[]> nameMatePair = new LinkedHashMap<String, SAMRecord[]>();
	private long numRecords = 0;
	private long numFastqPairsOut = 0;
	private long numFastqUnpairedOut = 0;
	private long numPassingBamOut = 0;
	private long numFailingBamOut = 0;
	private long numPropPairMatesNotFound = 0;
	private int numBasesSubsampled = 0;
	private int numPosInMem = 1000;
	private int maxNumNameMatePair = 50000;
	private int maxAlignmentsToCluster;
	private Random random = new Random(0);
	
	private File passingBamFile;
	private File failingBamFile;
	private File fastqFileReadOne;
	private File fastqFileReadTwo;
	private File fastqFileUnpaired;
	private Histogram famSizeHist = new Histogram(0, 100, 100);
	private boolean workFailed = false;
	
	//working objects
	private ArrayList<SAMRecord> workingRecords = null;
	private ArrayList<SAMRecord> workingMates = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> firstForward = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> firstReverse = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> secondForward = new ArrayList<SAMRecord>();
	private ArrayList<SAMRecord> secondReverse = new ArrayList<SAMRecord>();

	
	public OldBarcodeChromMerger(String chromosome, Consensus consensus) throws Exception {
		this.consensus = consensus;
		this.chromosome = chromosome;
		maxAlignmentsToCluster = consensus.getMaxAlignmentsToCluster();
		
		//make engines
		clusterEngine = new BarcodeClusterEngine(consensus.getMinBarcodeBaseQuality(), 
				consensus.getMinNumPassingBarcodeBases(), consensus.getMinFractionBarcodeIdentity());
		consensusEngine = new ConsensusEngine(consensus.getMinReadBaseQuality(), consensus.getMinReadBaseFractionIdentity());
		
		makeReadersAndWriters();	
	}
	
	public void run() {
		try {
			//load initial set of records to numPosInMem, if false then none found and work is complete
			if (loadRecords() == false) closeAndDeleteIO();
			else {
				processRecords();
				closeIO();
			}
			System.out.print(chromosome+" ");
		} catch (Exception e){
			System.err.println("\nError encountered when processing chromosome "+chromosome+"\n");
			e.printStackTrace();
			workFailed = true;
		}
	}
	


	private void processRecords() throws IOException {
//int counter = 0;
		while (orderedRecords.size() != 0){
			
//long start = System.nanoTime();  

//Integer position = orderedRecords.firstKey();

			//extract first list of records
			workingRecords = orderedRecords.remove(orderedRecords.firstKey());
			
			//check to see if over the limit and subsample
			checkSizeWorkingRecordsAndSubsample();
//workingRecords = orderedRecords.remove(position);

//System.out.println("\nchr"+chromosome+":"+position+"\t"+workingRecords.size());
			
			//process the records found at this particular unclipped start position
			processRecordsAtOnePosition();
			
//System.out.println((System.nanoTime()-start)+" time for processRecordsAtOnePosition()");
//start = System.nanoTime();
			
			//load more records to meet the numPosInMem, at some point there will be no more records for this chrom and the orderedRecords.size will drop to zero
			loadRecords();
//System.out.println((System.nanoTime()-start)+" time for loadRecords()");
//start = System.nanoTime();
			
			
			//if (++counter > 2000) break;
//System.out.println("Rec: "+workingRecords.size()+" FPO: "+numFastqPairsOut+" FUO: "+numFastqUnpairedOut+" BPO: "+numPassingBamOut+" BFO: "+numFailingBamOut);
			//if (position == 9413061){
				//for (SAMRecord s: workingRecords) System.out.println(s.getSAMString().trim());
			//}
//if (position == 9533516) break;

		}
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

	private void processRecordsAtOnePosition() throws IOException {
		//just one record? or a group that might need to be split for clustering
		if (workingRecords.size() == 1) {
			SAMRecord sam = workingRecords.get(0);
			processSolo(sam, sam.getFirstOfPairFlag());
		}
		else {
			processGroup();
		}
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
	
	/* SLOW
	private SAMRecord queryMate(SAMRecord sam){
		SAMRecord mate = null;
		//long start = System.nanoTime();	

		//add check if discordant, if so then fetch from head thread discordant pair.
		//StringBuilder sb = new StringBuilder();
		//sb.append("Mate lookup for\n\t"+sam.getSAMString().trim()+"\n");
		//first query local loaded hash, this might be null, or contain only one of the two pairs, the other being null
		SAMRecord[] pair = nameMatePair.get(sam.getReadName());
		if (pair != null){
			if (sam.getFirstOfPairFlag()) mate = pair[1];
			else mate = pair[0];

			//if (mate != null) sb.append("Pair found\n\t"+mate.getSAMString().trim()+"\n");
			//else sb.append("Pair found but mate was null\n");
		}
		//else sb.append("Pair NOT found\n");
		//still null?
		if (mate == null && sam.getMateUnmappedFlag() == false) {
			//this can take seconds!
			mate = bamReaderForMates.queryMate(sam); 

			//if (mate != null) sb.append("QueryMateSlow found\n\t"+mate.getSAMString().trim()+"\n");
			//else sb.append("QueryMateSlow failed mate null\n");
		}
		//long diff = System.nanoTime()-start;
		//if (diff > 300000) System.out.println("\n************"+diff+" SLOW time to queryMate()\n"+sb);		
		return mate;
	}*/
	
	private SAMRecord queryMate(SAMRecord sam){
		SAMRecord mate = null;
		
		//is it a proper pair? if so then query local hash
		if (sam.getProperPairFlag()){
StringBuilder sb = new StringBuilder();
sb.append("\nSAM: "+sam.getSAMString());
			SAMRecord[] pair = nameMatePair.get(sam.getReadName());
			if (pair != null){
				if (sam.getFirstOfPairFlag()) mate = pair[1];
				else mate = pair[0];
				sb.append("Pair found\n");
			}
			else {
				sb.append("Pair not found\n");
			}
			
			if (mate == null) {
				numPropPairMatesNotFound++;
				sb.append("Mate null\n");
				Misc.printErrAndExit(sb.toString());
			}
			else sb.append("Mate found\n"+mate.getSAMString());
			
			if (sam.getReadName().equals("NS500690:33:H5W2KBGXX:4:12608:14678:16772:BMF:ATAACGAGGGCCTCAG@CCCCFFFFFFFEFFF")) Misc.printErrAndExit(sb.toString());
			
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
		HashMap<String, SAMRecord>[] positionHash = null; //consensus.chromDiscordAlign.get(sam.getMateReferenceName());
		//it might not exist due to having no discordant alignments
		if (positionHash == null) return null;
		
		//fetch HashMap<String, SAMRecord> for a particular base position/ index
		HashMap<String, SAMRecord> recordsHash = positionHash[sam.getMateAlignmentStart()];
		//might not have any records at that position
		if (recordsHash == null) return null;
		
		//make name for the mate and see if in that hash
		String nameOrder;
		if (sam.getFirstOfPairFlag()) nameOrder = sam.getReadName()+ consensus.secondOfPair;
		else nameOrder = sam.getReadName()+ consensus.firstOfPair;
		
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

	private boolean loadRecords() {
		//any records to load?
		if (iterator.hasNext() == false) return false;

		//load records until numPosInMem is reached
		while (iterator.hasNext()){
			SAMRecord sam = iterator.next();
			numRecords++;
			
			//is it a failing alignment?
			if (sam.getReadUnmappedFlag() || sam.isSecondaryOrSupplementary()) writeOutFailingSam(sam);
				
			//nope, OK
			else{
				//get position and any existing ArrayList containing records
				Integer unclippedStart = sam.getUnclippedStart();
				ArrayList<SAMRecord> positionRecords = orderedRecords.get(unclippedStart);
				if (positionRecords == null) {
					positionRecords = new ArrayList<SAMRecord>();
					orderedRecords.put(unclippedStart, positionRecords);
				}
				//add record
				positionRecords.add(sam);
				
				//add record to mate hash
				SAMRecord[] pair = nameMatePair.get(sam.getReadName());
				if (pair == null){
					pair = new SAMRecord[2];
					nameMatePair.put(sam.getReadName(), pair);
				}
				if (sam.getFirstOfPairFlag()) pair[0] = sam;
				else pair[1] = sam;
				
				//reached number of positions required in tree?
				if (orderedRecords.size() >= numPosInMem) break;
			}
		}
		
		//check that nameMatePair hasn't gotten too big
		checkSizeOfNameMatePairHash();
		
		return true;
	}

	private void checkSizeOfNameMatePairHash() {
		//if too big, trims back the size of the hash.
		int diff = nameMatePair.size() - maxNumNameMatePair;
		if (diff > 0) {
			Iterator<Entry<String, SAMRecord[]>> it = nameMatePair.entrySet().iterator();
			for (int i=0; i< diff; i++){
				it.next();
				it.remove();
			}
		}
		
	}

	private void makeReadersAndWriters() throws Exception {
		String name = Misc.removeExtension( consensus.getBamFile().getName() );
		//bam readers
		bamReaderForLoading = consensus.getReaderfactory().open(consensus.getBamFile());
		bamReaderForMates = consensus.getReaderfactory().open(consensus.getBamFile());
		
		//bam writers
		passingBamFile = new File (consensus.getSaveDirectory(), "passingNotMerged_"+chromosome+"_"+name+".bam");
		passingBamWriter = consensus.getWriterFactory().makeBAMWriter(bamReaderForLoading.getFileHeader(), false, passingBamFile);
		failingBamFile = new File (consensus.getSaveDirectory(), "failingNotMerged_"+chromosome+"_"+name+".bam");
		failingBamWriter = consensus.getWriterFactory().makeBAMWriter(bamReaderForLoading.getFileHeader(), false, failingBamFile);
		
		iterator = bamReaderForLoading.queryOverlapping(chromosome, 0, 0);
		
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
		bamReaderForMates.close();
		iterator.close();
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

}
