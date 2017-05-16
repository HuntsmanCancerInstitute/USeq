package edu.utah.seq.barcodes;

import java.io.*;
import java.util.*;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.*;
import util.bio.annotation.Coordinate;
import util.gen.*;

/** This is a stand alone thread for loading alignments into a chunk, tracking all of the counters, and coordinating the writing of bam and fastq data.
 * @author Nix
 * */
public class BarcodeLoader extends Thread{

	//user defined fields
	private File bamFile;
	private File saveDirectory;

	//for chunking
	private int numRecordsInChunk;
	private int maxAlignmentsToCluster;

	//internal 
	private Consensus consensus;
	private SamReader samReader;
	protected SamReaderFactory readerfactory;
	private Gzipper failingSamWriter;
	private String samHeader;
	private ArrayList<Coordinate> chunks = new ArrayList<Coordinate>();
	private boolean completedInitialization = false;
	private HashMap<String, Integer> chromLength = new HashMap<String, Integer>();
	private File failingSamFile;
	private File passingSamFile;

	//stats collected from workers
	private long numRecords = 0;
	private long numBadBarcodes = 0;
	private long numFastqPairsOut = 0;
	private long numFastqUnpairedOut = 0;
	private long numPassingSamOut = 0;
	private long numFailingSamOut = 0;
	private long numModifiedBps = 0;
	private long totalPassingBps = 0;
	private Histogram famSizeHist = new Histogram(0, 100, 100);


	//constructor
	public BarcodeLoader(Consensus consensus) throws IOException{
		//set fields
		bamFile = consensus.getBamFile();
		saveDirectory = consensus.getSaveDirectory();
		numRecordsInChunk = consensus.getNumRecordsInChunk();
		maxAlignmentsToCluster = consensus.getMaxAlignmentsToCluster();
		this.consensus = consensus;
		
		makeFactories();
		
		start();
	}

	public void run(){
		
		System.out.print("Chunking genome...\n\t");
		openReaderAndFetchInfo();
		chunkGenome();
		
		//write out chunks
		Coordinate[] c = new Coordinate[chunks.size()];
		chunks.toArray(c);
		Coordinate.writeToFile(c, new File (saveDirectory,"chunks.bed"));
		
//chunks.clear();
//chunks.add(new Coordinate("21", 42876254, 48119869));
//System.out.println("Proc just chunk on 21");

		completedInitialization = true;
	}

	/**This shuts down the IO and prints the run stats.*/
	public void completeTasks() throws IOException {
		closeIO();
		printStats();
	}

	public void printStats() {
		//stats
		System.out.println("\nProcessing statistics:");
		System.out.println("Total MM records   : "+numRecords);
		System.out.println("Paired fastq       : "+numFastqPairsOut);
		System.out.println("Unpaired fastq     : "+numFastqUnpairedOut);
		System.out.println("Passing bam records: "+numPassingSamOut);
		System.out.println("Failing bam records: "+numFailingSamOut);
		System.out.println("Bad barcodes       : "+numBadBarcodes);
		System.out.println("Modified BPs       : "+numModifiedBps);
		System.out.println("Total BPs          : "+totalPassingBps);
		//histo
		System.out.println("\nFamily size histogram:");
		famSizeHist.printScaledHistogram();

	}
	
	public void loadChromLengths(){
		//fetch chroms to scan from bam header
		List<SAMSequenceRecord> chrList = samReader.getFileHeader().getSequenceDictionary().getSequences();
		//for each
		Iterator<SAMSequenceRecord> it = chrList.iterator();
		while (it.hasNext()){
			SAMSequenceRecord ssr = it.next();
			String chromName = ssr.getSequenceName();
			Integer chromLen = ssr.getSequenceLength() + 1000;
			this.chromLength.put(chromName, chromLen);
		}
	}

	public void openReaderAndFetchInfo(){
		try {
			//make reader and look for index
			samReader = readerfactory.open(bamFile);
			if ( samReader.hasIndex() == false) {
				samReader.close();
				Misc.printErrAndExit("\nError: cannot find the index for your bam file? "+bamFile);
			}
			
			//set text header for thread writers
			samHeader = samReader.getFileHeader().getTextHeader().trim();
			
			//create a writer for the failing records
			failingSamFile = new File (consensus.getSaveDirectory(), "0_failing.sam.gz");
			failingSamWriter = new Gzipper(failingSamFile);
			failingSamWriter.println(samHeader);
			failingSamFile.deleteOnExit();
			
			//write out a 0_passing.sam.gz for the header
			passingSamFile = new File (consensus.getSaveDirectory(), "0_passing.sam.gz");
			Gzipper passing = new Gzipper (passingSamFile);
			passing.println(samHeader);
			passing.close();
			passingSamFile.deleteOnExit();

			//leave samReader open!

		} catch (Exception e) {
			System.err.println("\nError loading chromsomes from "+bamFile.getName());
			e.printStackTrace();
			System.exit(0);
		} 
	}

	private void chunkGenome() {
		//get chromosomes
		List<SAMSequenceRecord> chrs = samReader.getFileHeader().getSequenceDictionary().getSequences();

		//for each chromosome
		for (SAMSequenceRecord r : chrs) {
			//create entry for the hash hash
			String chromName = r.getSequenceName();
			int length = r.getSequenceLength();
			
			//create counter for chunks
			int[] alignmentCounts = new int[length+150];
			
			//iterate through all records for this chromosome
			SAMRecordIterator it = samReader.queryOverlapping(chromName, 0, 0);

			while (it.hasNext()){
				SAMRecord sam = it.next();
//if (sam.getReadName().equals("NS500690:33:H5W2KBGXX:2:21212:18825:8799:BMF:GCAGTGCTTGCCGAGT@CCCCFFFFFFFFFFF")) {
	//System.out.println("Found it in discordant\n"+sam.getAlignmentStart()+" "+sam.getUnclippedStart());
//}
				//increment counter only if under max
				//remember sam uses 1 based numbering so need to subtract one to get to interbase
				//lastly, remember the position in the sam records from MatchMates is the unclipped start!
				int unclippedStart = sam.getAlignmentStart()-1;
				if (unclippedStart < 0) unclippedStart = 0;
				if (alignmentCounts[unclippedStart] < maxAlignmentsToCluster) alignmentCounts[unclippedStart]++;
			}
			it.close();
			
			chunkCounts(chromName, alignmentCounts);
		}
		System.out.println();
	}

	/**This walks the counts and creates chunks with numRecordsInChunk */
	private void chunkCounts(String chromName, int[] alignmentCounts) {
		int counter = 0;
		int start = -1;
		int lastNonZero = -1;
		for (int i=0; i< alignmentCounts.length; i++){
			//any counts?
			int currNum = alignmentCounts[i];
			if (currNum == 0) continue;
			lastNonZero = i;
			
			//first record?
			if (start == -1) start = i;
			
			//increment counter
			counter+= currNum;
			
			//make a new one? interbase coor so last base is excluded
			if (counter > numRecordsInChunk){
				chunks.add(new Coordinate(chromName, start, lastNonZero));
				//reset
				counter = currNum;
				start = lastNonZero;
			}
		}
		//make last?
		if (lastNonZero != -1) chunks.add(new Coordinate(chromName, start, lastNonZero+1));
	}

	public synchronized BarcodeChunk loadNextChunk() throws Exception {
		//any more chunks?
		if (chunks.size() == 0) return null;
		Coordinate coor = chunks.remove(0);
		
		//set fields
		String chromosome = coor.getChromosome();
		int start = coor.getStart();
		int stop = coor.getStop();
		
		SAMRecordIterator iterator = samReader.queryOverlapping(chromosome, start-1, stop+1);

		//any records to load?
		if (iterator.hasNext() == false) {
			iterator.close();
			return null;
		}

		//make array and hash to hold them
		ArrayList<SAMRecord>[] orderedRecords = new ArrayList[stop-start];

		//load all the records for this chunk
		while (iterator.hasNext()){
			SAMRecord sam = iterator.next();
			
			//remember that alignments from MatchMates have there position set to the unclipped start
			int unclippedStart = sam.getAlignmentStart()-1;
			
			//in range?
			if (unclippedStart< start || unclippedStart >= stop) continue;
			numRecords++;
			
			//get any existing ArrayList containing records
			int si = unclippedStart- start;
			ArrayList<SAMRecord> positionRecords = orderedRecords[si];
			if (positionRecords == null) orderedRecords[si] = new ArrayList<SAMRecord>();

			//add record? must watch out for spiking positions where alignment counts can exceed 100K
			if (orderedRecords[si].size() < maxAlignmentsToCluster) orderedRecords[si].add(sam);

			//would be nice to skip this for speed....
			else writeOutFailingSam(sam);
		}
		iterator.close();
		return new BarcodeChunk(chromosome, start, stop, orderedRecords);
	}

	private void makeFactories() throws IOException {
		readerfactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	}

	public void closeIO() throws IOException{
		samReader.close();
		failingSamWriter.close();
	}

	public SamReaderFactory getReaderfactory() {
		return readerfactory;
	}

	public synchronized void incrementFamSizeHist(Histogram h) throws Exception{
		famSizeHist.addCounts(h);
	}
	public synchronized void incrementNumBadBarcodes(long n) {
		numBadBarcodes += n;
	}
	public synchronized void incrementNumFastqPairsOut(long n) {
		numFastqPairsOut += n;
	}
	public synchronized void incrementNumModifiedBps(long n) {
		numModifiedBps += n;
	}
	public synchronized void incrementTotalPassingBps(long n) {
		totalPassingBps += n;
	}
	public synchronized void incrementNumFastqUnpairedOut(long n) {
		numFastqUnpairedOut += n;
	}
	public synchronized void incrementNumPassingSamOut(long n) {
		numPassingSamOut += n;
	}
	public synchronized void incrementNumFailingSamOut(long n) {
		numFailingSamOut += n;
	}
	public boolean isCompletedInitialization() {
		return completedInitialization;
	}

	public ArrayList<Coordinate> getChunks() {
		return chunks;
	}

	public String getSamHeader() {
		return samHeader;
	}
	public void writeOutFailingSam(SAMRecord sam) throws IOException{
		failingSamWriter.print(sam.getSAMString());
		numFailingSamOut++;
	}

	public File getFailingSamFile() {
		return failingSamFile;
	}

	public File getPassingSamFile() {
		return passingSamFile;
	}

	
}
