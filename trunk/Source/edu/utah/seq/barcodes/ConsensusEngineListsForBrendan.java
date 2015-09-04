package edu.utah.seq.barcodes;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**This calls consensus on List<SAMRecord> sams returning a modifed sam record for fastq extraction and realignment.
 * No checks are made on the input sams.  They are assumed to have the same strand and read first/ second order.
 * Bases with quality less than the minimumBaseQuality are N'ed and ignored.  
 * For each base index, the number of GATC's are tabulated. If a particular GATC base is >= minimumIdentity, it is taken
 * as the consensus and the maximum quality observed for that base assigned to represent it. If no base reaches 
 * the minimumIdentity, an N is assigned and the quality is set to 0.
 * 
 * @author David Nix*/
public class ConsensusEngineListsForBrendan {

	static final String FAMILY_SIZE_TAG = "FZ";
	
	//fields
	private int minimumBaseQuality = 20;
	private double minimumIdentity = 0.66;
	private char[][] seq;
	private byte[][] qual;
	
	//counters
	private double gCount = 0;
	private double aCount = 0;
	private double tCount = 0;
	private double cCount = 0;
	
	//max base quality observed
	private int gQual = 0;
	private int aQual = 0;
	private int tQual = 0;
	private int cQual = 0;
	
	//final merged seq and qual
	char[] consensusSeq = null;
	int[] consensusQual = null;
	String[] fastqRecord = new String[4];
	
	//constructors
	/**The minimumIdentity is assumed to be > 0.5, setting otherwise will lead to unexpected first in skewed behavior.*/
	public ConsensusEngineListsForBrendan(int minimumBaseQuality, double minimumIdentity){
		this.minimumBaseQuality = minimumBaseQuality;
		this.minimumIdentity = minimumIdentity;
		if (minimumIdentity <= 0.5) throw new IllegalArgumentException("ERROR: the minimumIdentity must be > 0.5");
	}
	public ConsensusEngineListsForBrendan(){}

	
	
	//methods
	/**Call this repeatedly to generate a consensus sequence and max observed quality for the majority base.  
	 * The first record with the modified qual and seq is written to output.
	*/
	public void callConsensus(List<SAMRecord> sams, SAMFileWriter output) throws IOException {
		
		loadSeqQualArrays(sams);
		nSeqQualArray();
		buildConsensus();
		
		//create consensus
		SAMRecord consensus = sams.get(0);
		consensus.setAttribute(FAMILY_SIZE_TAG, sams.size());
		addConsensusSeqQualTo(consensus);
		
		//would be good to add the barcode BC and QT tags in sam spec

		System.out.println(consensus.getSAMString().trim());
		//save it
		output.addAlignment(consensus);
	}


	private void addConsensusSeqQualTo(SAMRecord sam) {
		String seq = new String(consensusSeq);
		String qual = new String(convertSangerQualityScoresToAscii(consensusQual));
		sam.setReadString(seq);
		sam.setBaseQualityString(qual);
	}

	private void buildConsensus() {
		consensusSeq = new char[seq[0].length];
		consensusQual = new int[seq[0].length];
		
		//for each base
		for (int i=0; i< seq[0].length; i++){
			//clear counters and max base qual
			gCount= aCount= tCount= cCount= 0;
			gQual= aQual= tQual= cQual= 0;
			
			//increment counters for all reads at that position
			for (int j=0; j< seq.length; j++){
				if (seq[j][i] == 'G') {
					gCount++;
					if (qual[j][i] > gQual) gQual = qual[j][i];
				}
				else if (seq[j][i] == 'A') {
					aCount++;
					if (qual[j][i] > aQual) aQual = qual[j][i];
				}
				else if (seq[j][i] == 'T') {
					tCount++;
					if (qual[j][i] > tQual) tQual = qual[j][i];
				}
				else if (seq[j][i] == 'C') {
					cCount++;
					if (qual[j][i] > cQual) cQual = qual[j][i];
				}
			}
			//for this base index
			callConsensus(i);
		}
	}

	private void callConsensus(int i) {
		double total = gCount+ aCount+ tCount+ cCount;
		double frac = gCount/total;
		if (frac >= minimumIdentity) {
			consensusSeq[i] = 'G';
			consensusQual[i] = gQual;
			return;
		}
		frac = aCount/total;
		if (frac >= minimumIdentity) {
			consensusSeq[i] = 'A';
			consensusQual[i] = aQual;
			return;
		}
		frac = tCount/total;
		if (frac >= minimumIdentity) {
			consensusSeq[i] = 'T';
			consensusQual[i] = tQual;
			return;
		}
		frac = cCount/total;
		if (frac >= minimumIdentity) {
			consensusSeq[i] = 'C';
			consensusQual[i] = cQual;
			return;
		}
		consensusSeq[i] = 'N';
		consensusQual[i] = 0;
	}

	private void nSeqQualArray() {
		//for each sequence
		for (int i=0; i< seq.length; i++){
			//for each base
			for (int j=0; j< seq[i].length; j++){
				if (qual[i][j] < minimumBaseQuality) seq[i][j] = 'N';
			}
		}
	}

	private void loadSeqQualArrays(List<SAMRecord> sams) throws IOException {
		//get sizes
		int numRecords = sams.size();
		//read length
		int seqLength = sams.get(0).getReadLength();
		
		seq = new char[numRecords][seqLength];
		qual = new byte[numRecords][seqLength];
		for (int i=0; i< numRecords; i++){
			SAMRecord sam = sams.get(i);
			char[] s = sam.getReadString().toCharArray();
			if (s.length != seqLength) throw new IOException("Read seq lengths differ, all must be "+seqLength+", see "+sam.getSAMString().trim());
			seq[i] = s;
			byte[] q = sam.getBaseQualities();
			
			if (q.length != seqLength) throw new IOException("Read qual lengths differ, all must be "+seqLength+", see "+sam.getSAMString().trim());
			qual[i] = q;
		}
	}
	
	//from Seq class in USeq
	public static final char[] ORDERED_ASCII_CHAR = new char[]{'!', '\'', '#', '$', '%', '&', '\'', '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', 
		'6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 
		'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_', '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 
		'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '{', '|', '}', '~'};	
	public static char[] convertSangerQualityScoresToAscii(int[] qual) throws IndexOutOfBoundsException {
		char[] q = new char[qual.length];
		for (int i=0; i< qual.length; i++) q[i] = ORDERED_ASCII_CHAR[qual[i]];
		return q;
	}
	/**Takes a DNA seq and reverse comps it, only GATCN.*/	
	public static String reverseComplementDNA(String seq){
		int seqLen = seq.length();
		StringBuffer rcSeq = new StringBuffer(seqLen);
		char test;
		for (int i=seqLen-1; i>=0; i--){
			test = seq.charAt(i);
			switch (test){
				case 'A': rcSeq.append('T'); break;
				case 'C': rcSeq.append('G'); break;
				case 'G': rcSeq.append('C'); break;
				case 'T': rcSeq.append('A'); break;
				case 'N': rcSeq.append('N'); break;
				default: rcSeq.append(test); System.err.println("\nWarning: odd base in revComp-> '"+test+
				"' from "+seq+" Reverse Complement possibly incorrect!\n");
			}
		}
		return rcSeq.toString();
	}

	public static void main(String[] args) {
		try {
			//load some clustered alignments
			SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			SamReader bamReader = readerFactory.open(new File ("/Users/u0028003/Desktop/Consensus/test.sam"));
			SAMRecordIterator it = bamReader.iterator();
			List<SAMRecord> sams = new ArrayList<SAMRecord>();
			while (it.hasNext()) sams.add(it.next());
			bamReader.close();
			
			//make a reusable consensus engine
			ConsensusEngineListsForBrendan ce = new ConsensusEngineListsForBrendan();
			
			//call consensus returning a single fastq for realignment
			ce.callConsensus(sams, null);
			
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}



	}




	

}
