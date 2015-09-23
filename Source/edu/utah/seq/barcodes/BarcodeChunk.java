package edu.utah.seq.barcodes;

import htsjdk.samtools.SAMRecord;
import java.util.ArrayList;
import java.util.HashMap;

/**Just a container for a chunk of sequence to process with a BarcodeChromMerger*/
public class BarcodeChunk {
	private String chromosome;
	private int start;
	private int stop;
	private ArrayList<SAMRecord>[] orderedRecords = null;
	
	public BarcodeChunk(String chromosome, int start, int stop, ArrayList<SAMRecord>[] orderedRecords){
		this.chromosome = chromosome;
		this.start = start;
		this.stop = stop;
		this.orderedRecords = orderedRecords;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	public ArrayList<SAMRecord>[] getOrderedRecords() {
		return orderedRecords;
	}
}
