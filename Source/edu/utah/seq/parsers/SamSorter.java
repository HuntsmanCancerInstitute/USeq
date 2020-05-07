package edu.utah.seq.parsers;

import java.io.File;
import edu.utah.seq.data.sam.PicardSortSam;

public class SamSorter implements Runnable {

	private boolean failed = false;
	private File samToSort = null;
	private File bamOutput = null;

	public SamSorter (File samToSort, File bamOutput){
		this.samToSort = samToSort;
		this.bamOutput = bamOutput;
	}

	public void run() {	
		try {
			new PicardSortSam (samToSort, bamOutput, true);
		} catch (Exception e) {
			failed = true;
			bamOutput.delete();
			System.err.println("\nError: problem sorting "+samToSort );
			e.printStackTrace();
		}
	}

	public boolean isFailed() {
		return failed;
	}
}

