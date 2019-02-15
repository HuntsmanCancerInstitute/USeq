package edu.utah.seq.data;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import edu.utah.seq.data.sam.SamAlignment;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class SamAlignmentDepthLoader implements Runnable{
	
	private Bed[] regions;
	private Gzipper samOut;
	private File sam;
	private boolean failed = false;
	private String errorMessage = null;
	private SamReader toMatchReader;
	private SamReader toSubReader;
	private HashMap<String, ArrayList<SAMRecord>> nameReads = new HashMap<String, ArrayList<SAMRecord>>();
	private Random rng = new Random();
	private int loaderName;
	private double numPassing = 0;
	private double minFrac = 0.95;

	public SamAlignmentDepthLoader(int loaderName, Bed[] regions, SamAlignmentDepthMatcher srdm) throws FileNotFoundException, IOException{
		sam = new File(srdm.getTempDir(), loaderName+".sam.gz");
		sam.deleteOnExit();
		samOut = new Gzipper(sam);
		toMatchReader = srdm.getReaderFactory().open(srdm.getToMatchBam());
		toSubReader = srdm.getReaderFactory().open(srdm.getToSubSampleBam());
		this.loaderName = loaderName;
		this.regions = regions;
	}
	
	
	public void run() {
		try {
			//for each region
			int counter = 0;
			for (int i=0; i< regions.length; i++){
				if (++counter > 500) {
					IO.p(".");
					counter = 0;
				}
				//fetch count of alignments to match, number of alignments including first and second reads
				Bed bed = regions[i];
				double alignmentCountToMatch = fetchCount(bed);
				//if zero then no need to print anything
				if (alignmentCountToMatch == 0)  numPassing++;
				else {
					ArrayList<SAMRecord>[] rdmGrps = fetchReadGroups(bed);
					double numPrinted = printSetNumber(alignmentCountToMatch, rdmGrps);
					double frac = numPrinted/ alignmentCountToMatch;
					if (frac >= minFrac) numPassing++;
				}
			}
			
			//close IO
			samOut.close();
			toMatchReader.close();
			toSubReader.close();
		} catch (IOException e) {
			failed = true;
			errorMessage = e.toString();
			e.printStackTrace();
		}
		
	}

	private double printSetNumber(double alignmentCountToMatch, ArrayList<SAMRecord>[] rdmGrps) throws IOException {
		double numPrinted = 0;
		for (ArrayList<SAMRecord> al: rdmGrps){
			for (SAMRecord sam: al){
				samOut.print(sam.getSAMString());
				numPrinted++;
				if (numPrinted == alignmentCountToMatch) return numPrinted;
			}
		}
		return numPrinted;
	}


	private ArrayList<SAMRecord>[] fetchReadGroups(Bed bed) {
		SAMRecordIterator it = toSubReader.queryOverlapping(bed.getChromosome(), bed.getStart(), bed.getStop());
		nameReads.clear();
		while (it.hasNext()){
			SAMRecord sam = it.next();
			ArrayList<SAMRecord> al = nameReads.get(sam.getReadName());
			if (al == null){
				al = new ArrayList<SAMRecord>();
				nameReads.put(sam.getReadName(), al);
			}
			al.add(sam);
		}
		it.close();
		ArrayList<SAMRecord>[] rec = new ArrayList[nameReads.size()];
		nameReads.values().toArray(rec);
		this.randomize(rec, rng);
		return rec;
	}
	
	/**For multiple iterations.*/
	public static void randomize (ArrayList<SAMRecord>[] array, Random rng){     
	    // n is the number of items left to shuffle
	    for (int n = array.length; n > 1; n--) {
	        // Pick a random element to move to the end
	        int k = rng.nextInt(n);  // 0 <= k <= n - 1.
	        // Simple swap of variables
	        ArrayList<SAMRecord> tmp = array[k];
	        array[k] = array[n - 1];
	        array[n - 1] = tmp;
	    }
	}



	private int fetchCount(Bed bed) {
		SAMRecordIterator it = toMatchReader.queryOverlapping(bed.getChromosome(), bed.getStart(), bed.getStop());
		int counter = 0;
		while (it.hasNext()){
			it.next();
			counter++;
		}
		it.close();
		return counter;
	}
	
	public boolean isFailed() {
		return failed;
	}


	public String getErrorMessage() {
		return errorMessage;
	}


	public File getSam() {
		return sam;
	}


	public double getNumPassing() {
		return numPassing;
	}

}
