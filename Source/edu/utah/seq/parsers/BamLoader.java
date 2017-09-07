package edu.utah.seq.parsers;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import util.gen.Gzipper;
import util.gen.IO;

public class BamLoader implements Runnable {
	
	//fields
	private boolean failed = false;
	private int threadNumber = 0;
	private ArrayList<QueryInterval> al = new ArrayList<QueryInterval>();
	private QueryInterval[] toFetch = null;
	private SamReader samReader = null;
	private SamAlignmentLoader sal = null;
	private Gzipper[] writers = null;
	private File[] alignments = null;
	
	public BamLoader(SamAlignmentLoader sal, int threadNumber) throws Exception {
		this.sal = sal;
		this.threadNumber = threadNumber;
		SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		samReader = factory.open(sal.getBam());
		
		createWriters();
	}

	public void run() {	
		try {
			//get next chunk of work
			while (sal.loadRegions(al)){ 
				
				//pull regions to query
				toFetch = new QueryInterval[al.size()];
				al.toArray(toFetch);
				al.clear();
				
				//search for overlapping records
				SAMRecordIterator samIterator = samReader.queryOverlapping(toFetch);
				
				//for each record
				while (samIterator.hasNext()) {
					SAMRecord sam = samIterator.next();	
					parseSam(sam);
				}
				
				samIterator.close();
			}
			for (Gzipper g: writers) g.close();
			samReader.close();
			
			//check sizes and delete if empty
			IO.deleteZeroSizedFiles(alignments);
			
		} catch (Exception e) {
			failed = true;
			IO.deleteFiles(alignments);
			System.err.println("\nError: problem fetching alignments" );
			e.printStackTrace();
			try {
				samReader.close();
				for (Gzipper g: writers) g.close();
			} catch (IOException e1) {
			}
			
		}
	}
	
	public boolean isFailed() {
		return failed;
	}
	
	/*Extend and override these methods to extend the functionality*/
	
	private void parseSam(SAMRecord sam) throws Exception{
		//is it paired
		if (sam.getReadPairedFlag() && sam.getProperPairFlag()){
			if (sam.getFirstOfPairFlag()) writers[0].print(sam.getSAMString());
			else if (sam.getSecondOfPairFlag()) writers[1].print(sam.getSAMString());
		}
	}
	
	
	private void createWriters() throws Exception {
		alignments = new File[2];
		alignments[0] = new File (sal.getResultsDir(), threadNumber+"_R1.sam.gz");
		alignments[1] = new File (sal.getResultsDir(), threadNumber+"_R2.sam.gz");
		
		writers = new Gzipper[2];
		writers[0] = new Gzipper(alignments[0]);
		writers[1] = new Gzipper(alignments[1]);
	}
	
	
}