package edu.utah.seq.parsers.jpileup;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.data.sam.SamLayoutForMutation;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import util.bio.annotation.Bed;
import util.gen.Histogram;
import util.gen.IO;
import util.gen.Misc;

public class RCBamLoader implements Runnable { 

	//fields
	private boolean failed = false;
	private SamReader samReader = null;
	private int minBaseQuality = 0;
	private int minMappingQuality = 0;
	private int minimumReadDepth = 0;
	private PrintWriter bedOut= null;
	private Bed[] regions = null;
	private File coverageFile = null;
	private Histogram histogram = null;


	public RCBamLoader (UniObRC uniObRC, int loaderIndex, Bed[] regions) throws IOException{
		minBaseQuality = uniObRC.getMinBaseQuality();
		minMappingQuality = uniObRC.getMinMappingQuality();
		minimumReadDepth = uniObRC.getMinimumReadDepth();
		this.regions = regions;

		//create writers
		coverageFile = new File(uniObRC.getTempDir(), Misc.getRandomString(10)+"_"+loaderIndex+"_temp.bed");
		coverageFile.deleteOnExit();
		bedOut = new PrintWriter (new FileWriter(coverageFile));

		//create sam reader
		samReader = uniObRC.getSamFactory().open(uniObRC.getBamFile());
		if (samReader.hasIndex() == false) {
			failed = true;
			throw new IOException("Failed to find an index for "+uniObRC.getBamFile());
		}

		//create histogram to count read depth
		histogram = new Histogram(0, uniObRC.getMaximumCoverageCalculated(), (int)uniObRC.getMaximumCoverageCalculated());
	}

	public void run() {	
		try {
			//for each region
			int counter = 0;			
			for (Bed region: regions) {
				if (counter++ > 100) {
					counter = 0;
					IO.p(".");
				}
				String chr = region.getChromosome();

				//create a BaseCount for each base 
				RCBaseCount[] bamBC = pileup(chr, region.getStart(), region.getStop(), samReader);

				//for each base in the region
				for (int bp=0; bp< region.getLength(); bp++) {
					histogram.count(bamBC[bp].baseCount);
					if (bamBC[bp].baseCount >= minimumReadDepth) bedOut.println(chr+"\t"+bamBC[bp].bpPosition+"\t"+(bamBC[bp].bpPosition+1));
				}
			}
		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem processing\n" );
			e.printStackTrace();
		} finally {
			try {
				bedOut.close();
				samReader.close();
			} catch (IOException e) {}
		}
	}

	private RCBaseCount[] pileup(String chr, int start, int stop, SamReader samReader) throws Exception{

		RCBaseCount[] bc = new RCBaseCount[stop-start];
		int counter = 0;
		for (int i=start; i< stop; i++) bc[counter++] = new RCBaseCount(i);

		//fetch alignments
		SAMRecordIterator it = samReader.queryOverlapping(chr, start-1, stop+1);
		if (it.hasNext() == false) {
			it.close();
			return bc;
		}

		//for each record
		while (it.hasNext()){

			SAMRecord sam = it.next();
			if (sam.getMappingQuality() < minMappingQuality || sam.isSecondaryOrSupplementary() || sam.getDuplicateReadFlag() || sam.getReadFailsVendorQualityCheckFlag()) continue;

			//make a layout
			SamAlignment sa = new SamAlignment(sam.getSAMString().trim(), true);
			SamLayoutForMutation layout = new SamLayoutForMutation(sa);
			String readName = sa.getName();

			//for each base in the region see if it overlaps
			counter = 0;
			for (int i = start; i< stop; i++){
				int index = layout.findBaseIndex(i);
				//System.out.println(i+ " i \t index "+index);

				//present in the alignment?
				if (index == -1) {
					//System.out.println("\tNotFound");
					counter++;
					continue;
				}

				//watch out for -1, not set
				int qual = layout.getQual()[index];
				if (qual == -1) {
					counter++;
					continue;
				}

				//mate already processed
				if (bc[counter].readNames.contains(readName)) {
					counter++;
					continue;
				}
				else bc[counter].readNames.add(readName);

				char call = layout.getCall()[index];

				if (call == 'M') {
					char base = layout.getSeq()[index];
					if (base != 'N' && qual >= minBaseQuality) bc[counter].baseCount++;
				}
				//an deletion
				else if (call == 'D') {
					bc[counter].baseCount++;
				}
				//an insertion
				else if (call == 'I') {
					bc[counter].baseCount++;

					//advance till base pos changes
					int[] pos = layout.getPos();
					char[] calls = layout.getCall();
					int currPos = pos[index];

					for (int x=index+1; x<pos.length; x++) {
						//diff base? exit
						if (pos[x] != currPos) break;
						//same base, score M or dels
						else {
							call = calls[x];
							if (call == 'M') {
								char base = layout.getSeq()[index];
								if (base != 'N' && qual >= minBaseQuality) bc[counter].baseCount++;
							}
							//an deletion
							else if (call == 'D') bc[counter].baseCount++;
							//do nothing
						}
					}

				}

				//nope must be a masked base (S,H) or N so skip
				//else System.out.println("\tMask or N "+call);
				counter++;
			}
		}
		it.close();
		return bc;
	}

	public boolean isFailed() {
		return failed;
	}

	public File getCoverageFile() {
		return coverageFile;
	}

	public Histogram getHistogram() {
		return histogram;
	}
}
