package edu.utah.seq.parsers.jpileup;

import java.io.File;
import java.io.IOException;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.data.sam.SamLayoutForMutation;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.IO;

public class BamPileupLoader implements Runnable { 

	//fields
	private boolean failed = false;
	private SamReader[] samReaders = null;
	private int minBaseQuality = 0;
	private int minMappingQuality = 0;
	private Gzipper out = null;
	private IndexedFastaSequenceFile fasta = null;
	private Bed[] regions = null;
	private File resultsFile = null;
	private boolean includeOverlaps;


	public BamPileupLoader (BamPileup bamPileup, int loaderIndex, Bed[] regions) throws IOException{
		minBaseQuality = bamPileup.getMinBaseQuality();
		minMappingQuality = bamPileup.getMinMappingQuality();
		this.regions = regions;
		includeOverlaps = bamPileup.isIncludeOverlaps();

		//create gzipper for results
		resultsFile = new File(bamPileup.getTempDir(), loaderIndex+"_tempBamPileup.txt.gz");
		resultsFile.deleteOnExit();
		out = new Gzipper(resultsFile);
		
		//create sam readers
		File[] bamFiles = bamPileup.getBamFiles();
		samReaders = new SamReader[bamFiles.length];
		for (int i=0; i< samReaders.length; i++) samReaders[i] = bamPileup.getSamFactory().open(bamFiles[i]);

		//add header?
		if (loaderIndex == 0) {
			out.println("# MinMapQual\t"+minMappingQuality);
			out.println("# MinBaseQual\t"+minBaseQuality);
			out.println("# IncludeOverlappingBpCounts "+includeOverlaps);
			out.println("# Bed "+IO.getCanonicalPath(bamPileup.getBedFile()));
			for (int i=0; i<bamFiles.length; i++) out.println("# Bam "+i+" "+IO.getCanonicalPath(bamFiles[i]));
			out.println("# Chr\t1BasePos\tRef\tA,C,G,T,N,Del,Ins,FailBQ");
		}


		//create fasta sequence reader
		fasta = new IndexedFastaSequenceFile(bamPileup.getFastaFile());
		if (fasta.isIndexed() == false) throw new IOException("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n");


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
//IO.pl("\tReg "+region);
				//pileup each base from each bam
				BaseCount[][] bamBC = new BaseCount[samReaders.length][];
				for (int i=0; i< samReaders.length; i++ ) {
					bamBC[i]  = pileup(chr, region.getStart(), region.getStop(), samReaders[i]);
				}

				//for each base
				for (int bp=0; bp< region.getLength(); bp++) {
					out.print(chr+"\t"+ (bamBC[0][bp].bpPosition +1)+ "\t"+ bamBC[0][bp].ref);
//IO.p(chr+"\t"+ (bamBC[0][bp].bpPosition +1)+ "\t"+ bamBC[0][bp].ref);
					//for each bam
					for (int b = 0; b<samReaders.length; b++ ) {
						out.print("\t");
						out.print(bamBC[b][bp].toString());
//IO.p("\t"+bamBC[b][bp].toString());
					}
					out.println();
//IO.pl();
				}
			}
		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem processing\n" );
			e.printStackTrace();
		} finally {
			out.closeNoException();
			try {
				fasta.close();
				for (SamReader sr: samReaders) sr.close();
			} catch (IOException e) {}
		}
	}

	private BaseCount[] pileup(String chr, int start, int stop, SamReader samReader) throws Exception{
		//create container for counts
		ReferenceSequence p = fasta.getSubsequenceAt(chr, start+1, stop+1);
		char[] refSeq = new String(p.getBases()).toCharArray();
		BaseCount[] bc = new BaseCount[stop-start];
		int counter = 0;
		for (int i=start; i< stop; i++) bc[counter] = new BaseCount(i, refSeq[counter++]);

		//fetch alignments
		SAMRecordIterator it = samReader.queryOverlapping(chr, start-1, stop+1);
		if (it.hasNext() == false) {
			it.close();
			return bc;
		}

		//for each record
		while (it.hasNext()){

			SAMRecord sam = it.next();
			if (sam.getMappingQuality() < minMappingQuality) continue;

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
				if (includeOverlaps == false) {
					if (bc[counter].readNames.contains(readName)) continue;
					else bc[counter].readNames.add(readName);
				}


				char call = layout.getCall()[index];

				if (call == 'M') {
					char base = layout.getSeq()[index];
					//System.out.println("\tMatch "+layout.getSeq()[index]);
					if (base == 'N' ) bc[counter].n++;
					else if (qual < minBaseQuality) bc[counter].failQual++;
					else bc[counter].increment(base);
				}
				//an deletion
				else if (call == 'D') {
					bc[counter].del++;
					//System.out.println("\tdeletion");
				}
				//an insertion
				else if (call == 'I') {
					bc[counter].ins++;

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
								//System.out.println("\tMatch "+layout.getSeq()[index]);
								if (base == 'N' ) bc[counter].n++;
								else if (qual < minBaseQuality) bc[counter].failQual++;
								else bc[counter].increment(base);
							}
							//an deletion
							else if (call == 'D') {
								bc[counter].del++;
								//System.out.println("\tdeletion");
							}
							//do nothing
						}
					}
					//System.out.println("\tinsertion ");
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

	public File getResultsFile() {
		return resultsFile;
	}
}
