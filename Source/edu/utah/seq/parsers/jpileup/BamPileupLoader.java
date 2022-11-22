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
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class BamPileupLoader implements Runnable { 

	//fields
	private boolean failed = false;
	private SamReader[] samReaders = null;
	private int minBaseQuality = 0;
	private int minMappingQuality = 0;
	private PrintWriter pileupOut = null;
	private PrintWriter bedOut= null;
	private IndexedFastaSequenceFile fasta = null;
	private Bed[] regions = null;
	private File pileupFile = null;
	
	private boolean alignNoChr;
	private boolean includeOverlaps;
	private boolean printAll;
	private boolean verbose; 
	private int maxNumReads = 0;


	public BamPileupLoader (BamPileup bamPileup, int loaderIndex, Bed[] regions) throws IOException{
		minBaseQuality = bamPileup.getMinBaseQuality();
		minMappingQuality = bamPileup.getMinMappingQuality();
		this.regions = regions;
		includeOverlaps = bamPileup.isIncludeOverlaps();
		printAll = bamPileup.isPrintAll();
		verbose = bamPileup.isVerbose();
		alignNoChr = bamPileup.isAlignNoChr();

		//create writers
		pileupFile = new File(bamPileup.getTempDir(), loaderIndex+"_tempBamPileup.txt");
		pileupFile.deleteOnExit();
		pileupOut = new PrintWriter (new FileWriter(pileupFile));
		
		//create sam readers
		File[] bamFiles = bamPileup.getBamFiles();
		samReaders = new SamReader[bamFiles.length];
		for (int i=0; i< samReaders.length; i++) {
			samReaders[i] = bamPileup.getSamFactory().open(bamFiles[i]);
			if (samReaders[i].hasIndex() == false) {
				failed = true;
				throw new IOException("Failed to find an index for "+bamFiles[i]);
			}
		}

		//add header?
		if (loaderIndex == 0) {
			pileupOut.println("# MinMapQual\t"+minMappingQuality);
			pileupOut.println("# MinBaseQual\t"+minBaseQuality);
			pileupOut.println("# IncludeOverlappingBpCounts "+includeOverlaps);
			pileupOut.println("# PrintAll "+printAll);
			pileupOut.println("# Bed "+IO.getCanonicalPath(bamPileup.getBedFile()));
			for (int i=0; i<bamFiles.length; i++) pileupOut.println("# BamCram\t"+i+"\t"+IO.getCanonicalPath(bamFiles[i]));
			pileupOut.println("# Chr\t1BasePos\tRef\tA,C,G,T,N,Del,Ins,FailBQ");
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
					if (verbose) IO.p(".");
				}
				
				//need to account for Stubben's RNASeq alignments which are to GRCh38 where chroms are 1,2,3..X,Y,MT
				//this fix won't work for the Random and Alt contigs!
				String chr = region.getChromosome();
				String chrToSearchAli = chr;
				if (alignNoChr && chr.startsWith("chr")) {
					//watch out for mito
					if (chr.equals("chrM")) chrToSearchAli = "MT";
					else chrToSearchAli = chrToSearchAli.substring(3);
				}

				//pileup each base from each bam
				BaseCount[][] bamBC = new BaseCount[samReaders.length][];
				for (int i=0; i< samReaders.length; i++ ) {
					bamBC[i]  = pileup(chr, chrToSearchAli, region.getStart(), region.getStop(), samReaders[i]);
				}

				//for each base in the region
				StringBuilder sb  = null;
				boolean noCounts = true;
				for (int bp=0; bp< region.getLength(); bp++) {
					sb  = new StringBuilder();
					noCounts = true;
					
					sb.append(chr); sb.append("\t");
					sb.append((bamBC[0][bp].bpPosition +1)); sb.append("\t");
					sb.append(bamBC[0][bp].ref);
					
					//out.print(chr+"\t"+ (bamBC[0][bp].bpPosition +1)+ "\t"+ bamBC[0][bp].ref);

					//for each bam
					for (int b = 0; b<samReaders.length; b++ ) {
						sb.append("\t");
						if (bamBC[b][bp].loadStringBuilderWithCounts(sb)  == true) {
							noCounts = false;
						}
					}
					if (noCounts == false || printAll == true) pileupOut.println(sb.toString());
				}
			}
		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem processing\n" );
			e.printStackTrace();
		} finally {
			try {
				pileupOut.close();
				if(bedOut != null) bedOut.close();
				fasta.close();
				for (SamReader sr: samReaders) sr.close();
			} catch (IOException e) {}
		}
	}

	private BaseCount[] pileup(String chr, String chrToSearchAli, int start, int stop, SamReader samReader) throws Exception{
		//watch end
		int stopPlusOne = stop+1;
		int chromEnd = (int)fasta.getIndex().getIndexEntry(chrToSearchAli).getSize();
		if (stopPlusOne > chromEnd) stopPlusOne = chromEnd;
		
		//create container for counts
		ReferenceSequence p = fasta.getSubsequenceAt(chrToSearchAli, start+1, stopPlusOne);
		
		char[] refSeq = new String(p.getBases()).toCharArray();
		BaseCount[] bc = new BaseCount[stop-start];
		int counter = 0;
		for (int i=start; i< stop; i++) bc[counter] = new BaseCount(i, refSeq[counter++]);

		//fetch alignments, will throw exception if the chr isn't in the index
		SAMRecordIterator it = samReader.queryOverlapping(chrToSearchAli, start-1, stop+1);
		if (it.hasNext() == false) {
			it.close();
			return bc;
		}
		int numReadsProcessed = 0;
		//for each record
		while (it.hasNext()){
			
			SAMRecord sam = it.next();
			if (sam.getMappingQuality() < minMappingQuality || sam.isSecondaryOrSupplementary() || sam.getDuplicateReadFlag() || sam.getReadFailsVendorQualityCheckFlag()) continue;
			numReadsProcessed++;
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
					if (bc[counter].readNames.contains(readName)) {
						counter++;
						continue;
					}
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
		if (numReadsProcessed > maxNumReads) maxNumReads = numReadsProcessed;
		it.close();
		return bc;
	}

	public boolean isFailed() {
		return failed;
	}

	public File getPileupFile() {
		return pileupFile;
	}

	public int getMaxNumReads() {
		return maxNumReads;
	}
}
