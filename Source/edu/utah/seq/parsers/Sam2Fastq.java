package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import util.bio.seq.Seq;
import util.gen.*;

/**
 * App for extracting fastq data from a queryname sorted sam or bam file.  Doesn't have the memory leak in picard.  Writes as compressed gz. Error checks, unlike samtools.
 * @author Nix
 * */
public class Sam2Fastq {

	//user defined fields
	private File bamFile;
	private File unfilteredBamFile;
	private File saveDirectory;

	//internal fields
	private SamReader bamReader;
	private SamReader unfilteredBamReader = null;
	private SAMRecordIterator unfilteredIt = null;
	private SAMRecord lastSam = null;
	private Gzipper firstRead;
	private Gzipper secondRead;
	private Gzipper unpairedRead;
	private Gzipper failingSam;
	
	//counters
	private long numAlignments = 0;
	private long numPaired = 0;
	private long numUnpaired = 0;
	private long numFailingAlignments = 0;
	
	

	//constructors
	/**Stand alone.*/
	public Sam2Fastq(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//launch
		run();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}


	public void run(){
		try {
			walkAlignments();

			//close readers
			closeIO();

			//print stats
			System.out.println("Statistics:");
			System.out.println(numAlignments+"\tNum input alignments");
			System.out.println(numFailingAlignments+"\tNum non primary, secondary, supplemental or vendor failing QC alignments written to failed sam file.");
			System.out.println(numPaired+"\tNum alignment pairs exported as fastq");
			System.out.println(numUnpaired+"\tNum unpaired alignments exported as fastq");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem converting sams!\n");
		} 
	}


	private void walkAlignments() throws Exception {

		SAMRecordIterator it = bamReader.iterator();
		if (unfilteredBamFile != null) unfilteredIt = unfilteredBamReader.iterator();
		ArrayList<SAMRecord> sameNameSams = new ArrayList<SAMRecord>();
		String lastName = null;
		int counter = 0;
		//load all alignments with the same name (should be just two!)
		while (it.hasNext()) {
			SAMRecord sam = it.next();
			numAlignments++;
			if (counter++ > 1000000) {
				System.out.print(".");
				counter = 0;
			}
			//remove any that are secondary or not primary, want to keep unmapped in case their mate is mapped
			//calling all of the methods since these are not getting into the failed file for some reason
			if (sam.getSupplementaryAlignmentFlag() || sam.isSecondaryOrSupplementary() || sam.getNotPrimaryAlignmentFlag()){
				failingSam.println(sam.getSAMString().trim());
				numFailingAlignments++;
				continue;
			}
			
			//first one?
			if (lastName == null){
				lastName = sam.getReadName();
				sameNameSams.add(sam);
				continue;
			}

			//part of prior set?
			if (sam.getReadName().equals(lastName)) sameNameSams.add(sam);
			else {
				//no, so a new one so process old and...
				processSams(sameNameSams);
				//reset
				sameNameSams.clear();
				lastName = sam.getReadName();
				sameNameSams.add(sam);
			}
		}
		//process last
		processSams(sameNameSams);
		it.close();
		System.out.println();
	}
	
	
	/**Generates a fastq record from a sam record. Will use the original quality tag value if present, 
	 * e.g. OQ:Z:xxxxxxxxxxx , from GATK --emit_original_quals */
	public static String[] exportFastq(SAMRecord sam) {
		String[] fastq = new String[4];
		//names
		fastq[0] = "@"+sam.getReadName();
		fastq[2] = "+";
		
		String seq = sam.getReadString();
		String qual;
		
		//original qualities present?
		Object o = sam.getAttribute("OQ");
		if (o != null) qual = (String)o;
		else qual = sam.getBaseQualityString();

		//negative strand? reverse em
		if (sam.getReadNegativeStrandFlag()) {
			seq = Seq.reverseComplementDNA(seq);
			qual = new StringBuilder(qual).reverse().toString();
		}
		fastq[1] = seq;
		fastq[3] = qual;
		return fastq;
	}

	private void processSams(ArrayList<SAMRecord> sameNameSams) throws Exception {
		int num = sameNameSams.size();
		//a pair
		if (num == 2){
			//which to save?
			SAMRecord fop = sameNameSams.get(0);
			SAMRecord sop = sameNameSams.get(1);
			if (fop.getFirstOfPairFlag() == false){
				SAMRecord s = fop;
				fop = sop;
				sop = s;
			}
			for (String x : exportFastq(fop)) firstRead.println(x);
			for (String x : exportFastq(sop)) secondRead.println(x);
			numPaired++;
		}
		//just one
		else if (num == 1){
			SAMRecord s = sameNameSams.get(0);
			SAMRecord m = null;
			
			//look for mate? 
			if (unfilteredBamFile != null){
				//attempt to find from raw
				ArrayList<SAMRecord> rawSams = fetchUnfilteredSams(s.getReadName());
				if (rawSams == null){
					//this can sometimes happen because the read names are different and thus the order
					//reset the iterator and try again
					unfilteredIt.close();
					unfilteredIt = unfilteredBamReader.iterator();
					rawSams = fetchUnfilteredSams(s.getReadName());
					if (rawSams == null){
						//still null so fail it and reset iterator
						for (String x : exportFastq(s)) unpairedRead.println(x);
						numUnpaired++;
						unfilteredIt.close();
						unfilteredIt = unfilteredBamReader.iterator();
						return;
					}
				}
				//paired?
				boolean firstOfPair = false;
				if (s.getReadPairedFlag()){
					firstOfPair = s.getFirstOfPairFlag();
					for (SAMRecord sam: rawSams){
						if (sam.getSupplementaryAlignmentFlag() || sam.isSecondaryOrSupplementary() || sam.getNotPrimaryAlignmentFlag()) continue;
						else {
							if (sam.getFirstOfPairFlag() != firstOfPair) {
								m=sam;
								break;
							}
						}
					}
				}
				//single
				else {
					int pos = s.getUnclippedStart();
					boolean strand = s.getReadNegativeStrandFlag();
					for (SAMRecord sam: rawSams){
						if (sam.getSupplementaryAlignmentFlag() || sam.isSecondaryOrSupplementary() || sam.getNotPrimaryAlignmentFlag()) continue;
						else {
							if (sam.getUnclippedStart() != pos || sam.getReadNegativeStrandFlag() != strand) {
								m = sam;
								//use mate to define firstOfPair for the unpaired read
								if (m.getFirstOfPairFlag() == true) firstOfPair = false;
								else firstOfPair = true;
								break;
							}
						}
					}
					
				}

				if (m == null) {
					for (String x : exportFastq(s)) unpairedRead.println(x);
					numUnpaired++;
				}
				else {
					if (firstOfPair){
						for (String x : exportFastq(s)) firstRead.println(x);
						for (String x : exportFastq(m)) secondRead.println(x);
					}
					else {
						for (String x : exportFastq(m)) firstRead.println(x);
						for (String x : exportFastq(s)) secondRead.println(x);
					}
					numPaired++;
				}
			}
			//nope, just print unpaired
			else {
				for (String x : exportFastq(s)) unpairedRead.println(x);
				numUnpaired++;
			}
		}
		
		//more than one found, this shouldn't happen.
		else {
			//must be more than 2, throw error!
			StringBuilder sb = new StringBuilder("\nMore than 2 alignments found with the same fragment name?\n");
			for (SAMRecord s: sameNameSams) sb.append(s.getSAMString());
			throw new Exception(sb.toString());
		}
	}




	private ArrayList<SAMRecord> fetchUnfilteredSams(String readName) {
		ArrayList<SAMRecord> matches = new ArrayList<SAMRecord>();
		boolean found = false;
		
		//check last
		if (lastSam != null && lastSam.getReadName().equals(readName)) {
			matches.add(lastSam);
			found = true;
		}
	
		while (unfilteredIt.hasNext()) {
			lastSam = unfilteredIt.next();
			if (lastSam.getReadName().equals(readName)){
				matches.add(lastSam);
				found = true;
			}
			else if (found) return matches;
		}
		return null;
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new Sam2Fastq(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': bamFile = new File(args[++i]); break;
					case 'u': unfilteredBamFile = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group()+", Try -h");
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for bam files
		if (bamFile != null && bamFile.canRead() == false) Misc.printErrAndExit("\nError: cannot read your xxx.bam file?\n"+bamFile);

		//look for save file, can be null
		if (saveDirectory == null ) Misc.printErrAndExit("\nError: please provide a directory to save the passing and failing bam alignments.\n");
		saveDirectory.mkdirs();

		//make IO
		makeIO();

	}	

	private void makeIO() {
		try {
			SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			bamReader = readerFactory.open(bamFile);

			
			//check that it isn't sorted by coordinate
			if (bamReader.getFileHeader().getSortOrder().equals(SortOrder.coordinate) || bamReader.hasIndex()) {
				bamReader.close();
				throw new IOException ("\nError: you bam file cannot be sorted by coordinate, must be by query name.\n");
			};

			String name = Misc.removeExtension(bamFile.getName());
			File fail = new File(saveDirectory, name+"_Fail.sam.gz");
			File first = new File(saveDirectory, name+"_1.fastq.gz");
			File second = new File(saveDirectory, name+"_2.fastq.gz");
			File single = new File(saveDirectory, name+"_unpaired.fastq.gz");

			firstRead = new Gzipper(first);
			secondRead = new Gzipper(second);
			unpairedRead = new Gzipper(single);
			failingSam = new Gzipper(fail);

			//write header
			String header = bamReader.getFileHeader().getTextHeader().trim();
			failingSam.println(header);
			
			//raw reader?
			if (unfilteredBamFile !=null) unfilteredBamReader = readerFactory.open(unfilteredBamFile);
			
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: problem opening the bam IO readers or writers \n");
		}
	}

	private void closeIO() {
		try {
			bamReader.close();
			if (unfilteredBamFile !=null) unfilteredBamReader.close();
			firstRead.close();
			secondRead.close();
			unpairedRead.close();
			failingSam.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: critical, problem closing the bam IO \n");
		}
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Sam 2 Fastq: March 2018                           **\n" +
				"**************************************************************************************\n" +
				"Given a query name sorted bam/sam alignment file, writes out fastq data for paired and \n"+
				"unpaired alignemnts. Any non-primary, secondary, or supplemental \n"+
				"alignments are written to a failed sam file.  This app doesn't have the memory leak\n"+
				"found in Picard, writes gzipped fastq, and error checks the reads. Provide an\n"+
				"unfiltered bam to use in retrieving missing mate alignments. S2F will remove pre and\n"+
				"post naming info from the BamBlaster restoring the original fragment names.\n"+

				"\nOptions:\n"+
				"-s Path to a directory for saving parsed data..\n"+
				"-a Path to a query name sorted bam/sam alignment file. \n"+
				"-u (Optional) Path to a query name sorted unfiltered bam/sam alignment file for use\n"+
				"      in fetching missing mates of the first bam.\n"+

				"\n"+

				"Example: java -Xmx2G -jar pathTo/USeq/Apps/Sam2Fastq -a myQNSorted.bam -s S2F/\n\n"+

				"**************************************************************************************\n");

	}


}
