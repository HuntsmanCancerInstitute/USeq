package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import edu.utah.seq.barcodes.ConsensusEngine;
import edu.utah.seq.data.sam.PicardSortSam;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import util.gen.*;

/**
 * App for extracting fastq data from a queryname sorted sam or bam file.  Doesn't have the memory leak in picard.  Writes as compressed gz. Error checks, unlike samtools.
 * @author Nix
 * */
public class Sam2Fastq {

	//user defined fields
	private File bamFile;
	private File saveDirectory;

	//internal fields
	private SamReader bamReader;
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
		} 
	}


	private void walkAlignments() throws Exception {

		SAMRecordIterator it = bamReader.iterator();
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
			//remove any that are secondary or fail vendor QC, want to keep unmapped incase their mate is mapped
			//calling all of the methods since these are not getting into the failed file for some reason
			if (sam.getSupplementaryAlignmentFlag() || sam.isSecondaryOrSupplementary() || sam.getReadFailsVendorQualityCheckFlag() || sam.getNotPrimaryAlignmentFlag()){
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


	private void processSams(ArrayList<SAMRecord> sameNameSams) throws Exception {
		int num = sameNameSams.size();
		if (num == 2){
			//which to save?
			SAMRecord fop = sameNameSams.get(0);
			SAMRecord sop = sameNameSams.get(1);
			if (fop.getFirstOfPairFlag() == false){
				SAMRecord s = fop;
				fop = sop;
				sop = s;
			}
			for (String x : ConsensusEngine.exportFastq(fop)) firstRead.println(x);
			for (String x : ConsensusEngine.exportFastq(sop)) secondRead.println(x);
		}
		else if (num == 1){
			//no mate 
			SAMRecord s = sameNameSams.get(0);
			for (String x : ConsensusEngine.exportFastq(s)) unpairedRead.println(x);
		}
		else {
			//must be more than 2, throw error!
			StringBuilder sb = new StringBuilder("\nMore than 2 alignments found with the same fragment name?\n");
			for (SAMRecord s: sameNameSams) sb.append(s.getSAMString());
			throw new Exception(sb.toString());
		}
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
			
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: problem opening the bam IO readers or writers \n");
		}
	}

	private void closeIO() {
		try {
			bamReader.close();
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
				"**                                 Sam 2 Fastq: July 2016                           **\n" +
				"**************************************************************************************\n" +
				"Given a queryname sorted bam/sam alignment file, writes out fastq data for paired and \n"+
				"unpaired alignemnts. Any non primary, secondary, supplemental or vendor failing QC\n"+
				"alignments are written to a failed sam file.  This app doesn't have the memory leak\n"+
				"found in Picard, writes gzipped fastq, and error checks the reads.\n"+

				"\nOptions:\n"+
				"-s Path to a directory for saving parsed data..\n"+
				"-a Path to a query name sorted bam/sam alignment file. \n"+

				"\n"+

				"Example: java -Xmx2G -jar pathTo/USeq/Apps/Sam2Fastq -a myQNSorted.bam -s S2F/\n\n"+

				"**************************************************************************************\n");

	}


}
