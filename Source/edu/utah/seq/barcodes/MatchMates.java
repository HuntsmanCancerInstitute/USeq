package edu.utah.seq.barcodes;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import edu.utah.seq.data.sam.PicardSortSam;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import util.gen.*;

/**
 * App for attaching mates to alignments and modifying the start to enable consensus calling on barcoded alignments.
 * Do not use the output without running it through the Consensus app first!
 * Writes summary stats to file in json format.
 * @author Nix
 * */
public class MatchMates {

	//user defined fields
	private File bamFile;
	private File saveDirectory;

	//internal fields
	private SamReader bamReader;
	private SAMFileWriter passingBamWriter;
	private SAMFileWriter failingBamWriter;
	private File jsonOutputFile;

	//counters
	private int numRawAlignments = 0;
	private int numPaired = 0;
	private int numUnpaired = 0;
	private int numUnalignedPairs = 0;
	private int numUnalignedSingles = 0;
	private int numFailingAlignments;
	private File passingBamFile;

	//constructors
	/**Stand alone.*/
	public MatchMates(String[] args){
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
			System.out.println("MM statistics:");
			System.out.println(numRawAlignments+"\tNum input alignments");
			System.out.println(numFailingAlignments+"\tNum non primary or vendor failing QC alignments");
			System.out.println(numPaired+"\tNum alignment pairs joined and saved");
			System.out.println(numUnpaired+"\tNum unpaired alignments passing filters and saved");
			System.out.println(numUnalignedPairs+"\tNum unaligned pairs failed");
			System.out.println(numUnalignedSingles+"\tNum unaligned unpaired alignments failed");

			if (jsonOutputFile !=null) saveJson();
			
			//sort by coordinate
			System.out.println("\nSorting passing bam file by coordinate...");
			File sortedBam = new File(saveDirectory, "passingMM.sorted.bam");
			new PicardSortSam(passingBamFile, sortedBam);

		} catch (Exception e) {
			e.printStackTrace();
		} 
	}

	
	@SuppressWarnings("unchecked")
	private void saveJson() {
		try {			
			//output simple json, DO NOT change the key names without updated downstream apps that read this file!
			Gzipper gz = new Gzipper(jsonOutputFile);
			gz.println("{");
			gz.printJson("numberInputAlignments", numRawAlignments, true);
			gz.printJson("numberFailingAlignments", numFailingAlignments, true);
			gz.printJson("numberPairedAlignments", numPaired, true);
			gz.printJson("numberUnpairedAlignments", numUnpaired, true);
			gz.printJson("numberUnalignedPairs", numUnalignedPairs, true);
			gz.printJson("numberUnalignedSingles", numUnalignedSingles, false);
			gz.println("}");
			gz.close();
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem writing json file! "+jsonOutputFile);
		}
	}

	private void walkAlignments() throws Exception {
		//load first
		SAMRecordIterator it = bamReader.iterator();
		ArrayList<SAMRecord> sameNameSams = new ArrayList<SAMRecord>();
		SAMRecord s = it.next();
		String lastName = s.getReadName();
		sameNameSams.add(s);
		numRawAlignments++;

		//load rest collecting all alignments with the same name (should be just two!)
		while (it.hasNext()) {
			SAMRecord sam = it.next();
			numRawAlignments++;

			//remove any that are secondary or fail vendor QC, want to keep unmapped incase their mate is mapped
			//calling all of the methods since these are not getting into the failed file for some reason
			if (sam.getSupplementaryAlignmentFlag() || sam.isSecondaryOrSupplementary() || sam.getReadFailsVendorQualityCheckFlag() || sam.getNotPrimaryAlignmentFlag()){
				failingBamWriter.addAlignment(sam);
				numFailingAlignments++;
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
			//is fop aligned?
			if (fop.getReadUnmappedFlag() == false) printPair (fop, sop);
			//is sop aligned
			else if (sop.getReadUnmappedFlag() == false) printPair(sop, fop);
			//ugg both unaligned so fail em
			else {
				failingBamWriter.addAlignment(sop);
				failingBamWriter.addAlignment(fop);
				numUnalignedPairs++;
			}
		}
		else if (num == 1){
			//no mate 
			SAMRecord s = sameNameSams.get(0);
			if (s.getReadUnmappedFlag() == false) printUnpaired(sameNameSams.get(0));
			else {
				failingBamWriter.addAlignment(s);
				numUnalignedSingles++;
			}
		}
		else {
			//must be more than 2, throw error!
			StringBuilder sb = new StringBuilder("More than 2 alignments found with the same fragment name?\n");
			for (SAMRecord s: sameNameSams) sb.append(s.getSAMString());
			throw new Exception(sb.toString());
		}
	}


	private void printPair(SAMRecord a, SAMRecord b) throws Exception {
		//clean args 
		filterAttributesForRG_AS(a);
		filterAttributesForAS(b);

		//add String rep of b to a
		String bString = b.getSAMString().trim();
		int nameLength = a.getReadNameLength();
		//remove the name and tab
		bString = bString.substring(nameLength +1);
		//replace tabs with spaces
		bString = Misc.TAB.matcher(bString).replaceAll(" ");
		//set
		a.setAttribute("MT", bString);

		//reset start position to unclipped start so it sorts on unclipped
		a.setAlignmentStart(a.getUnclippedStart());

		//save to file
		passingBamWriter.addAlignment(a);
		numPaired++;
	}

	private void printUnpaired(SAMRecord a) throws Exception {
		//clean args 
		filterAttributesForRG_AS(a);

		//reset start position to unclipped start so it sorts on unclipped
		a.setAlignmentStart(a.getUnclippedStart());

		//save to file
		passingBamWriter.addAlignment(a);
		numUnpaired++;
	}



	/**Removes everything but RG: AS:*/
	private void filterAttributesForRG_AS(SAMRecord sam) throws Exception{
		String rg = sam.getStringAttribute("RG");
		Integer as = sam.getIntegerAttribute("AS");
		if (rg == null || as == null) throw new Error("Failed to find an RG or AS attribute in "+sam.getSAMString());
		sam.clearAttributes();
		sam.setAttribute("RG", rg);
		sam.setAttribute("AS", as);
	}
	
	/**Removes everything but AS:*/
	private void filterAttributesForAS(SAMRecord sam) throws Exception{
		Integer as = sam.getIntegerAttribute("AS");
		if (as == null) throw new Error("Failed to find an AS attribute in "+sam.getSAMString());
		sam.clearAttributes();
		sam.setAttribute("AS", as);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MatchMates(args);
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
					case 'b': bamFile = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'j': jsonOutputFile = new File(args[++i]); break;
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
		
		if (jsonOutputFile == null) jsonOutputFile = new File (saveDirectory, "matchMates.json.gz");

		//make IO
		makeIO();

	}	

	private void makeIO() {
		try {
			SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			if (bamFile != null) bamReader = readerFactory.open(bamFile);
			else bamReader = readerFactory.open(SamInputResource.of(System.in));
			
			//check that it isn't sorted by coordinate
			if (bamReader.getFileHeader().getSortOrder().equals(SortOrder.coordinate) || bamReader.hasIndex()) {
				bamReader.close();
				throw new IOException ("\nError: you bam file cannot be sorted by coordinate, must be by query name.\n");
			};

			passingBamFile = new File(saveDirectory, "passMM.bam");
			passingBamFile.deleteOnExit();
			File fail = new File(saveDirectory, "failMM.bam");
			SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
			writerFactory.setCreateIndex(false);
			writerFactory.setTempDirectory(saveDirectory);

			//must explicit set into the header that it is sorted for samtools proc alignments
			passingBamWriter = writerFactory.makeBAMWriter(bamReader.getFileHeader(), false, passingBamFile);
			failingBamWriter = writerFactory.makeBAMWriter(bamReader.getFileHeader(), false, fail);

		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: problem opening the bam IO readers or writers \n");
		}
	}

	private void closeIO() {
		try {
			bamReader.close();
			passingBamWriter.close();
			failingBamWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: critical, problem closing the bam IO \n");
		}
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                 MatchMates: Sept 2015                            **\n" +
				"**************************************************************************************\n" +
				"This app attaches mates of aligned first of pair reads to the attributes and modifies\n"+
				"the start position to enable sorting by unclipped start. Call Consensus to cluster and\n"+
				"collapse alignments with related molecular barcodes.\n"+

				"\nOptions:\n"+
				"-s (Required) Provide a directory path for saving the modified alignments.\n"+
				"-b Path to a query name sorted bam/sam alignment file, defaults to reading from STDIN. \n"+
				"-j Write summary stats in json format to this file, defaults to save directory.\n"+

				"\n"+

				"Example: myAligner | java -Xmx2G -jar pathTo/USeq/Apps/MatchMates -s ReadyForConsensus\n\n"+

				"**************************************************************************************\n");

	}


}
