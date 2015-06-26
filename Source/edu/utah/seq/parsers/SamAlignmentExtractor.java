package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import htsjdk.samtools.*;
import util.bio.annotation.Bed;
import util.bio.seq.Seq;
import util.gen.*;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.useq.data.RegionScoreText;


/**
 * @author Nix
 * */
public class SamAlignmentExtractor {

	//user defined fields
	private File[] bamFiles;
	private File bedFile;
	private File saveFile;
	private int minimumReadDepth = 1;
	private int maximumReadDepth = -1;
	private int minimumMappingQuality = -1;
	private int alignmentScoreThreshold = -1;
	private boolean biggerASIsBetter = true;
	private boolean calcExtraStats = false;
	private int minimumBaseQuality = 20;

	//internal fields
	private SamReader[] samReaders;
	private HashMap<String,RegionScoreText[]> chromRegions;
	private String chromosome;
	private static Pattern CIGAR_SUB = Pattern.compile("(\\d+)([MSDHN])");
	private Gzipper samOut;
	private HashSet<String> hits = null;
	private boolean printCoverageStats = false;
	private SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	private int numDroppedAlignments = 0;
	private int numFetchedAlignments = 0;
	private int numPrintedAlignments = 0;
	private String useqArguments;
	
	//only for extra stats
	private int numberPassingAlignments = 0;
	private double numberPassingBaseScores = 0.0;
	private double numberTotalBaseScores = 0.0;
	
	//constructors
	/**Stand alone.*/
	public SamAlignmentExtractor(String[] args){
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
			//make readers on each bam file
			makeSamReaders();

			//make output writer
			File tempSam = new File(saveFile + ".temp.sam.gz");
			tempSam.deleteOnExit();
			samOut = new Gzipper(tempSam);
			
			//add header from first reader
			String header = samReaders[0].getFileHeader().getTextHeader();
			if (header != null){
				samOut.println(header.trim());
				//add program args
				samOut.println("@PG\tID:USeq SamAlignmentExtractor "+Misc.getDateTime()+"\tCL:"+useqArguments);
			}
			
			//calc extra stats?
			if (calcExtraStats){
				System.out.println("Calculating extra coverage statistics...");
				for (int i=0; i< samReaders.length; i++) numberPassingAlignments += statAllAlignments(samReaders[i]);
			}

			//for each chromosome of gene models
			System.out.println("Scanning regions by chromosome...\n");
			if (printCoverageStats) System.out.println("RegionCoordinates\t#Alignments");
			Iterator<String> it = chromRegions.keySet().iterator();
			while (it.hasNext()){
				chromosome = it.next();
				scanRegions();
			}
			if (printCoverageStats) System.out.println();
			
			//look for no hits?
			if (hits != null){
				System.out.println("Scanning for alignments that don't intersect...");
				System.out.println("#Ints\t"+hits.size());
				printNoHits();
				System.out.println("#NonInts\t"+numPrintedAlignments);
			}
			else {
				String passed = Num.formatNumber((double)numPrintedAlignments/(double)numFetchedAlignments, 4);
				System.out.println(numFetchedAlignments +"\t#Intersecting alignments");
				System.out.println(numDroppedAlignments +"\t#Failed flags or filters (MQ:"+minimumMappingQuality+" AS:"+alignmentScoreThreshold+" BIB:"+biggerASIsBetter+")");
				System.out.println(numPrintedAlignments +"\t#Saved to sam file ("+passed+")");
				if (calcExtraStats){
					String onTarget = Num.formatNumber((double)numPrintedAlignments/(double)numberPassingAlignments, 4);
					System.out.println(onTarget +"\tFraction on target");
					String passBases = Num.formatNumber(numberPassingBaseScores/numberTotalBaseScores, 4);
					System.out.println(passBases +"\tFraction on target base quality scores >= "+minimumBaseQuality);
				}
			}

			//close readers
			closeSamReaders();
			samOut.close();
			
			//sort?
			System.out.println("\nSorting...");
			new PicardSortSam(samOut.getGzipFile(), saveFile);

		} catch (Exception e) {
			e.printStackTrace();
		} 
	}

	private void printNoHits() throws Exception{
		//for each sam reader
		for (int i=0; i< samReaders.length; i++) {
			//make gzipper
			Gzipper samOutNoHit = new Gzipper(new File (bamFiles[i].getParentFile(), Misc.removeExtension(bamFiles[i].getName())+"NoInt.sam.gz"));
			//add header
			String header = samReaders[i].getFileHeader().getTextHeader();
			samOutNoHit.println(header.trim());
			//get iterator
			SAMRecordIterator it = samReaders[i].iterator();
			while (it.hasNext()) {
				SAMRecord sam = it.next();
				if (passThresholds(sam, false) == false) continue;
				String samString = sam.getSAMString().trim();
				if (hits.contains(samString) == false) {
					samOutNoHit.println(samString);
					numPrintedAlignments++;
				}
			}
			it.close();
			samOutNoHit.close();
		}
		
	}

	/**Checks that the alignment actually touches down on at least one base of the region to avoid spanners.*/
	public ArrayList<String> fetchAlignments (RegionScoreText ei, SamReader reader){
		ArrayList<String> al = new ArrayList<String>();
		SAMRecordIterator i = reader.queryOverlapping(chromosome, (ei.getStart()+1), ei.getStop());
		while (i.hasNext()) {
			SAMRecord sam = i.next();
			if (passThresholds(sam, true) == false) continue;
			//fetch blocks of actual alignment
			ArrayList<int[]> blocks = fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);
			//check to see if any intersect the exon
			for (int[] b : blocks){
				if (ei.intersects(b[0], b[1])){
					al.add(sam.getSAMString().trim());
					break;
				}
			}
		}
		i.close();
		i = null;
		return al;
	}
	
	/**Counts the number of alignments that pass thresholds.*/
	public int statAllAlignments (SamReader reader){
		int numPassingAlignments = 0;
		SAMRecordIterator i = reader.iterator();
		while (i.hasNext()) {
			SAMRecord sam = i.next();
			if (passThresholds(sam, false) == false) continue;
			numPassingAlignments++;
		}
		i.close();
		i = null;
		return numPassingAlignments;
	}
	
	public boolean passThresholds(SAMRecord sam, boolean incrementCounters){
		if (incrementCounters) numFetchedAlignments++;
		//any problematic flags?
		if (sam.getReadUnmappedFlag() || sam.getNotPrimaryAlignmentFlag() || sam.getReadFailsVendorQualityCheckFlag()) {
			if (incrementCounters) numDroppedAlignments++;
			return false;
		}
		//pass mapping quality?
		if (minimumMappingQuality != -1 && sam.getMappingQuality() < minimumMappingQuality) {		
			if (incrementCounters) numDroppedAlignments++;
			return false;
		}
		//pass alignment score? with novoalignments (~30 pt penalty per mismatch) smaller AS is better 
		//with bwa bigger scores are better (readLength - #SNV*5 + #INDEL*7)
		if (alignmentScoreThreshold != -1){
			Object obj = sam.getAttribute("AS");
			if (obj != null){
				Integer as = (Integer)obj;
				if (biggerASIsBetter){ 
					if (as.intValue() < alignmentScoreThreshold) {
						if (incrementCounters) numDroppedAlignments++;
						return false;
					}
				}
				else {
					if (as.intValue() > alignmentScoreThreshold) {
						if (incrementCounters) numDroppedAlignments++;
						return false;
					}
				}
			}
		}
		return true;
	}

	/**Assumes interbase coordinates for start and returned blocks.*/
	public static ArrayList<int[]> fetchAlignmentBlocks(String cigar, int start){
		//for each cigar block
		Matcher mat = CIGAR_SUB.matcher(cigar);
		ArrayList<int[]> blocks = new ArrayList<int[]>();
		while (mat.find()){
			String call = mat.group(2);
			int numberBases = Integer.parseInt(mat.group(1));
			//a match
			if (call.equals("M")) {
				blocks.add(new int[]{start, start+numberBases});
			}
			//just advance for all but insertions which should be skipped via the failure to match
			start += numberBases;
		}
		return blocks;
	}

	/**For each sam file, fetches all sam alignments.*/
	public ArrayList<String> fetchOverlappingAlignments (RegionScoreText ei){
		ArrayList<String> ohMy = new ArrayList<String>();
		for (int i=0; i< samReaders.length; i++) {
			ohMy.addAll(fetchAlignments(ei, samReaders[i]));
		}
		return ohMy;
	}

	private void scanRegions() throws IOException{
		RegionScoreText[] regions = chromRegions.get(chromosome);
		HashSet<String> uniqueAlignments = new HashSet<String>();
		//for each region in a particular chromosome
		for (int i=0; i< regions.length; i++){
			//fetch the overlapping alignments
			ArrayList<String> alignments = fetchOverlappingAlignments (regions[i]);
			
			//pass min and max?
			int numAlignments = alignments.size();
			if (numAlignments >= minimumReadDepth){
				//check max?
				if (maximumReadDepth != -1){
					if (numAlignments > maximumReadDepth) continue;
				}
				//print em
				if (printCoverageStats) System.out.println(chromosome+"\t"+regions[i].toString()+"\t"+alignments.size());
				//strip those already fetched to avoid repeat extraction!
				ArrayList<String> toPrint = new ArrayList<String>(alignments.size());
				for (String sam: alignments){
					if (uniqueAlignments.contains(sam) == false){
						uniqueAlignments.add(sam);
						toPrint.add(sam);
					}
				}
				samOut.println(toPrint);
				numPrintedAlignments+=toPrint.size();
				if (calcExtraStats) incrementBaseQualities(toPrint);
			}
		}
		//save hits?
		if (hits != null) hits.addAll(uniqueAlignments);
	}

	private void incrementBaseQualities(ArrayList<String> toPrint) {
		//parse base qualities
		for (String sam : toPrint){
			String[] fields = Misc.TAB.split(sam);
			int passingBases = Seq.countNumberOfPassingQualityScores(fields[10], minimumBaseQuality);
			if (passingBases == -1) Misc.printErrAndExit("\nERROR: problem fetching base qualities for "+sam+"\n");
			numberPassingBaseScores+= passingBases;
			numberTotalBaseScores+= fields[10].length();
		}
	}


	public void makeSamReaders(){
		samReaders = new SamReader[bamFiles.length];
		for (int i=0; i< samReaders.length; i++) {
			samReaders[i] = readerFactory.open(bamFiles[i]);
		}
	}

	public void closeSamReaders(){
		try {
			for (int i=0; i< samReaders.length; i++) samReaders[i].close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}





	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SamAlignmentExtractor(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		useqArguments = IO.fetchUSeqVersion()+" "+Misc.stringArrayToString(args, " ");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': bamFiles = IO.extractFiles(args[++i], ".bam"); break;
					case 'r': bedFile = new File(args[++i]); break;
					case 'b': saveFile = new File(args[++i]); break;
					case 'p': printCoverageStats = true; break;
					case 'q': minimumMappingQuality = Integer.parseInt(args[++i]); break;
					case 'l': alignmentScoreThreshold = Integer.parseInt(args[++i]); break;
					case 'm': biggerASIsBetter = false; break;
					case 'c': calcExtraStats = true; break;
					case 'u': minimumBaseQuality = Integer.parseInt(args[++i]); calcExtraStats = true; break;
					case 'n': hits = new HashSet<String>(); break;
					case 'i': minimumReadDepth = Integer.parseInt(args[++i]); break;
					case 'x': maximumReadDepth = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for bam files
		if (bamFiles == null || bamFiles.length == 0) Misc.printErrAndExit("\nError: cannot find any treatment xxx.bam files?\n");

		//look for bai index files
		lookForBaiIndexes(bamFiles, false);

		//look for bed
		if (bedFile == null || bedFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your bed file?\n");
		chromRegions = Bed.parseBedFile(bedFile, true, false);
		
		//look for save file, can be null
		if (saveFile == null || saveFile.getName().endsWith(".bam") == false) Misc.printErrAndExit("\nError: Cannot find your save file or it doesn't end in .bam\n");

	}	

	/**Looks for xxx.bam.bai and xxx.bai for each bamFile, prints error and exits if missing.*/
	public static void lookForBaiIndexes (File[] bamFiles, boolean onlyResetLastModifiedDate){
		for (File f: bamFiles){
			File index = new File (f+".bai");
			if (index.exists() == false){
				int len = f.toString().length() - 3;
				index = new File(f.toString().substring(0, len) + "bai");
				if (onlyResetLastModifiedDate == false && index.exists() == false) Misc.printErrAndExit("\nError: failed to find a xxx.bai index file for -> "+f);
			}
			//reset date?
			if (index.exists() && index.lastModified() < f.lastModified()) index.setLastModified(f.lastModified()+1);
		}
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Sam Alignment Extractor: June 2015                    **\n" +
				"**************************************************************************************\n" +
				"Given a bed file containing regions of interest, parses the intersecting alignments.\n"+
				"Provides a variety of options for filtering the alignments and calculating QC\n"+
				"statistics.\n"+

				"Required Options:\n"+
				"-a Alignment file or directory containing xxx.bam files with their associated\n" +
				"       xxx.bai indexs sorted by coordinate. Multiple files are merged.\n" +
				"-r A bed file (chr, start, stop,...) of regions to intersect, full path, see,\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1\n"+
				"-b Provide a bam file path for saving extracted alignments, must end in .bam\n"+
				
				"\nDefault Options:\n"+
				"-p Print per region coverage stats to stdout.\n"+
				"-i Minimum read depth, defaults to 1\n"+
				"-x Maximum read depth, defaults to unlimited\n"+
				"-q Miminum mapping quality, defaults to 0, no filtering, recommend 13.\n"+
				"-l Alignment score threshold, defaults to no filtering. \n"+
				"-m Smaller alignment scores are better (novo), defaults to bigger are better (bwa).\n"+
				"-n Save alignments that don't intersect the regions of interest. Memory intensive!\n"+
				"-c Calculate extra aligment statistics (frac on target, QS>=20), time intensive.\n"+
				"-u Minimum base quality score for extra alignment calculations, defaults to 20.\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/SamAlignmentExtractor -q 13 -l 128 -a\n" +
				"      /Data/ExonCaptureAlignmentsX1/ -r /Data/SNPCalls/9484X1Calls.bed.gz \n"+
				"      -b /Data/Extracted/onTarget.bam -c -u 30\n" +

		"**************************************************************************************\n");

	}


}
