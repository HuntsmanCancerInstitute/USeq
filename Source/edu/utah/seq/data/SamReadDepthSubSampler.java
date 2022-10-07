package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;

import util.bio.annotation.Bed;
import util.gen.*;
import java.util.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.ValidationStringency;
import edu.utah.seq.analysis.DefinedRegionDifferentialSeq;
import edu.utah.seq.data.sam.ComparatorSamAlignmentName;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.useq.data.RegionScoreText;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class SamReadDepthSubSampler{
	//user defined fields
	private File bamFile;
	private File bedRegionsFile;
	private float minimumPosteriorProbability = 13;
	private float maximumAlignmentScore = 300;
	private String adapter = "chrAdapt";
	private String phiX = "chrPhiX";
	private int targetReadDepth = -1;
	private boolean keepPairsTogether = false;
	
	//internal fields
	private Gzipper out = null;
	private File workingGzippedFile = null;
	private File workingBamFile = null;
	private HashMap<String, RegionScoreText[]> chrRegions;
	private String[] toSkip = {"Un", "random", "decoy", "alt", "HLA"};
	
	//alignment counts for sam files
	private long numberAlignmentsFailingQualityScore = 0;
	private long numberAlignmentsFailingAlignmentScore = 0;
	private long numberControlAlignments = 0;
	private long numberAlignmentsFailingQC = 0;
	private long numberAlignmentsUnmapped = 0;
	private long numberPassingAlignments = 0;
	private long numberAlignmentsSaved = 0;
	


	//constructors
	public SamReadDepthSubSampler(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		
		System.out.println("Parsing alignments...");
		parseWorkingSAMFile();
		
		System.out.println("\n\nSorting...");
		new PicardSortSam(workingGzippedFile, workingBamFile, true);

		printStats();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}
	

	private void printStats(){
		//Alignment filtering stats
		double total = numberAlignmentsFailingQualityScore + numberAlignmentsFailingAlignmentScore + numberControlAlignments + numberAlignmentsFailingQC + numberAlignmentsUnmapped + numberPassingAlignments;
		System.out.println("\nFiltering statistics for "+(int)total+" alignments:");
		System.out.println(numberAlignmentsFailingQualityScore +"\tFailed mapping quality score ("+minimumPosteriorProbability+")");
		System.out.println(numberAlignmentsFailingAlignmentScore +"\tFailed alignment score ("+maximumAlignmentScore+")");
		System.out.println(numberControlAlignments +"\tAligned to chrPhiX* or chrAdapt*");
		System.out.println(numberAlignmentsFailingQC +"\tFailed vendor QC");
		System.out.println(numberAlignmentsUnmapped +"\tAre unmapped");
		System.out.println(numberPassingAlignments +"\tPassed filters ("+Num.formatPercentOneFraction(((double)numberPassingAlignments)/ total)+")");
		System.out.println(numberAlignmentsSaved +"\tPassed read depth filter ("+Num.formatPercentOneFraction(((double)numberAlignmentsSaved)/ (double)numberPassingAlignments)+") and saved");
		System.out.println();
	}


	public void parseWorkingSAMFile(){
		try{
			//make reader
			SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
			SAMFileHeader sfh = reader.getFileHeader();
			SortOrder so = sfh.getSortOrder();
			if (so.equals(SortOrder.coordinate) == false) {
				throw new Exception("\n\tSAM file does not appear to be sorted by coordinate? Aborting.\n");
			}
			
			//make gzipper and write out header
			String fileName = Misc.removeExtension(bamFile.getName()) + "RDSub"+targetReadDepth;
			workingGzippedFile = new File(bamFile.getParentFile(), fileName+".sam.gz");
			workingBamFile = new File(bamFile.getParentFile(), fileName+".bam");
			workingGzippedFile.deleteOnExit();
			out = new Gzipper(workingGzippedFile);
			SAMFileHeader h = reader.getFileHeader();
			out.println(h.getSAMString().trim());
			
			//for each chrom
			for (String chr: chrRegions.keySet()){
				if (passChrCheck(chr)){
					//for each region
					RegionScoreText[] regions = chrRegions.get(chr);
					IO.p(chr+" ");
					for (RegionScoreText region: regions){
						SAMRecordIterator iterator = reader.queryOverlapping(chr, region.getStart(), region.getStop());
						if (iterator.hasNext()) loadParse(iterator);
						iterator.close();
					}
				}
			}
			reader.close();
			out.close();
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing alignment file\n");
		}
	}

	private boolean passChrCheck(String chr) {
		for (int i=0; i< toSkip.length; i++) {
			if (chr.contains(toSkip[i])){
				IO.p("Skipping-"+chr+" ");
				return false;
			}
		}
		return true;
	}


	private void loadParse(SAMRecordIterator iterator) throws IOException {
		int numBadLines = 0;
		int min = Integer.MAX_VALUE; 
		int max = Integer.MIN_VALUE;
		ArrayList<SamAlignment> samAL = new ArrayList<SamAlignment>();
		
		//check and load records
		while (iterator.hasNext()){
			SAMRecord samRecord = iterator.next();
			//this is a bit inefficient but gives absolute control on the sam data
			SamAlignment sa;
			String samLine = samRecord.getSAMString().trim();
			try {
				sa = new SamAlignment(samLine, false);
			} catch (Exception e) {
				System.out.println("\nSkipping malformed sam alignment ->\n"+samRecord.getSAMString()+"\n"+e.getMessage());
				if (numBadLines++ > 1000) Misc.printErrAndExit("\nAboring: too many malformed SAM alignments.\n");
				continue;
			}
			//is it aligned?
			if (sa.isUnmapped()) {
				numberAlignmentsUnmapped++;
				continue;
			}
			//does it pass the vendor qc?
			if (sa.failedQC()) {
				numberAlignmentsFailingQC++;
				continue;
			}
			//skip phiX, adapter, lambda
			String chr = sa.getReferenceSequence();
			if (chr.startsWith(phiX) || chr.startsWith(adapter)) {
				numberControlAlignments++;
				continue;
			}
			//does it pass the scores threshold?
			if (sa.getAlignmentScore() > maximumAlignmentScore) {
				numberAlignmentsFailingAlignmentScore++;
				continue;
			}
			if (sa.getMappingQuality() < minimumPosteriorProbability) {
				numberAlignmentsFailingQualityScore++;
				continue;
			}
			//increment counter, save, and check start and stop
			numberPassingAlignments++;
			samAL.add(sa);
			if (sa.getUnclippedStart() < min) min = sa.getUnclippedStart();
			int end = sa.countLengthOfAlignment() + sa.getPosition();
			if (end > max) max = end;	
		}
		processChromAlignments(samAL, min, max);
	}


	private void processChromAlignments(ArrayList<SamAlignment> samAL, int min, int max) throws IOException {
		if (samAL.size() == 0) return;
		SamAlignment[] sams = new SamAlignment[samAL.size()];
		samAL.toArray(sams);
		samAL = null;
		
		//sort by name
		ComparatorSamAlignmentName comp = new ComparatorSamAlignmentName();
		Arrays.sort(sams, comp);
		
		//build read groups
		ArrayList<ReadGroup> rgsAL = new ArrayList<ReadGroup>();
		String name = sams[0].getName();
		ReadGroup rg = new ReadGroup(sams[0]);
		rgsAL.add(rg);
		for (int i=1; i< sams.length; i++){
			String testName = sams[i].getName();
			if (testName.equals(name)) rg.alignments.add(sams[i]);
			else {
				name = testName;
				rg = new ReadGroup(sams[i]);
				rgsAL.add(rg);
			}
		}
		sams = null;
		
		//randomize the read groups
		ReadGroup[] rgs = new ReadGroup[rgsAL.size()];
		rgsAL.toArray(rgs);
		Misc.randomize(rgs, System.currentTimeMillis());
		rgsAL = null;
		
		//make short[] to hold depth info
		short[] rd = new short[1000+max-min];
		
		//for each read group assess whether any of the alignments cover a base with < requested depth
		if (keepPairsTogether){
			for (int i=0; i< rgs.length; i++){
				ArrayList<SamAlignment> sas = rgs[i].alignments;
				//check first alignment in the group, print all if pass
				for (SamAlignment sa : sas){
					if (checkSetDepth(sa, rd, min)){
						for (SamAlignment s : sas) out.println(s.getUnmodifiedSamRecord());
						numberAlignmentsSaved += sas.size();
					}
					break;
				}
			}
		}
		else {
			for (int i=0; i< rgs.length; i++){
				ArrayList<SamAlignment> sas = rgs[i].alignments;
				//check each alignment 
				for (SamAlignment sa : sas){
					if (checkSetDepth(sa, rd, min)){
						out.println(sa.getUnmodifiedSamRecord());
						numberAlignmentsSaved++;
					}
				}
			}
		}
		//call garbage collection
		rd = null;
		rgs = null;
		System.gc(); System.gc(); System.gc();
	}

	
	private boolean checkSetDepth(SamAlignment sa, short[] rd, int minPos) {
		int zeroPos = sa.getPosition()- minPos;
		//fetch alignment blocks
		ArrayList<int[]> blocks = DefinedRegionDifferentialSeq.fetchAlignmentBlocks(sa.getCigar(), zeroPos);
		if (incrementCounters(blocks, rd)){
			//increment read depth counters
			for (int[] startStop: blocks){
				for (int i= startStop[0]; i< startStop[1]; i++) rd[i]++;
			}
			return true;
		}
		return false;
	}
	
	private boolean incrementCounters(ArrayList<int[]> blocks, short[] rd){
		for (int[] startStop: blocks){
			for (int i= startStop[0]; i< startStop[1]; i++){
				if (rd[i] >= targetReadDepth) {
					return false;
				}
			}
		}
		return true;
	}

	private class ReadGroup{
		ArrayList<SamAlignment> alignments = new ArrayList<SamAlignment>();
		public ReadGroup(SamAlignment sam){
			alignments.add(sam);
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SamReadDepthSubSampler(args);
	}		

	/**This method will process each argument and assign new varibles*/
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
					case 'x': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'q': minimumPosteriorProbability = Float.parseFloat(args[++i]); break;
					case 't': targetReadDepth = Integer.parseInt(args[++i]); break;
					case 'p': keepPairsTogether = true; break;
					case 'b': bedRegionsFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (bamFile == null || bamFile.exists() == false) Misc.printExit("\nError: cannot find your xxx.bam file!\n");
		if (targetReadDepth < 0) Misc.printErrAndExit("\nPlease provide a target read depth.\n");
		if (bedRegionsFile == null || bedRegionsFile.exists() == false) Misc.printErrAndExit("\nError: cannot find your bed regions file? "+bedRegionsFile);
		
		chrRegions = Bed.parseBedFile(bedRegionsFile, true, false);
	
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Sam Read Depth Sub Sampler: Oct 2022                     **\n" +
				"**************************************************************************************\n" +
				"Filters, randomizes, and subsamples a coordinate sorted bam alignment file to a target\n"+
				"base level read depth over each of the provided regions. Depending on the gaps between\n"+
				"your regions, you may need to remove duplicate lines, e.g. 'sort -u body.sam > uni.sam'\n"+

				"\nOptions:\n"+
				"-a Alignment xxx.bam file, coordinate sorted with index.\n" +
				"-b Bed file of regions to subsample (e.g. use Sam2USeq -c 1 -b hg38StdChrms.bed)\n"+
				"-t Target read depth.\n"+
				
				"\nDefault Options:\n"+
				"-p Keep read groups together.  Causes greater variation in depth.\n"+
				"-x Maximum alignment score. Defaults to 300, smaller numbers are more stringent.\n"+
				"-q Minimum mapping quality score. Defaults to 13, bigger numbers are more stringent.\n" +
				"      For RNASeq data, set this to 0.\n" +

				"\nExample: java -Xmx25G -jar pathToUSeq/Apps/SamReadDepthSubSampler -x 240 -q 20 -a\n" +
				"      /Novo/Run7/full.bam -n 100 -b regionsWith1PlusAlignment.bed.gz \n\n" +


		"**************************************************************************************\n");

	}

}
