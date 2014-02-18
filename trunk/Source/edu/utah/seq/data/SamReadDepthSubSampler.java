package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import edu.utah.seq.analysis.DefinedRegionDifferentialSeq;
import edu.utah.seq.data.sam.ComparatorSamAlignmentName;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.data.sam.SamAlignment;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class SamReadDepthSubSampler{
	//user defined fields
	private File[] samFiles;
	private float minimumPosteriorProbability = 13;
	private float maximumAlignmentScore = 300;
	private String adapter = "chrAdapt";
	private String phiX = "chrPhiX";
	private int targetReadDepth = -1;
	private boolean keepPairsTogether = false;
	
	//internal fields
	private String samHeader = null;
	private Gzipper out = null;
	private File workingGzippedFile = null;
	private File workingBamFile = null;
	
	//alignment counts for sam files
	private long numberAlignmentsFailingQualityScore = 0;
	private long numberAlignmentsFailingAlignmentScore = 0;
	private long numberControlAlignments = 0;
	private long numberAlignmentsFailingQC = 0;
	private long numberAlignmentsUnmapped = 0;
	private long numberPassingAlignments = 0;
	private File workingFile = null;
	private long numberAlignmentsSaved = 0;
	


	//constructors
	public SamReadDepthSubSampler(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		System.out.println("Parsing alignments...");
		
			//for each file, filter, split into numberChunks chunks
			for (int i=0; i< samFiles.length; i++){
				workingFile = samFiles[i];
				System.out.print("\t"+workingFile.getName());
				//parse it
				if (parseWorkingSAMFile() == false) Misc.printExit("\n\tError: failed to parse, aborting.\n");
				//sort results by coordinate
				System.out.println(" Sorting");
				new PicardSortSam(workingGzippedFile, workingBamFile);
				//print stats
				printStats();
			}
				
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
		
		//reset counters
		samHeader = null;
		numberAlignmentsFailingQualityScore = 0;
		numberAlignmentsFailingAlignmentScore = 0;
		numberControlAlignments = 0;
		numberAlignmentsFailingQC = 0;
		numberAlignmentsUnmapped = 0;
		numberPassingAlignments = 0;
		numberAlignmentsSaved = 0;
	}


	public boolean parseWorkingSAMFile(){
		try{
			//make reader
			SAMFileReader reader = new SAMFileReader(workingFile);	
			reader.setValidationStringency(ValidationStringency.SILENT);
			SAMFileHeader sfh = reader.getFileHeader();
			SortOrder so = sfh.getSortOrder();
			if (so.equals(SortOrder.coordinate) == false) {
				System.err.println("\n\tSAM file does not appear to be sorted by coordinate? Aborting.\n");
				reader.close();
				System.exit(1);
			}
			
			//make gzipper and write out header
			String fileName = Misc.removeExtension(workingFile.getName()) + "RDSub"+targetReadDepth;
			workingGzippedFile = new File(workingFile.getParentFile(), fileName+".sam.gz");
			workingBamFile = new File(workingFile.getParentFile(), fileName+".bam");
			workingGzippedFile.deleteOnExit();
			out = new Gzipper(workingGzippedFile);
			samHeader = sfh.getTextHeader().trim();
			out.println(samHeader);
			
			SAMRecordIterator iterator = reader.iterator();
			String chrom = null;
			int numBadLines = 0;
			int min = Integer.MAX_VALUE; 
			int max = Integer.MIN_VALUE;
			ArrayList<SamAlignment> samAL = new ArrayList<SamAlignment>();
			
			//for each record, load all of a particular chrom and process
			while (iterator.hasNext()){
				SAMRecord samRecord = iterator.next();
				//print status blip
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
				
				//new chrom?
				if (chrom == null || chrom.equals(sa.getReferenceSequence()) == false){
					//process old
					processChromAlignments(samAL, min, max);
					//reset
					min = Integer.MAX_VALUE; 
					max = Integer.MIN_VALUE;
					chrom = sa.getReferenceSequence();
					samAL = new ArrayList<SamAlignment>();
					System.out.print(" "+chrom);
					
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
			//process last
			System.out.print(" "+chrom);
			processChromAlignments(samAL, min, max);
			
			//cleanup
			out.close();
			reader.close();
			
		} catch (Exception e){
			System.err.println("\nError parsing Novoalign file or writing split binary chromosome files.\nToo many open files? Too many chromosomes? " +
			"If so then login as root and set the default higher using the ulimit command (e.g. ulimit -n 10000)\n");
			e.printStackTrace();
			return false;
		}
		return true;
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
				//check each alignment in the group, print all if pass
				for (SamAlignment sa : sas){
					if (checkSetDepth(sa, rd, min)){
						for (SamAlignment s : sas) out.println(s.getUnmodifiedSamRecord());
						numberAlignmentsSaved += sas.size();
						break;
					}
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
		boolean increment = false;
		//scan bases
		for (int[] startStop: blocks){
			for (int i= startStop[0]; i< startStop[1]; i++){
				if (rd[i] < targetReadDepth) {
					increment = true;
					break;
				}
			}
		}
		if (increment){
			//increment read depth counters
			for (int[] startStop: blocks){
				for (int i= startStop[0]; i< startStop[1]; i++) rd[i]++;
			}
			return true;
		}
		
		return false;
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
					case 'a': samFiles = IO.extractFiles(new File(args[++i]),".bam"); break;
					case 'x': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'q': minimumPosteriorProbability = Float.parseFloat(args[++i]); break;
					case 't': targetReadDepth = Integer.parseInt(args[++i]); break;
					case 'p': keepPairsTogether = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (samFiles == null || samFiles.length ==0 || samFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.bam file(s)!\n");
		if (targetReadDepth < 0) Misc.printErrAndExit("\nPlease provide a target read depth.\n");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              SamReadDepthSubSampler: Feb 2014                    **\n" +
				"**************************************************************************************\n" +
				"Filters, randomizes, subsamples each coordinate sorted bam alignment file to a target\n"+
				"base level read depth. Useful for reducing extreem read depths over localized areas.\n"+

				"\nOptions:\n"+
				"-a Alignment file or directory containing coordinate sorted xxx.bam files. Each is \n" +
				"      processed independently.\n" +
				"-t Target read depth.\n"+
				
				"\nDefault Options:\n"+
				"-p Keep read groups together.  Causes greater variation in depth.\n"+
				"-x Maximum alignment score. Defaults to 300, smaller numbers are more stringent.\n"+
				"-q Minimum mapping quality score. Defaults to 13, bigger numbers are more stringent.\n" +
				"      For RNASeq data, set this to 0.\n" +

				"\nExample: java -Xmx25G -jar pathToUSeq/Apps/SamReadDepthSubSampler -x 240 -q 20 -a\n" +
				"      /Novo/Run7/ -n 100 \n\n" +


		"**************************************************************************************\n");

	}

}
