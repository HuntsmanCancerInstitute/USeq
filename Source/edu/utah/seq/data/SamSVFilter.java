package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;

import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.*;

import java.util.*;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.ValidationStringency;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.data.sam.SamAlignmentFlags;
import edu.utah.seq.useq.data.IntersectingRegions;
import edu.utah.seq.useq.data.Region;

/**
 * @author david.nix@hci.utah.edu/elainegee 
 **/
public class SamSVFilter{
	//user defined fields
	private File[] samFiles;
	private File saveDirectory;
	private File bedFile;
	private float minimumPosteriorProbability = 5;
	private float maximumAlignmentScore = 1000;
	private int minimumSoftMaskedBases = 10;
	private String[] chromToSkip = new String[]{"chrAdap", "chrPhi", "chrM", "random", "chrUn"};
	private boolean sortFinal = true;
	private boolean deSoftMask = false;

	//internal fields
	private String samHeader;
	private Gzipper gzipOutPassSpan;
	private Gzipper gzipOutPassSoft;
	private Gzipper gzipOutPassSingle;
	private Gzipper gzipOutFail;
	private HashMap<String, IntersectingRegions> regions;
	private TreeMap<String, Integer> links = new TreeMap<String, Integer>();

	//alignment counts for sam files
	private long numberAlignments = 0;
	private long numberAlignmentsFailingQualityScore = 0;
	private long numberAlignmentsFailingAlignmentScore = 0;
	private long numberChrSkippedAlignments = 0;
	private long numberAlignmentsFailingQC = 0;
	private long numberAlignmentsUnmapped = 0;
	
	//pass the filters
	//private long numberPassingButNoIntersection = 0;
	//private long numberPassingButFailingSoftMasking = 0;
	//private long numberPassingButAllInternalToOneRegion = 0;
	//private long numberPassingAlignments = 0;
	private long numberFailingNotPaired = 0;
	private long numberFailingBothMissingRegions = 0;
	private long numberFailingSoftMaskThreshold = 0;
	private long numberPassingSameRegionSoftMask = 0;
	private long numberPassingDiffRegions = 0;
	private long numberPassingSingleRegion = 0;

	//constructors
	public SamSVFilter(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);

		//load region booleans
		System.out.println("Loading regions into fast lookup arrays...");
		regions = IntersectingRegions.parseRegions(bedFile);

		//for each file, filter, and sort
		System.out.println("Processing alignments...");
		try {
			for (int i=0; i< samFiles.length; i++){
				System.out.print("\n"+samFiles[i]);

				//make gzippers
				String name = Misc.removeExtension(samFiles[i].getName());
				File span = new File (saveDirectory, "passSpan_"+name+".sam.gz");
				File soft = new File (saveDirectory, "passSoft_"+name+".sam.gz");
				File single = new File (saveDirectory, "passSingle_"+name+".sam.gz");
				File fail = new File (saveDirectory, "fail_"+name+".sam.gz");
				gzipOutPassSpan = new Gzipper(span);
				gzipOutPassSoft = new Gzipper(soft);
				gzipOutFail = new Gzipper(fail);
				gzipOutPassSingle = new Gzipper(single);

				//parse it
				if (parseWorkingSAMFile(samFiles[i]) == false) {
					gzipOutPassSpan.close();
					gzipOutPassSoft.close();
					gzipOutPassSingle.close();
					gzipOutFail.close();
					span.delete();
					soft.delete();
					fail.delete();
					Misc.printErrAndExit("\n\tERROR: failed to parse, skipping.\n");
				}
				
				//close gzippers
				gzipOutPassSpan.close();
				gzipOutPassSoft.close();
				gzipOutPassSingle.close();
				gzipOutFail.close();
				
				//sort
				if (sortFinal && samHeader != null){
					System.out.println("Sorting...");
					File bam = new File (saveDirectory, "passSpan_"+name+".bam");
					//sort and convert to BAM
					new PicardSortSam (span, bam);
					span.delete();
					bam = new File (saveDirectory, "passSoft_"+name+".bam");
					//sort and convert to BAM
					new PicardSortSam (soft, bam);
					soft.delete();
					bam = new File (saveDirectory, "passSingle_"+name+".bam");
					//sort and convert to BAM
					new PicardSortSam (single, bam);
					single.delete();
					bam = new File (saveDirectory, "fail_"+name+".bam");
					//sort and convert to BAM
					new PicardSortSam (fail, bam);
					fail.delete();
				}
				
				//write out stats
				printStats();
				
				//print hits?
				printLinks();
			}
		} catch (Exception e) {
			e.printStackTrace();
		} 

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}

	private void printLinks() {
		System.out.println("\nTarget region link counts:");
		for (String coor : links.keySet()){
			System.out.println(coor+"\t"+links.get(coor));
		}
		
	}

	public boolean parseWorkingSAMFile(File workingFile){
		try{
			zeroStats();
			
			//make reader, check sort order
			SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(workingFile);
			SAMRecordIterator iterator = reader.iterator();
			SAMFileHeader sfh = reader.getFileHeader();
			SortOrder so = sfh.getSortOrder();
			if (so.equals(SortOrder.coordinate)) {
				System.err.println("\n\tSAM file appears to be sorted by coordinate. Analysis requires file to be sorted by name.");
				reader.close();
				return false;
			}
			
			//write header to gzippers
			samHeader = sfh.getTextHeader().trim();
			gzipOutFail.println(samHeader);
			gzipOutPassSpan.println(samHeader);
			gzipOutPassSingle.println(samHeader);
			gzipOutPassSoft.println(samHeader);
			
			int counter =0;
			int numBadLines = 0;
			ArrayList<SamAlignment> samAL = new ArrayList<SamAlignment>();
			String name = "";
			
			//for each record
			while (iterator.hasNext()){
				SAMRecord samRecord = iterator.next();

				//print status blip
				if (++counter == 1000000){
					System.out.print(".");
					counter = 0;
				}

				//this is a bit inefficient but gives absolute control on the sam data
				SamAlignment sa;
				String samLine = samRecord.getSAMString().trim();
				try {
					sa = new SamAlignment(samLine, false);
				} catch (Exception e) {
					System.out.println("\nSkipping malformed sam alignment ->\n"+samRecord.getSAMString()+"\n"+e.getMessage());
					if (numBadLines++ > 1000) Misc.printErrAndExit("\nAboring: too many malformed SAM alignments.\n");
					gzipOutFail.println(samLine);
					continue;
				}
				numberAlignments++; 
				boolean pass = true;
				
				//is it aligned?
				if (sa.isUnmapped()) {
					numberAlignmentsUnmapped++;
					pass = false;
				}
				//does it pass the vendor qc?
				if (sa.failedQC()) {
					numberAlignmentsFailingQC++;
					pass = false;
				}
				//does it pass the scores threshold?
				if (sa.getAlignmentScore() > maximumAlignmentScore) {
					numberAlignmentsFailingAlignmentScore++;
					pass = false;
				}
				if (sa.getMappingQuality() < minimumPosteriorProbability) {
					numberAlignmentsFailingQualityScore++;
					pass = false;
				}
				//skip phiX, adapter, lambda, chrM, etc...
				String chr = sa.getReferenceSequence();
				if (isChromToSkip(chr)){
					numberChrSkippedAlignments++;
					pass = false;
				}
				
				if (pass == false){
					gzipOutFail.println(samLine);
					continue;
				}
				
				//process if the current SAM record is a new alignment pair as now both pairs are loaded into samAL
				if (name.equals(sa.getName()) == false){
					processBlock(samAL); 
					name = sa.getName();
					samAL.clear();
				}
				samAL.add(sa);
			}
			reader.close();
			System.out.println();
			
			//process last block
			if (samAL.size() !=0) processBlock(samAL);
			
		} catch (Exception e){
			System.err.println("\nError parsing alignment file.\n");
			e.printStackTrace();
			return false;
		}
		return true;
	}
	

	
	/**This saves alignment pairs where each map to a different region or where both do and one has significant external soft masking*/
	private void processBlock(ArrayList<SamAlignment> samAL) throws IOException {
		boolean print = true;
		SamAlignment one = null;
		SamAlignment two = null;
		String regionOne = null;
		String regionTwo = null;
		Gzipper out  = null;
		//just one?
		int num = samAL.size();
		if (num != 2) {
			print = false;
			numberFailingNotPaired+=num;
		}
		else {
			one = samAL.get(0);
			two = samAL.get(1);
			regionOne =  fetchIntersectingRegion(one);
			regionTwo =  fetchIntersectingRegion(two);
			//both miss?
			if (regionOne == null && regionTwo == null) {
				print = false;
				numberFailingBothMissingRegions+=2;
			}
			//both hit?
			else if (regionOne != null && regionTwo != null) {
				//same region? check if pass soft masking
				if (regionOne.equals(regionTwo)) {
					if (failsMinSoftMasking (one, two)) {
						print = false;
						numberFailingSoftMaskThreshold+=2;
					}
					else {
						numberPassingSameRegionSoftMask+=2;
						out = gzipOutPassSoft;
					}
				}
				//nope diff regions
				else {
					numberPassingDiffRegions+=2;
					out = gzipOutPassSpan;
				}
			}
			//must be just one hit
			else {
				numberPassingSingleRegion+=2;
				out = gzipOutPassSingle;
			}
		}
		
		//print them
		if (print){
			if(deSoftMask) {
				one.deSoftMaskCigar();
				two.deSoftMaskCigar();
			}
			
			// Write to file only if a minimum of one read is in the bed file
			if (!(regionOne == null && regionTwo == null)) {
				out.println(one.getUnmodifiedSamRecord());
				out.println(two.getUnmodifiedSamRecord());
			}
			
			
			//make link
			String key;
			if (regionOne == null || regionTwo == null){
				if (regionOne != null) key = regionOne;
				else key = regionTwo;
			}
			else {
				int comp = regionOne.compareTo(regionTwo);
				if (comp < 0) key = regionOne+"_"+regionTwo;
				else if (comp > 0) key = regionTwo+"_"+regionOne;
				else key = regionOne;
			}
			int count = 1;
			
			if (links.containsKey(key)) count = links.get(key) + 1;
			links.put(key, count); 
		}
		else {
			for (SamAlignment sam : samAL) gzipOutFail.println(sam.getUnmodifiedSamRecord());
		}
		
		
	}

	private boolean failsMinSoftMasking(SamAlignment one, SamAlignment two) {
		//which is left? (lower chr position)
		SamAlignment left;
		SamAlignment right;
		if (one.getPosition() <= two.getPosition()) {
			left = one;
			right = two;
		}
		else {
			left = two;
			right = one;
		}
		//Check softclipping on left of pair
		int num_left = left.countLengthOfSidedSoftMaskedBases(true); //if true, count soft clipping on left (else right)
		int num_right = right.countLengthOfSidedSoftMaskedBases(false);
		if (num_left >= minimumSoftMaskedBases || num_right >= minimumSoftMaskedBases) return false;
		// Check softclipping on right of pair
		num_left = left.countLengthOfSidedSoftMaskedBases(false); 
		num_right = right.countLengthOfSidedSoftMaskedBases(true);
		if (num_left >= minimumSoftMaskedBases || num_right >= minimumSoftMaskedBases) return false;
		return true;
	}

	/**Returns null if unmapped, or 0 or more than one intersecting regions are found.*/
	private String fetchIntersectingRegion(SamAlignment sam){
		IntersectingRegions ir = null;
		ir = regions.get(sam.getReferenceSequence());
		if (ir == null) return null;
		int start = sam.getUnclippedStart();
		HashSet<Integer> hits = ir.fetchIntersectingIndexes(start, sam.countLengthOfAlignment()+ start);
		if (hits== null || hits.size() != 1) return null;
		Region region = ir.getRegions()[hits.iterator().next()];
		return sam.getReferenceSequence() +":"+ region.getStart()+"-"+region.getStop();
	}
		

	public void printStats(){
		//Alignment filtering stats
		System.out.println("\nProcessing statistics for number of alignments...");
		System.out.println(numberAlignments+"\tTotal SAM Records");
		System.out.println("\t"+numberAlignmentsFailingQualityScore +"\tFailed mapping quality score ("+minimumPosteriorProbability+").");
		System.out.println("\t"+numberAlignmentsFailingAlignmentScore +"\tFailed alignment score ("+maximumAlignmentScore+").");
		System.out.println("\t"+numberChrSkippedAlignments +"\tAligned to control chromosomes.");
		System.out.println("\t"+numberAlignmentsFailingQC +"\tFailed vendor QC.");
		System.out.println("\t"+numberAlignmentsUnmapped +"\tAre unmapped.");
		System.out.println("\t"+numberFailingNotPaired +"\tNot paired.");
		System.out.println("\t"+numberFailingBothMissingRegions +"\tNeither of pair align to a target region.");
		System.out.println("\t"+numberPassingSingleRegion +"\tOne of pair aligns to a target region, the other somewhere else. Single.");
		System.out.println("\t"+numberFailingSoftMaskThreshold +"\tBoth align to the same region but fail soft masking threshold  ("+minimumSoftMaskedBases+").");
		System.out.println("\t"+numberPassingSameRegionSoftMask +"\tBoth align to the same region yet pass soft masking threshold ("+minimumSoftMaskedBases+"). Soft.");
		System.out.println("\t"+numberPassingDiffRegions +"\tBoth align to different target regions. Span.");
	}
	
	private void zeroStats() {
		numberAlignments = 0;
		numberAlignmentsFailingQualityScore = 0;
		numberAlignmentsFailingAlignmentScore = 0;
		numberChrSkippedAlignments = 0;
		numberAlignmentsFailingQC = 0;
		numberAlignmentsUnmapped = 0;
		numberFailingNotPaired = 0;
		numberFailingBothMissingRegions = 0;
		numberFailingSoftMaskThreshold = 0;
		numberPassingSameRegionSoftMask = 0;
		numberPassingDiffRegions = 0;
		numberPassingSingleRegion = 0;
		samHeader = null;
		links.clear();
	}

	private boolean isChromToSkip(String chr) {
		for (int i=0; i< chromToSkip.length; i++){
			if (chromToSkip[i].contains(chr)) return true;
		}
		return false;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		/*
		System.out.println("##### JVM Heap statistics [GB] #####");
		int gb = 1024*1024*1024;
		Runtime runtime = Runtime.getRuntime(); 
        //Print used memory
        System.out.println("Used Memory: "
            + (runtime.totalMemory() - runtime.freeMemory()) / gb);
	 
        //Print free memory
        System.out.println("Free Memory: "
            + runtime.freeMemory() / gb);
	         
        //Print total available memory
        System.out.println("Total Available Memory: " + runtime.totalMemory() / gb);
	 
        //Print Maximum available memory
        System.out.println("Max Available Memory: " + runtime.maxMemory() / gb);
		System.out.println("####################################");
		*/
		
		new SamSVFilter(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		String controlChroms = null;
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': forExtraction = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'x': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'q': minimumPosteriorProbability = Float.parseFloat(args[++i]); break;
					case 'm': minimumSoftMaskedBases = Integer.parseInt(args[++i]); break;
					case 'd': sortFinal = false; break;
					case 'e': deSoftMask = true; break;
					case 'c': controlChroms = args[++i]; break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		File[][] tot = new File[4][];
		tot[0] = IO.extractFiles(forExtraction,".sam");
		tot[1] = IO.extractFiles(forExtraction,".sam.gz");
		tot[2] = IO.extractFiles(forExtraction,".sam.zip");
		tot[3] = IO.extractFiles(forExtraction,".bam");
		samFiles = IO.collapseFileArray(tot);
		if (samFiles == null || samFiles.length ==0 || samFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.sam(.zip/.gz) file(s)!\n");
		if (saveDirectory != null) saveDirectory.mkdirs();
		if (saveDirectory == null || saveDirectory.isDirectory() == false) Misc.printErrAndExit("\nPlease enter a directory to use in saving your results.\n");
		if (bedFile == null || bedFile.canRead() == false) Misc.printErrAndExit("\nPlease enter bed file to use in filtering your alignments.\n");
		
		//reset chromToSkip?
		if (controlChroms != null) chromToSkip = controlChroms.split(",");
		
		
	}	

	public static void printDocs(){
		String version = IO.fetchUSeqVersion();
		if (version.length() == 0) version = "USeq_xxx";
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Sam SV Filter: Oct 2015                            **\n" +
				"**************************************************************************************\n" +
				"Filters SAM records based on their intersection with a list of target regions for\n" +
				"structural variation analysis. Both mates of a paired alignment are kept if they align to\n" + 
				"at least one target region. These are split into those that align to different targets,\n" +
				"(span) the same target with sufficient softmasking on either the left or right side of the\n" +
				"mate pair (soft), or to one target and somewhere else outside of the bed file (single).\n"+

				"\nOptions:\n"+
				"-a Alignment file or directory containing NAME sorted SAM/BAM files. Multiple files\n"+
				"       are processed independantly. Xxx.sam(.gz/.zip) or xxx.bam are OK. Assumes only\n"+
				"       uniquely aligned reads. Remove duplicates with Picard's MarkDuplicates app.\n" +
				"-s Save directory for the results.\n"+
				"-b Bed file (tab delim: chr, start, stop, ...) of target regions interbase coordinates.\n"+

				"\nDefault Options:\n"+
				"-n Mark passing alignments as secondary. Needed for Delly with -n 30 novoalignments.\n"+
				"-d Don't coordinate sort and index alignments.\n"+
				"-x Maximum alignment score. Defaults to 1000, smaller numbers are more stringent.\n"+
				"-q Minimum mapping quality score. Defaults to 5, bigger numbers are more stringent.\n" +
				"-c Chromosomes to skip, defaults to 'chrAdap,chrPhi,chrM,random,chrUn'. Any SAM\n"+
				"       record chromosome name that contains one will be failed.\n"+ 
				"-m Minimum number of soft masked bases needed to keep paired alignments. Both must intersect\n"+
				"       a target region in the bed file, defaults to 10\n"+
				//"-e Replace S in cigar with M, only for visualization in IGB, not for downstream use.\n"+

				"\nExample: java -Xmx25G -jar pathTo/"+version+"/Apps/SamSVFilter -x 150 -q 13 -a\n" +
				"      /Novo/Run7/ -s /Novo/Run7/SSVF/ -c 'chrPhi,_random,chrUn_' \n\n" +


				"**************************************************************************************\n");

	}

}
