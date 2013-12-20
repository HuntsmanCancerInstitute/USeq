package edu.utah.seq.data;

import java.io.*;
import java.util.regex.*;

import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.*;

import java.util.*;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.useq.data.IntersectingRegions;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class SamSVFilter{
	//user defined fields
	private File[] samFiles;
	private File saveDirectory;
	private File bedFile;
	private float minimumPosteriorProbability = 10;
	private float maximumAlignmentScore = 300;
	private int minimumSoftMaskedBases = 5;
	private int maximumInsertSize = 2000;
	private String[] chromToSkip = new String[]{"chrAdap", "chrPhi", "chrM", "_random", "chrUn_"};
	private boolean sortFinal = true;
	private boolean removePairsInSameRegion = true;

	//internal fields
	private String samHeader;
	private Gzipper gzipOutPass;
	private Gzipper gzipOutFail;
	HashMap<String, IntersectingRegions> regions;

	//alignment counts for sam files
	private long numberAlignments = 0;
	private long numberAlignmentsFailingQualityScore = 0;
	private long numberAlignmentsFailingAlignmentScore = 0;
	private long numberChrSkippedAlignments = 0;
	private long numberAlignmentsFailingQC = 0;
	private long numberAlignmentsUnmapped = 0;
	private long numberProperPairedAlignments = 0;
	
	//pass the filters
	private long numberPassingButNoIntersection = 0;
	private long numberPassingButFailingSoftMasking = 0;
	private long numberPassingButAllInternalToOneRegion = 0;
	private long numberPassingAlignments = 0;

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
				File pass = new File (saveDirectory, "pass_"+name+".sam.gz");
				File fail = new File (saveDirectory, "fail_"+name+".sam.gz");
				gzipOutPass = new Gzipper(pass);
				gzipOutFail = new Gzipper(fail);

				//parse it
				if (parseWorkingSAMFile(samFiles[i]) == false) {
					gzipOutPass.close();
					gzipOutFail.close();
					pass.delete();
					fail.delete();
					Misc.printErrAndExit("\n\tERROR: failed to parse, skipping.\n");
				}
				
				//close gzippers
				gzipOutPass.close();
				gzipOutFail.close();
				
				//sort
				if (sortFinal && samHeader != null){
					System.out.println("Sorting...");
					File bam = new File (saveDirectory, "pass_"+name+".bam");
					//sort and convert to BAM
					new PicardSortSam (pass, bam);
					pass.delete();
					bam = new File (saveDirectory, "fail_"+name+".bam");
					//sort and convert to BAM
					new PicardSortSam (fail, bam);
					fail.delete();
				}
				
				//write out stats
				printStats();
			}
		} catch (Exception e) {
			e.printStackTrace();
		} 

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}

	public boolean parseWorkingSAMFile(File workingFile){
		try{
			zeroStats();
			
			//make reader, check sort order
			SAMFileReader reader = new SAMFileReader(workingFile);	
			reader.setValidationStringency(ValidationStringency.SILENT);
			SAMRecordIterator iterator = reader.iterator();
			SAMFileHeader sfh = reader.getFileHeader();
			SortOrder so = sfh.getSortOrder();
			if (so.equals(SortOrder.coordinate)) {
				System.err.println("\n\tSAM file appears to be sorted by coordinate.");
				reader.close();
				return false;
			}
			
			//write header to gzippers
			samHeader = sfh.getTextHeader().trim();
			gzipOutFail.println(samHeader);
			gzipOutPass.println(samHeader);

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
				//proper pair in close prox?
				int insertSize = Math.abs(sa.getInferredInsertSize());
				if (sa.isAProperPairedAlignment() && insertSize <= maximumInsertSize && insertSize !=0){
					numberProperPairedAlignments++;
					pass = false;
				}
				
				if (pass == false){
					gzipOutFail.println(samLine);
					continue;
				}
				
				//is it an old alignment?
				if (name.equals(sa.getName()) == false){
					processBlock (samAL);
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

	private void processBlock(ArrayList<SamAlignment> samAL) throws IOException {
		//do any of these alignments fall within a region?
		int num = samAL.size();
		IntersectingRegions ir = null;
		String chrom = "";
		
		//one alignment must intersect
		boolean intersection = false;
		HashSet<Integer> indexHits = null;
		int i;
		//scan records
		for (i=0; i< num; i++){
			SamAlignment sam = samAL.get(i);
			//fetch IR?
			if (sam.getReferenceSequence().equals(chrom) == false){
				chrom = sam.getReferenceSequence();
				ir = regions.get(chrom);
			}
			if (ir == null) continue;
			
			//any intersection?
			int start = sam.getUnclippedStart();
			indexHits = ir.fetchIntersectingIndexes(start, sam.countLengthOfAlignment()+ start);
			if (indexHits != null){
				intersection = true;
				break;
			}
		}
		
		//no intersection write all to fail
		if (intersection == false){
			for (SamAlignment sam : samAL){
				gzipOutFail.println(sam.getUnmodifiedSamRecord());
				numberPassingButNoIntersection++;
			}
			return;
		}
		
		//just one alignment, if so then must have suff soft masking to keep
		boolean softMaskGood = false;
		if (samAL.size() == 1){
			//must have soft masking to be good
			SamAlignment sam = samAL.get(0);
			int numS = sam.countLengthOfSoftMaskedBases();		
			if (numS >= minimumSoftMaskedBases) softMaskGood = true;
			else {
				gzipOutFail.println(sam.getUnmodifiedSamRecord());
				numberPassingButFailingSoftMasking++;
				return;
			}
		}
		
		//OK at least one intersects, its hits are in indexHits, its index is i
		
		//did prior alignments not hit a region?
		//do they want to not filter for same region reads?
		//only one alignment with a passing softMask?
		//already hitting multiple regions?
		if (i!=0  || removePairsInSameRegion == false || softMaskGood || indexHits.size() > 1){
			for (SamAlignment sam : samAL){
				gzipOutPass.println(sam.getUnmodifiedSamRecord());
				numberPassingAlignments++;
			}
			return;
		}
		
		//scan subsequent records for different regions
		boolean printAll = false;
		i++;
		for (; i< num; i++){
			SamAlignment sam = samAL.get(i);
			//fetch IR?
			if (sam.getReferenceSequence().equals(chrom) == false){
				//different chrom so must fall out of other region
				printAll = true;
				break;
			}
			//any intersection?
			HashSet<Integer> newHits = ir.fetchIntersectingIndexes(sam.getUnclippedStart(), sam.countLengthOfAlignment()+ sam.getUnclippedStart());
			if (newHits == null){
				//outside of any region thus print all
				printAll = true;
				break;
			}
			//ok more hits, add to original and see if change
			indexHits.addAll(newHits);
			//diff number?
			if (indexHits.size() != 1){
				//yes, added new region so print all
				printAll = true;
				break;
			}
		}
		
		if (printAll){
			for (SamAlignment sam : samAL){
				gzipOutPass.println(sam.getUnmodifiedSamRecord());
				numberPassingAlignments++;
			}
		}
		else{
			for (SamAlignment sam : samAL){
				gzipOutFail.println(sam.getUnmodifiedSamRecord());
				numberPassingButAllInternalToOneRegion++;  
			}
		}
	}
	
	public void printStats(){
		//Alignment filtering stats
		System.out.println("Processing statistics...");
		System.out.println(numberAlignments+"\tSAM Records");
		System.out.println("\t"+numberAlignmentsFailingQualityScore +"\tFailed mapping quality score ("+minimumPosteriorProbability+")");
		System.out.println("\t"+numberAlignmentsFailingAlignmentScore +"\tFailed alignment score ("+maximumAlignmentScore+")");
		System.out.println("\t"+numberChrSkippedAlignments +"\tAligned to control chromosomes");
		System.out.println("\t"+numberAlignmentsFailingQC +"\tFailed vendor QC");
		System.out.println("\t"+numberAlignmentsUnmapped +"\tAre unmapped");
		System.out.println("\t"+numberProperPairedAlignments +"\tProper paired");
		
		long numPassingFilters = numberPassingButNoIntersection+ numberPassingButFailingSoftMasking+ numberPassingButAllInternalToOneRegion+ numberPassingAlignments;
		System.out.println(numPassingFilters+"\tSAM Records remaining");
		System.out.println("\t"+numberPassingButNoIntersection +"\tNo intersection with regions");
		System.out.println("\t"+numberPassingButFailingSoftMasking +"\tNon paired failing minimum softmask threshold");
		System.out.println("\t"+numberPassingButAllInternalToOneRegion +"\tInternal to one region");
		System.out.println(numberPassingAlignments +"\tPassing all filters");
	}
	
	private void zeroStats() {
		numberAlignments = 0;
		numberAlignmentsFailingQualityScore = 0;
		numberAlignmentsFailingAlignmentScore = 0;
		numberChrSkippedAlignments = 0;
		numberAlignmentsFailingQC = 0;
		numberAlignmentsUnmapped = 0;
		numberPassingButNoIntersection = 0;
		numberPassingButFailingSoftMasking = 0;
		numberPassingButAllInternalToOneRegion = 0;
		numberPassingAlignments = 0;
		samHeader = null;
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
					case 'i': maximumInsertSize = Integer.parseInt(args[++i]); break;
					case 'd': sortFinal = false; break;
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
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Sam SV Filter: Dec 2013                          **\n" +
				"**************************************************************************************\n" +
				"Filters SAM records based on their intersection with a list of target regions for\n"+
				"structural variation analysis. Paired alignments are kept if one pair lands in a target"+
				"region and the other aligns somewhere outside that target.  Non paired alignments\n"+
				"must intersect a target region and pass the minimum soft masked base filter. Proper\n"+
				"pairs are removed.\n" +

				"\nOptions:\n"+
				"-a Alignment file or directory containing NAME sorted SAM/BAM files. Multiple files\n"+
				"       are processed independantly. Xxx.sam(.gz/.zip) or xxx.bam are OK. Assumes only\n"+
				"       uniquely aligned reads.\n" +
				"-s Save directory for the results.\n"+
				"-b Bed file (tab delim: chr, start, stop, ...) interbase coordinates.\n"+

				"\nDefault Options:\n"+
				"-d Don't coordinate sort and index alignments.\n"+
				"-x Maximum alignment score. Defaults to 300, smaller numbers are more stringent.\n"+
				"-q Minimum mapping quality score. Defaults to 10, bigger numbers are more stringent.\n" +
				"-c Chromosomes to skip, defaults to 'chrAdap,chrPhi,chrM,_random,chrUn_'. Any SAM\n"+
				"       record chromosome name that contains one will be failed.\n"+ 
				"-m Minimum soft masked bases for keeping non paired alignments, defaults to 5\n"+
				"-i Maximum insert size for scoring proper paired alignments, defaults to 2000\n"+

				"\nExample: java -Xmx25G -jar pathToUSeq/Apps/SamSVFilter -x 150 -q 13 -a\n" +
				"      /Novo/Run7/ -s /Novo/Run7/Filt/ -c 'chrPhi,_random,chrUn_' \n\n" +


				"**************************************************************************************\n");

	}

}
