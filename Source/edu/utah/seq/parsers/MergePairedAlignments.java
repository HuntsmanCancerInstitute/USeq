
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import edu.utah.seq.data.sam.*;

/**
 * @author david.nix@hci.utah.edu 
 * Similar to MergePairedSamAlignments but cut down in options for speed. Now works on coordinate and query name sorted files.
 **/
public class MergePairedAlignments{

	//user defined fields
	private File alignmentFile;
	private File saveFile;
	private float maximumAlignmentScore = 120;
	private float minimumMappingQualityScore = 13;
	private boolean secondPairReverseStrand = false;
	private boolean skipMergingPairs = false;
	private int minimumDiffQualScore = 3;
	private double minimumFractionInFrameMismatch = 0.05;
	private int maximumProperPairDistanceForMerging = 5000;
	private String namePhiXChromosome = "chrPhiX";
	private String nameAdapterChromosome = "chrAdap";
	private boolean printPairedStatsAndHistogram = false;
	private int numberAlignmentsToLoad = 100000;

	//counters for initial filtering
	private int numberAlignments = 0;
	private int numberUnmapped = 0;
	private int numberFailingVendorQC = 0;
	private int numberFailingAlignmentScore = 0;
	private int numberFailingMappingQualityScore = 0;
	private int numberAdapter = 0;
	private int numberPhiX = 0;
	private int numberDuplicates = 0;
	private int numberPassingAlignments = 0;

	//for trimming/ merging paired data
	ArrayList<SamAlignment> firstPairs = new ArrayList<SamAlignment>();
	ArrayList<SamAlignment> secondPairs = new ArrayList<SamAlignment>();
	private int numberPrintedAlignments = 0;
	private double numberOverlappingBases = 0;
	private double numberNonOverlappingBases = 0;
	private double numberMergedPairs = 0;
	private double numberFailedMergedPairs = 0;
	private int numberNonPairedAlignments = 0;
	private int numberNonProperPairedAlignments = 0;
	private int numberUnmappedMatePairedAlignments = 0;
	private int numberAlignmentsMissingPair;
	private int numberPairsFailingChrDistStrand;
	private boolean crossCheckMateCoordinates = true;
	private int numberPairsFailingMateCrossCoordinateCheck;
	private int numberRepeatAlignmentsLackingMate;

	//internal fields
	private Gzipper samOut = null;
	LinkedHashMap<String, ArrayList<SamAlignment>> alignments = new LinkedHashMap<String, ArrayList<SamAlignment>>();
	public static Pattern CIGAR_SUB = Pattern.compile("(\\d+)([MDIN])");
	public static Pattern CIGAR_BAD = Pattern.compile(".*[^\\dMDIN].*");
	private Histogram insertSize = new Histogram(0,2001,400);
	private String programArguments;
	private SAMRecordIterator samIterator;
	
private HashSet<String> processed = new HashSet<String>();

	//constructors
	public MergePairedAlignments(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void doWork() throws IOException{
		//make gzipper for writing good alignments
		samOut = new Gzipper(saveFile);

		System.out.print("\nParsing");
		parseAlignmentFile();
			
		//close the gzipper
		samOut.close();

		//stats
		System.out.println("\n\nStats (some flags aren't set so be suspicious of zero read catagories):");
		System.out.println(numberAlignments+"\tTotal # alignments from sam/bam file");
		System.out.println(format(numberUnmapped, numberAlignments)+"Unmapped reads");
		System.out.println(format(numberFailingVendorQC,numberAlignments)+"Alignments failing vendor/ platform QC");
		System.out.println(format(numberFailingAlignmentScore,numberAlignments)+"Alignments failing alignment score ("+(int)maximumAlignmentScore+")");
		System.out.println(format(numberFailingMappingQualityScore,numberAlignments)+"Alignments failing mapping quality score ("+(int)minimumMappingQualityScore+")");
		if (numberAdapter!=0) System.out.println(format(numberAdapter,numberAlignments)+"Adapter alignments");
		System.out.println(format(numberPhiX,numberAlignments)+"PhiX alignments");
		System.out.println(format(numberDuplicates,numberAlignments)+"Duplicate alignments");
		System.out.println(format(numberPassingAlignments,numberAlignments)+"Alignments passing all filters");
		System.out.println();
		System.out.println("Paired alignment stats:");  
		if (printPairedStatsAndHistogram){
			System.out.println(numberNonPairedAlignments+"\t# Non paired alignments");
			System.out.println(numberNonProperPairedAlignments+"\t# Non proper paired alignments");
			System.out.println(numberUnmappedMatePairedAlignments+"\t# Non mapped mate paired alignments");
			System.out.println(numberAlignmentsMissingPair+"\t# Alignments missing mate paired alignment");
			System.out.println(numberPairsFailingChrDistStrand+"\t# Paired alignments failing chromosome, distance, or strand check");
			System.out.println(numberPairsFailingMateCrossCoordinateCheck+"\t# Paired alignments failing mate pair cross coordinate check");
			System.out.println(numberRepeatAlignmentsLackingMate+"\t# Repeat alignments lacking a mate");
			System.out.println((int)numberFailedMergedPairs+"\t# Proper paired alignments that could not be unambiguously merged");
		}
		System.out.println(format(numberMergedPairs, (numberMergedPairs+numberFailedMergedPairs))+"Proper paired alignments that were merged");
		double totalBases = numberNonOverlappingBases + numberOverlappingBases;
		double fractionOverlap = numberOverlappingBases/totalBases;
		String fractionString = Num.formatNumber(fractionOverlap, 4);
		System.out.println(fractionString+"\tFraction overlapping bases in paired alignments");		
		//histogram
		if (printPairedStatsAndHistogram && insertSize.getTotalBinCounts() !=0){
			System.out.println("\nMapped genomic insert length distribution for merged paired alignments:\n");
			insertSize.setSkipZeroBins(true);
			insertSize.setTrimLabelsToSingleInteger(true);
			insertSize.printScaledHistogram();
		}
		System.out.println("\n"+format(numberPrintedAlignments,numberAlignments)+"Alignments written to SAM/BAM file.");
	}

	public String format (double stat, double total){
		double frac = stat/total;
		return (int)stat +"\t"+Num.formatNumber(frac, 4)+"\t";
	}

	public void parseAlignmentFile(){
		String line = null;
		try {
			//get reader 
			SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			SamReader samReader = factory.open(alignmentFile);

			//add header to results file
			String[] header = samReader.getFileHeader().getTextHeader().split("\\n");
			for (String h: header) samOut.println(h);
			samOut.println("@PG\tID:MergePairedSamAlignments\tCL: "+programArguments);

			samIterator = samReader.iterator();
			ArrayList<SamAlignment> alignmentsToMerge = null;
			int dotCounter = 0;
			
			//load up hash with lots of records
			if (loadAlignments(numberAlignmentsToLoad) == false){
				//no more to add so process all
				processAlignmentHash();
				return;
			}
		
			//more to work with so....
			while (true) {
				if (++dotCounter > 100000){
					System.out.print(".");
					dotCounter = 0;
				}
				
				//grab first set
				alignmentsToMerge = fetchRemoveAlignmentBlock();
				
				//process
				filterPrintAlignments(alignmentsToMerge);
				
				//add more to hash
				if (loadAlignments(alignmentsToMerge.size()) == false){
					//no more to add so process all
					processAlignmentHash();
					break;
				}
			}
			
			samIterator.close();
			samReader.close();
		} catch (Exception e) {
			System.err.println("\nError processing -> "+line);
			e.printStackTrace();	
		} 
	}
	
	/**Runs through the hash, then nulls it.*/
	private void processAlignmentHash() {
		for (String readName: alignments.keySet()){
			ArrayList<SamAlignment> alignmentsToMerge = alignments.get(readName);
			filterPrintAlignments(alignmentsToMerge);
		}
		alignments = null;
	}

	/**Fetches ArrayList of alignments, then removes it from the hash.
	 * Returns null if no more.*/
	private ArrayList<SamAlignment> fetchRemoveAlignmentBlock(){
		try {
			String readName = alignments.keySet().iterator().next();
			ArrayList<SamAlignment> alignmentsToMerge = alignments.get(readName);
			alignments.remove(readName);			
			return alignmentsToMerge;
		} catch (Exception e){
			return null;
		}
	}
	
	/**Adds alignments to the hash, returns whether they were all loaded.*/
	private boolean loadAlignments(int numberToLoad) {
		String line = null;
		int counter = 0;
		try {
			while (samIterator.hasNext()) {
				//fetch and parse sam record
				SAMRecord sam = samIterator.next();
				line = sam.getSAMString().trim();
				SamAlignment sa = new SamAlignment(line, false);
				numberAlignments++;
			
				//check it
				if (checkSamAlignment(sa, line) == false) continue;
				
				//fetch or make ArrayList and add it
				ArrayList<SamAlignment> al = alignments.get(sa.getName());
				if (al == null){
					al = new ArrayList<SamAlignment>();
					alignments.put(sa.getName(), al);	
				}
				al.add(sa);		

				//check to see if the numberToLoad has been met.
				counter++;
				if (numberToLoad == counter) return true;
			}
		} catch (Exception e){
			System.err.println("\nError processing -> "+line);
			e.printStackTrace();
		}
		return false;
	}

	/**Checks a bunch of flags and scores to see if alignment should be saved for attempted merging.*/
	public boolean checkSamAlignment(SamAlignment sa, String line) throws IOException{
		boolean isGood = true;
		//is it aligned?
		if (sa.isUnmapped()) {
			numberUnmapped++;
			isGood = false;
		}
		//does it pass the vendor qc?
		if (sa.failedQC()) {
			numberFailingVendorQC++;
			isGood = false;
		}
		//increment phiX
		boolean firstIsPhiX = sa.getReferenceSequence().startsWith(namePhiXChromosome);
		if (firstIsPhiX) {
			numberPhiX++;
			isGood = false;
		}
		//skip adapter
		boolean firstIsAdapt = sa.getReferenceSequence().startsWith(nameAdapterChromosome);
		if (firstIsAdapt){
			numberAdapter++;
			isGood = false;
		}
		//does it pass the scores threshold?
		int alignmentScore = sa.getAlignmentScore();
		if (alignmentScore != Integer.MIN_VALUE){
			if (alignmentScore > maximumAlignmentScore){
				numberFailingAlignmentScore++;
				isGood = false;
			}
		}
		//check mapping quality for genomic match reads?
		if (minimumMappingQualityScore !=0){
			if (sa.getMappingQuality() < minimumMappingQualityScore){
				numberFailingMappingQualityScore++;
				isGood = false;
			}
		}
		//duplicate?
		if (sa.isADuplicate()){
			numberDuplicates++;
			isGood = false;
		}
		
		//did any fail?
		if (isGood == false) return false;
		
		//OK, it passes individual checks, increment counter, now check for pairing.
		numberPassingAlignments++;

		//modify second read if phiX or adapter; some paired reads have a mate hitting the control chroms
		SamAlignmentFlags saf = null;
		if (sa.isPartOfAPairedAlignment() && sa.isMateUnMapped() == false){
			if ( sa.getMateReferenceSequence().startsWith(namePhiXChromosome) || sa.getMateReferenceSequence().startsWith(nameAdapterChromosome)) {
				saf = new SamAlignmentFlags(sa.getFlags());
				saf.setMateUnMapped(true);
				saf.setaProperPairedAlignment(false);
				sa.setUnMappedMate();
				sa.setFlags(saf.getFlags());
			}
		}
		
		//is it part of a pair?
		if (sa.isPartOfAPairedAlignment()== false){
			numberNonPairedAlignments++;
			samOut.println(sa);
			numberPrintedAlignments++;
			return false;
		}

		//is it not part of a proper pair?
		if (sa.isAProperPairedAlignment() == false){
			numberNonProperPairedAlignments++;
			samOut.println(sa);
			numberPrintedAlignments++;
			return false;
		}
		//is mate unmapped
		if (sa.isMateUnMapped()){
			numberUnmappedMatePairedAlignments++;
			samOut.println(sa);
			numberPrintedAlignments++;
			return false;
		}
		//OK looks good to add for potential pairing
		return true;
	}

	/**Takes a block of alignments all originating from the same fragment.  Attempts to identify pairs and merge them.*/
	public void filterPrintAlignments(ArrayList<SamAlignment> al){
		try {

			//Split by first and second pairs
			firstPairs.clear();
			secondPairs.clear();
			for (SamAlignment sam : al) {
				if (sam.isFirstPair()) firstPairs.add(sam);
				else if (sam.isSecondPair()) secondPairs.add(sam);
				else Misc.printErrAndExit("\nError: seeing non first second pairing! Aborting.\n");
			}
			int numFirstPairs = firstPairs.size();
			int numSecondPairs = secondPairs.size();

			//missing partners due to filtering?
			if (numFirstPairs == 0 || numSecondPairs == 0){
				for (SamAlignment sam : al) {
					numberAlignmentsMissingPair++;
					samOut.println(sam);
					numberPrintedAlignments++;
				}
				return;
			}

			//just one of each, these should be pairs, so do detail diagnostics
			if (numFirstPairs == 1 && numSecondPairs == 1){
				SamAlignment first = firstPairs.get(0);
				SamAlignment second = secondPairs.get(0);
				processPair(first, second);
				return;
			}

			//nope looks like we've got repeat matches. Is it OK to check mate coordinates.
			if (crossCheckMateCoordinates == false){
				Misc.printErrAndExit("\n\nError, aborting, found repeat match pairs yet you've indicated you don't " +
						"want to cross check mate coordinates.  This is the only way to correctly pair repeat " +
						"matches. Restart after resolving -> "+firstPairs.get(0).getName()+"\n");
			}
			//walk through each first looking at seconds to attempt to merge
			SamAlignment[] firsts = fetchSamAlignmentArray(firstPairs);
			SamAlignment[] seconds = fetchSamAlignmentArray(secondPairs);
			//for each first
			for (int i=0; i< firsts.length; i++){
				//for each second
				for (int j=0; j< seconds.length; j++){
					//null? already processed?
					if (seconds[j] == null) continue;
					//check mate info
					if (testMateChrPosition(firsts[i], seconds[j])){
						//OK looks like they are a pair
						//now test if it is OK to attempt a merge
						if (testChrDistStrnd(firsts[i], seconds[j])){
							//attempt a merge
							mergeAndScorePair(firsts[i], seconds[j]);
						}
						//failed distance or strand so not a proper pair
						else {
							samOut.println(firsts[i]);
							samOut.println(seconds[j]);
							numberPrintedAlignments+=2;
							numberPairsFailingChrDistStrand++;
						}
						//regardless of out come set to null so they aren't processed again
						firsts[i] = null;
						seconds[j] = null;
						break;
					}
					//nope, not a mate, check next second
				}


			}
			//OK now walk through arrays and print any that are not null
			for (SamAlignment s: firsts){
				if (s !=null){
					numberRepeatAlignmentsLackingMate++;
					samOut.println(s);
					numberPrintedAlignments++;
				}
			}
			for (SamAlignment s: seconds){
				if (s !=null){
					numberRepeatAlignmentsLackingMate++;
					samOut.println(s);
					numberPrintedAlignments++;
				}
			}




		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem printing alignment block!?\n");
		}
	}

	private SamAlignment[] fetchSamAlignmentArray(ArrayList<SamAlignment> al){
		SamAlignment[] sam = new SamAlignment[al.size()];
		al.toArray(sam);
		return sam;
	}


	private void processPair(SamAlignment first, SamAlignment second) throws IOException {
		//check if the pair is from the same chromosome, acceptable distance, and proper strand
		if (testChrDistStrnd(first, second)){
			//cross validate mate information? this is often messed up especially if from the STP app.
			if (crossCheckMateCoordinates){
				if (testMateChrPosition(first, second)){
					//attempt a merge
					mergeAndScorePair(first, second);
				}
				else {
					samOut.println(first);
					samOut.println(second);
					numberPrintedAlignments+=2;
					numberPairsFailingMateCrossCoordinateCheck++;
				}
			}
			else {
				//attempt a merge
				mergeAndScorePair(first, second);
			}
		}
		else {
			samOut.println(first);
			samOut.println(second);
			numberPrintedAlignments+=2;
			numberPairsFailingChrDistStrand++;
		}

	}

	/**Attempts to merge a proper paired alignment.  Increments counters and sends examples to the gzipper.*/
	private void mergeAndScorePair(SamAlignment first, SamAlignment second) throws IOException{
		//collect string rep since the merge method will modify the SamAlignment while merging so if it fails you can output the unmodified SamAlignment
		String firstSamString = first.toString();
		String secondSamString = second.toString();
		
		//do they want to skip
		if (skipMergingPairs){
			samOut.println(firstSamString);
			samOut.println(secondSamString);
			numberPrintedAlignments+=2;
			return;
		}
		
		SamAlignment mergedSam = mergePairedAlignments(first, second);
		
		//failed to merge?
		if (mergedSam!=null) {
			numberMergedPairs++;
			samOut.println(mergedSam);
			numberPrintedAlignments++;
			//if (first.getCigar().contains("D") || second.getCigar().contains("D") || first.getCigar().contains("I") || second.getCigar().contains("I")){
			//System.out.println(first.getName()+"\t"+first.getReferenceSequence()+"\t"+first.getPosition()+"\t"+first.getCigar()+"\t"+second.getCigar());
			//}
		}
		else {
			samOut.println(firstSamString);
			samOut.println(secondSamString);
			numberPrintedAlignments+=2;
			numberFailedMergedPairs++;
		}
	}

	/**Checks chrom and position of mate info.  This info is often messed up so it probably is a good idea to watch this outcome.*/
	private boolean testMateChrPosition(SamAlignment first, SamAlignment second) {
		//check mate position
		if (first.getMatePosition() != second.getPosition()) return false;
		if (second.getMatePosition() != first.getPosition()) return false;

		//check chromosome
		String mateChrom = first.getMateReferenceSequence();

		if (mateChrom.equals("=") == false && mateChrom.equals(first.getReferenceSequence()) == false) return false;
		mateChrom = second.getMateReferenceSequence();
		if (mateChrom.equals("=") == false && mateChrom.equals(second.getReferenceSequence()) == false) return false;
		return true;
	}


	/**Checks chrom, max distance, correct strand.  Required before attempting a merge!  Otherwise you'll be merging improper pairs.*/
	private boolean testChrDistStrnd(SamAlignment first, SamAlignment second) {
		if (first.getReferenceSequence().equals(second.getReferenceSequence())) {
			//within acceptable distance
			int diff = Math.abs(first.getPosition()- second.getPosition());
			if (diff < maximumProperPairDistanceForMerging){
				//correct strand?
				if (secondPairReverseStrand){
					if (first.isReverseStrand() != second.isReverseStrand()) return false;
				}
				else if (first.isReverseStrand() == second.isReverseStrand()) return false;
			}
			else return false;
		}
		else return false;

		return true;
	}

	/**Attempts to merge alignments. Doesn't check if proper pairs!  Returns null if it cannot. This modifies the input SamAlignments so print first before calling*/
	private SamAlignment mergePairedAlignments(SamAlignment first, SamAlignment second) {
		
		//trim them of soft clipped info
		first.trimMaskingOfReadToFitAlignment();
		second.trimMaskingOfReadToFitAlignment();

		//look for bad CIGARs
		if (CIGAR_BAD.matcher(first.getCigar()).matches()) Misc.printErrAndExit("\nError: unsupported cigar string! See -> "+first.toString()+"\n");
		if (CIGAR_BAD.matcher(second.getCigar()).matches()) Misc.printErrAndExit("\nError: unsupported cigar string! See -> "+second.toString()+"\n");

		//order left and right
		SamAlignment left = first;
		SamAlignment right = second;
		if (first.getPosition() > second.getPosition()) {
			right = first;
			left = second;
		}

		//System.out.println("Name "+left.getName());

		//fetch genomic space coordinates
		int startLeft = left.getPosition();
		int stopLeft = startLeft + countLengthOfCigar(left.getCigar());
		int startRight = right.getPosition();
		int stopRight = startRight + countLengthOfCigar(right.getCigar());
		int stop = stopRight;
		if (stopLeft > stop) stop = stopLeft;

		//any Is in left that precede the start of right?
		int numAdders = countIs(left.getCigar(), startRight-startLeft);

		//make arrays to hold sequence and qualities in cigar space
		int size = numAdders + stop-startLeft;

		SamLayout leftLayout = new SamLayout(size);
		SamLayout rightLayout = new SamLayout(size);

		//layout data
		leftLayout.layoutCigar(startLeft, left);
		rightLayout.layoutCigar(startLeft-numAdders, right);
		
		//increment overlap and insert size
		int[] overNonOver = SamLayout.countOverlappingBases(leftLayout, rightLayout);
		numberOverlappingBases+= overNonOver[0];
		numberNonOverlappingBases+= overNonOver[1];
		//set insert length
		insertSize.count(overNonOver[2]);

		//merge layouts, modifies original layouts so print first if you want to see em before mods.
		SamLayout mergedSamLayout = SamLayout.mergeLayouts(leftLayout, rightLayout, minimumDiffQualScore, minimumFractionInFrameMismatch);

		//System.out.println("MergedLayout");
		//mergedSamLayout.print();

		if (mergedSamLayout == null) {
			//add failed merge tag
			left.addMergeTag(false);
			right.addMergeTag(false);
			return null;
		}
		else {
			//make merged
			SamAlignment mergedSam = makeSamAlignment(first.isReverseStrand(), left, right, mergedSamLayout, startLeft);
			return mergedSam;
		}

	}

	public static SamAlignment makeSamAlignment(boolean isReverseStrand, SamAlignment first, SamAlignment second, SamLayout merged, int position){
		SamAlignment mergedSam = new SamAlignment();
		//<QNAME>
		mergedSam.setName(first.getName());
		//<FLAG>
		SamAlignmentFlags saf = new SamAlignmentFlags();
		saf.setReverseStrand(isReverseStrand);
		mergedSam.setFlags(saf.getFlags());
		//<RNAME>
		mergedSam.setReferenceSequence(first.getReferenceSequence());
		//<POS>
		mergedSam.setPosition(position);
		//<MAPQ>, bigger better
		int mqF = first.getMappingQuality();
		int mqS = second.getMappingQuality();
		if (mqF > mqS) mergedSam.setMappingQuality(mqF);
		else mergedSam.setMappingQuality(mqS);
		//<CIGAR>
		mergedSam.setCigar(merged.fetchCigar());
		//<MRNM> <MPOS> <ISIZE>
		mergedSam.setUnMappedMate();
		//<SEQ> <QUAL>
		merged.setSequenceAndQualities(mergedSam);
		/////tags, setting to read with better as score.
		//alternative score, smaller better
		int asF = first.getAlignmentScore();
		int asS = second.getAlignmentScore();
		if (asF != Integer.MIN_VALUE && asS != Integer.MIN_VALUE){
			if (asF < asS) mergedSam.setTags(first.getTags());
			else mergedSam.setTags(second.getTags());
		}
		else mergedSam.setTags(first.getTags());
		//add merged tag
		mergedSam.addMergeTag(true);
		return mergedSam;
	}

	/**Counts the number Is relative to the stop.*/
	public static int countIs (String cigar, int stop){
		if (stop == 0) return 0;
		int length = 0;
		int numIs = 0;
		//for each M D I or N block  MDIN
		Matcher mat = CIGAR_SUB.matcher(cigar);
		while (mat.find()){
			//pass
			int num = Integer.parseInt(mat.group(1));
			//is call I 
			if (mat.group(2).equals("I")){			
				for (int i=0; i< num; i++){
					length++;	
					numIs++;
					if (length >= stop) break;
					//numIs++;
				}
			}
			//nope just addit
			else {
				length += num;
				if (length >= stop) break;
			}

		}	
		return numIs;
	}

	/**Counts the number bases in the cigar string. Only counts M D I and N.*/
	public static int countLengthOfCigar (String cigar){
		int length = 0;
		//for each M D I or N block  MDIN
		Matcher mat = CIGAR_SUB.matcher(cigar);
		while (mat.find()){
			length += Integer.parseInt(mat.group(1));
		}
		return length;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MergePairedAlignments(args);
	}		


	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		programArguments = useqVersion+" "+Misc.stringArrayToString(args, " ");
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': alignmentFile = new File(args[++i]); break;
					case 's': saveFile = new File(args[++i]); break;
					case 'm': crossCheckMateCoordinates = false; break;
					case 'r': secondPairReverseStrand = true; break;
					case 'k': skipMergingPairs = true; break;
					case 'p': printPairedStatsAndHistogram = true; break;
					case 'd': maximumProperPairDistanceForMerging = Integer.parseInt(args[++i]); break;
					case 'a': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'q': minimumMappingQualityScore = Float.parseFloat(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (alignmentFile == null || alignmentFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find your xxx.bam/.sam(.gz OK) results file!\n");

		//check save file
		if (saveFile != null){
			//end in sam.gz?
			if (saveFile.getName().endsWith(".sam.gz") == false) Misc.printErrAndExit("\nError: your save file must end with the xxx.sam.gz extension.\n");
			try {
				saveFile.createNewFile();
			} catch (IOException e) {
				e.printStackTrace();
			}
			if (saveFile.canWrite() == false) Misc.printErrAndExit("\nError: cannot create or modify your indicated save file -> "+saveFile);
		}
		else {
			String saveFileString = Misc.removeExtension(alignmentFile.getName());
			String mq = ((int)minimumMappingQualityScore)+"MQ";
			saveFile = new File (alignmentFile.getParentFile(), saveFileString+"_"+mq+(int)maximumAlignmentScore+"AS.sam.gz");
		}
		//print info
		System.out.println(maximumAlignmentScore+ "\tMaximum alignment score (AS).");
		System.out.println(minimumMappingQualityScore+ "\tMinimum mapping quality score (MQ).");
		System.out.println(secondPairReverseStrand +"\tSecond read pair's strand has been reversed.");
		System.out.println(crossCheckMateCoordinates +"\tCross check read mate coordinates.");
		System.out.println(maximumProperPairDistanceForMerging +"\tMaximum bp distance for merging paired alignments.");
		if (skipMergingPairs) System.out.println(skipMergingPairs +"\tSkip merging paired alignments, useful for testing effect of merging on downstream analysis.");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Merge Paired Alignments: Nov 2014                     **\n" +
				"**************************************************************************************\n" +
				"Merges proper paired alignments that pass a variety of checks and thresholds. Only\n" +
				"unambiguous pairs will be merged. Increases base calling accuracy in overlap and helps\n" +
				"avoid non-independent variant observations and other double counting issues. Identical\n" +
				"overlapping bases are assigned the higher quality scores. Disagreements are resolved\n" +
				"toward the higher quality base. If too close in quality, then the quality is set to 0.\n" +

				"\nOptions:\n"+
				"-f Path to a xxx.bam/sam.gz file containing paired alignments, sort order irrelevant.\n" +

				"\nDefault Options:\n"+
				"-s Path to a xxx.sam.gz file for saving results, defaults to that inferred by -f\n"+
				"-a Maximum alignment score (AS:i: tag). Defaults to 120, smaller numbers are more\n" +
				"      stringent. Approx 30pts per mismatch for novoalignments.\n"+
				"-q Minimum mapping quality score, defaults to 13, larger numbers are more stringent.\n" +
				"      Set to 0 if processing splice junction indexed RNASeq data.\n"+
				"-r The second paired alignment's strand has been reversed. Defaults to not reversed.\n" +
				"-d Maximum acceptible base pair distance for merging, defaults to 5000.\n"+
				"-m Don't cross check read mate coordinates, needed for merging repeat matches. Defaults\n" +
				"      to checking.\n"+
				"-k Skip merging paired alignments. Defaults to merging. Useful for testing effect of\n" +
				"      merging on downstream analysis.\n"+
				"-p Print paired alignment statistics and insert size histogram.\n"+

				"\nExample: java -Xmx20G -jar pathToUSeq/Apps/MergePairedBamAlignments -f /Bams/ms.bam\n" +
				"     -p -s /Bams/MergedPairs/ms.mergedPairs.sam.gz -d 10000 \n\n" +

		"**************************************************************************************\n");

	}

	public int getNumberPassingAlignments() {
		return numberPassingAlignments;
	}

	public int getMinimumDiffQualScore() {
		return minimumDiffQualScore;
	}

	public double getMinimumFractionInFrameMismatch() {
		return minimumFractionInFrameMismatch;
	}	

}
