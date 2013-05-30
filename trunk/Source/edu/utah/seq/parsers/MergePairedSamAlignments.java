
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;

import util.gen.*;
import java.util.*;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;
import edu.utah.seq.data.sam.*;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class MergePairedSamAlignments{

	//user defined fields
	private File[] dataFiles;
	private File saveFile;
	private File outputFile;
	private float maximumAlignmentScore = 120;
	private float minimumMappingQualityScore = 13;
	private boolean secondPairReverseStrand = false;
	private boolean removeControlAlignments = false;
	private boolean skipMergingPairs = false;
	private boolean onlyMergeOverlappingAlignments = true;
	private int minimumDiffQualScore = 3;
	private double minimumFractionInFrameMismatch = 0.05;
	private int maximumProperPairDistanceForMerging = 5000;

	//counters for initial filtering
	private int numberAlignments = 0;
	private int numberUnmapped = 0;
	private int numberFailingVendorQC = 0;
	private int numberPassingAlignments = 0;
	private int numberFailingAlignmentScore = 0;
	private int numberFailingMappingQualityScore = 0;
	private int numberAdapter = 0;
	private int numberPhiX = 0;

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
	private PrintWriter samOut;
	private Gzipper failedSamOut = null;
	private LinkedHashSet<String> samHeader = new LinkedHashSet<String>();
	public static Pattern CIGAR_SUB = Pattern.compile("(\\d+)([MDIN])");
	public static Pattern CIGAR_BAD = Pattern.compile(".*[^\\dMDIN].*");
	private Histogram insertSize = new Histogram(0,2001,400);
	private String programArguments;
	private boolean alignmentsOverlap;

	//constructors
	public MergePairedSamAlignments(String[] args){
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
		//make print writer
		outputFile = new File(saveFile+"_temp");
		samOut = new PrintWriter( new FileWriter (outputFile));

		//make gzipper for failed examples
		String name = Misc.removeExtension(saveFile.getName());
		File failedReadOutputFile = new File(saveFile.getParentFile(), name+"_UnMappedPoorScore.sam.gz");
		failedSamOut = new Gzipper(failedReadOutputFile);

		//for each file, parse and save to disk	
		System.out.println("\nParsing, filtering, and merging SAM/BAM files...");
		for (int i=0; i< dataFiles.length; i++){
			System.out.print("\t"+dataFiles[i].getName());
			if (dataFiles[i].getName().endsWith(".bam")) parseBamFile(dataFiles[i]);
			else parseTextFile(dataFiles[i]); 
			System.out.println();
		}

		//close the writers
		samOut.close();
		failedSamOut.close();

		//add header and output results, this deletes the outputFile too
		if (saveFile.getName().endsWith(".sam")){
			System.out.println("\nAdding SAM header and gzip compressing xxx.sam file...");
			addHeaderAndCompress();
		}
		else {
			System.out.println("\nAdding SAM header, sorting, and writing bam output with Picard's SortSam...");
			addHeaderAndSort();
		}

		//stats
		double fractionPassing = ((double)numberPassingAlignments)/((double)numberAlignments);
		System.out.println("\nStats (some flags aren't set so be suspicious of zero read catagories):\n");
		System.out.println("\t"+numberAlignments+"\tTotal # alignments from sam/bam file");
		System.out.println("\t"+numberPassingAlignments+"\tAlignments passing individual read filters ("+Num.formatPercentOneFraction(fractionPassing)+")");
		System.out.println("\t\t"+numberUnmapped+"\t# Unmapped reads");
		System.out.println("\t\t"+numberFailingVendorQC+"\t# Alignments failing vendor/ platform QC");
		System.out.println("\t\t"+numberFailingAlignmentScore+"\t# Alignments failing alignment score ("+(int)maximumAlignmentScore+")");
		System.out.println("\t\t"+numberFailingMappingQualityScore+"\t# Alignments failing mapping quality score ("+(int)minimumMappingQualityScore+")");
		System.out.println("\t\t"+numberAdapter+"\t# Adapter alignments");
		System.out.println("\t\t"+numberPhiX+"\t# PhiX alignments");
		System.out.println();
		System.out.println("Paired alignment stats:\n");
		System.out.println("\t\t"+numberNonPairedAlignments+"\t# Non paired alignments");
		System.out.println("\t\t"+numberNonProperPairedAlignments+"\t# Non proper paired alignments");
		System.out.println("\t\t"+numberUnmappedMatePairedAlignments+"\t# Non mapped mate paired alignments");
		System.out.println("\t\t"+numberAlignmentsMissingPair+"\t# Alignments missing mate paired alignment");
		System.out.println("\t\t"+numberPairsFailingChrDistStrand+"\t# Paired alignments failing chromosome, distance, or strand check");
		System.out.println("\t\t"+numberPairsFailingMateCrossCoordinateCheck+"\t# Paired alignments failing mate pair cross coordinate check");
		System.out.println("\t\t"+numberRepeatAlignmentsLackingMate+"\t# Repeat alignments lacking a mate");
		System.out.println("\t\t"+(int)numberFailedMergedPairs+"\t# Proper paired alignments that could not be unambiguously merged");
		System.out.println("\t\t"+(int)numberMergedPairs+"\t# Proper paired alignments that were merged");
		double totalBases = numberNonOverlappingBases + numberOverlappingBases;
		double fractionOverlap = numberOverlappingBases/totalBases;
		String fractionString = Num.formatNumber(fractionOverlap, 4);
		System.out.println("\t\t\t"+fractionString+"\tFraction overlapping bases in paired alignments");		
		//histogram
		if (insertSize.getTotalBinCounts() !=0){
			System.out.println("\nMapped genomic insert length distribution for merged paired alignments:\n");
			insertSize.setSkipZeroBins(true);
			insertSize.setTrimLabelsToSingleInteger(true);
			insertSize.printScaledHistogram();
		}

		System.out.println("\n\t"+numberPrintedAlignments+"\t# Alignments written to SAM/BAM file.");


	}

	public boolean parseTextFile(File samFile){
		BufferedReader in = null;
		int numBadLines = 0;
		try {
			in = IO.fetchBufferedReader(samFile);
			String line;
			String priorReadName = "";
			boolean priorSet = false;
			ArrayList<SamAlignment> alignmentsToSave = new ArrayList<SamAlignment>();
			int dotCounter = 0;
			while ((line=in.readLine())!= null) {
				if (++dotCounter > 1000000){
					System.out.print(".");
					dotCounter = 0;
				}
				line = line.trim();

				//skip blank lines
				if (line.length() == 0) continue;

				//header line?
				if (line.startsWith("@")){
					samHeader.add(line);
					continue;
				}

				SamAlignment sa;
				try {
					sa = new SamAlignment(line, false);
				} catch (MalformedSamAlignmentException e) {
					System.out.println("\nSkipping malformed sam alignment -> "+e.getMessage());
					if (numBadLines++ > 100) Misc.printErrAndExit("\nAboring: too many malformed SAM alignments.\n");
					continue;
				}
				numberAlignments++;

				if (checkSamAlignment(sa, line) == false) continue;
				
				//phiX? don't merge since this throws the per cycle error estimation
				if (sa.getReferenceSequence().startsWith("chrPhiX")){
					samOut.println(sa);
					continue;
				}

				String readName = sa.getName();

				//prior set
				if (priorSet == false){
					priorSet = true;
					priorReadName = readName;
					alignmentsToSave.add(sa);
				}
				//is it an old read name?
				else if (readName.equals(priorReadName)) alignmentsToSave.add(sa);
				//nope new read so process read alignment block
				else {
					filterPrintAlignments(alignmentsToSave);
					//set prior 
					priorReadName = readName;
					//clear 
					alignmentsToSave.clear();
					//add new
					alignmentsToSave.add(sa);
				}
			}
			//process last read alignment block
			filterPrintAlignments(alignmentsToSave);

		} catch (Exception e) {
			e.printStackTrace();
			return false;
		} finally {
			try {
				if (in != null) in.close();
			} catch (IOException e) {}
		}
		return true;
	}



	public boolean parseBamFile(File bamFile){
		SAMFileReader samReader = null;
		int numBadLines = 0;
		try {
			String line;
			String priorReadName = "";
			boolean priorSet = false;
			ArrayList<SamAlignment> alignmentsToSave = new ArrayList<SamAlignment>();
			int dotCounter = 0;

			samReader = new SAMFileReader(bamFile);

			//check sort order
			if (samReader.getFileHeader().getSortOrder().compareTo(SortOrder.coordinate) == 0){
				Misc.printErrAndExit("\nError, your bam file appears sorted by coordinate. Sort by query name and restart.\n");
			}

			//load header 
			String[] header = samReader.getFileHeader().getTextHeader().split("\\n");
			for (String h: header) samHeader.add(h);


			SAMRecordIterator it = samReader.iterator();

			while (it.hasNext()) {
				SAMRecord sam = it.next();
				line = sam.getSAMString().trim();

				if (++dotCounter > 1000000){
					System.out.print(".");
					dotCounter = 0;
				}

				SamAlignment sa;
				try {
					sa = new SamAlignment(line, false);
				} catch (MalformedSamAlignmentException e) {
					System.out.println("\nSkipping malformed sam alignment -> "+e.getMessage());
					if (numBadLines++ > 100) Misc.printErrAndExit("\nAboring: too many malformed SAM alignments.\n");
					continue;
				}
				numberAlignments++;

				if (checkSamAlignment(sa, line) == false) continue;

				String readName = sa.getName();

				//prior set
				if (priorSet == false){
					priorSet = true;
					priorReadName = readName;
					alignmentsToSave.add(sa);
				}
				//is it an old read name?
				else if (readName.equals(priorReadName)) alignmentsToSave.add(sa);
				//nope new read so process read alignment block
				else {
					filterPrintAlignments(alignmentsToSave);
					//set prior 
					priorReadName = readName;
					//clear 
					alignmentsToSave.clear();
					//add new
					alignmentsToSave.add(sa);
				}
			}
			//process last read alignment block
			filterPrintAlignments(alignmentsToSave);

		} catch (Exception e) {
			e.printStackTrace();
			return false;
		} finally {
			if (samReader != null) samReader.close();
		}
		return true;
	}

	/**Checks a bunch of flags and scores to see if alignment should be saved.*/
	public boolean checkSamAlignment(SamAlignment sa, String line) throws IOException{
		//is it aligned?
		if (sa.isUnmapped()) {
			numberUnmapped++;
			failedSamOut.println(line);
			return false;
		}

		//does it pass the vendor qc?
		if (sa.failedQC()) {
			numberFailingVendorQC++;
			failedSamOut.println(line);
			return false;
		}

		//increment phiX
		boolean firstIsPhiX = sa.getReferenceSequence().startsWith("chrPhiX");
		if (firstIsPhiX){
			numberPhiX++;
		}
		
		//skip adapter
		boolean firstIsAdapt = sa.getReferenceSequence().startsWith("chrAdapt");
		if (firstIsAdapt){
			numberAdapter++;
			return false;
		}

		//does it pass the scores threshold?
		int alignmentScore = sa.getAlignmentScore();
		if (alignmentScore != Integer.MIN_VALUE){
			if (alignmentScore > maximumAlignmentScore){
				numberFailingAlignmentScore++;
				failedSamOut.println(line);
				return false;
			}
		}

		//check mapping quality for genomic match reads?
		if (minimumMappingQualityScore !=0){
			if (sa.getMappingQuality() < minimumMappingQualityScore){
				numberFailingMappingQualityScore++;
				failedSamOut.println(line);
				return false;
			}
		}

		//modify second read if phiX or adapter; some paired reads have a mate hitting the control chroms
		SamAlignmentFlags saf = null;
		if (removeControlAlignments && sa.isPartOfAPairedAlignment() && sa.isMateUnMapped() == false){
			if ( sa.getMateReferenceSequence().startsWith("chrPhiX") || sa.getMateReferenceSequence().startsWith("chrAdapt")) {
				saf = new SamAlignmentFlags(sa.getFlags());
				saf.setMateUnMapped(true);
				saf.setaProperPairedAlignment(false);
				sa.setUnMappedMate();
				sa.setFlags(saf.getFlags());
			}
		}

		//OK, it passes individual checks, increment counter, now check for pairing.
		numberPassingAlignments++;

		//is it not part of a pair?
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
		
		//do they want to not merge non overlapping alignemnts
		if (onlyMergeOverlappingAlignments == true && alignmentsOverlap == false){
			samOut.println(firstSamString);
			samOut.println(secondSamString);
			numberPrintedAlignments+=2;
		}
		
		else {
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
		alignmentsOverlap = (overNonOver[0] !=0);
		//set insert length
		insertSize.count(overNonOver[2]);
		
		/*System.out.println("\nNumAdders "+numAdders);
		System.out.println("PreFirstLayout");
		leftLayout.print();
		System.out.println("PreSecondLayout");
		rightLayout.print();
		Misc.printArray(overNonOver);*/
		
		//skip merging non overlapping alignments
		if (onlyMergeOverlappingAlignments && alignmentsOverlap == false) return null;

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
					if (length >= stop) break;
					numIs++;
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


	public void addHeaderAndCompress() throws IOException{
		Gzipper gz = new Gzipper(saveFile);

		//add header 
		//add this program info
		samHeader.add("@PG\tID:MergePairedSamAlignments\tCL: "+programArguments);
		Iterator<String> it = samHeader.iterator();
		while (it.hasNext()) gz.println(it.next());

		//add file contents
		gz.print(outputFile);

		//close it
		gz.close();

		//delete old files
		saveFile.delete();
		outputFile.delete();
	}

	public void addHeaderAndSort() throws IOException{
		//check header size
		if (samHeader.size() == 0){
			System.err.println("\nError! no sam header was found in your file(s), cannot make a sorted bam.  Saving as an unsorted sam.");
			//replace .bam with .sam
			String name = saveFile.getName();
			name = name.replace(".bam", ".sam");
			saveFile = new File(saveFile.getParentFile(), name);
			addHeaderAndCompress();
			return;
		}

		File headerFile = new File (saveFile+"_temp.sam");
		PrintWriter out = new PrintWriter(new FileWriter(headerFile));

		//add this program info
		samHeader.add("@PG\tID:MergePairedSamAlignments\tCL: "+programArguments);
		Iterator<String> it = samHeader.iterator();
		while (it.hasNext()) out.println(it.next());

		//add file contents
		BufferedReader in = new BufferedReader (new FileReader(outputFile));
		String line;
		while ((line = in.readLine()) != null) out.println(line);

		//close 
		in.close();
		out.close();

		//sort and convert to BAM
		new PicardSortSam (headerFile, saveFile);

		//delete old files
		headerFile.delete();
		outputFile.delete();
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MergePairedSamAlignments(args);
	}		


	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
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
					case 'f': forExtraction = new File(args[++i]); break;
					case 's': saveFile = new File(args[++i]); break;
					case 'm': crossCheckMateCoordinates = false; break;
					case 'o': onlyMergeOverlappingAlignments = false; break;
					case 'r': secondPairReverseStrand = true; break;
					case 'k': skipMergingPairs = true; break;
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

		//pull files
		if (forExtraction == null) Misc.printExit("\nError: please indicate which xxx.sam/bam(.zip/.gz) file(s) to process.\n");
		File[][] tot = new File[4][];
		tot[0] = IO.extractFiles(forExtraction,".sam");
		tot[1] = IO.extractFiles(forExtraction,".sam.gz");
		tot[2] = IO.extractFiles(forExtraction,".sam.zip");
		tot[3] = IO.extractFiles(forExtraction,".bam");

		dataFiles = IO.collapseFileArray(tot);
		if (dataFiles == null || dataFiles.length==0) dataFiles = IO.extractFiles(forExtraction);
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.sam/bam(.zip/.gz) file(s)!\n");

		//check save file
		if (saveFile != null){
			try {
				saveFile.createNewFile();
			} catch (IOException e) {
				e.printStackTrace();
			}
			if (saveFile.canWrite() == false) Misc.printErrAndExit("\nError: cannot create or modify your indicated save file -> "+saveFile);
			if (saveFile.getName().endsWith(".sam") == false && saveFile.getName().endsWith(".bam") == false)  Misc.printErrAndExit("\nError: your indicated save file must end with xxx.sam or xxx.bam -> "+saveFile);
		}
		else {
			String saveFileString;
			if (dataFiles.length == 1) saveFileString = Misc.removeExtension(dataFiles[0].toString());
			else saveFileString = dataFiles[0].getParent();

			String mq = ((int)minimumMappingQualityScore)+"MQ";
			saveFile = new File (saveFileString+"_MPSA"+mq+(int)maximumAlignmentScore+"AS.bam");
		}

		//print info
		System.out.println(maximumAlignmentScore+ "\tMaximum alignment score.");
		System.out.println(minimumMappingQualityScore+ "\tMinimum mapping quality score.");
		System.out.println(removeControlAlignments +"\tRemove control chrPhiX and chrAdapter alignments.");
		System.out.println(secondPairReverseStrand +"\tSecond read pair's strand has been reversed.");
		System.out.println(crossCheckMateCoordinates +"\tCross check read mate coordinates.");
		System.out.println(maximumProperPairDistanceForMerging +"\tMaximum bp distance for merging paired alignments.");
		System.out.println(onlyMergeOverlappingAlignments +"\tOnly merge overlapping alignments.");
		if (skipMergingPairs) System.out.println(skipMergingPairs +"\tSkip merging paired alignments, useful for testing effect of merging on downstream analysis.");



	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          MergePairedSamAlignments: Dec 2012                      **\n" +
				"**************************************************************************************\n" +
				"Merges proper paired alignments that pass a variety of checks and thresholds. Only\n" +
				"unambiguous pairs will be merged. Increases base calling accuracy in overlap and helps\n" +
				"avoid non-independent variant observations and other double counting issues. Identical\n" +
				"overlapping bases are assigned the higher quality scores. Disagreements are resolved\n" +
				"toward the higher quality base. If too close in quality, then the quality is set to 0.\n" +
				"Be certain your input bam/sam file(s) are sorted by query name, NOT coordinate. \n" +

				"\nOptions:\n"+
				"-f The full path file or directory containing raw xxx.sam(.gz/.zip OK)/.bam file(s)\n" +
				"      paired alignments that are sorted by query name (standard novoalign output).\n" +
				"      Multiple files will be merged.\n" +

				"\nDefault Options:\n"+
				"-s Save file, defaults to that inferred by -f. If an xxx.sam extension is provided,\n" +
				"      the alignments won't be sorted by coordinate and saved as a bam file.\n"+
				"-a Maximum alignment score (AS:i: tag). Defaults to 120, smaller numbers are more\n" +
				"      stringent. Approx 30pts per mismatch for novoalignments.\n"+
				"-q Minimum mapping quality score, defaults to 13, larger numbers are more stringent.\n" +
				"      Set to 0 if processing splice junction indexed RNASeq data.\n"+
				"-r The second paired alignment's strand is reversed. Defaults to not reversed.\n" +
				"-d Maximum acceptible base pair distance for merging, defaults to 5000.\n"+
				"-m Don't cross check read mate coordinates, needed for merging repeat matches. Defaults\n" +
				"      to checking.\n"+
				"-o Merge all proper paired alignments. Defaults to only merging those that overlap.\n"+
				"-k Skip merging paired alignments. Defaults to merging. Useful for testing effect of\n" +
				"      merging on downstream analysis.\n"+

				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/MergePairedSamAlignments -f /Novo/Run7/\n" +
				"     -c -s /Novo/STPParsedBams/run7.bam -d 10000 \n\n" +

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
