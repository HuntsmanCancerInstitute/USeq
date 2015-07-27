package edu.utah.seq.parsers;

import java.io.*;

import util.bio.seq.Seq;
import util.gen.*;

import java.util.*;
import java.util.regex.Matcher;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import edu.utah.seq.data.ChromData;
import edu.utah.seq.data.sam.*;

/**Used to parse one chromosome worth of alignments in a thread*/
public class PairedAlignmentChrParser extends Thread{

	//fields
	private String chromosome = null;
	private LinkedHashMap<String, ArrayList<SamAlignment>> alignments = new LinkedHashMap<String, ArrayList<SamAlignment>>();
	private Gzipper samOut = null;
	private SAMRecordIterator samIterator;
	private ArrayList<SamAlignment> firstPairs = new ArrayList<SamAlignment>();
	private ArrayList<SamAlignment> secondPairs = new ArrayList<SamAlignment>();
	private MergePairedAlignments mpa = null;
	private boolean crossCheckMateCoordinates;
	private int maximumProperPairDistanceForMerging;
	private boolean secondPairReverseStrand;
	private String namePhiXChromosome = "chrPhiX";
	private String nameAdapterChromosome = "chrAdap";
	private int minimumDiffQualScore = 3;
	private double minimumFractionInFrameMismatch = 0.05;
	private float maximumAlignmentScore;
	private float minimumMappingQualityScore;
	private boolean complete = false;
	private boolean started = false;
	private boolean saveSams = false;
	private boolean removeDuplicates = false;
	private boolean onlyMergeOverlappingAlignments = true;

	//local counters
	private int numberAlignments = 0;
	private int numberAlignmentsMissingPair = 0;
	private int numberPrintedAlignments = 0;
	private int numberPairsFailingChrDistStrand = 0;
	private int numberRepeatAlignmentsLackingMate = 0;
	private int numberPairsFailingMateCrossCoordinateCheck = 0;
	private int numberMergedPairs = 0;
	private int numberFailedMergedPairs = 0;
	private double numberOverlappingBases = 0;
	private double numberNonOverlappingBases = 0;
	private int numberUnmapped = 0;
	private int numberFailingVendorQC = 0;
	private int numberPhiX = 0;
	private int numberAdapter = 0;
	private int numberFailingAlignmentScore = 0;
	private int numberFailingMappingQualityScore = 0;
	private int numberDuplicates = 0;
	private int numberPassingAlignments = 0;
	private int numberNonPairedAlignments = 0;
	private int numberNonProperPairedAlignments = 0;
	private int numberUnmappedMatePairedAlignments = 0;
	private Histogram insertSize = new Histogram(0,2001,400);
	private ChromData chromData;
	private DataOutputStream chromDataOut;
	private boolean queryUnMapped = false;
	private boolean calculateBaseQualities = false;
	private double numberPassingBases = 0;
	private double numberPassingQ20Bases = 0;
	private double numberPassingQ30Bases = 0;

	public PairedAlignmentChrParser (String chromosome, MergePairedAlignments mpa){
		//set params
		this.chromosome = chromosome;
		this.mpa = mpa;
		queryUnMapped = chromosome.equals("unmapped");
		crossCheckMateCoordinates = mpa.isCrossCheckMateCoordinates();
		maximumProperPairDistanceForMerging = mpa.getMaximumProperPairDistanceForMerging();
		secondPairReverseStrand = mpa.isSecondPairReverseStrand();
		maximumAlignmentScore = mpa.getMaximumAlignmentScore();
		minimumMappingQualityScore = mpa.getMinimumMappingQualityScore();
		saveSams = mpa.isSaveSams();
		removeDuplicates = mpa.isRemoveDuplicates();
		onlyMergeOverlappingAlignments = mpa.isOnlyMergeOverlappingAlignments();
		calculateBaseQualities = mpa.isCalculateBaseQualities();
	}

	public void run(){
		try {
			started = true;
			
			File chromDataFile = null;
			if (saveSams == false){ 
				//make binary data container, nonstranded
				chromDataFile = new File(mpa.getSaveDirectory(), chromosome+".ChromData");
				chromDataOut = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(chromDataFile)));
				chromData = new ChromData (1000000, 0, chromosome, ".", chromDataFile, chromDataOut);
			}
			
			//parse
			parseAlignmentFile();
			
			//clean up
			if (chromDataFile != null){
				chromDataOut.close();
				if (numberPrintedAlignments == 0) chromDataFile.deleteOnExit();
			}
			alignments = null;
			firstPairs = null;
			secondPairs = null;
			
			//finish
			complete = true;
		} catch (Exception e) {
			Misc.printErrAndExit("Error parsing "+chromosome+" data from "+mpa.getBamFile());
			e.printStackTrace();
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
		Matcher mat = MergePairedAlignments.CIGAR_SUB.matcher(cigar);
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
		Matcher mat = MergePairedAlignments.CIGAR_SUB.matcher(cigar);
		while (mat.find()){
			length += Integer.parseInt(mat.group(1));
		}
		return length;
	}

	/**Attempts to merge alignments. Doesn't check if proper pairs!  Returns null if it cannot. This modifies the input SamAlignments so print first before calling*/
	public SamAlignment mergePairedAlignments(SamAlignment first, SamAlignment second) {

		//trim them of soft clipped info
		first.trimMaskingOfReadToFitAlignment();
		second.trimMaskingOfReadToFitAlignment();

		//look for bad CIGARs
		if (MergePairedAlignments.CIGAR_BAD.matcher(first.getCigar()).matches()) Misc.printErrAndExit("\nError: unsupported cigar string! See -> "+first.toString()+"\n");
		if (MergePairedAlignments.CIGAR_BAD.matcher(second.getCigar()).matches()) Misc.printErrAndExit("\nError: unsupported cigar string! See -> "+second.toString()+"\n");

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
		numberOverlappingBases+= (double)overNonOver[0];
		numberNonOverlappingBases+= (double)overNonOver[1];
		//set insert length
		insertSize.count(overNonOver[2]);
		
		//skip merging non overlapping alignments?
		if (onlyMergeOverlappingAlignments && overNonOver[0] == 0) return null;

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
			if (removeDuplicates) isGood = false;
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
			save(sa);
			return false;
		}

		//is it not part of a proper pair?
		if (sa.isAProperPairedAlignment() == false){
			numberNonProperPairedAlignments++;
			save(sa);
			
			return false;
		}
		//is mate unmapped
		if (sa.isMateUnMapped()){
			numberUnmappedMatePairedAlignments++;
			save(sa);
			return false;
		}
		//OK looks good to add for potential pairing
		return true;
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

	/**Runs through the hash, then nulls it.*/
	private void processAlignmentHash() {
		for (String readName: alignments.keySet()){
			ArrayList<SamAlignment> alignmentsToMerge = alignments.get(readName);
			filterPrintAlignments(alignmentsToMerge);
		}
		alignments = null;
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
					save(sam);
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
							save(firsts[i]);
							save(seconds[j]);
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
					save(s);
				}
			}
			for (SamAlignment s: seconds){
				if (s !=null){
					numberRepeatAlignmentsLackingMate++;
					save(s);
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

	private void processPair(SamAlignment first, SamAlignment second) throws Exception {
		//check if the pair is from the same chromosome, acceptable distance, and proper strand
		if (testChrDistStrnd(first, second)){
			//cross validate mate information? this is often messed up especially if from the STP app.
			if (crossCheckMateCoordinates){
				if (testMateChrPosition(first, second)){
					//attempt a merge
					mergeAndScorePair(first, second);
				}
				else {
					save(first);
					save(second);
					numberPairsFailingMateCrossCoordinateCheck++;
				}
			}
			else {
				//attempt a merge
				mergeAndScorePair(first, second);
			}
		}
		else {
			save(first);
			save(second);
			numberPairsFailingChrDistStrand++;
		}

	}

	/**Attempts to merge a proper paired alignment.  Increments counters and sends examples to the gzipper.*/
	private void mergeAndScorePair(SamAlignment first, SamAlignment second) throws Exception{
		//collect string rep since the merge method will modify the SamAlignment while merging 
		//so if it fails you can output the unmodified SamAlignment
		String firstSamString = first.toString();
		String secondSamString = second.toString();

		SamAlignment mergedSam = mergePairedAlignments(first, second);
		//failed to merge?
		if (mergedSam!=null) {
			numberMergedPairs++;
			save(mergedSam);
		}
		else {
			save(new SamAlignment(firstSamString, false));
			save(new SamAlignment(secondSamString, false));
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

	public void parseAlignmentFile(){
		String line = null;
		try {
			//get reader and interator
			SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			SamReader samReader = factory.open(mpa.getBamFile());
			if (queryUnMapped) samIterator = samReader.queryUnmapped();
			else samIterator = samReader.queryOverlapping(chromosome, 0, 0);
			
			//any records?
			if (samIterator.hasNext() == false) return;

			//create gzipper?
			if (saveSams){
				String name = Misc.removeExtension(mpa.getBamFile().getName());
				samOut = new Gzipper (new File (mpa.getSaveDirectory(), chromosome+"_"+name+".sam.gz"));
				//add header to results file
				samOut.println(mpa.getSamHeader());
			}

			ArrayList<SamAlignment> alignmentsToMerge = null;
			int dotCounter = 0;

			//load up hash with lots of records
			if (loadAlignments(mpa.getNumberAlignmentsToLoad()) == false) processAlignmentHash();
			
			//more to work with so....
			else {
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
			}
			samIterator.close();
			samReader.close();
			if (saveSams) samOut.close();
		} catch (Exception e) {
			System.err.println("\nError processing -> "+line+"\nFrom "+mpa.getBamFile());
			e.printStackTrace();	
		} 
	}
	/**Fetches ArrayList of alignments, then removes it from the hash.
	 * Returns null if no more.*/
	private ArrayList<SamAlignment> fetchRemoveAlignmentBlock(){
		try {
			Iterator<String> it = alignments.keySet().iterator();
			if (it.hasNext() == false) return null;
			String readName = it.next();
			ArrayList<SamAlignment> alignmentsToMerge = alignments.get(readName);
			alignments.remove(readName);			
			return alignmentsToMerge;
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
	}
	private void save(SamAlignment sam) throws IOException{
		//save string?
		if (saveSams) samOut.println(sam);
		else {
			//get start and end
			int start = sam.getUnclippedStart(); //is minus one needed?
			int end = SamAlignment.countLengthOfCIGAR(sam.getCigar()) + start +1;
			//set first and last?
			if (start < chromData.firstBase) chromData.firstBase = start;
			if (end > chromData.lastBase) chromData.lastBase = end;
			//save data start and cigar
			chromDataOut.writeInt(start);
			chromDataOut.writeUTF(sam.getCigar());
		}
		numberPrintedAlignments++;
		
		//calculate base level quality stats? this is slow!
		if (calculateBaseQualities) incrementBaseQualities(sam);
	}
	
	private void incrementBaseQualities(SamAlignment sam){
		SamLayoutLite sll = new SamLayoutLite(sam);
		String qual = sll.getMIQualities();
		int[] values = Seq.convertSangerQualityScores(qual);
		numberPassingBases+= values.length;
		for (int i=0; i< values.length; i++){
			if (values[i] > 19) {
				numberPassingQ20Bases++;
				if (values[i] > 29) numberPassingQ30Bases++;
			}
		}
	}

	public boolean isComplete() {
		return complete;
	}

	public int getNumberAlignments() {
		return numberAlignments;
	}

	public int getNumberAlignmentsMissingPair() {
		return numberAlignmentsMissingPair;
	}

	public int getNumberPrintedAlignments() {
		return numberPrintedAlignments;
	}

	public int getNumberPairsFailingChrDistStrand() {
		return numberPairsFailingChrDistStrand;
	}

	public int getNumberRepeatAlignmentsLackingMate() {
		return numberRepeatAlignmentsLackingMate;
	}

	public int getNumberPairsFailingMateCrossCoordinateCheck() {
		return numberPairsFailingMateCrossCoordinateCheck;
	}

	public int getNumberMergedPairs() {
		return numberMergedPairs;
	}

	public int getNumberFailedMergedPairs() {
		return numberFailedMergedPairs;
	}

	public double getNumberOverlappingBases() {
		return numberOverlappingBases;
	}

	public double getNumberNonOverlappingBases() {
		return numberNonOverlappingBases;
	}

	public int getNumberUnmapped() {
		return numberUnmapped;
	}

	public int getNumberFailingVendorQC() {
		return numberFailingVendorQC;
	}

	public int getNumberPhiX() {
		return numberPhiX;
	}

	public int getNumberAdapter() {
		return numberAdapter;
	}

	public int getNumberFailingAlignmentScore() {
		return numberFailingAlignmentScore;
	}

	public int getNumberFailingMappingQualityScore() {
		return numberFailingMappingQualityScore;
	}

	public int getNumberDuplicates() {
		return numberDuplicates;
	}

	public int getNumberPassingAlignments() {
		return numberPassingAlignments;
	}

	public int getNumberNonPairedAlignments() {
		return numberNonPairedAlignments;
	}

	public int getNumberNonProperPairedAlignments() {
		return numberNonProperPairedAlignments;
	}

	public int getNumberUnmappedMatePairedAlignments() {
		return numberUnmappedMatePairedAlignments;
	}

	public Histogram getInsertSize() {
		return insertSize;
	}

	public String getChromosome() {
		return chromosome;
	}

	public boolean isStarted() {
		return started;
	}

	public ChromData getChromData() {
		return chromData;
	}

	public double getNumberPassingBases() {
		return numberPassingBases;
	}

	public double getNumberPassingQ20Bases() {
		return numberPassingQ20Bases;
	}

	public double getNumberPassingQ30Bases() {
		return numberPassingQ30Bases;
	}


}
