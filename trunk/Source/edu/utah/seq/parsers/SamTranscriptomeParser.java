
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;
import edu.utah.seq.data.sam.*;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class SamTranscriptomeParser{
	//fields
	private File[] dataFiles;
	private File saveFile;
	private File outputFile;
	private String replacementHeader = null;
	private float maximumAlignmentScore = 90;
	private float minimumMappingQualityScore = 0;
	private int numberAlignments = 0;
	private int numberUnmapped = 0;
	private int numberFailingVendorQC = 0;
	private int numberPassingAlignments = 0;
	private int numberFailingAlignmentScore = 0;
	private int numberFailingMappingQualityScore = 0;
	private int numberAdapter = 0;
	private int numberPhiX = 0;
	private int numberPrintedAlignments = 0;
	private double numberOverlappingBases = 0;
	private double numberNonOverlappingBases = 0;
	private int maxMatches = 1;
	private Gzipper samOut;
	private Gzipper failedSamOut = null;
	private boolean saveUnmappedAndFailedScore = false;
	private HashMap <String, Integer> chromLength = new HashMap <String, Integer>();
	private String programArguments;
	private String genomeVersion = null;
	private static final Pattern TAB = Pattern.compile("\\t");
	private Pattern CIGAR_BAD = Pattern.compile(".*[^\\dMDIN].*");
	private boolean reverseStrand = false;
	private boolean reverseBoth = false;
	private boolean removeControlAlignments = true;
	private boolean randomPickAlignment = false;
	private HashSet<SamAlignment> uniques = new HashSet<SamAlignment>();
	private ArrayList<SamAlignment> firstPair = new ArrayList<SamAlignment>();
	private ArrayList<SamAlignment> secondPair = new ArrayList<SamAlignment>();
	private Random random = new Random();
	private boolean verbose = true;
	//for trimming/ merging paired data
	private boolean mergePairedAlignments = false;
	private int minimumDiffQualScore = 3;
	private double minimumFractionInFrameMismatch = 0.01;
	private int maximumProperPairDistanceForMerging = 300000;
	private double numberMergedPairs = 0;
	private double numberFailedMergedPairs = 0;
	private int maximumMappingQuality = 0;

	//constructors
	public SamTranscriptomeParser(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			if (verbose) System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	//for integration with RNASeq app
	public SamTranscriptomeParser(File[] samFiles, File saveFile, float maximumAlignmentScore, String genomeVersion, int maxMatches, boolean verbose, boolean flipped, boolean bothFlipped) throws IOException{
		this.dataFiles = samFiles;
		this.saveFile = saveFile;
		this.maximumAlignmentScore = maximumAlignmentScore;
		this.genomeVersion = genomeVersion;
		this.maxMatches = maxMatches;
		this.randomPickAlignment = false;
		this.verbose = verbose;
		this.reverseStrand = flipped;
		this.reverseBoth = bothFlipped;
		doWork();
	}

	public void doWork() throws IOException{
		//make print writer
		outputFile = new File(saveFile+"_temp.sam.gz");
		//samOut = new PrintWriter( new FileWriter (outputFile));
		samOut = new Gzipper(outputFile);
		if (replacementHeader != null) samOut.println(replacementHeader);

		if (saveUnmappedAndFailedScore) {
			String name = Misc.removeExtension(saveFile.getName());
			File failedReadOutputFile = new File(saveFile.getParentFile(), name+"_UnMappedPoorScore.sam.gz");
			failedSamOut = new Gzipper(failedReadOutputFile);
			if (replacementHeader != null) failedSamOut.println(replacementHeader);
		}

		//for each file, parse and save to disk	
		if (verbose) System.out.println("\nParsing, filtering, and merging SAM files...");
		for (int i=0; i< dataFiles.length; i++){
			if (verbose) System.out.print("\t"+dataFiles[i].getName());
			parseFile(dataFiles[i]); 
			if (verbose) System.out.println();
		}

		//close the writers
		samOut.close();
		if (saveUnmappedAndFailedScore) failedSamOut.close();

		//add header and output results, this deletes the outputFile too
		if (saveFile.getName().endsWith(".sam.gz")){
			//replacement header present? just change name
			if (replacementHeader != null){
				outputFile.renameTo(saveFile);
			}
			//nope, have to append header
			else {
				if (verbose) System.out.println("\nAdding SAM header and gzip compressing xxx.sam file...");
				addHeaderAndCompress();
			}
		}
		else {
			if (verbose) System.out.println("\nAdding SAM header, sorting, and writing bam output with Picard's SortSam...");
			addHeaderAndSort();
		}

		//stats
		double fractionPassing = ((double)numberPassingAlignments)/((double)numberAlignments);
		if (verbose) System.out.println("\nStats (some flags aren't set so be suspicious of zero read catagories):\n");
		System.out.println("\t"+numberAlignments+"\tTotal # Alignments from raw sam file");
		System.out.println("\t"+numberPassingAlignments+"\tAlignments passing filters ("+Num.formatPercentOneFraction(fractionPassing)+")");
		System.out.println("\t\t"+numberUnmapped+"\t# Unmapped Reads");
		System.out.println("\t\t"+numberFailingVendorQC+"\t# Alignments failing vendor/ platform QC");
		System.out.println("\t\t"+numberFailingAlignmentScore+"\t# Alignments failing alignment score");
		System.out.println("\t\t"+numberFailingMappingQualityScore+"\t# Alignments failing mapping quality score");
		System.out.println("\t\t"+numberAdapter+"\t# Adapter alignments");
		System.out.println("\t\t"+numberPhiX+"\t# PhiX alignments");
		System.out.println();
		System.out.println("\t"+numberPrintedAlignments+"\t# Alignments written to SAM/BAM file. These passed the maxMatch, collapsed coordinate, and possibly merge pairs filters.");
		//pair overlap stats?
		if (mergePairedAlignments) {
			double fractionFailed = numberFailedMergedPairs/ (numberFailedMergedPairs+ numberMergedPairs);
			System.out.println("\t"+Num.formatNumber(fractionFailed, 4)+"\tFraction proper paired alignments that could not be merged.");
			double totalBases = numberNonOverlappingBases + numberOverlappingBases;
			double fractionOverlap = numberOverlappingBases/totalBases;
			String fractionString = Num.formatNumber(fractionOverlap, 4);
			System.out.println("\t"+fractionString+"\tFraction overlapping bases in proper paired alignments.");
		}
	}

	public boolean parseFile(File samFile){
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
					if (verbose) System.out.print(".");
					dotCounter = 0;
				}
				line = line.trim();

				//skip blank lines
				if (line.length() == 0) continue;

				//header line?
				if (line.startsWith("@")){
					//parse genome version?
					if (genomeVersion == null) {
						String[] tokens = TAB.split(line);
						for (String t : tokens){
							if (t.startsWith("AS:")) {
								genomeVersion = t;
								break;
							}
						}
					}
					continue;
				}

				SamAlignment sa;
				try {
					sa = new SamAlignment(line, true);
				} catch (MalformedSamAlignmentException e) {
					if (verbose) System.out.println("\nSkipping malformed sam alignment -> "+e.getMessage());
					if (numBadLines++ > 1000) Misc.printErrAndExit("\nAboring: too many malformed SAM alignments.\n");
					continue;
				}
				numberAlignments++;

				//is it aligned?
				if (sa.isUnmapped()) {
					numberUnmapped++;
					if (saveUnmappedAndFailedScore) failedSamOut.println(line);
					continue;
				}

				//does it pass the vendor qc?
				if (sa.failedQC()) {
					numberFailingVendorQC++;
					continue;
				}

				//skip phiX and adapter
				boolean firstIsPhiX = sa.getReferenceSequence().startsWith("chrPhiX");
				if (firstIsPhiX){
					numberPhiX++;
					if (removeControlAlignments) continue;
				}
				boolean firstIsAdapt = sa.getReferenceSequence().startsWith("chrAdapt");
				if (firstIsAdapt){
					numberAdapter++;
					if (removeControlAlignments) continue;
				}

				//does it pass the scores threshold?
				int alignmentScore = sa.getAlignmentScore();
				if (alignmentScore != Integer.MIN_VALUE){
					if (alignmentScore > maximumAlignmentScore){
						numberFailingAlignmentScore++;
						if (saveUnmappedAndFailedScore) failedSamOut.println(line);
						continue;
					}
				}

				//check mapping quality for genomic match reads?
				if (minimumMappingQualityScore !=0 && sa.isSpliceJunction() == false){
					if (sa.getMappingQuality() < minimumMappingQualityScore){
						numberFailingMappingQualityScore++;
						if (saveUnmappedAndFailedScore) failedSamOut.println(line);
						continue;
					}
				}
				
				//reset maximumMappingQuality?
				if (maximumMappingQuality !=0){
					if (sa.getMappingQuality() > maximumMappingQuality) sa.setMappingQuality(maximumMappingQuality);
				}

				//modify second read if phiX or adapter; some paired reads have a mate hitting the control chroms
				SamAlignmentFlags saf = null;
				if (removeControlAlignments && sa.isPartOfAPairedAlignment() && sa.isMateUnMapped() == false){
					if ( sa.getMateReferenceSequence().startsWith("chrPhiX") || sa.getMateReferenceSequence().startsWith("chrAdapt")) {
						saf = new SamAlignmentFlags(sa.getFlags());
						saf.setMateUnMapped(true);
						sa.setUnMappedMate();
						sa.setFlags(saf.getFlags());
					}
				}

				//OK, it passes, increment counter
				numberPassingAlignments++;

				//set inferred insert size and mate position to zero
				sa.setInferredInsertSize(0);

				//reverse strands of both alignments
				if (reverseBoth) {
					if (saf == null) {
						saf = new SamAlignmentFlags(sa.getFlags());
					}
					if (saf.isReverseStrand()) {
						saf.setReverseStrand(false);
					} else {
						saf.setReverseStrand(true);
					}
					sa.setFlags(saf.getFlags());
				}

				//reverse second alignment?
				if (reverseStrand && sa.isSecondPair()) {
					//sa.printFlags();
					if (saf == null) saf = new SamAlignmentFlags(sa.getFlags());
					if (saf.isReverseStrand()) saf.setReverseStrand(false);
					else saf.setReverseStrand(true);
					sa.setFlags(saf.getFlags());	
				}

				//convert possible splice junctions to genomic coordinates, toss MD and RG tags
				if (sa.convertTranscriptomeAlignment(true) == false) {
					System.err.println("Failed to convert, skippping ->\n"+sa);
					numberPassingAlignments--;
					continue;
				}

				String readName = sa.getName();

				//prior set
				if (priorSet == false){
					priorSet = true;
					priorReadName = readName;
					alignmentsToSave.add(sa);
				}
				//is it an old read?
				else if (readName.equals(priorReadName)){
					alignmentsToSave.add(sa);
				}
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

	/**Takes a block of alignments all originating from the same fragment.  Collapses those with the same coordinates and CIGAR.
	 * Saves them if they pass a whole set of filters.*/
	public void filterPrintAlignments(ArrayList<SamAlignment> al){
		try {
			//collapse alignments with same coordinates, CIGAR, and from same read pair 
			uniques.clear();
			for (SamAlignment sam : al) uniques.add(sam); 

			//set booleans for different cases
			firstPair.clear();
			secondPair.clear();
			boolean firstPairPresent = false;
			boolean secondPairPresent = false;
			boolean nonPairedPresent = false;
			for (SamAlignment sam : uniques) {
				if (sam.isFirstPair()) {
					firstPairPresent = true;
					firstPair.add(sam);
				}
				else if (sam.isSecondPair()) {
					secondPairPresent = true;
					secondPair.add(sam);
				}
				else {
					nonPairedPresent = true;
					firstPair.add(sam);
				}
			}

			int numberFirstPair = firstPair.size();
			int numberSecondPair = secondPair.size();

			//fix mate info in pairs? Can only do this if one first and one second.  Don't know how to join up repeat matches?
			//merge ?
			if (secondPairPresent && numberFirstPair == 1 && numberSecondPair == 1){
				SamAlignment first = firstPair.get(0);
				SamAlignment second = secondPair.get(0);

				//any junctions?
				if (first.isConvertedJunctionCoordinates()){
					second.setMateReferenceSequence(first.getReferenceSequence());
					second.setMatePosition(first.getPosition());
				}
				if (second.isConvertedJunctionCoordinates()){
					first.setMateReferenceSequence(second.getReferenceSequence());
					first.setMatePosition(second.getPosition());
				}

				//merge pairs?
				if (mergePairedAlignments) {
					//same chromosome?
					if (first.getReferenceSequence().equals(second.getReferenceSequence())) {
						//within acceptable distance
						int diff = Math.abs(first.getPosition()- second.getPosition());
						if (diff < maximumProperPairDistanceForMerging){
							//correct strand?
							boolean merge = false;
							if (reverseStrand){
								if (first.isReverseStrand() == second.isReverseStrand()) merge = true;
							}
							else if (first.isReverseStrand() != second.isReverseStrand()) merge = true;
							//attempt a merge?
							if (merge){								
								SamAlignment mergedSam = mergePairedAlignments (first, second);

								//success?
								if (mergedSam != null) {
									printSam(mergedSam,1);
									numberMergedPairs++;
									return;
								}
								else numberFailedMergedPairs++;
							}
						}
					}
				}

			}

			//print em?
			if (numberFirstPair <= maxMatches && numberSecondPair <=maxMatches){
				//print reads
				for (SamAlignment sam : firstPair) printSam(sam, numberFirstPair);
				if (numberSecondPair!=0) for (SamAlignment sam : secondPair) printSam(sam, numberSecondPair);
			}
			//don't print all, maybe just random pick?
			else {
				//pick and print random alignments?
				if (randomPickAlignment){
					//first pair or non paired
					if (firstPairPresent || nonPairedPresent){
						int index = random.nextInt(numberFirstPair);
						printSam(firstPair.get(index), numberFirstPair);
					}
					//second pair?
					if (secondPairPresent){
						int index = random.nextInt(numberSecondPair);
						printSam(secondPair.get(index), numberSecondPair);
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem printing alignment block!?\n");
		}
	}
	/*Attempts to merge alignments. Doesn't check if proper pairs!  Returns null if it cannot. This modifies the input SamAlignments so print first before calling
	private SamAlignment mergePairedAlignments(SamAlignment first, SamAlignment second) {
		//trim them of soft clipped info
		first.trimMaskingOfReadToFitAlignment();
		second.trimMaskingOfReadToFitAlignment();

		//look for bad CIGARs
		if (CIGAR_BAD.matcher(first.getCigar()).matches()) Misc.printErrAndExit("\nError: unsupported cigar string! See -> "+first.toString()+"\n");
		if (CIGAR_BAD.matcher(second.getCigar()).matches()) Misc.printErrAndExit("\nError: unsupported cigar string! See -> "+second.toString()+"\n");

		//fetch coordinates
		int startBaseFirst = first.getPosition();
		int stopBaseFirst = startBaseFirst + countLengthOfCigar(first.getCigar());
		int startBaseSecond = second.getPosition();
		int stopBaseSecond = startBaseSecond + countLengthOfCigar(second.getCigar());

		//make arrays to hold sequence and qualities
		int start = startBaseFirst;
		if (startBaseSecond < start) start = startBaseSecond;
		int stop = stopBaseFirst;
		if (stopBaseSecond > stop) stop = stopBaseSecond;
		int size = stop-start;

		SamLayout firstLayout = new SamLayout(size);
		SamLayout secondLayout = new SamLayout(size);

		//layout data
		firstLayout.layoutCigar(start, first);
		secondLayout.layoutCigar(start, second);



		//merge layouts, modifies original layouts so print first if you want to see em before mods.
		SamLayout mergedSamLayout = SamLayout.mergeLayouts(firstLayout, secondLayout, minimumDiffQualScore, minimumFractionInFrameMismatch);

		if (mergedSamLayout == null) {
			//if (true){


				//add failed merge tag
				first.addMergeTag(false);
				second.addMergeTag(false);

			return null;
		}

		else {
			//calculate overlap
			int[] overNonOver = SamLayout.countOverlappingBases(firstLayout, secondLayout);
			numberOverlappingBases+= overNonOver[0];
			numberNonOverlappingBases+= overNonOver[1];

			//make merged
			SamAlignment mergedSam = makeSamAlignment(first, second, mergedSamLayout, start);
			return mergedSam;
		}

	}*/

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
		int stopLeft = startLeft + MergePairedSamAlignments.countLengthOfCigar(left.getCigar());
		int startRight = right.getPosition();
		int stopRight = startRight + MergePairedSamAlignments.countLengthOfCigar(right.getCigar());
		int stop = stopRight;
		if (stopLeft > stop) stop = stopLeft;

		//any Is in left that preceed the start of right?
		int numAdders = MergePairedSamAlignments.countIs(left.getCigar(), startRight-startLeft);

		//make arrays to hold sequence and qualities in cigar space
		int size = numAdders + stop-startLeft;

		SamLayout leftLayout = new SamLayout(size);
		SamLayout rightLayout = new SamLayout(size);

		//layout data
		leftLayout.layoutCigar(startLeft, left);
		rightLayout.layoutCigar(startLeft-numAdders, right);

		//System.out.println("\nNumAdders "+numAdders);
		//System.out.println("PreFirstLayout");
		//leftLayout.print();
		//System.out.println("PreSecondLayout");
		//rightLayout.print();


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
			//calculate overlap
			int[] overNonOver = SamLayout.countOverlappingBases(leftLayout, rightLayout);
			numberOverlappingBases+= overNonOver[0];
			numberNonOverlappingBases+= overNonOver[1];

			//make merged
			SamAlignment mergedSam = MergePairedSamAlignments.makeSamAlignment(first.isReverseStrand(), left, right, mergedSamLayout, startLeft);

			return mergedSam;
		}

	}

	public void printSam(SamAlignment sam, int numberRepeats) throws IOException{
		numberPrintedAlignments++;
		//add/ replace IH tag for number of "Number of stored alignments in SAM that contains the query in the current record"
		sam.addRepeatTag(numberRepeats);
		samOut.println(sam);

		//get chromosome
		String chrom = sam.getReferenceSequence();
		//get current max length
		int maxLength =0;
		if (chromLength.containsKey(chrom)) maxLength = (chromLength.get(chrom)).intValue();

		//calc end position
		int endPosition = sam.countLengthOfAlignment() + sam.getPosition() + 10000;

		//reset length?
		if (endPosition > maxLength) {
			chromLength.put(chrom, new Integer(endPosition));
			//System.out.println(chrom+" "+endPosition+" Setting max from curr.");
			maxLength = endPosition;
		}

		//how about the mate?  even if it is missing this is checked  and can cause an error
		//increment chrom length for second pair
		if (sam.isPartOfAPairedAlignment() && sam.isMateUnMapped()==false){

			int mateMaxLength = 0;

			//get mat chromosome name
			String mateChrom = sam.getMateReferenceSequence();

			//is mate the same chrom? if so then set as above
			if (mateChrom.equals("=")) {
				mateChrom = chrom;
				mateMaxLength = maxLength;
			}
			//nope different chromosome so set
			else if (chromLength.containsKey(mateChrom)) mateMaxLength = chromLength.get(mateChrom);

			//don't know the end so guess. SAM stinks!
			int mateEndPosition = sam.getMatePosition() + 100000;

			if (mateEndPosition > mateMaxLength) {
				chromLength.put(mateChrom, new Integer(mateEndPosition));
				//System.out.println(mateChrom+" "+mateEndPosition+" Setting max from mate.");
			}
		}
	}

	public ArrayList<String> fetchSamHeader() {
		ArrayList<String> al = new ArrayList<String>();
		//add unsorted
		al.add("@HD\tVN:1.0\tSO:unsorted");
		//add program
		al.add("@PG\tID:SamTranscriptomeParser\tCL: args "+programArguments);
		//add readgroup
		al.add("@RG\tID:unknownReadGroup\tSM:unknownSample");
		//as sq lines for each chromosome @SQ	SN:chr10	AS:mm9	LN:129993255
		String gv = "";
		if (genomeVersion != null) gv = "\tAS:" +genomeVersion;
		//remove = chromosomes
		chromLength.remove("=");
		for (String chromosome: chromLength.keySet()){
			int length = chromLength.get(chromosome);
			al.add("@SQ\tSN:"+chromosome+ gv+ "\tLN:"+length);
		}
		return al;
	}

	public void addHeaderAndCompress() throws IOException{
		Gzipper gz = new Gzipper(saveFile);

		//add header lines
		ArrayList<String> header = fetchSamHeader();
		gz.println(header);

		//add file contents
		gz.print(outputFile);

		//close it
		gz.close();

		//delete old files
		outputFile.delete();
	}

	public void addHeaderAndSort() throws IOException{

		//if no header need to add header and copy over sam data
		File toSortFile;
		if (replacementHeader == null){

			toSortFile = new File (saveFile+"_temp.sam");
			PrintWriter out = new PrintWriter(new FileWriter(toSortFile));

			//add header lines
			ArrayList<String> header = fetchSamHeader();
			for (String s : header) out.println(s);

			//add file contents
			BufferedReader in = IO.fetchBufferedReader(outputFile);
			String line;
			while ((line = in.readLine()) != null) out.println(line);

			//close 
			in.close();
			out.close();

		}
		else toSortFile = outputFile;

		//sort and convert to BAM
		new PicardSortSam (toSortFile, saveFile);

		//delete old files
		toSortFile.delete();
		outputFile.delete();
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SamTranscriptomeParser(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		File replacementHeaderFile = null;
		String useqVersion = IO.fetchUSeqVersion();
		programArguments = useqVersion+" "+Misc.stringArrayToString(args, " ");
		if (verbose) System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': forExtraction = new File(args[++i]); break;
					case 'n': maxMatches = Integer.parseInt(args[++i]); break;
					case 'r': reverseStrand = true; break;
					case 'u': saveUnmappedAndFailedScore = true; break;
					case 's': saveFile = new File(args[++i]); break;
					case 'h': replacementHeaderFile = new File(args[++i]); break;
					case 'c': removeControlAlignments = false; break;
					case 'b': reverseBoth = true; break;
					case 'd': randomPickAlignment = true; break;
					case 'p': mergePairedAlignments = true; break;
					case 'q': maximumProperPairDistanceForMerging = Integer.parseInt(args[++i]); mergePairedAlignments = true; break;
					case 'a': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'x': maximumMappingQuality = Integer.parseInt(args[++i]); break;
					case 'm': minimumMappingQualityScore = Float.parseFloat(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//incompatible scores?
		if (minimumMappingQualityScore > 0){
			if (maxMatches > 1 || randomPickAlignment == true) Misc.printExit("Error: One cannot have a positive minimum mapping quality score and ALSO" +
			" enable more than one maximum matches or random picks! Either set the minimum mapping quality to 0, the default, OR set the maximum matches to 1 and not select the pick random alignment option.\n");
		}

		//pull files
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".sam");
		tot[1] = IO.extractFiles(forExtraction,".sam.gz");
		tot[2] = IO.extractFiles(forExtraction,".sam.zip");

		dataFiles = IO.collapseFileArray(tot);
		if (dataFiles == null || dataFiles.length==0) dataFiles = IO.extractFiles(forExtraction);
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.txt(.zip/.gz) file(s)!\n");

		//check save file
		if (saveFile != null){
			if (saveFile.getName().endsWith(".sam") == false && saveFile.getName().endsWith(".bam") == false)  Misc.printErrAndExit("\nError: your indicated save file must end with xxx.sam or xxx.bam -> "+saveFile);
			if (saveFile.getName().endsWith(".sam") ){
				saveFile = new File (saveFile.getParentFile(), saveFile.getName()+".gz");
			}
		}
		else {
			String saveFileString;
			if (dataFiles.length == 1) saveFileString = Misc.removeExtension(dataFiles[0].toString());
			else saveFileString = dataFiles[0].getParent();
			String randomSelected = "";
			if (randomPickAlignment) randomSelected = "Rnd";
			String mq = "";
			if (minimumMappingQualityScore !=0) mq = ((int)minimumMappingQualityScore)+"MQ";
			saveFile = new File (saveFileString+"_STP"+mq+maxMatches+"N"+randomSelected+(int)maximumAlignmentScore+"A.bam");
		}

		//load header?
		if (replacementHeaderFile != null){
			StringBuilder sb = new StringBuilder();
			String[] lines = IO.loadFile(replacementHeaderFile);
			if (lines.length ==0) Misc.printErrAndExit("\nError: replacement header contains no comment lines?\n");
			for (String l: lines) {
				sb.append(l);
				sb.append("\n");
			}
			//add program
			sb.append("@PG\tID:SamTranscriptomeParser\tCL: args "+programArguments);
			replacementHeader = sb.toString();
		}


		//print info
		if (verbose) {
			System.out.println(maximumAlignmentScore+ "\tMaximum alignment score.");
			System.out.println(minimumMappingQualityScore+ "\tMinimum mapping quality score.");
			System.out.println(maxMatches+"\tMaximum locations each read may align.");
			System.out.println(reverseStrand +"\tReverse the strand of the second paired alignemnt.");
			System.out.println(reverseBoth + "\tReverse the strand of both alignments.");
			System.out.println(saveUnmappedAndFailedScore +"\tSave unmapped and low score reads.");
			System.out.println(removeControlAlignments +"\tRemove control chrPhiX and chrAdapter alignments.");
			System.out.println(randomPickAlignment +"\tRandomly choose an alignment from read blocks that fail the max locations threshold.");
			System.out.println(mergePairedAlignments +"\tMerge proper paired alignments.");
			if (mergePairedAlignments) System.out.println(maximumProperPairDistanceForMerging +"\tMaximum bp distance for merging paired alignments.\n");

		}

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Sam Transcriptome Parser: June 2013                     **\n" +
				"**************************************************************************************\n" +
				"STP takes SAM alignment files that were aligned against chromosomes and extended\n" +
				"splice junctions (see MakeTranscriptome app), converts the coordinates to genomic\n" +
				"space and sorts and saves the alignments in BAM format. Although alignments don't need\n" +
				"to be sorted by chromosome and position, it is assumed all the alignments for a given\n" +
				"fragment are grouped together. \n" +

				"\nOptions:\n"+
				"-f The full path file or directory containing raw xxx.sam(.gz/.zip OK) file(s).\n" +
				"      Multiple files will be merged.\n" +

				"\nDefault Options:\n"+
				"-s Save file, defaults to that inferred by -f. If an xxx.sam extension is provided,\n" +
				"      the alignments won't be sorted by coordinate or saved as a bam file.\n"+
				"-a Maximum alignment score. Defaults to 90, smaller numbers are more stringent.\n" +
				"      Approx 30pts per mismatch.\n"+
				"-m Minimum mapping quality score, defaults to 0 (no filtering), larger numbers are\n" +
				"      more stringent. Only applies to genomic matches, not splice junctions. Set to 13\n" +
				"      or more to require near unique alignments.\n"+
				"-x Maximum mapping quality, reset reads with a mapping quality greater than the max to\n"+
				"      this max.\n"+
				"-n Maximum number of locations each read may align, defaults to 1 (unique matches).\n"+
				"-d If the maximum number of locations threshold fails, save one randomly picked repeat\n" +
				"      alignment per read.\n"+
				"-r Reverse the strand of the second paired alignment. Reversing the strand is\n" +
				"      needed for proper same strand visualization of paired stranded Illumina data.\n"+
				"-b Reverse the strand of both pairs.  Use this option if you would like the orientation\n" +
				"      of the alignments to match the orientation of the annotation in Illumina stranded \n" +
				"      UTP sequencing.\n" +
				"-u Save unmapped reads and those that fail the alignment score.\n"+
				"-c Don't remove chrAdapt and chrPhiX alignments.\n"+
				"-p Merge proper paired unique alignments. Those that cannot be unambiguously merged\n" +
				"      are left as pairs. Recommended to avoid double counting errors and increase\n" +
				"      base calling accuracy. For paired Illumina UTP data, use -p -r -b .\n"+
				"-q Maximum acceptable  base pair distance for merging, defaults to 300000.\n"+
				"-h Full path to a txt file containing a sam header, defaults to autogenerating the\n"+
				"      header from the read data.\n"+

				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/SamTranscriptomeParser -f /Novo/Run7/\n" +
				"     -m 20 -s /Novo/STPParsedBams/run7.bam -p -r \n\n" +

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
