package edu.utah.diagnostics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;
import util.gen.Gzipper;
import util.gen.Histogram;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;
import edu.utah.seq.data.sam.MalformedSamAlignmentException;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.data.sam.SamAlignmentFlags;
import edu.utah.seq.data.sam.SamLayout;
import edu.utah.seq.parsers.MergePairedSamAlignments;

public class GenerateOverlapStats {
	//user defined fields
		private File[] dataFiles;
		private File saveFile;
		private File outputFile;
		private float maximumAlignmentScore = 120;
		private float minimumMappingQualityScore = 0;
		private boolean secondPairReverseStrand = false;
		private boolean removeControlAlignments = false;
		private boolean skipMergingPairs = false;
		private boolean onlyMergeOverlappingAlignments = true;
		private int minimumDiffQualScore = 3;
		private double minimumFractionInFrameMismatch = 0.05;
		private int maximumProperPairDistanceForMerging = 5000;
		private File logFile=null;
		private boolean writeSam = true;

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
		public GenerateOverlapStats(String[] args){
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
			
			//for each file, parse and save to disk	
			System.out.println("\nParsing, filtering, and merging BAM files...");
			for (int i=0; i< dataFiles.length; i++){
				System.out.print("\t"+dataFiles[i].getName());
				parseBamFile(dataFiles[i]);
				System.out.println();
			}

			//stats
			PrintStream oldStream = System.out;
			try {
				if (logFile != null) {
					PrintStream customStream =new PrintStream(new FileOutputStream(logFile));
					System.setOut(customStream);
				}
				
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
				
				System.out.close();
				System.setOut(oldStream);
			} catch (FileNotFoundException ioex) {
				System.out.println("Could not write to log file: " + ioex.getMessage());
				ioex.printStackTrace();
				System.exit(1);
			}
		
			


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
				int timCounter = 0;
				int matched = 0;
				int unmatched = 0;

				samReader = new SAMFileReader(bamFile);
				
				HashMap<String,SamAlignment> storedReads = new HashMap<String,SamAlignment>();
				HashSet<String> skipped = new HashSet<String>();

				SAMRecordIterator it = samReader.iterator();

				while (it.hasNext()) {
					SAMRecord sam = it.next();
					line = sam.getSAMString().trim();
					
					if (++timCounter > 100000) {
						System.out.println(matched + "  " + unmatched + " " + storedReads.size() + " " + skipped.size());
//						for (SamAlignment s: storedReads.values()) {
//							if (sam.getAlignmentStart() - s.getPosition() > 1000) {
//								System.out.println(s.toString());
//							}
//						}
						timCounter = 0;
					}

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
					
					String readName = sa.getName();

					if (checkSamAlignment(sa, line) == false) {
						if (storedReads.containsKey(readName)) {
							storedReads.remove(readName);
						} else {
							if (skipped.contains(readName)) {
								skipped.remove(readName);
							} else {
								skipped.add(readName);
							}
							
						}
						continue;
					}
					
					if (skipped.contains(readName)) {
						skipped.remove(readName);
						continue;
					}

					
					if (Math.abs(sam.getInferredInsertSize()) > 10000 || sam.getReferenceIndex() != sam.getMateReferenceIndex() || sa.isMateUnMapped() || sa.isUnmapped()) {
						continue;
					}
					
					if (!storedReads.containsKey(readName)) {
						storedReads.put(readName,sa);
						unmatched += 1;
						continue;
					}
					
					matched += 1;
					
					//Grab first observed
					SamAlignment firstObserved = storedReads.get(readName);
					
					//determine first read
					if (firstObserved.isFirstPair()) {
						alignmentsToSave.add(firstObserved);
						alignmentsToSave.add(sa);
					} else {
						alignmentsToSave.add(sa);
						alignmentsToSave.add(firstObserved);
					}
				
					//merge stuff
					filterPrintAlignments(alignmentsToSave);
					
					//Clear array
					alignmentsToSave.clear();
					
					//Clean hashmap
					storedReads.remove(readName);
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
				return false;
			}

			//does it pass the vendor qc?
			if (sa.failedQC()) {
				numberFailingVendorQC++;
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
					return false;
				}
			}

			//check mapping quality for genomic match reads?
			if (minimumMappingQualityScore !=0){
				if (sa.getMappingQuality() < minimumMappingQualityScore){
					numberFailingMappingQualityScore++;
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
				numberPrintedAlignments++;
				return false;
			}

			//is it not part of a proper pair?
			if (!sa.isAProperPairedAlignment() && !sa.isUnmapped() && !sa.isMateUnMapped()){
				numberNonProperPairedAlignments++;
				numberPrintedAlignments++;
				return false;
			}
			//is mate unmapped
			if (sa.isMateUnMapped()){
				numberUnmappedMatePairedAlignments++;
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
			
			mergePairedAlignments(first, second);

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
		private void mergePairedAlignments(SamAlignment first, SamAlignment second) {
			
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

		public static void main(String[] args) {
			if (args.length ==0){
				printDocs();
				System.exit(0);
			}
			new GenerateOverlapStats(args);
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
						case 'm': crossCheckMateCoordinates = false; break;
						case 'r': secondPairReverseStrand = true; break;
						case 'd': maximumProperPairDistanceForMerging = Integer.parseInt(args[++i]); break;
						case 'a': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
						case 'q': minimumMappingQualityScore = Float.parseFloat(args[++i]); break;
						case 'l': logFile = new File(args[++i]); break;
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
			tot[0] = IO.extractFiles(forExtraction,".bam");

			dataFiles = IO.collapseFileArray(tot);
			if (dataFiles == null || dataFiles.length==0) dataFiles = IO.extractFiles(forExtraction);
			if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.bam file(s)!\n");

			//print info
			System.out.println(maximumAlignmentScore+ "\tMaximum alignment score.");
			System.out.println(minimumMappingQualityScore+ "\tMinimum mapping quality score.");
			System.out.println(removeControlAlignments +"\tRemove control chrPhiX and chrAdapter alignments.");
			System.out.println(secondPairReverseStrand +"\tSecond read pair's strand has been reversed.");
			System.out.println(crossCheckMateCoordinates +"\tCross check read mate coordinates.");
			System.out.println(maximumProperPairDistanceForMerging +"\tMaximum bp distance for merging paired alignments.");
		}	

		public static void printDocs(){
			System.out.println("\n" +
					"**************************************************************************************\n" +
					"**                          Generate Overlas: Dec 2012                      **\n" +
					"**************************************************************************************\n" +
					"Merges proper paired alignments that pass a variety of checks and thresholds. Only\n" +
					"unambiguous pairs will be merged. Increases base calling accuracy in overlap and helps\n" +
					"avoid non-independent variant observations and other double counting issues. Identical\n" +
					"overlapping bases are assigned the higher quality scores. Disagreements are resolved\n" +
					"toward the higher quality base. If too close in quality, then the quality is set to 0.\n" +
					"Be certain your input bam/sam file(s) are sorted by query name, NOT coordinate. \n" +

					"\nOptions:\n"+
					"-f The full path file or directory containing raw xxx.sam(.gz/.zip OK)/.bam file(s)\n" +
					"      paired alignments. \n" +
					"      Multiple files will be merged.\n" +

					"\nDefault Options:\n"+
					"-a Maximum alignment score (AS:i: tag). Defaults to 120, smaller numbers are more\n" +
					"      stringent. Approx 30pts per mismatch for novoalignments.\n"+
					"-q Minimum mapping quality score, defaults to 13, larger numbers are more stringent.\n" +
					"      Set to 0 if processing splice junction indexed RNASeq data.\n"+
					"-r The second paired alignment's strand is reversed. Defaults to not reversed.\n" +
					"-d Maximum acceptible base pair distance for merging, defaults to 5000.\n"+
					"-m Don't cross check read mate coordinates, needed for merging repeat matches. Defaults\n" +
					"      to checking.\n"+
					"-l Output file name.  Write merging statitics to file instead of standard output.\n" +

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
