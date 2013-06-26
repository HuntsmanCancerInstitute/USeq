package edu.utah.seq.data;

import java.io.*;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import edu.utah.seq.data.sam.MalformedSamAlignmentException;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.data.sam.SamAlignment;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import util.bio.parsers.MultiFastaParser;
import util.gen.*;

/** Calculates per cycle error rates.
 * @author Nix
 * */
public class CalculatePerCycleErrorRate {

	//user fields
	private File[] alignmentFiles;
	private File fastaFile;
	private File logFile;
	private String readNamePrefix = null;

	//internal
	private String chromName;
	private int chromLength;
	private Pattern cigarSub = Pattern.compile("(\\d+)([MSDHIN])");
	private char[] chromSeq;
	private int[] numberCycles;
	private double[] numberAlignments;
	private double[] numberAlignmentsWithInsertions;
	private double[] numberAlignmentsWithDeletions;
	private double[] numberAlignmentsWithSoftMasking;
	private double[] numberAlignmentsWithHardMasking;
	private double[][] correctBases;
	private double[][] incorrectBases;
	private int fileIndex;
	private boolean isBam;
	private boolean deleteTempFiles = true;
	private Pattern headerLine = Pattern.compile("^@[HSRPC][DQGO]\\s.+");
	

	/**For stand alone app.*/
	public CalculatePerCycleErrorRate(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		//process args
		processArgs(args);

		doWork();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	public void doWork(){

		//fetch seq
		MultiFastaParser p = new MultiFastaParser(fastaFile);
		chromSeq = p.getSeqs()[0].toUpperCase().toCharArray();
		if (p.isFastaFound() == false) Misc.printErrAndExit("\nError: something is wrong with your fasta sequence.  Check that it is correct.\n");
		chromLength = chromSeq.length+1;
		chromName = p.getNames()[0];

		//make arrays to hold data
		correctBases = new double[alignmentFiles.length][];
		incorrectBases = new double[alignmentFiles.length][];
		numberAlignments = new double[alignmentFiles.length];
		numberAlignments = new double[alignmentFiles.length];
		numberAlignmentsWithInsertions = new double[alignmentFiles.length];
		numberAlignmentsWithDeletions = new double[alignmentFiles.length];
		numberAlignmentsWithSoftMasking = new double[alignmentFiles.length];
		numberAlignmentsWithHardMasking = new double[alignmentFiles.length];

		//process each file
		System.out.println("Processing...");
		for (int i=0; i< alignmentFiles.length; i++){
			fileIndex = i;
			correctBases[fileIndex] = new double[1000];
			incorrectBases[fileIndex] = new double[1000];
			System.out.println("\t"+ alignmentFiles[fileIndex].getName());
			isBam = alignmentFiles[fileIndex].getName().endsWith(".bam");
			if (isBam) {
				scanBamFile(alignmentFiles[fileIndex]);
			}
			else {
				File tempSam = parseSamFile();
				if (tempSam == null) {
					System.out.println("\t\tError! No "+chromName+" matching alignments? Skipping.");
				}
				else {
					try {
						String name = Misc.removeExtension(tempSam.getCanonicalPath());
						File tempBam = new File (name+".bam");
						File tempBai = new File (name+".bai");
						if (deleteTempFiles) {
							tempBam.deleteOnExit();
							tempBai.deleteOnExit();
						}
						//sort and convert to BAM
						new PicardSortSam (tempSam, tempBam);
						//scan it!
						scanBamFile(tempBam);

					} catch (IOException e) {
						System.out.println("\nProblem sorting temp sam file "+tempSam);
						e.printStackTrace();
					}
				}
			}
		}
		printReport();
	}

	private void printReport() {
		PrintStream oldStream = System.out;
		try {
			if (logFile != null) {
				System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream(logFile))));
			}
	
			System.out.println("Per cycle error rates for aligned bases:\n");
			//print header
			System.out.print("Cycle#");
			for (int i=0; i< alignmentFiles.length; i++) System.out.print("\t"+Misc.removeExtension(alignmentFiles[i].getName()));
			System.out.println();
	
			//make array lists to hold errors
			ArrayList<Double>[] fractionError = new ArrayList[alignmentFiles.length];
			for (int i=0; i< alignmentFiles.length; i++) fractionError[i] = new ArrayList<Double>();
	
			//for each row that contains any data correctBases[fileIndex][0-1000]
			for (int i=0; i< correctBases[0].length; i++){
				//skip?
				double total = 0;
				for (int j=0; j< alignmentFiles.length; j++) {
					total+= correctBases[j][i];
					total+= incorrectBases[j][i];
				}
				if (total == 0.0) continue;
	
				StringBuilder sb = new StringBuilder();
				sb.append(i+1);
				//for each dataset
				for (int j=0; j< alignmentFiles.length; j++){
					sb.append("\t");
					double error = incorrectBases[j][i]/ (incorrectBases[j][i] + correctBases[j][i]);
					sb.append(error);
					fractionError[j].add(error);
					/*sb.append("\t");
					sb.append(incorrectBases[j][i]);
					sb.append("\t");
					sb.append(correctBases[j][i]);
					 */
				}
				System.out.println(sb);
			}
	
			//calculate average error
			double[] averageError = new double[alignmentFiles.length];
			for (int i=0; i< alignmentFiles.length; i++) {
				double[] err = Num.arrayListOfDoubleToArray(fractionError[i]);
				averageError[i] = Num.mean(err);
			}
			System.out.println("\nMean\t"+Num.doubleArrayToString(averageError, "\t"));
	
			//print alignment stats
			for (int i=0; i<numberAlignmentsWithInsertions.length; i++ ) numberAlignmentsWithInsertions[i] = 100* numberAlignmentsWithInsertions[i]/numberAlignments[i];
			for (int i=0; i<numberAlignmentsWithDeletions.length; i++ ) numberAlignmentsWithDeletions[i] = 100* numberAlignmentsWithDeletions[i]/numberAlignments[i];
			for (int i=0; i<numberAlignmentsWithSoftMasking.length; i++ ) numberAlignmentsWithSoftMasking[i] = 100* numberAlignmentsWithSoftMasking[i]/numberAlignments[i];
			for (int i=0; i<numberAlignmentsWithHardMasking.length; i++ ) numberAlignmentsWithHardMasking[i] = 100* numberAlignmentsWithHardMasking[i]/numberAlignments[i];
	
			System.out.println("\n# Alignments\t"+Num.doubleArrayToString(numberAlignments, 0, "\t"));
			System.out.println("% WithInsertions\t"+Num.doubleArrayToString(numberAlignmentsWithInsertions, 3, "\t"));
			System.out.println("% WithDeletions\t"+Num.doubleArrayToString(numberAlignmentsWithDeletions, 3, "\t"));
			System.out.println("% WithSoftMasking\t"+Num.doubleArrayToString(numberAlignmentsWithSoftMasking, 3, "\t"));
			System.out.println("% WithHardMasking\t"+Num.doubleArrayToString(numberAlignmentsWithHardMasking, 3, "\t"));
			
			System.out.close();
			System.setOut(oldStream);
		} catch (FileNotFoundException ex) {
			System.out.println("Could not create the log file: " + ex.getMessage());
			ex.printStackTrace();
			System.exit(1);
		}


	}


	public void scanBamFile(File bamFile){
		numberCycles = estimateNumberOfCyclesBam(bamFile);

		SAMFileReader samReader = null;

		try {
			samReader = new SAMFileReader(bamFile);

			SAMRecordIterator it;
			if (samReader.hasIndex()) it = samReader.queryOverlapping(chromName, 0, chromLength);
			else it = samReader.iterator();

			int counter = 0;
			while (it.hasNext()) {
				counter++;
				SAMRecord sam = it.next();

				//is it aligned?
				if (sam.getReadUnmappedFlag()) continue;

				//does it pass the vendor qc?
				if (sam.getReadFailsVendorQualityCheckFlag()) continue;
				
				//prefix?
				if (readNamePrefix != null) {
					if (sam.getReadName().startsWith(readNamePrefix) == false) continue;
				}

				
				scoreAlignment(sam);
					
			}
			samReader.close();
		} catch (Exception e){
			System.err.println("\nError parsing sam file or writing split binary chromosome files.\n\nToo many open files exception? Too many chromosomes? " +
			"If so then login as root and set the default higher using the ulimit command (e.g. ulimit -n 10000)\n");
			e.printStackTrace();
			System.exit(1);
		}
	}

	public File parseSamFile(){
		BufferedReader in = null;
		PrintWriter samOut = null;
		int numBadLines = 0;
		File samOutFile = null;
		int numberAlignments = 0;
		try {
			//make print writer
			samOutFile = new File(Misc.removeExtension(alignmentFiles[fileIndex].getCanonicalPath())+"_temp.sam");
			if (deleteTempFiles) samOutFile.deleteOnExit();
			samOut = new PrintWriter( new FileWriter (samOutFile));
			in = IO.fetchBufferedReader(alignmentFiles[fileIndex]);
			String line;
			String sqMatcher = "SN:"+chromName;
			while ((line=in.readLine())!= null) {
				line = line.trim();
				if (line.length() == 0) continue;
				//is it a header line?
				if (headerLine.matcher(line).matches()){
					//an SQ?
					if (line.startsWith("@SQ")){
						if (line.contains(sqMatcher)) samOut.println(line);
					}
					else samOut.println(line);
					continue;
				}

				SamAlignment sa;
				try {
					sa = new SamAlignment(line, true);
				} catch (Exception e) {
					System.err.println("Skipping malformed sam alignment -> "+e.getMessage());
					if (numBadLines++ > 1000) Misc.printErrAndExit("Aboring: too many malformed SAM alignments");
					continue;
				}

				//is it aligned? failed QC?
				if (sa.isUnmapped() || sa.failedQC()) continue;

				//map to fasta
				if (sa.getReferenceSequence().equals(chromName) == false) continue;
				
				//prefix?
				if (readNamePrefix != null) {
					if (sa.getName().startsWith(readNamePrefix) == false) continue;
				}

				samOut.println(line);
				numberAlignments++;

			}
		} catch (Exception e) {
			System.err.println("\nError parsing SAM file.\n");
			e.printStackTrace();
		} finally {
			try {
				if (in != null) in.close();
				if (samOut != null) samOut.close();
			} catch (IOException e) {}
		}
		//any alignments?
		if (numberAlignments == 0) return null;
		return samOutFile;
	}


	private int[] estimateNumberOfCyclesBam(File bamFile) {
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader samReader = new SAMFileReader(bamFile);
		//samReader.setValidationStringency(ValidationStringency.SILENT);
		SAMRecordIterator it;
		if (samReader.hasIndex()) it = samReader.queryOverlapping(chromName, 0, chromLength);
		else it = samReader.iterator();
		int counter1 = 0;
		int counter2 = 0;
		int maxFirstPair = 0;
		int maxSecondPair = 0;
		Pattern pat = Pattern.compile("(\\d+)M");
		while (it.hasNext()) {
			SAMRecord sam = it.next();
		
			//is it aligned?
			if (sam.getReadUnmappedFlag()) continue;

			//does it pass the vendor qc?
			if (sam.getReadFailsVendorQualityCheckFlag()) continue;

			Matcher mat = pat.matcher(sam.getCigarString());
			if (mat.matches()){
				int cycles = Integer.parseInt(mat.group(1));
				if (sam.getReadNegativeStrandFlag()) {
					if (maxSecondPair < cycles) {
						maxSecondPair = cycles;
						counter1++;
					}

				}
				else if (maxFirstPair < cycles) {
					maxFirstPair = cycles;
					counter2++;
				}
			}
			//Removing exit condition.  There are instances were the first 10000 reads are read1 and when a read2 is encountered, the program crashes.
			if (counter1 > 10000 && counter2 > 10000) break;
			

		}
		it.close();
		//need to subtract one to keep in frame
		return new int[]{maxFirstPair-1, maxSecondPair-1};
	}

	private int scoreAlignment(SAMRecord sam) {
		//get start and end, watch out for those that are off ends
		int start = sam.getUnclippedStart() -1;
		if (start < 0) return -1;
		int end = sam.getAlignmentEnd();
		if (end >= chromLength) return -2;

		numberAlignments[fileIndex]++;

		//need to flip index for negative strand reads!
		int cycleSubtractor = 0;
		if (sam.getReadNegativeStrandFlag()){
			if (sam.getReadPairedFlag() == false || sam.getFirstOfPairFlag()) cycleSubtractor = numberCycles[0];
			else cycleSubtractor = numberCycles[1];
		}

		//walk through Ms
		//for each cigar block
		Matcher mat = cigarSub.matcher(sam.getCigarString());
		char[] seq = sam.getReadString().toUpperCase().toCharArray();
		int seqIndex =0;
		boolean insertionFound = false;
		boolean deletionFound = false;
		boolean softMaskingFound = false;
		boolean hardMaskingFound = false;
		int numberMismatches = 0;
		
		

		while (mat.find()){
			String call = mat.group(2);
			int numberBases = Integer.parseInt(mat.group(1));
			
			
			//a match
			if (call.equals("M")) {
				for (int i = 0; i< numberBases; i++) {
					if (chromSeq[start] == seq[seqIndex]) {
						if (cycleSubtractor !=0) correctBases[fileIndex][cycleSubtractor - seqIndex]++;
						else correctBases[fileIndex][seqIndex]++;
					}
					else if (seq[seqIndex] != 'N') {
						if (cycleSubtractor !=0) incorrectBases[fileIndex][cycleSubtractor - seqIndex]++;
						else incorrectBases[fileIndex][seqIndex]++;
						numberMismatches++;
					}
					start++;
					seqIndex++;
				}
			}
			//a hard mask
			else if (call.equals("H")){
				//just increment start
				start+= numberBases;
				hardMaskingFound = true;
			}
			//a soft mask
			else if (call.equals("S")){
				//advance both
				start += numberBases;
				seqIndex += numberBases;
				softMaskingFound = true;
			}
			// a deletion
			else if (call.equals("D")){
				//just increment start
				start+= numberBases;
				deletionFound = true;
			}
			else if (call.equals("I")){
				//just advance seqIndex
				seqIndex += numberBases;
				insertionFound = true;
			}
			else Misc.printErrAndExit("\nError: unsupported CIGAR string see -> \n"+sam.getSAMString()+"\n");

		}

		//increment counters
		if (insertionFound) numberAlignmentsWithInsertions[fileIndex]++;
		if (deletionFound) numberAlignmentsWithDeletions[fileIndex]++;
		if (hardMaskingFound) numberAlignmentsWithHardMasking[fileIndex]++;
		if (softMaskingFound) numberAlignmentsWithSoftMasking[fileIndex]++;

		return numberMismatches;

	}

	private int scoreAlignmentDebug(SAMRecord sam) {
		//get start and end
		int start = sam.getUnclippedStart() -1;
		if (start < 0) return -1;
		int end = sam.getAlignmentEnd();
		if (end >= chromLength) return -2;
		System.out.println(sam.getReadName() + " "+ sam.getCigarString() +" "+start);
		System.out.println(sam.getReadString());
		for (int x=start; x< end; x++) System.out.print(chromSeq[x]);
		System.out.println();

		//walk through Ms
		//for each cigar block
		Matcher mat = cigarSub.matcher(sam.getCigarString());
		char[] seq = sam.getReadString().toUpperCase().toCharArray();
		int seqIndex =0;
		int numberMismatches = 0;
		while (mat.find()){
			String call = mat.group(2);
			int numberBases = Integer.parseInt(mat.group(1));
			//a match
			if (call.equals("M")) {
				for (int i = 0; i< numberBases; i++) {
					int match = -1;
					if (chromSeq[start] == seq[seqIndex]) {
						match = 0;
					}
					else if (seq[seqIndex] != 'N') {
						numberMismatches++;
						match = 1;
					}
					System.out.println("\t"+ seqIndex+"\t"+chromSeq[start]+"\t"+seq[seqIndex]+"\t"+match);
					start++;
					seqIndex++;
				}
			}
			//a hard mask
			else if (call.equals("H")){
				//just increment start
				start+= numberBases;
			}
			//a soft mask
			else if (call.equals("S")){
				//advance both
				start += numberBases;
				seqIndex += numberBases;
			}
			// a deletion
			else if (call.equals("D")){
				//just increment start
				start+= numberBases;
			}
			else if (call.equals("I")){
				//just advance seqIndex
				seqIndex += numberBases;
			}
			else Misc.printErrAndExit("\nError: unsupported CIGAR string see -> \n"+sam.getSAMString()+"\n");

		}

		return numberMismatches;

	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CalculatePerCycleErrorRate(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': forExtraction = new File(args[i+1]); i++; break;
					case 'f': fastaFile = new File (args[++i]); break;
					case 't': deleteTempFiles = false; break;
					case 'n': readNamePrefix = args[++i]; break;
					case 'o': logFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//pull files
		if (forExtraction == null ) Misc.printExit("\nError: cannot find your xxx.bam file(s) or xxx.sam(.zip/.gz)!\n");
		File[][] tot = new File[4][];
		tot[0] = IO.extractFiles(forExtraction,".sam");
		tot[1] = IO.extractFiles(forExtraction,".sam.gz");
		tot[2] = IO.extractFiles(forExtraction,".sam.zip");
		tot[3] = IO.extractFiles(forExtraction,".bam");

		alignmentFiles = IO.collapseFileArray(tot);
		if (alignmentFiles == null || alignmentFiles.length==0) alignmentFiles = IO.extractFiles(forExtraction);
		if (alignmentFiles == null || alignmentFiles.length ==0 || alignmentFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.bam file(s) or xxx.sam(.zip/.gz)!\n");

		//check fasta
		if (fastaFile == null || fastaFile.canRead() == false) Misc.printExit("\nError: cannot find or read your fasta file!\n");
		if (fastaFile.isDirectory()) Misc.printExit("\nError: please provide a single fasta file for scoring. Not a directory.\n");

		//check results
		if (fastaFile == null ) Misc.printExit("\nError: please enter a file for saving the results!\n");

	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Calculate Per Cycle Error Rate : Feb 2013                 **\n" +
				"**************************************************************************************\n" +
				"Calculates per cycle error rates provided a sorted indexed bam file and a fasta\n" +
				"sequence file. Only checks CIGAR M bases not masked or INDEL bases.\n\n" +

				"Required Options:\n"+
				"-b Full path to a coordinate sorted bam file (xxx.bam) with its associated (xxx.bai)\n" +
				"      index or directory containing such. Multiple files are processed independently.\n" +
				"      Unsorted xxx.sam(.gz/.zip OK) files also work but are processed rather slowly.\n"+
				"-f Full path to the single fasta file you wish to use in calculating the error rate.\n" +
				"-n Require read names to begin with indicated text, defaults to accepting everything.\n"+
				"-o Path to log file.  Write coverage statistics to a log file instead of stdout.\n" +

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/CalculatePerCycleErrorRate -b /Data/Bam/\n"+
				"     -f /Fastas/chrPhiX_Illumina.fasta.gz -n HWI\n\n"+

		"**************************************************************************************\n");

	}


}
