package edu.utah.seq.data;

import java.io.*;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import edu.utah.seq.data.sam.MalformedSamAlignmentException;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.data.sam.SamAlignment;
import htsjdk.samtools.*;
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
	private double firstReadMaximumError = 0.01;
	private double secondReadMaximumError = 0.0175;
	private double maximumFractionFailingCycles = 0.1;

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
	private int index;
	private boolean isBam;
	private boolean deleteTempFiles = true;
	private Pattern headerLine = Pattern.compile("^@[HSRPC][DQGO]\\s.+");
	private boolean mergeStrands = true;
	private double[] numberFailingCycles;
	private int totalCycles;

	/**For stand alone app.*/
	public CalculatePerCycleErrorRate(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		//process args
		processArgs(args);

		doWork();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Num.formatNumber(diffTime, 2)+" min\n");
	}

	public void doWork(){

		//fetch seq
		MultiFastaParser p = new MultiFastaParser(fastaFile);
		chromSeq = p.getSeqs()[0].toUpperCase().toCharArray();
		if (p.isFastaFound() == false) Misc.printErrAndExit("\nError: something is wrong with your fasta sequence.  Check that it is correct.\n");
		chromLength = chromSeq.length+1;
		chromName = p.getNames()[0];

		int numberDatasets = alignmentFiles.length;
		if (mergeStrands== false) numberDatasets = numberDatasets*2;
		
		//make arrays to hold data
		correctBases = new double[numberDatasets][];
		incorrectBases = new double[numberDatasets][];
		for (int i=0; i< numberDatasets; i++){
			correctBases[i] = new double[1000];
			incorrectBases[i] = new double[1000];
		}
		numberAlignments = new double[numberDatasets];
		numberAlignments = new double[numberDatasets];
		numberAlignmentsWithInsertions = new double[numberDatasets];
		numberAlignmentsWithDeletions = new double[numberDatasets];
		numberAlignmentsWithSoftMasking = new double[numberDatasets];
		numberAlignmentsWithHardMasking = new double[numberDatasets];

		//process each file
		System.out.println("Processing...");
		for (int i=0; i< alignmentFiles.length; i++){
			fileIndex = i;
			if (mergeStrands == false) index = i*2;
			else index = fileIndex;
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
			int numberDatasets = alignmentFiles.length;
			if (mergeStrands == false) numberDatasets = numberDatasets*2;
	
			System.out.println("Per cycle error rates for aligned bases:\n");
			//print header
			String[] names = new String[numberDatasets];
			int nameIndex = 0;
			System.out.print("Cycle#");
			for (int i=0; i< alignmentFiles.length; i++) {
				String name = Misc.removeExtension(alignmentFiles[i].getName());
				if (mergeStrands) {
					System.out.print("\t"+name);
					names[nameIndex++] = name;
				}
				else {
					System.out.print("\t"+name+"_1\t"+name+"_2");
					names[nameIndex++] = name;
					names[nameIndex++] = name;
				}
			}
			System.out.println();
	
			//make array lists to hold errors
			ArrayList<Double>[] fractionError = new ArrayList[numberDatasets];
			for (int i=0; i< numberDatasets; i++) fractionError[i] = new ArrayList<Double>();
			numberFailingCycles = new double[numberDatasets];
			totalCycles = 0;
			//for each row that contains any data correctBases[fileIndex][0-1000]
			for (int i=0; i< correctBases[0].length; i++){
				//for each dataset sum correct and in correct bases
				double total = 0;
				for (int j=0; j< numberDatasets; j++) {
					total+= correctBases[j][i];
					total+= incorrectBases[j][i];
				}
				if (total == 0.0) continue;
				totalCycles++;
				
				StringBuilder sb = new StringBuilder();
				sb.append(i+1);
				//for each dataset
				for (int j=0; j< numberDatasets; j++){
					sb.append("\t");
					double error = incorrectBases[j][i]/ (incorrectBases[j][i] + correctBases[j][i]);
					sb.append(error);
					fractionError[j].add(error);
					/*sb.append("\t");
					sb.append(incorrectBases[j][i]);
					sb.append("\t");
					sb.append(correctBases[j][i]);*/
					
					//increment failed cycles?
					boolean isEven = Num.isEven(j);
					if (mergeStrands || isEven){
						//use first strand error rate
						if (error>= firstReadMaximumError) numberFailingCycles[j]++;
					}
					//use second
					else if (error>= secondReadMaximumError) numberFailingCycles[j]++;						
				}
				System.out.println(sb);
			}
	
			//calculate average error
			double[] averageError = new double[numberDatasets];
			for (int i=0; i< numberDatasets; i++) {
				double[] err = Num.arrayListOfDoubleToArray(fractionError[i]);
				averageError[i] = Num.mean(err);
			}
			System.out.println("\nMean\t"+Num.doubleArrayToString(averageError, "\t"));
			
			//calc failing 
			LinkedHashSet<String> badGuys = new LinkedHashSet<String>();
			for (int i=0; i< numberDatasets; i++){
				numberFailingCycles[i] = numberFailingCycles[i]/totalCycles;
				if (numberFailingCycles[i] > maximumFractionFailingCycles) badGuys.add(names[i]);
			}
			System.out.println("\nFractionFailingCycles\t"+Num.doubleArrayToString(numberFailingCycles, 4, "\t"));
	
			//print alignment stats
			for (int i=0; i<numberAlignmentsWithInsertions.length; i++ ) numberAlignmentsWithInsertions[i] = numberAlignmentsWithInsertions[i]/numberAlignments[i];
			for (int i=0; i<numberAlignmentsWithDeletions.length; i++ ) numberAlignmentsWithDeletions[i] = numberAlignmentsWithDeletions[i]/numberAlignments[i];
			for (int i=0; i<numberAlignmentsWithSoftMasking.length; i++ ) numberAlignmentsWithSoftMasking[i] = numberAlignmentsWithSoftMasking[i]/numberAlignments[i];
			for (int i=0; i<numberAlignmentsWithHardMasking.length; i++ ) numberAlignmentsWithHardMasking[i] = numberAlignmentsWithHardMasking[i]/numberAlignments[i];
	
			System.out.println("\n# Alignments\t"+Num.doubleArrayToString(numberAlignments, 0, "\t"));
			System.out.println("FractionWithInsertions\t"+Num.doubleArrayToString(numberAlignmentsWithInsertions, 4, "\t"));
			System.out.println("FractionWithDeletions\t"+Num.doubleArrayToString(numberAlignmentsWithDeletions, 4, "\t"));
			System.out.println("FractionWithSoftMasking\t"+Num.doubleArrayToString(numberAlignmentsWithSoftMasking, 4, "\t"));
			System.out.println("FractionWithHardMasking\t"+Num.doubleArrayToString(numberAlignmentsWithHardMasking, 4, "\t"));
			
			
			System.out.println("\nDatasets failing thresholds ("+Num.formatNumber(maximumFractionFailingCycles, 4)+
					" maxFracFailingCycles, "+Num.formatNumber(firstReadMaximumError, 4)+
					" maxFirstReadOrMergeError, "+Num.formatNumber(secondReadMaximumError, 4)+
					" maxSecondReadError):");
			if (badGuys.size() != 0) System.out.println(Misc.hashSetToString(badGuys, ", "));
			else System.out.println("None");
			
			if (logFile != null) {
				System.out.close();
				System.setOut(oldStream);
			}
			
		} catch (FileNotFoundException ex) {
			System.err.println("Could not create the log file: " + ex.getMessage());
			ex.printStackTrace();
			System.exit(1);
		}


	}


	public void scanBamFile(File bamFile){
		
		SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SamReader samReader = factory.open(bamFile);

		try {
			//look for chrom to scan
			if (samReader.getFileHeader().getSequence(chromName) == null) {
				System.err.println ("\nError: failed to find any alignments for '"+chromName+"' in "+bamFile.getName()+", aborting!\n");
				System.exit(0); //don't want to exit with 1 or tomato will crash out.
			}
			numberCycles = estimateNumberOfCyclesBam(bamFile);
			
			SAMRecordIterator it = samReader.queryOverlapping(chromName, 0, chromLength);

			while (it.hasNext()) {
				SAMRecord sam = it.next();

				//is it aligned?
				if (sam.getReadUnmappedFlag()) continue;
				
				if (sam.getNotPrimaryAlignmentFlag()) continue;

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

	/**Returns null if the chromName doesn't exist or no records were found.
	 * @throws Exception */
	private int[] estimateNumberOfCyclesBam(File bamFile) throws Exception {
		SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SamReader samReader = factory.open(bamFile);
		SAMRecordIterator it = samReader.queryOverlapping(chromName, 0, chromLength);
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
		samReader.close();
		//need to subtract one to keep in frame
		return new int[]{maxFirstPair-1, maxSecondPair-1};
	}

	private int scoreAlignment(SAMRecord sam) {
		int localIndex = index;
		if (mergeStrands == false && sam.getReadPairedFlag() == true && sam.getSecondOfPairFlag() == true) localIndex++;

		//get start and end, watch out for those that are off ends
		int start = sam.getUnclippedStart() -1;
		if (start < 0) return -1;
		int end = sam.getAlignmentEnd();
		if (end >= chromLength) return -2;

		numberAlignments[localIndex]++;

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
						if (cycleSubtractor !=0) correctBases[localIndex][cycleSubtractor - seqIndex]++;
						else correctBases[localIndex][seqIndex]++;
					}
					else if (seq[seqIndex] != 'N') {
						if (cycleSubtractor !=0) incorrectBases[localIndex][cycleSubtractor - seqIndex]++;
						else incorrectBases[localIndex][seqIndex]++;
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
		if (insertionFound) numberAlignmentsWithInsertions[localIndex]++;
		if (deletionFound) numberAlignmentsWithDeletions[localIndex]++;
		if (hardMaskingFound) numberAlignmentsWithHardMasking[localIndex]++;
		if (softMaskingFound) numberAlignmentsWithSoftMasking[localIndex]++;

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
		Pattern pat = Pattern.compile("-[a-z12]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': alignmentFiles = IO.extractFiles(new File(args[++i]),".bam"); break;
					case 'f': fastaFile = new File (args[++i]); break;
					case 't': deleteTempFiles = false; break;
					case 's': mergeStrands = false; break;
					case 'n': readNamePrefix = args[++i]; break;
					case 'o': logFile = new File(args[++i]); break;
					case '1': firstReadMaximumError = Double.parseDouble(args[++i]); break;
					case '2': secondReadMaximumError = Double.parseDouble(args[++i]); break;
					case 'c': maximumFractionFailingCycles = Double.parseDouble(args[++i]); break;
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
		if (alignmentFiles == null || alignmentFiles.length == 0) Misc.printExit("\nError: cannot find your xxx.bam file(s)!\n");
		//look for indexes
		OverdispersedRegionScanSeqs.lookForBaiIndexes(alignmentFiles, false);
		
		//check fasta
		if (fastaFile == null || fastaFile.canRead() == false) Misc.printExit("\nError: cannot find or read your fasta file!\n");
		if (fastaFile.isDirectory()) Misc.printExit("\nError: please provide a single fasta file for scoring. Not a directory.\n");

		//check results
		if (fastaFile == null ) Misc.printExit("\nError: please enter a file for saving the results!\n");

	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Calculate Per Cycle Error Rate : March 2015                **\n" +
				"**************************************************************************************\n" +
				"Calculates per cycle error rates provided a sorted indexed bam file and a fasta\n" +
				"sequence file. Only checks CIGAR M bases not masked or INDEL bases.\n\n" +

				"Required Options:\n"+
				"-b Full path to a coordinate sorted bam file (xxx.bam) with its associated (xxx.bai)\n" +
				"      index or directory containing such. Multiple files are processed independently.\n" +
				"-f Full path to the single fasta file you wish to use in calculating the error rate.\n" +
				"-n Require read names to begin with indicated text, defaults to accepting everything.\n"+
				"-o Path to log file.  Write coverage statistics to a log file instead of stdout.\n" +
				"-s Perform separate first read second read analysis, defaults to merging.\n"+
				"-c Maximum fraction failing cycles, defaults to 0.1\n"+
				"-1 Maximum first read or merged read error rate, defaults to 0.01\n"+
				"-2 Maximum second read error rate, defaults to 0.0175\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/CalculatePerCycleErrorRate -b /Data/Bam/\n"+
				"     -f /Fastas/chrPhiX_Illumina.fasta.gz -n HWI\n\n"+

		"**************************************************************************************\n");

	}
}
