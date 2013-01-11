package edu.utah.seq.data;

import java.io.*;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import net.sf.samtools.*;
import util.bio.parsers.MultiFastaParser;
import util.gen.*;

/** Calculates per cycle error rates.
 * @author Nix
 * */
public class CalculatePerCycleErrorRate {

	//user fields
	private File[] bamFiles;
	private File fastaFile;

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
		correctBases = new double[bamFiles.length][];
		incorrectBases = new double[bamFiles.length][];
		numberAlignments = new double[bamFiles.length];
		numberAlignmentsWithInsertions = new double[bamFiles.length];
		numberAlignmentsWithDeletions = new double[bamFiles.length];
		numberAlignmentsWithSoftMasking = new double[bamFiles.length];
		numberAlignmentsWithHardMasking = new double[bamFiles.length];
		
		//process each bam file
		for (int i=0; i< bamFiles.length; i++){
			fileIndex = i;
			scanBamFile();
		}
		
		printReport();

	}

	private void printReport() {
		
		System.out.println("Per cycle error rates for aligned bases:\n");
		//print header
		System.out.print("Cycle#");
		for (int i=0; i< bamFiles.length; i++) System.out.print("\t"+Misc.removeExtension(bamFiles[i].getName()));
		System.out.println();
		
		//make array lists to hold errors
		ArrayList<Double>[] fractionError = new ArrayList[bamFiles.length];
		for (int i=0; i< bamFiles.length; i++) fractionError[i] = new ArrayList<Double>();
		
		//for each row that contains any data correctBases[fileIndex][0-1000]
		for (int i=0; i< correctBases[0].length; i++){
			//skip?
			double total = 0;
			for (int j=0; j< bamFiles.length; j++) {
				total+= correctBases[j][i];
				total+= incorrectBases[j][i];
			}
			if (total == 0.0) continue;
			
			StringBuilder sb = new StringBuilder();
			sb.append(i+1);
			//for each dataset
			for (int j=0; j< bamFiles.length; j++){
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
		double[] averageError = new double[bamFiles.length];
		for (int i=0; i< bamFiles.length; i++) {
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
		
		
	}
	

	public void scanBamFile(){
		
		correctBases[fileIndex] = new double[1000];
		incorrectBases[fileIndex] = new double[1000];
		
		SAMFileReader samReader = null;

			try {
				samReader = new SAMFileReader(bamFiles[fileIndex]);
				//calc number cycles for first and second reads (second read might be null)
				numberCycles = estimateNumberOfCycles(samReader);
				
				SAMRecordIterator it = samReader.queryOverlapping(chromName, 0, chromLength);

				//int counter = 0;
				while (it.hasNext()) {
					SAMRecord sam = it.next();

					//is it aligned?
					if (sam.getReadUnmappedFlag()) continue;

					//does it pass the vendor qc?
					if (sam.getReadFailsVendorQualityCheckFlag()) continue;
					
					//if (sam.getReadNegativeStrandFlag() == false) continue;
					//if (sam.getFirstOfPairFlag()) continue;
					//second pair, reverse all
					
					//if (sam.getCigarString().equals("101M") == false) continue;
					//if (sam.getReadString().contains("N") == false) continue;
					
					int numberMismatches = scoreAlignment(sam);
					
					//scoreAlignmentDebug(sam);
					//if (counter++ == 9) break;
					
					
					
					//if (numberMismatches > 10) scoreAlignmentDebug(sam);
					
					

				}
				samReader.close();
			} catch (Exception e){
				System.err.println("\nError parsing sam file or writing split binary chromosome files.\n\nToo many open files exception? Too many chromosomes? " +
				"If so then login as root and set the default higher using the ulimit command (e.g. ulimit -n 10000)\n");
				e.printStackTrace();
				System.exit(1);
			}
	}


	private int[] estimateNumberOfCycles(SAMFileReader samReader) {
		SAMRecordIterator it = samReader.queryOverlapping(chromName, 0, chromLength);
		int counter = 0;
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
					if (maxSecondPair < cycles) maxSecondPair = cycles;
				}
				else if (maxFirstPair < cycles) maxFirstPair = cycles;
			}
			//exit
			if (counter++ > 10000) break;

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
			if (sam.getFirstOfPairFlag()) cycleSubtractor = numberCycles[0];
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
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bamFiles = IO.extractFiles(new File (args[++i]), ".bam"); break;
					case 'f': fastaFile = new File (args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//check bam files
		if (bamFiles == null || bamFiles.length == 0 || bamFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.bam file(s)!\n");
		OverdispersedRegionScanSeqs.lookForBaiIndexes(bamFiles, false);
		
		//check fasta
		if (fastaFile == null || fastaFile.canRead() == false) Misc.printExit("\nError: cannot find or read your fasta file!\n");
		
		//check results
		if (fastaFile == null ) Misc.printExit("\nError: please enter a file for saving the results!\n");

	}	
	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Calculate Per Cycle Error Rate : Jan 2013                 **\n" +
				"**************************************************************************************\n" +
				"Calculates per cycle error rates provided a sorted indexed bam file and a fasta\n" +
				"sequence file. Only checks CIGAR M bases not masked or INDEL bases.\n\n" +

				"Required Options:\n"+
				"-b Full path to a coordinate sorted bam file (xxx.bam) or directory containing such.\n" +
				"      Multiple files are processed independently.\n"+
				"-f Full path to the fasta file you wish to use in calculating the per cycle error rate.\n" +

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/CalculatePerCycleErrorRate -b /Data/Bam/\n"+
				"     -f /Fastas/chrPhiX_Illumina.fasta.gz \n\n"+

		"**************************************************************************************\n");

	}
		

}
