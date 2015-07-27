package edu.utah.seq.data;

import java.io.*;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
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
	private double firstReadMaximumError = 0.01;
	private double secondReadMaximumError = 0.0175;
	private double maximumFractionFailingCycles = 0.1;
	private File jsonOutputFile = null;

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
	private double[][] perCycleErrorRatePerDataset;
	private int fileIndex;
	private int index;
	private boolean mergeStrands = true;
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
			scanBamFile(alignmentFiles[fileIndex]);
		}
		double totalAlignments = Num.sumArray(numberAlignments);
		if (totalAlignments == 0.0) System.out.println("\nNo "+chromName+" alignments found for analysis?!");
		else printReport();
		
		if (jsonOutputFile != null) saveJson();
	}

	private void printReport() {
		PrintStream oldStream = System.out;
		try {
			if (logFile != null)  System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream(logFile))));
			
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
			double[] numberFailingCycles = new double[numberDatasets];
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
				//cycle #
				sb.append(i+1);
				//for each dataset
				for (int j=0; j< numberDatasets; j++){
					sb.append("\t");
					double error = incorrectBases[j][i]/ (incorrectBases[j][i] + correctBases[j][i]);
					sb.append(error);
					fractionError[j].add(error);
					/*
					sb.append("\t");
					sb.append(incorrectBases[j][i]);
					sb.append("\t");
					sb.append(correctBases[j][i]);
					*/
					
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
			double[] averageErrorPerDataset = new double[numberDatasets];
			perCycleErrorRatePerDataset = new double[numberDatasets][];
			for (int i=0; i< numberDatasets; i++) {
				perCycleErrorRatePerDataset[i] = Num.arrayListOfDoubleToArray(fractionError[i]);
				averageErrorPerDataset[i] = Num.mean(perCycleErrorRatePerDataset[i]);
			}
			System.out.println("\nMean\t"+Num.doubleArrayToString(averageErrorPerDataset, "\t"));
			
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

	@SuppressWarnings("unchecked")
	private void saveJson() {
		try {
			//calc stats just for first two datasets, read one and read two
			double[] perCycleErrorRatesReadOne = perCycleErrorRatePerDataset[0];
			double[] perCycleErrorRatesReadTwo = perCycleErrorRatePerDataset[1];
			int stop = perCycleErrorRatesReadOne.length-1;
			
			//zero first and last to skip first and last cycles, these are likely artifacts of the alignment process
			perCycleErrorRatesReadOne[0] = 0;
			perCycleErrorRatesReadOne[stop] = 0;
			perCycleErrorRatesReadTwo[0] = 0;
			perCycleErrorRatesReadTwo[stop] = 0;
			
			double totalOne = 0;
			double totalTwo = 0;
			double onePerOne = 0;
			double onePerTwo = 0;
			for (int i=1; i< stop; i++){
				totalOne+= perCycleErrorRatesReadOne[i];
				totalTwo+= perCycleErrorRatesReadTwo[i];
				if (perCycleErrorRatesReadOne[i]>=0.01) onePerOne++;
				if (perCycleErrorRatesReadTwo[i]>=0.01) onePerTwo++;
			}
			double numCyclesMinFirstLast = stop-1;
			double meanPerCycleErrorReadOne = totalOne/numCyclesMinFirstLast;
			double fractionCyclesWithOnePercentErrorReadOne = onePerOne/numCyclesMinFirstLast;
			double meanPerCycleErrorReadTwo = totalTwo/numCyclesMinFirstLast;
			double fractionCyclesWithOnePercentErrorReadTwo = onePerTwo/numCyclesMinFirstLast;

			int numberAlignmentsReadOne = (int)numberAlignments[0];
			double fractionAlignmentsWithINDELsReadOne = (numberAlignmentsWithInsertions[0]+numberAlignmentsWithDeletions[0])/ numberAlignments[0];
			int numberAlignmentsReadTwo = (int)numberAlignments[1];
			double fractionAlignmentsWithINDELsReadTwo = (numberAlignmentsWithInsertions[1]+numberAlignmentsWithDeletions[1])/ numberAlignments[1];
			
			//save json, DO NOT change the key names without updated downstream apps that read this file!
			Gzipper gz = new Gzipper(jsonOutputFile);
			gz.println("{");
			gz.printJson("numberAlignmentsReadOne", numberAlignmentsReadOne, true);
			gz.printJson("numberAlignmentsReadTwo", numberAlignmentsReadTwo, true);
			gz.printJson("fractionAlignmentsWithINDELsReadOne", fractionAlignmentsWithINDELsReadOne, true);
			gz.printJson("fractionAlignmentsWithINDELsReadTwo", fractionAlignmentsWithINDELsReadTwo, true);
			gz.printJson("meanPerCycleErrorReadOne", meanPerCycleErrorReadOne, true);
			gz.printJson("meanPerCycleErrorReadTwo", meanPerCycleErrorReadTwo, true);
			gz.printJson("fractionCyclesWithOnePercentErrorReadOne", fractionCyclesWithOnePercentErrorReadOne, true);
			gz.printJson("fractionCyclesWithOnePercentErrorReadTwo", fractionCyclesWithOnePercentErrorReadTwo, true);
			gz.printJson("perCycleErrorRatesReadOne", perCycleErrorRatesReadOne, true);
			gz.printJson("perCycleErrorRatesReadTwo", perCycleErrorRatesReadTwo, false);
			gz.println("}");
			gz.close();
			
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem writing json file! "+jsonOutputFile);
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
		//Pattern pat = Pattern.compile("(\\d+)M");
		while (it.hasNext()) {
			SAMRecord sam = it.next();
		
			//is it aligned?
			if (sam.getReadUnmappedFlag()) continue;

			//does it pass the vendor qc?
			if (sam.getReadFailsVendorQualityCheckFlag()) continue;

			//Matcher mat = pat.matcher(sam.getCigarString());
			int cycles = sam.getReadString().length();
			if (sam.getReadNegativeStrandFlag()) {
				if (maxSecondPair < cycles) {
					maxSecondPair = cycles;
					counter1++;
				}
			} else if (maxFirstPair < cycles) {
				maxFirstPair = cycles;
				counter2++;
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
		if (start < 0) {
			return -1;
		}
		int end = sam.getAlignmentEnd();
		if (end >= chromLength) {
			return -2;
		}

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
					//an N in reference?
					if (chromSeq[start] != 'N'){
						if (chromSeq[start] == seq[seqIndex]) {
							if (cycleSubtractor !=0) correctBases[localIndex][cycleSubtractor - seqIndex]++;
							else correctBases[localIndex][seqIndex]++;
						}
						//else if (seq[seqIndex] != 'N') {, count N's as an error.
						else {
							if (cycleSubtractor !=0) incorrectBases[localIndex][cycleSubtractor - seqIndex]++;
							else incorrectBases[localIndex][seqIndex]++;
							numberMismatches++;
						}
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
					case 's': mergeStrands = false; break;
					case 'o': logFile = new File(args[++i]); break;
					case '1': firstReadMaximumError = Double.parseDouble(args[++i]); break;
					case '2': secondReadMaximumError = Double.parseDouble(args[++i]); break;
					case 'c': maximumFractionFailingCycles = Double.parseDouble(args[++i]); break;
					case 'j': jsonOutputFile = new File(args[++i]); mergeStrands = false; break;
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
		
		if (mergeStrands && jsonOutputFile != null) Misc.printErrAndExit("\nError: a merged strand analysis isn't permitted when exporting json statistics. Include -s or drop -j.\n");

	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Calculate Per Cycle Error Rate : July 2015                **\n" +
				"**************************************************************************************\n" +
				"Calculates per cycle snv error rate provided a sorted indexed bam file and a fasta\n" +
				"sequence file. Only checks CIGAR M bases not masked or INDEL bases.\n\n" +

				"Required Options:\n"+
				"-b Full path to a coordinate sorted bam file (xxx.bam) with its associated (xxx.bai)\n" +
				"      index or directory containing such. Multiple files are processed independently.\n" +
				"-f Full path to the single fasta file you wish to use in calculating the error rate.\n" +
				
				"\nDefault Options:\n"+
				"-s Perform separate first read second read analysis, defaults to merging.\n"+
				"-c Maximum fraction failing cycles, defaults to 0.1\n"+
				"-1 Maximum first read or merged read error rate, defaults to 0.01\n"+
				"-2 Maximum second read error rate, defaults to 0.0175\n"+
				"-o Write coverage statistics to this log file instead of stdout.\n" +
				"-j Write summary stats in json format to this file. Only stats for the first bam file\n"+ 
				"      are saved. Only separate strand analysis permitted.\n"+

				"\nExample: java -Xmx1500M -jar pathTo/USeq/Apps/CalculatePerCycleErrorRate\n"+
				"     -b /Data/Bam/ -f /Fastas/chrPhiX_Illumina.fasta.gz \n\n"+

		"**************************************************************************************\n");

	}
}
