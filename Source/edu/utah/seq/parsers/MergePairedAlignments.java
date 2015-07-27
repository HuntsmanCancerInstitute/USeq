package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;
import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import edu.utah.seq.data.ChromDataSave;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * @author david.nix@hci.utah.edu 
 * Similar to MergePairedSamAlignments but cut down in options for speed. 
 * Now works on coordinate and query name sorted files. Multi threaded.
 **/
public class MergePairedAlignments {

	//user defined fields
	private File bamFile;
	private File saveDirectory;
	private float maximumAlignmentScore = 300;
	private float minimumMappingQualityScore = 0;
	private boolean secondPairReverseStrand = false;
	private boolean skipMergingPairs = false;
	private boolean removeDuplicates = false;
	private int maximumProperPairDistanceForMerging = 5000;
	private boolean printPairedStatsAndHistogram = true;
	private boolean crossCheckMateCoordinates = true;
	private int numberConcurrentThreads = 0;
	private boolean saveSams = false;
	private boolean onlyMergeOverlappingAlignments = true;
	private File jsonOutputFile = null;

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
	private int numberPrintedAlignments = 0;
	private double numberOverlappingBases = 0;
	private double numberNonOverlappingBases = 0;
	private double numberMergedPairs = 0;
	private double numberFailedMergedPairs = 0;
	private int numberNonPairedAlignments = 0;
	private int numberNonProperPairedAlignments = 0;
	private int numberUnmappedMatePairedAlignments = 0;
	private int numberAlignmentsMissingPair = 0;
	private int numberPairsFailingChrDistStrand = 0;
	private int numberPairsFailingMateCrossCoordinateCheck = 0;
	private int numberRepeatAlignmentsLackingMate = 0;
	
	//base scores
	private boolean calculateBaseQualities = true;
	private double numberPassingBases = 0;
	private double numberPassingQ20Bases = 0;
	private double numberPassingQ30Bases = 0;

	//internal fields
	private int numberAlignmentsToLoad = 100000;
	public static Pattern CIGAR_SUB = Pattern.compile("(\\d+)([MDIN])");
	public static Pattern CIGAR_BAD = Pattern.compile(".*[^\\dMDIN].*");
	private Histogram insertSize = new Histogram(0,2001,400);
	private String programArguments;
	private ArrayList<String> chromosomes = null;
	private String samHeader;
	private ArrayList<PairedAlignmentChrParser> parsers;

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

		//fetch file header and chromosomes
		System.out.println("\nFetching header and chroms...");
		fetchChromosomesAndHeader();

		//for each chromosome, create a thread object but don't start
		parsers = new ArrayList<PairedAlignmentChrParser>();
		//add parser for unmapped reads, note this doesn't fetch all of them, need to scan rest of chromosomes
		parsers.add(new PairedAlignmentChrParser("unmapped", this));
		//add chrom parsers
		for (int i=0; i< chromosomes.size(); i++) parsers.add(new PairedAlignmentChrParser(chromosomes.get(i), this));
		//for (int i=23; i< 26; i++) parsers.add(new PairedAlignmentChrParser(chromosomes.get(i), this));

		System.out.print("\nRunning threads...");
		runThreads();
		System.out.println();
		
		//do it before removing parsers
		collectStats();
		
		//save chrom data info
		if (saveSams == false){
			//strip out parsers with no printed alignments
			for (int i=0; i< parsers.size(); i++){
				if (parsers.get(i).getNumberPrintedAlignments() == 0){
					parsers.remove(i);
					i--;
				}
			}
			//save chrom data
			ChromDataSave[] cd = new ChromDataSave[parsers.size()];
			for (int i=0; i< parsers.size(); i++) cd[i] = new ChromDataSave(parsers.get(i).getChromData());
			File al = new File (saveDirectory, "chromData.obj");
			IO.saveObject(al, cd);
			//save count for scaling if they need it
			File countsFile = new File (saveDirectory, "count.obj");
			IO.saveObject(countsFile, new Long(numberPrintedAlignments));
		}
		printStats();
		if (jsonOutputFile != null) saveJson();
	}

	/**Runs each parser in its own thread to max threads defined by user.*/
	private void runThreads() {
		while (true){ 
			try {
				//how many running
				int numRunning = 0;
				ArrayList<PairedAlignmentChrParser> toStart = new ArrayList<PairedAlignmentChrParser>();
				for (PairedAlignmentChrParser p : parsers){
					if (p.isStarted() && p.isComplete()==false) numRunning++;
					if (p.isStarted() == false) toStart.add(p);
				}
				int numToStart = numberConcurrentThreads - numRunning;
				if (numToStart > 0){
					int stop = numToStart;
					if (stop > toStart.size()) stop = toStart.size();
					for (int i=0; i< stop; i++){
						PairedAlignmentChrParser p = toStart.get(i);
						p.start();
						//System.out.println("Starting "+p.getChromosome());
					}
				}
				boolean allComplete = true;
				for (PairedAlignmentChrParser p : parsers){
						if (p.isComplete() == false) {
						allComplete = false;
						break;
					}
				}
				if (allComplete) break;
				//sleep
				Thread.sleep(60000);
			} catch (InterruptedException ie) {
				ie.printStackTrace();
				Misc.printErrAndExit("\nProblem running threads!");
			}
		}
	}

	/**Sums the stats from the individual parsers into this.*/
	private void collectStats() {
		try {			
			//from each parser collect
			for (int i=0; i< parsers.size(); i++){
				PairedAlignmentChrParser parser = parsers.get(i);
				numberAlignments += parser.getNumberAlignments();
				numberUnmapped += parser.getNumberUnmapped();
				numberFailingVendorQC += parser.getNumberFailingVendorQC() ;
				numberFailingAlignmentScore += parser.getNumberFailingAlignmentScore() ;
				numberFailingMappingQualityScore += parser.getNumberFailingMappingQualityScore() ;
				numberAdapter += parser.getNumberAdapter() ;
				numberPhiX += parser.getNumberPhiX() ;
				numberDuplicates += parser.getNumberDuplicates() ;
				numberPassingAlignments += parser.getNumberPassingAlignments() ;
				numberPrintedAlignments += parser.getNumberPrintedAlignments() ;
				numberOverlappingBases += parser.getNumberOverlappingBases() ;
				numberNonOverlappingBases += parser.getNumberNonOverlappingBases() ;
				numberMergedPairs += parser.getNumberMergedPairs() ;
				numberFailedMergedPairs += parser.getNumberFailedMergedPairs() ;
				numberNonPairedAlignments += parser.getNumberNonPairedAlignments() ;
				numberNonProperPairedAlignments += parser.getNumberNonProperPairedAlignments() ;
				numberUnmappedMatePairedAlignments += parser.getNumberUnmappedMatePairedAlignments() ;
				numberAlignmentsMissingPair += parser.getNumberAlignmentsMissingPair() ;
				numberPairsFailingChrDistStrand += parser.getNumberPairsFailingChrDistStrand() ;
				numberPairsFailingMateCrossCoordinateCheck += parser.getNumberPairsFailingMateCrossCoordinateCheck() ;
				numberRepeatAlignmentsLackingMate += parser.getNumberRepeatAlignmentsLackingMate() ;
				insertSize.addCounts(parser.getInsertSize());
				numberPassingBases += parser.getNumberPassingBases();
				numberPassingQ20Bases += parser.getNumberPassingQ20Bases();
				numberPassingQ30Bases += parser.getNumberPassingQ30Bases();
			};
		} catch (Exception e) {
			System.err.println("Problem collecting stats!");
			e.printStackTrace();
		}
	}

	public void printStats(){
		//stats, don't modify these descriptions without updating the CollectBamStats app!
		System.out.println("\nStats (some flags aren't set so be suspicious of zero read catagories):");
		System.out.println("Total # alignments\t"+numberAlignments);
		System.out.println("Unmapped reads\t"+format(numberUnmapped, numberAlignments));
		System.out.println("Alignments failing vendor or platform QC\t"+format(numberFailingVendorQC,numberAlignments));
		System.out.println("Alignments failing alignment score\t"+format(numberFailingAlignmentScore,numberAlignments));
		System.out.println("Alignments failing mapping quality score\t"+format(numberFailingMappingQualityScore,numberAlignments));
		if (numberAdapter!=0) System.out.println("Adapter alignments\t"+format(numberAdapter,numberAlignments));
		System.out.println("PhiX alignments\t"+format(numberPhiX,numberAlignments));
		System.out.println("Duplicate alignments\t"+format(numberDuplicates,numberAlignments));
		System.out.println("Alignments passing all filters\t"+format(numberPassingAlignments,numberAlignments));
		System.out.println();
		if (calculateBaseQualities){
			System.out.println("Base stats for passing alignments:");
			System.out.println("Total # aligned non-overlapping bps\t"+(long)numberPassingBases);
			System.out.println("Q20 bps\t"+format(numberPassingQ20Bases, numberPassingBases));
			System.out.println("Q30 bps\t"+format(numberPassingQ30Bases, numberPassingBases));
			System.out.println();
		}
		System.out.println("Paired alignment stats:");  
		if (printPairedStatsAndHistogram){
			System.out.println("Non paired alignments\t"+numberNonPairedAlignments);
			System.out.println("Non proper paired alignments\t"+numberNonProperPairedAlignments);
			System.out.println("Non mapped mate paired alignments\t"+numberUnmappedMatePairedAlignments);
			System.out.println("Alignments missing mate paired alignment\t"+numberAlignmentsMissingPair);
			System.out.println("Paired alignments failing chromosome, distance, or strand check\t"+numberPairsFailingChrDistStrand);
			System.out.println("Paired alignments failing mate pair cross coordinate check\t"+numberPairsFailingMateCrossCoordinateCheck);
			System.out.println("Repeat alignments lacking a mate\t"+numberRepeatAlignmentsLackingMate);
			System.out.println("Proper paired alignments that could not be unambiguously merged\t"+(int)numberFailedMergedPairs);
		}
		System.out.println("Proper paired alignments that were merged\t"+format(numberMergedPairs, (numberMergedPairs+numberFailedMergedPairs)));
		double totalBases = numberNonOverlappingBases + numberOverlappingBases;
		double fractionOverlap = numberOverlappingBases/totalBases;
		String fractionString = Num.formatNumber(fractionOverlap, 4);
		System.out.println("Fraction overlapping bases in paired alignments\t"+fractionString);		
		//histogram
		if (printPairedStatsAndHistogram && insertSize.getTotalBinCounts() !=0){
			System.out.println("\nMapped genomic insert length distribution for merged paired alignments:\n");
			insertSize.setSkipZeroBins(true);
			insertSize.setTrimLabelsToSingleInteger(true);
			insertSize.printScaledHistogram();
			System.out.println();
		}
		double mean = insertSize.getStandardDeviation().getMean();
		System.out.println("Mean insert size\t"+Num.formatNumber(mean, 1));
		System.out.println("Alignments written to file\t"+format(numberPrintedAlignments,numberAlignments));
	}
	
	@SuppressWarnings("unchecked")
	private void saveJson() {
		try {
			//save json obj, DO NOT change the key names without updated downstream apps that read this file!
			Gzipper gz = new Gzipper(jsonOutputFile);
			gz.println("{");
			gz.printJson("fractionOverlappingBpsInPairedReads", numberOverlappingBases/(numberNonOverlappingBases + numberOverlappingBases), true);
			gz.printJson("meanInsertSize", insertSize.getStandardDeviation().getMean(), true);
			gz.printJson("numberPassingBps", (long)numberPassingBases, true);
			gz.printJson("fractionPassingQ20bps", numberPassingQ20Bases/numberPassingBases, true);
			gz.printJson("fractionPassingQ30bps", numberPassingQ30Bases/numberPassingBases, false);
			gz.println("}");
			gz.close();
			
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem writing json file! "+jsonOutputFile);
		}
	}

	public static String format (double stat, double total){
		double frac = stat/total;
		return (long)stat +"\t"+Num.formatNumber(frac, 4);
	}

	public void fetchChromosomesAndHeader(){
		try {
			//get reader 
			SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			SamReader samReader = factory.open(bamFile);

			//fetch header
			samHeader = samReader.getFileHeader().getTextHeader().trim();
			samHeader = samHeader+"\n@PG\tID:MergePairedAlignments\tCL: "+programArguments;

			//get chromosomes
			List<SAMSequenceRecord> chrs = samReader.getFileHeader().getSequenceDictionary().getSequences();
			chromosomes = new ArrayList<String>();
			for (SAMSequenceRecord r : chrs) chromosomes.add(r.getSequenceName());
			//chromosomes.add("22");
			//chromosomes.add("21");

			samReader.close();
		} catch (Exception e) {
			System.err.println("\nError loading header and chromsomes from "+bamFile.getName());
			e.printStackTrace();
			System.exit(0);
		} 
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
					case 'b': bamFile = new File(args[++i]); break;
					case 'd': saveDirectory = new File(args[++i]); break;
					case 's': saveSams = true; break;
					case 'm': crossCheckMateCoordinates = false; break;
					case 'r': secondPairReverseStrand = true; break;
					case 'k': skipMergingPairs = true; break;
					case 'u': removeDuplicates = true; break;
					case 'p': printPairedStatsAndHistogram = false; break;
					case 'i': maximumProperPairDistanceForMerging = Integer.parseInt(args[++i]); break;
					case 't': numberConcurrentThreads = Integer.parseInt(args[++i]); break;
					case 'o': onlyMergeOverlappingAlignments = false; break;
					case 'j': jsonOutputFile = new File(args[++i]); break;
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

		if (bamFile == null || bamFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find your xxx.bam/.sam(.gz OK) results file!\n");
		OverdispersedRegionScanSeqs.lookForBaiIndexes(new File[]{bamFile}, true);

		//check save dir
		if (saveDirectory == null) Misc.printErrAndExit("\nError: cannot find your save output directory!\n");
		saveDirectory.mkdirs();
		if (saveDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: your save output directory does not appear to be a directory?\n");

		//number of threads to use?
		if (numberConcurrentThreads == 0) numberConcurrentThreads = Runtime.getRuntime().availableProcessors();
		
		//print info
		System.out.println("Settings:");
		System.out.println(maximumAlignmentScore+ "\tMaximum alignment score (AS).");
		System.out.println(minimumMappingQualityScore+ "\tMinimum mapping quality score (MQ).");
		System.out.println(removeDuplicates +"\tRemove all alignments marked duplicate.");
		System.out.println(secondPairReverseStrand +"\tSecond read pair's strand has been reversed.");
		System.out.println(crossCheckMateCoordinates +"\tCross check read mate coordinates.");
		System.out.println(maximumProperPairDistanceForMerging +"\tMaximum bp distance for merging paired alignments.");
		if (skipMergingPairs) System.out.println(skipMergingPairs +"\tSkip merging paired alignments, useful for testing effect of merging on downstream analysis.");
		System.out.println(numberConcurrentThreads +"\tNumber concurrent threads.");
		System.out.println(saveSams +"\tSave sam alignments.");
		System.out.println(onlyMergeOverlappingAlignments +"\tOnly merge overlapping alignments.");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Merge Paired Alignments: July 2015                    **\n" +
				"**************************************************************************************\n" +
				"Merges proper paired alignments that pass a variety of checks and thresholds. Only\n" +
				"unambiguous pairs will be merged. Increases base calling accuracy in overlap and helps\n" +
				"avoid non-independent variant observations and other double counting issues. Identical\n" +
				"overlapping bases are assigned the higher quality scores. Disagreements are resolved\n" +
				"toward the higher quality base. If too close in quality, then the quality is set to 0.\n" +

				"\nOptions:\n"+
				"-b Path to a coordinate sorted xxx.bam file containing paired alignments.\n" +
				"-d Path to a directory for saving the results.\n"+

				"\nDefault Options:\n"+
				"-s Save merged xxx.sam.gz alignments instead of binary ChromData. Either can be used\n"+
				"      in Sam2USeq for read coverage analysis, the ChromData is much faster.\n"+
				"-u Remove all alignments marked as duplicates, defaults to keeping.\n"+
				"-a Maximum alignment score (AS:i: tag). Defaults to 300, smaller numbers are more\n" +
				"      stringent for novoalign where each mismatch is ~30pts.\n"+
				"-q Minimum mapping quality score, defaults to 0, larger numbers are more stringent.\n" + 
				"-r The second paired alignment's strand has been reversed. Defaults to not reversed.\n" +
				"-i Maximum acceptible base pair distance for merging, defaults to 5000.\n"+
				"-m Don't cross check read mate coordinates, needed for merging repeat matches. Defaults\n" +
				"      to checking.\n"+
				"-o Merge all proper paired alignments. Defaults to only merging those that overlap.\n"+
				"-p Don't print detailed paired alignment statistics and insert size histogram.\n"+
				"-t Number concurrent threads to run, defaults to the max available to the jvm.\n"+
				"-j Write summary stats in json format to this file.\n"+

				"\nExample: java -Xmx20G -jar pathToUSeq/Apps/MergePairedBamAlignments -f /Bams/ms.bam\n" +
				"     -p -s /Bams/MergedPairs/ms.mergedPairs.sam.gz -d 10000 \n\n" +

				"**************************************************************************************\n");

	}

	public File getBamFile() {
		return bamFile;
	}

	public File getSaveDirectory() {
		return saveDirectory;
	}

	public int getNumberAlignmentsToLoad() {
		return numberAlignmentsToLoad;
	}

	public String getSamHeader() {
		return samHeader;
	}

	public boolean isCrossCheckMateCoordinates() {
		return crossCheckMateCoordinates;
	}

	public boolean isSkipMergingPairs() {
		return skipMergingPairs;
	}

	public int getMaximumProperPairDistanceForMerging() {
		return maximumProperPairDistanceForMerging;
	}

	public boolean isSecondPairReverseStrand() {
		return secondPairReverseStrand;
	}

	public float getMaximumAlignmentScore() {
		return maximumAlignmentScore;
	}

	public float getMinimumMappingQualityScore() {
		return minimumMappingQualityScore;
	}

	public boolean isSaveSams() {
		return saveSams;
	}

	public boolean isRemoveDuplicates() {
		return removeDuplicates;
	}

	public boolean isOnlyMergeOverlappingAlignments() {
		return onlyMergeOverlappingAlignments;
	}

	public boolean isCalculateBaseQualities() {
		return calculateBaseQualities;
	}
}
