package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import util.bio.annotation.Bed;
import util.gen.*;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.useq.data.RegionScoreText;


/**
 * @author Nix
 * App for extracting sam alignments that are on target and pass basic QC flags and MQ and AS scores.
 * Writes summary stats to file in json format as well as detailed QC to stdout
 * */
public class SamAlignmentExtractor {

	//user defined fields
	private File bamFile;
	private File bedFile;
	private File saveDirectory;
	private int minimumMappingQuality = -1;
	private int minimumFamilySize = 0;
	private double alignmentScoreThreshold = -1.0;
	private boolean biggerASIsBetter = true;
	private boolean divideAlignmentScoreByCigarM = false;
	private boolean writeOffTargetToPass = false;
	private boolean removeSecSupNotPrim = true;
	private boolean skipWritingFail = false;

	//internal fields
	private HashMap<String,RegionScoreText[]> chromRegions;
	private static Pattern CIGAR_SUB = Pattern.compile("(\\d+)([MSDHN])");
	private static Pattern FAMILY_SIZE = Pattern.compile(".+:FS:(\\d+).*");
	private SamReader bamReader;
	private SAMFileWriter passingBamWriter;
	private SAMFileWriter failingBamWriter;
	private String workingChromosome;
	private boolean[] workingCoveredBases;
	private RegionScoreText[] workingRegions;
	private int numRawAlignments = 0;
	private int numPassingBasicOnTargetYetFailingMQ = 0;
	private int numPassingBasicOnTargetYetFailingAS = 0;
	private int numPassingBasicAndOnTarget = 0;
	private int numFailingBasic = 0;
	private int numFailingFS = 0;
	private int numPassingBasicOnTargetAndScores = 0;
	private int numPassingBasicOnTargetAndScoresYetMarkedAsADuplicate = 0;
	private File jsonOutputFile;
	private Histogram histogram = new Histogram(1, 250, 249);
	
	
	//constructors
	/**Stand alone.*/
	public SamAlignmentExtractor(String[] args){
		long startTime = System.currentTimeMillis();
		
		//set fields
		processArgs(args);
		
		printSettings();
		
		//launch
		run();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}

	public void printSettings(){
		System.out.println("Settings:");
		System.out.println("Bam\t"+bamFile);
		System.out.println("Bed\t"+bedFile);
		System.out.println("Save dir\t"+saveDirectory);
		System.out.println("Min Map Qual\t"+minimumMappingQuality);
		System.out.println("Align Score\t"+alignmentScoreThreshold);
		System.out.println("Bigger AS is better\t"+biggerASIsBetter);
		System.out.println("Divide AS by cigar M\t"+divideAlignmentScoreByCigarM);
		System.out.println("Min Family Size\t"+minimumFamilySize);
		System.out.println("Pass off target align\t"+writeOffTargetToPass);
		System.out.println("Fail Sec, Sup, Non Prim align\t"+removeSecSupNotPrim);
		System.out.println("Skip writing failing align\t"+skipWritingFail);
	}
	
	public void run(){
		try {
			
			//fetch chroms to scan from bam header
			List<SAMSequenceRecord> chrList = bamReader.getFileHeader().getSequenceDictionary().getSequences();

			//for each
			//TODO: thread this! would need to write out individual them merge... might be slow.. the bottleneck is the io not the compute so parallel won't gain much
			System.out.print("\nWalking bam alignments");
			Iterator<SAMSequenceRecord> it = chrList.iterator();
			while (it.hasNext()){
				SAMSequenceRecord ssr = it.next();
				workingChromosome = ssr.getSequenceName();
				//any regions?
				workingRegions = chromRegions.get(workingChromosome);
				if (workingRegions == null) printAllToFail();
				else {
					makeMask();
					walkChromAlignments();
				}
				System.out.print(".");
			}
			
			//save all unmapped to fail
			printUnmapped();
			
			//close readers
			closeIO();
			
			//print stats
			System.out.println("\n\nSAE statistics:");
			printStatLine(numFailingBasic, numRawAlignments, "Unmapped, failing vendor QC, possibly secondary/ non prim");
			printStatLine(numPassingBasicAndOnTarget, (numRawAlignments-numFailingBasic), "On target regions");
			printStatLine(numPassingBasicOnTargetYetFailingMQ, numPassingBasicAndOnTarget, "Failing MQ score ("+minimumMappingQuality+")");
			String dir = "<";
			if (biggerASIsBetter == false) dir = ">";
			printStatLine(numPassingBasicOnTargetYetFailingAS, numPassingBasicAndOnTarget, "Failing AS score ("+dir+alignmentScoreThreshold+")");
			if (minimumFamilySize !=0) printStatLine(numFailingFS, (numRawAlignments-numFailingBasic), "Failing minimum molecular barcode family size ("+minimumFamilySize+")");
			printStatLine(numPassingBasicOnTargetAndScoresYetMarkedAsADuplicate, numPassingBasicOnTargetAndScores, "Duplicates (not filtered)");
			printStatLine(numPassingBasicOnTargetAndScores, numRawAlignments, "Passing all filters");
			
			System.out.println("\nFamily Size Histogram for alignments passing filters:");
			if (minimumFamilySize !=0) histogram.printScaledHistogram();
			
			if (jsonOutputFile !=null) saveJson();

		} catch (Exception e) {
			e.printStackTrace();
		} 
	}
	
	@SuppressWarnings("unchecked")
	private void saveJson() {
		try {
			//calc stats
			Integer numberUnfilteredAlignments = new Integer(numRawAlignments);
			Double fractionAlignmentsOnTarget = new Double((double)numPassingBasicAndOnTarget/(double)(numRawAlignments-numFailingBasic));
			Double fractionOnTargetAndPassQCScoreFilters = new Double((double)numPassingBasicOnTargetAndScores/(double)numRawAlignments);
			Double estimatedFractionDuplicateAlignments = new Double( (double)numPassingBasicOnTargetAndScoresYetMarkedAsADuplicate/ (double) numPassingBasicOnTargetAndScores  );
			Double fractionAlignmentsPassQCScoreFilters = new Double ( ((double)(numRawAlignments - numFailingBasic))/ (double)numRawAlignments);
			
			//output simple json, DO NOT change the key names without updated downstream apps that read this file!
			Gzipper gz = new Gzipper(jsonOutputFile);
			gz.println("{");
			gz.printJson("numberUnfilteredAlignments", numberUnfilteredAlignments, true);
			gz.printJson("fractionAlignmentsPassQCScoreFilters", fractionAlignmentsPassQCScoreFilters, true);
			gz.printJson("fractionAlignmentsOnTarget", fractionAlignmentsOnTarget, true);
			gz.printJson("fractionOnTargetAndPassQCScoreFilters", fractionOnTargetAndPassQCScoreFilters, true);
			gz.printJson("estimatedFractionDuplicateAlignments", estimatedFractionDuplicateAlignments, true);
			gz.printJson("targetRegionsFileName", bedFile.getName(), true);
			gz.printJson("mappingQualityThreshold", minimumMappingQuality, true);
			gz.printJson("alignmentScoreThreshold", alignmentScoreThreshold, true);
			gz.printJson("divideAlignmentScoreByCigarM", divideAlignmentScoreByCigarM, false);
			gz.println("}");
			gz.close();
			
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem writing json file! "+jsonOutputFile);
		}
	}


	public static void printStatLine(double numerator, double denomenator, String name){
		double fraction = numerator/ denomenator;
		if (denomenator == 0) fraction = 0;
		System.out.println(Num.formatNumber(fraction, 3)+"\t("+(long)numerator+"/"+(long)denomenator+")\t"+name);
	}
	
	private void printAllToFail() {
		if (skipWritingFail) return;
		SAMRecordIterator it = bamReader.queryOverlapping(workingChromosome, 0, 0);
		while (it.hasNext()) {
			SAMRecord sam = it.next();
			failingBamWriter.addAlignment(sam);
			numRawAlignments++;
			//fail basic?
			if (passBasic(sam) == false) numFailingBasic++;
			
			
		}
		it.close();
	}
	
	private void printUnmapped() {
		if (skipWritingFail) return;
		SAMRecordIterator it = bamReader.queryUnmapped();
		while (it.hasNext()) {
			failingBamWriter.addAlignment(it.next());
			numRawAlignments++;
			numFailingBasic++;
		}
		it.close();
	}


	public void makeMask (){
		//find max base
		int maxBase = findMaxBase();
		//make boolean array to hold whether it's flagged, initially they are all false
		workingCoveredBases = new boolean[maxBase+10000];
		//for each RegionScoreText[] scan and throw booleans to true
		for (int i=0; i<workingRegions.length; i++){
			int start = workingRegions[i].getStart();
			int end = workingRegions[i].getStop();
			for (int j=start; j< end; j++) workingCoveredBases[j] = true;
		}
	}

	public int findMaxBase (){
		int max = workingRegions[workingRegions.length-1].getStop();
		for (int i=0; i< workingRegions.length; i++){
			if (workingRegions[i].getStop()> max) max = workingRegions[i].getStop();
		}
		return max;
	}
	
	
	/**Assumes interbase coordinates, looking to see if just one M base hits a region.*/
	public boolean intersect(SAMRecord sam) throws Exception{
		//past the boolean[]?
		if (sam.getAlignmentEnd() >= workingCoveredBases.length) return false;
		int start = sam.getUnclippedStart()-1;
		//for each cigar block
		Matcher mat = CIGAR_SUB.matcher(sam.getCigarString());
		while (mat.find()){
			String call = mat.group(2);
			int numberBases = Integer.parseInt(mat.group(1));
			//a match
			if (call.equals("M")) {
				int end = start+ numberBases;
				for (int i=start; i< end; i++){
					if (workingCoveredBases[i]) return true;
				}
			}
			//advance for all but insertions which should be skipped via the failure to match
			start += numberBases;
		}
		return false;
	}
	
	private boolean passBasic(SAMRecord sam){
		if (sam.getReadUnmappedFlag() || sam.getReadFailsVendorQualityCheckFlag()) return false;
		else if (removeSecSupNotPrim ){
			if (sam.isSecondaryOrSupplementary() || sam.getNotPrimaryAlignmentFlag()) return false;
		}
		return true;
	}

	
	private void walkChromAlignments() {
		try {
			SAMRecordIterator it = bamReader.queryOverlapping(workingChromosome, 0, 0);
			while (it.hasNext()) {
				SAMRecord sam = it.next();
				numRawAlignments++;

				//fail basic flags?
				if (passBasic(sam) == false){
					if (skipWritingFail == false) failingBamWriter.addAlignment(sam);
					numFailingBasic++;
					continue;
				}
				
				//does it intersect? 
				boolean onTarget = false;
				if (workingRegions != null) onTarget = intersect(sam);
				if (onTarget == false) {
					if (writeOffTargetToPass == false){
						if (skipWritingFail == false) failingBamWriter.addAlignment(sam);
						continue;
					}
				}
				else numPassingBasicAndOnTarget++;

				//check both scores
				boolean passScores = true;
				//pass mapping quality?
				if (minimumMappingQuality != -1) {		
					if (sam.getMappingQuality() < minimumMappingQuality) {
						if (onTarget) numPassingBasicOnTargetYetFailingMQ++;
						passScores = false;
					}
				}
				//pass alignment score? with novoalignments (~30 pt penalty per mismatch) smaller AS is better 
				//with bwa bigger scores are better (readLength - #SNV*5 + #INDEL*7)
				if (alignmentScoreThreshold != -1.0){
					Object obj = sam.getAttribute("AS");
					if (obj != null){
						Integer as = (Integer)obj;
						double asScore = as.doubleValue();
						
						//scale the score?
						if (divideAlignmentScoreByCigarM){
							double numM = SamAlignment.countLengthOfM(sam.getCigarString());
							asScore = asScore/numM;
						}
						if (biggerASIsBetter){ 
							if (asScore < alignmentScoreThreshold) {
								if (onTarget) numPassingBasicOnTargetYetFailingAS++;
								passScores = false;
							}
						}
						else {
							if (asScore > alignmentScoreThreshold) {
								if (onTarget) numPassingBasicOnTargetYetFailingAS++;
								passScores = false;
							}
						}
					}
				}
				//family size
				if (minimumFamilySize !=0 && passScores){
					int size = 1;
					Matcher fs = FAMILY_SIZE.matcher(sam.getReadName());
					if (fs.matches()) size = Integer.parseInt(fs.group(1));
					histogram.count(size);
					if (size < minimumFamilySize) {
						passScores = false;
						numFailingFS++;
					}
				}
				
				//pass? or fail? scores
				if (passScores) {
					passingBamWriter.addAlignment(sam);
					if (onTarget) {
						numPassingBasicOnTargetAndScores++;
						if (sam.getDuplicateReadFlag()) numPassingBasicOnTargetAndScoresYetMarkedAsADuplicate++;
					}
				}
				else if (skipWritingFail == false) failingBamWriter.addAlignment(sam);
			}
			it.close();
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nError: problem walking chromosome "+workingChromosome);
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SamAlignmentExtractor(args);
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
					case 'b': bamFile = new File(args[++i]); break;
					case 'r': bedFile = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'q': minimumMappingQuality = Integer.parseInt(args[++i]); break;
					case 'm': minimumFamilySize = Integer.parseInt(args[++i]); break;
					case 'a': alignmentScoreThreshold = Double.parseDouble(args[++i]); break;
					case 'n': biggerASIsBetter = false; break;
					case 'd': divideAlignmentScoreByCigarM = true; break;
					case 'f': writeOffTargetToPass = true; break;
					case 'w': skipWritingFail = true; break;
					case 'x': removeSecSupNotPrim = false; break;
					case 'j': jsonOutputFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group()+", Try -h");
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for bam files
		if (bamFile == null || bamFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your xxx.bam file?\n"+bamFile);
		lookForBaiIndexes(new File[]{bamFile}, false);
		
		//look for bed
		if (bedFile == null || bedFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your region bed file?\n"+bedFile);
		chromRegions = Bed.parseBedFile(bedFile, true, false);
		
		//look for save file, can be null
		if (saveDirectory == null ) Misc.printErrAndExit("\nError: please provide a directory to save the passing and failing bam alignments.\n");
		saveDirectory.mkdirs();

		//make IO
		makeIO();
		
	}	

	private void makeIO() {
		try {
			SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			bamReader = readerFactory.open(bamFile);
			if (bamReader.hasIndex() == false) {
				bamReader.close();
				Misc.printErrAndExit("\nError: cannot find an index for your bam file?\n"+bamFile);
			}
			String rootName = Misc.removeExtension(bamFile.getName());
			File pass = new File(saveDirectory, rootName+"_passSAE.bam");
			File fail = new File(saveDirectory, rootName+"_failSAE.bam");
			SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
			writerFactory.setCreateIndex(true);
			writerFactory.setTempDirectory(saveDirectory);
			//must explicit set into the header that it is sorted for samtools proc alignments
			bamReader.getFileHeader().setSortOrder(SortOrder.coordinate);
			passingBamWriter = writerFactory.makeBAMWriter(bamReader.getFileHeader(), true, pass);
			if (skipWritingFail == false) failingBamWriter = writerFactory.makeBAMWriter(bamReader.getFileHeader(), true, fail);
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: problem closing the bam reader for "+bamFile);
		}
	}
	
	private void closeIO() {
		try {
			bamReader.close();
			passingBamWriter.close();
			if (skipWritingFail == false) failingBamWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError: critical, problem closing the bam IO \n");
		}
	}
	
	/**Looks for xxx.bam.bai and xxx.bai for each bamFile, prints error and exits if missing.*/
	public static void lookForBaiIndexes (File[] bamFiles, boolean onlyResetLastModifiedDate){
		for (File f: bamFiles){
			File index = new File (f+".bai");
			if (index.exists() == false){
				int len = f.toString().length() - 3;
				index = new File(f.toString().substring(0, len) + "bai");
				if (onlyResetLastModifiedDate == false && index.exists() == false) Misc.printErrAndExit("\nError: failed to find a xxx.bai index file for -> "+f);
			}
			//reset date?
			if (index.exists() && index.lastModified() < f.lastModified()) index.setLastModified(f.lastModified()+1);
		}
	}



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Sam Alignment Extractor: March 2017                      **\n" +
				"**************************************************************************************\n" +
				"Splits an alignment file into those that pass or fail thresholds and intersects\n"+
				"regions of interest. Calculates a variety of QC statistics.\n"+

				"\nRequired Options:\n"+
				"-b Bam alignment file with its associated xxx.bai index, sorted by coordinate. \n"+
				"-r A regions bed file (chr, start, stop,...) to intersect, full path, see,\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1 , gz/zip OK.\n"+
				"-s Provide a directory path for saving the filtered alignments\n"+
				
				"\nDefault Options:\n"+
				"-q Miminum mapping quality, defaults to no filtering, recommend 13.\n"+
				"-a Alignment score threshold, defaults to no filtering. \n"+
				"-n Smaller alignment scores are better (novo), defaults to bigger are better (bwa).\n"+
				"-d Divide alignment score by the number of CIGAR M bases.\n"+
				"-m Minimum molecular barcode family size, defaults to 0. Requires :FS:# in read name.\n"+
				"-j Write summary stats in json format to this file.\n"+
				"-f Save off target alignments that meet thresholds to the pass file, defaults to fail.\n"+
				"-x Save secondary, supplemental, and non primary alignments, that pass the thresholds\n"+
				"       defaults to fail.\n"+
				"-w Skip writing failing alignments. Speeds up processing but kills the stats.\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/SamAlignmentExtractor -q 20 -a 0.75 -d -b\n" +
				"      /Data/raw.bwaMem.bam -r /Data/targetsPad25bp.bed.gz -s/Data/SAE/ \n\n"+

		"**************************************************************************************\n");

	}


}
