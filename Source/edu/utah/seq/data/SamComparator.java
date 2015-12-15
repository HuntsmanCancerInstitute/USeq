package edu.utah.seq.data;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.data.sam.PicardSortSam;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

/**Compares two sam/ bam files sorted by coordinate, tosses alignments that aren't exact.*/
public class SamComparator {

	//fields
	private File firstSamFile;
	private File secondSamFile;
	private File saveDirectory;
	private boolean checkSequence = false;
	private Gzipper goodSamFirst;
	private Gzipper badSamFirst;
	private Gzipper goodSamSecond;
	private Gzipper badSamSecond;
	private File[] gzippedFiles = new File[4];
	private int numberMatches = 0;
	private int numberMisMatchesFirst = 0;
	private int numberMisMatchesSecond = 0;
	private int numberChangedAlignments = 0;
	private SamReader samReaderFirst;
	private SamReader samReaderSecond;
	private String firstFlag="f";
	private String secondFlag ="s";
	private boolean processFirstChrom = false;

	//current chrom
	private String workingChrom;
	private LinkedHashMap<String, SAMRecord> workingNameRecord = new LinkedHashMap<String, SAMRecord>();
	private HashSet<String> checked = new HashSet<String>();
	private boolean printMisMatches = false;

	public SamComparator (String[] args){
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);

			makeGzippers();

			walk();

			closeGzippersAndSort();
			
			printStats();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");

		} catch (Exception e) {
			System.err.println("\nProblem running app, discard results! Aborting.\n");
			e.printStackTrace();
		}
	}

	private void printStats() {
		System.out.println("\nParsing stats:");
		System.out.println(numberMatches+"\tNumber matches");
		System.out.println(numberMisMatchesFirst+"\tNumber mismatches first");
		System.out.println(numberMisMatchesSecond+"\tNumber mismatches second");
		System.out.println(numberChangedAlignments+"\tNumber changed alignments present in files");
	}

	private void makeGzippers() {
		try {
			String name = Misc.removeExtension(firstSamFile.getName());
			gzippedFiles[0] = new File (saveDirectory, name+"_Match.sam.gz");
			goodSamFirst = new Gzipper(gzippedFiles[0]);
			gzippedFiles[1] = new File (saveDirectory, name+"_MisMatch.sam.gz");
			badSamFirst = new Gzipper(gzippedFiles[1]);
			name = Misc.removeExtension(secondSamFile.getName());
			gzippedFiles[2] = new File (saveDirectory, name+"_Match.sam.gz");
			goodSamSecond = new Gzipper(gzippedFiles[2]);
			gzippedFiles[3] = new File (saveDirectory, name+"_MisMatch.sam.gz");
			badSamSecond = new Gzipper(gzippedFiles[3]);
			for (File f : gzippedFiles) f.deleteOnExit();

		} catch (Exception e) {
			System.err.println("\nProblem making gzippers! ");
			e.printStackTrace();
		} 
	}

	private void closeGzippersAndSort(){
		try {
			goodSamFirst.close();
			badSamFirst.close();
			goodSamSecond.close();
			badSamSecond.close();
			
			for (int i=0; i<gzippedFiles.length; i++){
				File bam = new File (saveDirectory, Misc.removeExtension(gzippedFiles[i].getName())+".bam");
				new PicardSortSam(gzippedFiles[i], bam);
			}
			
		} catch (IOException e) {
			System.err.println("\nProblem closing gzippers! ");
			e.printStackTrace();
		}
	}

	private void walk() throws Exception{
		//get readers and iterators 
		SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		samReaderFirst = factory.open(firstSamFile);
		samReaderSecond = factory.open(secondSamFile);
		
		//get and check chromosomes
		List<SAMSequenceRecord> firstChroms = samReaderFirst.getFileHeader().getSequenceDictionary().getSequences();
		List<SAMSequenceRecord> secondChroms = samReaderSecond.getFileHeader().getSequenceDictionary().getSequences();
		if (checkChromNames(firstChroms, secondChroms) == false) {
			samReaderFirst.close();
			samReaderSecond.close();
			Misc.printErrAndExit("\nDifferent chromosomes are present in the two sam/bam files, aborting.\n");
		}
		
		//add headers to output files
		addHeaders();
		
		//for each chrom
		for (SAMSequenceRecord r : firstChroms){
		//for ( int i=0; i< 1; i++){
			workingChrom = r.getSequenceName();
			//workingChrom="chr20";
			
			System.out.println(workingChrom);
			
			//load first records into a hash
			loadFirstHash();
			
			//walk through second records
			SAMRecordIterator it = samReaderSecond.queryOverlapping(workingChrom, 0, 0);
			while (it.hasNext()){
				SAMRecord second = it.next();
				//make name
				String name = second.getReadName();
				if (second.getReadPairedFlag()){
					if (second.getFirstOfPairFlag()) name =  name + firstFlag;
					else name = name + secondFlag;
				}
				
				//look for it in the hash
				SAMRecord first = workingNameRecord.get(name);
				if (first == null) {
					badSamSecond.println(second.getSAMString().trim());
					numberMisMatchesSecond++;
				}
				
				//match
				else matchAndSave(first, second);
			}
			it.close();
			
			//walk over hash of firsts
			for (String n: workingNameRecord.keySet()){
				if (checked.contains(n)) continue;
				SAMRecord s = workingNameRecord.get(n);
				badSamFirst.println(s.getSAMString().trim());
				numberMisMatchesFirst++;
			}
		}
		
		//stop?
		if (processFirstChrom) return;
	}
	
	

	private void addHeaders() throws IOException {
		String firstHeader = samReaderFirst.getFileHeader().getTextHeader().trim();
		goodSamFirst.println(firstHeader);
		badSamFirst.println(firstHeader);
		String secondHeader = samReaderSecond.getFileHeader().getTextHeader().trim();
		goodSamSecond.println(secondHeader);
		badSamSecond.println(secondHeader);
	}

	private boolean checkChromNames(List<SAMSequenceRecord> firstChroms, List<SAMSequenceRecord> secondChroms) {
		//diff lengths?
		if (firstChroms.size() != secondChroms.size()) return false;
		//check names
		HashSet<String> firstNames = new HashSet<String>();
		for (SAMSequenceRecord r : firstChroms) firstNames.add(r.getSequenceName());
		for (SAMSequenceRecord r : secondChroms) {
			if (firstNames.contains(r.getSequenceName()) == false) return false;
		}
		return true;
	}

	private void loadFirstHash() {
		workingNameRecord.clear();
		checked.clear();
		SAMRecordIterator it = samReaderFirst.queryOverlapping(workingChrom, 0, 0);
		//add remaining until chrom name changes
		while (it.hasNext()){
			SAMRecord next = it.next();
			String name = next.getReadName();
			if (next.getReadPairedFlag()){
				if (next.getFirstOfPairFlag()) name =  name + firstFlag;
				else name = name + secondFlag;
			}
			workingNameRecord.put(name, next);
		}
		it.close();
	}

	/*Just checks chrom and alignment position*/
	private void matchAndSave(SAMRecord samFirst, SAMRecord samSecond) throws IOException {
		boolean match = false;
		//check chrom
		if (samFirst.getReferenceName().equals(samSecond.getReferenceName())) {
			//check position
			if (samFirst.getUnclippedStart() == samSecond.getUnclippedStart()) {
				//check sequence?
				if (checkSequence){
					if (samFirst.getReadString().equals(samSecond.getReadString())) match = true;
				}
				else match = true;
			}
			
		}
		//match?
		if (match){
			numberMatches++;
			goodSamFirst.println(samFirst.getSAMString().trim());
			goodSamSecond.println(samSecond.getSAMString().trim());
		}
		else {
			numberChangedAlignments++;
			numberMisMatchesFirst++;
			numberMisMatchesSecond++;
			badSamFirst.println(samFirst.getSAMString().trim());
			badSamSecond.println(samSecond.getSAMString().trim());
			
			if (printMisMatches) {
				System.out.println("First\t"+samFirst.getSAMString().trim());
				System.out.println("Second\t"+samSecond.getSAMString().trim());
			}
			
		}
		
		String name = samFirst.getReadName();
		if (samFirst.getReadPairedFlag()){
			if (samFirst.getFirstOfPairFlag()) name =  name + firstFlag;
			else name = name + secondFlag;
		}
		checked.add(name);
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SamComparator(args);
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
					case 'a': firstSamFile = new File(args[++i]); break;
					case 'b': secondSamFile = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'p': printMisMatches = true; break;
					case 'f': processFirstChrom = true; break;
					case 'e': checkSequence = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (firstSamFile == null || firstSamFile.canRead() == false) Misc.printErrAndExit("\nCan't find your first sam/bam file? Aborting!\n");
		if (secondSamFile == null || secondSamFile.canRead() == false) Misc.printErrAndExit("\nCan't find your second sam/bam file? Aborting!\n");
		if (saveDirectory == null) Misc.printErrAndExit("\nCan't find your save directory? Aborting!\n");
		saveDirectory.mkdirs();

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Sam Comparator  : Dec 2015                           **\n" +
				"**************************************************************************************\n" +
				"Compares coordinate sorted, unique, alignment sam/bam files.  Splits alignments into\n"+
				"those that match or mismatch chrom and position (or sequence).\n\n"+

				"Required:\n"+
				"-a Full path sam/bam file name. zip/gz OK.\n"+
				"-b Full path sam/bam file name. zip/gz OK.\n"+
				"-s Full path to a directory to save the results.\n" +
				"-p Print paired mismatches to screen.\n" +
				"-f Only process first chrom, defaults to all.\n"+
				"-e Check sequence of pairs.\n\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/SamComparator -a /hg19/ref.sam.gz\n" +
				"       -b /hg19/alt.sam.gz -s /hg19/SplitAlignments/\n\n" +

				"**************************************************************************************\n");

	}
}
