package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;

import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.*;

import java.util.*;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.ValidationStringency;
import edu.utah.seq.data.*;
import edu.utah.seq.data.sam.PicardSortSam;
import edu.utah.seq.data.sam.SamAlignment;

/**Parses a Novoalign bisulfite alignment sam/bam files. PointData scores are set to 1.
 * @author david.nix@hci.utah.edu 
 **/
public class NovoalignBisulfiteParser{
	//fields
	private File[] samFiles;
	private String versionedGenome;
	private HashMap<String, File> chromosomeFastaFiles = new HashMap<String, File>();
	private File workingFile;
	private File saveDirectory;
	private File nonConvertedCPointDataDirectory;
	private File convertedCPointDataDirectory;
	private float minimumPosteriorProbability = 13;
	private float maximumAlignmentScore = 300;
	private int minimumBaseScore = 13;
	private int minimumStandAloneBaseScore = 13;
	private boolean printBed = false;
	private boolean removeDuplicateReads = false;
	public static final Pattern COMMENT = Pattern.compile("^#.*");
	public static final Pattern CTStart = Pattern.compile("CT.*");
	public static final Pattern MD_BIS_SUB_MATCHER= Pattern.compile("(\\d+)([^0-9])");
	private boolean samFormat = true;
	private boolean useParsedFiles = false;
	private String adapter = "chrAdapt";
	private String phiX = "chrPhiX";
	private String alt = "_alt";
	private PrintWriter nonConvertedCs = null;
	private PrintWriter convertedCs = null;
	private HashMap<String, Long> nonConvertedCTypes = new HashMap<String, Long>();
	private HashMap<String, Long> convertedCTypes = new HashMap<String, Long>();
	private HashMap <String, DataOutputStream> chromOut = new HashMap <String, DataOutputStream>();
	private ArrayList<File> parsedBinaryDataFiles = new ArrayList<File>();
	private ArrayList<Point> nonConvertedCPoints = new ArrayList<Point>();
	private ArrayList<Point> convertedCPoints = new ArrayList<Point>();
	
	private String strand;
	private String chromosome = "";
	private String genomicSequence = null;
	private int genomicSequenceLengthMinus3 = 0;
	private long totalConvertedCsSequenced = 0;
	private long totalNonConvertedCsSequenced = 0;
	private long bpPairedSequence = 0;
	private long bpPairedOverlappingSequence = 0;
	private double bpPassQual = 0;
	private double bpFailQual = 0;
	private ArrayList<BaseObservation> boAL = new ArrayList<BaseObservation>();
	private HashMap<Integer, BaseObservation> firstHM = new HashMap<Integer, BaseObservation>();
	
	//alignment counts for sam files
	private long numberAlignmentsFailingQualityScore = 0;
	private long numberAlignmentsFailingAlignmentScore = 0;
	private long numberControlAlignments = 0;
	private long numberAlignmentsFailingQC = 0;
	private long numberAlignmentsUnmapped = 0;
	private long numberPassingAlignments = 0;

	//constructors
	public NovoalignBisulfiteParser(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		
		try {
			if (removeDuplicateReads) {
				//for each input alignment file, sort sams and remove duplicate reads
				System.out.println("\nSorting alignment data and removing duplicate reads...\n");
				for (int i=0; i < samFiles.length; i++) {
					File sortFile = new File(samFiles[i].toString().replaceAll(".sam.gz", "_sort.bam"));

					//call Picard's SortSam to coordinate-sort input sam files
					new PicardSortSam(samFiles[i], sortFile, true);
					File dupeFile = new File(samFiles[i].toString().replaceAll(".sam.gz", "_dup.bam"));
					File metricsFile = new File(samFiles[i].toString().replaceAll(".sam.gz", "_metrics.txt"));

					//call Picard's MarkDuplicates to remove duplicate reads
					//new PicardMarkDuplicates(sortFile, dupeFile, metricsFile);
					//File sortFileIndex = new File (sortFile.toString().replaceAll("_sort.bam", "_sort.bai")); 

					//delete intermediate bam/bai files
					//sortFile.delete();
					//sortFileIndex.delete();
				}
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		}
				
		System.out.println("Splitting text alignment data by chromosome and filtering...");

		//look for parsed files 1st
		if (useParsedFiles && fetchParsedFiles()){
			System.err.println("\nWARNING: Parsed files found! Using these instead.  Be sure prior txt parsing completed. Delete the chrXXXX+/- files from "+saveDirectory+" if you would like to reparse the txt alignment files.");
		}
		else {
			//for each file, parse, filter, split by chrom and strand and save to disk	
			for (int i=0; i< samFiles.length; i++) {
				if (removeDuplicateReads) {
					//set working objects and parse tag file text
					workingFile = new File(samFiles[i].toString().replaceAll(".sam.gz", "_dup.bam"));
				}
				else {
					workingFile = samFiles[i];
				}
				System.out.print("\t"+workingFile);
				//parse
				if (parseWorkingSAMFile() == false) Misc.printErrAndExit("\nERROR: failed to parse, aborting!");
			}
			closeWriters();
			
			//Alignment filtering stats
			double total = numberAlignmentsFailingQualityScore + numberAlignmentsFailingAlignmentScore + numberControlAlignments + numberAlignmentsFailingQC + numberAlignmentsUnmapped + numberPassingAlignments;
			System.out.println("\nFiltering statistics for "+(int)total+" alignments:");
			System.out.println(numberAlignmentsFailingQualityScore +"\tFailed mapping quality score ("+minimumPosteriorProbability+")");
			System.out.println(numberAlignmentsFailingAlignmentScore +"\tFailed alignment score ("+maximumAlignmentScore+")");
			System.out.println(numberControlAlignments +"\tAligned to phiX, adapters, or alt");
			System.out.println(numberAlignmentsFailingQC +"\tFailed vendor QC");
			System.out.println(numberAlignmentsUnmapped +"\tAre unmapped\n");
			
			System.out.println(numberPassingAlignments +"\tPassed filters ("+Num.formatPercentOneFraction(((double)numberPassingAlignments)/ total)+")");

		}
		//for each binary data file, parse methylC's and nonMethylC's
		System.out.print("\nParsing binary data...\n\t");
		parseBinaryData();
		System.out.println();

		//print stats
		System.out.println("\nCounts for non-converted C contexts:\n"+nonConvertedCTypes+"\n");
		System.out.println("\nCounts for converted C contexts:\n"+convertedCTypes+"\n");
		System.out.println(totalNonConvertedCsSequenced +"\tTotal non-converted Cs sequenced");
		System.out.println(totalConvertedCsSequenced +"\tTotal converted Cs sequenced");
		double fractionNonConverted = (double)(totalNonConvertedCsSequenced)/(double)(totalNonConvertedCsSequenced+totalConvertedCsSequenced);
		System.out.println(Num.formatNumber(fractionNonConverted,3) +"\tFraction non converted C's.");
		
		System.out.println();
		double fractionPassQual = bpPassQual / (bpFailQual+ bpPassQual);
		System.out.println(Num.formatNumber(fractionPassQual,3)+ "\tFraction bp passing quality ("+minimumBaseScore+")"); 


		if (bpPairedOverlappingSequence != 0){
			System.out.println();
			System.out.println(bpPairedOverlappingSequence+"\tBPs overlapping paired sequence");
			System.out.println(bpPairedSequence+"\tBPs paired sequence");
			double fractionOverlap = (double)(bpPairedOverlappingSequence)/(double)bpPairedSequence;
			System.out.println(Num.formatNumber(fractionOverlap,3) +"\tFraction overlapping bps from paired reads.");
		}

		//delete temp files
		for (int i=0; i< parsedBinaryDataFiles.size(); i++) parsedBinaryDataFiles.get(i).delete();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}
	
	/**For processing blocks of SAMRecords*/
	public NovoalignBisulfiteParser (String chromosome, String genomicSequence){
		this.chromosome = chromosome;
		this.genomicSequence = genomicSequence;
		genomicSequenceLengthMinus3 = genomicSequence.length()-3;
	}

	public boolean fetchParsedFiles(){
		Pattern p = Pattern.compile(".+[-+]$");
		File[] files = IO.extractFiles(saveDirectory);
		for (int i=0; i< files.length; i++){
			if (p.matcher(files[i].getName()).matches()) parsedBinaryDataFiles.add(files[i]);
		}
		if (parsedBinaryDataFiles.size()!=0) return true;
		return false;
	}

	public void parseBinaryData(){
		//make PrintWriters
		try {
			if (printBed) {
				nonConvertedCs = new PrintWriter( new FileWriter( new File( saveDirectory, "nonConvertedCs.bed")));
				convertedCs = new PrintWriter( new FileWriter( new File( saveDirectory, "convertedCs.bed")));
			}

		} catch (IOException e1) {
			e1.printStackTrace();
		}

		//for each binary file
		File[] files = new File[parsedBinaryDataFiles.size()];
		parsedBinaryDataFiles.toArray(files);
		Arrays.sort(files);

		//for each chromosome
		String oldChrom = "";
		for (int i=0; i< files.length; i++){
			chromosome = files[i].getName();
			strand = chromosome.substring(chromosome.length()-1);
			chromosome = chromosome.substring(0,chromosome.length()-1);
			int readLength =0;
			//fetch the sequence
			if (oldChrom.equals(chromosome) ==false){
				File seqFile = chromosomeFastaFiles.get(chromosome);
				if (seqFile == null) {
					System.err.println("\n\tWarning, could not find a fasta file for "+chromosome+", skipping!");
					continue;
				}
				else System.out.print(chromosome+" ");				
				MultiFastaParser mfp = new MultiFastaParser(seqFile);
				genomicSequence = mfp.getSeqs()[0];
				genomicSequenceLengthMinus3 = genomicSequence.length()-3;
				oldChrom = chromosome;
			}

			//parse binary data: readID, position, sequence, baseScores, baseObservations
			DataInputStream dis = null;
			ParsedAlignment oldPA = null;
			try {
				dis = new DataInputStream(new BufferedInputStream(new FileInputStream(files[i])));
				while (true){
					//System.out.println("Reading "+files[i].getName());
					//read in data and parse
					ParsedAlignment pa = new ParsedAlignment(dis, samFormat, this);
					if (pa.isPassesParsing() == false) continue;
					else if (oldPA == null) oldPA = pa;
					//paired alignment
					else if (oldPA.getReadID().equals(pa.getReadID())){					
						//calc overlap
						int leftPos = oldPA.getPosition();
						int rightPos = pa.getPosition();
						int l;
						int r;
						BaseObservation[] first;
						BaseObservation[] second;

						if (leftPos > rightPos) {
							l = rightPos;
							r = leftPos;
							first = pa.getBo();
							second = oldPA.getBo();
						}
						else {
							l = leftPos;
							r = rightPos;
							second = pa.getBo();
							first = oldPA.getBo();
						} 
						int stopL = l+ oldPA.getBases().length;
						long diff = stopL - r;							
						//any overlap?
						if (diff > 0) {
							addBaseObservations( processPairs(r, stopL, first, second) );
						}
						else {
							addBaseObservations(first);
							addBaseObservations(second);
						}
						oldPA = null;
					}
					//non paired alignment
					else {
						addBaseObservations(oldPA.getBo());
						oldPA = pa;
					}
				}
			} catch (EOFException eof){				
				//add last
				if (oldPA != null) addBaseObservations(oldPA.getBo());
				//write out PointData?
				Point[] p = null;
				PointData pd = null;
				Info info = null;
				//methylC
				if (nonConvertedCPoints.size() !=0){
					p = new Point[nonConvertedCPoints.size()];
					nonConvertedCPoints.toArray(p);
					nonConvertedCPoints.clear();
					Arrays.sort(p, new ComparatorPointPosition());
					p= Point.sumIdenticalPositionScores(p);
					pd = Point.extractPositionScores(p);
					//make an Info object public Info (String name, String versionedGenome, String chromosome, String strand, int readLength, HashMap<String,String> notes){
					info = new Info("NonConvertedC", versionedGenome, chromosome, strand, readLength, null);
					pd.setInfo(info);
					//save
					pd.writePointData(nonConvertedCPointDataDirectory);
				}
				//non methyl C
				if (convertedCPoints.size() !=0){
					p = new Point[convertedCPoints.size()];
					convertedCPoints.toArray(p);
					convertedCPoints.clear();
					Arrays.sort(p, new ComparatorPointPosition());
					p= Point.sumIdenticalPositionScores(p);
					pd = Point.extractPositionScores(p);
					info = new Info("ConvertedC", versionedGenome, chromosome, strand, readLength, null);
					pd.setInfo(info);
					//save
					pd.writePointData(convertedCPointDataDirectory);
				}
			}
			catch (Exception e){
				e.printStackTrace();
				Misc.printErrAndExit("\nError encountered in parsing this binary file? "+files[i]);
			} finally {
				if (dis != null) {
					try {
						dis.close();
					} catch (IOException ignored) {}
				}
			}
		}

		//close print writers
		if (printBed) {
			nonConvertedCs.close();
			convertedCs.close();
		}
	}

	/*This is called by the AllelicMethylationDetector and not NBP*/
	public ArrayList<AlignmentPair> parseSamAlignments(ArrayList<SAMRecord> samAL){
		
				//make ParsedAlignments from SAM records
				ParsedAlignment[] pas = new ParsedAlignment[samAL.size()];
				for (int i=0; i< pas.length; i++) pas[i] = new ParsedAlignment(samAL.get(i), this);
		
				//sort by readID
				Arrays.sort(pas);
				
				//merge pairs if present
				ArrayList<AlignmentPair> pairs = new ArrayList<AlignmentPair>();
				ParsedAlignment oldPA = null;
				for (ParsedAlignment pa : pas){ 
					if (pa.isPassesParsing() == false) continue;
					else if (oldPA == null) oldPA = pa;
					//paired alignment
					else if (oldPA.getReadID().equals(pa.getReadID())){					
						//calc overlap
						int leftPos = oldPA.getPosition();
						int rightPos = pa.getPosition();
						int l;
						int r;
						BaseObservation[] first;
						BaseObservation[] second;
						ParsedAlignment firstPA;
						ParsedAlignment secondPA;

						if (leftPos > rightPos) {
							l = rightPos;
							r = leftPos;
							first = pa.getBo();
							second = oldPA.getBo();
							firstPA = pa;
							secondPA = oldPA;
						}
						else {
							l = leftPos;
							r = rightPos;
							second = pa.getBo();
							first = oldPA.getBo();
							secondPA = pa;
							firstPA = oldPA;
						} 
						int stopL = l+ oldPA.getBases().length;
						long diff = stopL - r;							
						//any overlap?
						if (diff > 0) {
							pairs.add(new AlignmentPair (firstPA, secondPA, processPairs(r, stopL, first, second)));
						}
						else {
							pairs.add(new AlignmentPair (firstPA, secondPA));
						}
						oldPA = null;
					}
					//non paired alignment
					else {
						oldPA = pa;
						pairs.add(new AlignmentPair (oldPA));
					}
				}
						
				//add last
				if (oldPA != null) pairs.add(new AlignmentPair (oldPA));
				
			return pairs;
			
	}

	/**Closes writers.*/
	public void closeWriters(){
		try{
			Iterator<DataOutputStream> it = chromOut.values().iterator();
			while (it.hasNext()) it.next().close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public boolean parseWorkingSAMFile(){
		try{
			//make reader
			SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(workingFile);
			SAMRecordIterator iterator = reader.iterator();
			
			int counter =0;
			String currentChromStrand = "";
			DataOutputStream dos = null;
			int numBadLines = 0;
			//for each record
			while (iterator.hasNext()){
				SAMRecord samRecord = iterator.next();
				
				//print status blip
				if (++counter == 2500000){
					System.out.print(".");
					counter = 0;
				}

				//this is a bit inefficient but gives absolute control on the sam data
				SamAlignment sa;
				try {
					sa = new SamAlignment(samRecord.getSAMString().trim(), false);				
				} catch (Exception e) {
					System.out.println("\nSkipping malformed sam alignment ->\n"+samRecord.getSAMString()+"\n"+e.getMessage());
					if (numBadLines++ > 1000) Misc.printErrAndExit("\nAboring: too many malformed SAM alignments.\n");
					continue;
				}
				
				//is it aligned?
				if (sa.isUnmapped()) {
					numberAlignmentsUnmapped++;
					continue;
				}
				
				//does it pass the vendor qc?
				if (sa.failedQC()) {
					numberAlignmentsFailingQC++;
					continue;
				}
				
				//skip phiX and adapter
				String ref = sa.getReferenceSequence();
				if (ref.startsWith(phiX) || ref.startsWith(adapter) || ref.endsWith(alt)) {
					numberControlAlignments++;
					continue;
				}

				//does it pass the scores threshold?
				if (sa.getAlignmentScore() > maximumAlignmentScore) {
					numberAlignmentsFailingAlignmentScore++;
					continue;
				}
				if (sa.getMappingQuality() < minimumPosteriorProbability) {
					numberAlignmentsFailingQualityScore++;
					continue;
				}

				//increment counter
				numberPassingAlignments++;
				
				//make readID
				String readID = sa.getName();
				
				//set position
				int position = sa.getPosition();
				
				//clear masking references
				sa.trimMaskingOfReadToFitAlignment();

				//set isGA and strand based the ZB:Z:CT or GA flag
				String isGA; 
				String chromosomeStrand = null;
				int za = sa.getCtGaTag();
				if (za == 1) {
					isGA = "t";
					chromosomeStrand = sa.getReferenceSequence() + "-";
				}
				else if (za == 2) {
					isGA = "f";
					chromosomeStrand = sa.getReferenceSequence() + "+";
				}
				else throw new Exception("\nError: Failed to find a ZB:Z:GA or ZB:Z:CT tag? Is this an old novoalignment? Update and realign?\n");
				
				//get PrintWriter
				if (currentChromStrand.equals(chromosomeStrand) == false){
					currentChromStrand = chromosomeStrand;
					if (chromOut.containsKey(currentChromStrand)) dos = chromOut.get(currentChromStrand);
					else {
						//make and set file
						File f = new File(saveDirectory, currentChromStrand);
						parsedBinaryDataFiles.add(f);
						dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(f)));
						chromOut.put(currentChromStrand, dos);
					}
				}
				
				//save data: readID, position, sequence, baseScores, t or f + cigar
				dos.writeUTF(readID);
				dos.writeInt(position);
				dos.writeUTF(sa.getSequence());
				dos.writeUTF(sa.getQualities());
				dos.writeUTF(isGA + sa.getCigar());

				//calculate fraction overlap for paired reads
				//is it a L read?
				if (sa.isFirstPair() && sa.isPartOfAPairedAlignment()){
					//check that chrom is same
					if (sa.getMateReferenceSequence().equals("=") || sa.getMateReferenceSequence().equals(sa.getReferenceSequence())) {
						//calc diff
						int leftPos = position;
						int rightPos = sa.getMatePosition();
						int l;
						int r;
						if (leftPos > rightPos) {
							l = rightPos;
							r = leftPos;
						}
						else {
							l = leftPos;
							r = rightPos;
						} 
						int sequenceLength = sa.getSequence().length();
						int stopL = l+ sequenceLength;
						long diff = stopL - r;
						if (diff >= 0) bpPairedOverlappingSequence += diff;
						bpPairedSequence += (2 * sequenceLength);
					}
				}
			}
			reader.close();
			System.out.println();
		} catch (Exception e){
			System.err.println("\nError parsing Novoalign file or writing split binary chromosome files.\n");
			e.printStackTrace();
			return false;
		}
		return true;
	}



	private void addBaseObservations (BaseObservation[] bos){
		for (int i=0; i< bos.length; i++) {
			if (bos[i].isConverted()){
				if (printBed) convertedCs.println(bos[i].getBedLine());
				convertedCPoints.add(bos[i].getPoint());
				totalConvertedCsSequenced++;
				long count = 1;
				if (convertedCTypes.containsKey(bos[i].getGenSeq())){
					count = convertedCTypes.get(bos[i].getGenSeq()).longValue() + 1;
				}
				convertedCTypes.put(bos[i].getGenSeq(), new Long(count));
			}
			else {
				if (printBed) nonConvertedCs.println(bos[i].getBedLine());
				nonConvertedCPoints.add(bos[i].getPoint());
				totalNonConvertedCsSequenced++;
				long count = 1;
				if (nonConvertedCTypes.containsKey(bos[i].getGenSeq())){
					count = nonConvertedCTypes.get(bos[i].getGenSeq()).longValue() + 1;
				}
				nonConvertedCTypes.put(bos[i].getGenSeq(), new Long(count));
			}
		}
	}

	/**Flattens two overlapping BaseObservations.  Only saves consensus BOs in the range defined by start and stop(excluded) or where one score is much better than the other.*/
	private BaseObservation[] processPairs(int start, int stop, BaseObservation[] first, BaseObservation[] second){
		boAL.clear();
		firstHM.clear();

		//add BO's from first that preceed the range, add others to a hash for fast look up
		for (int i=0; i< first.length; i++) {
			if (first[i].getPosition().intValue() < start) boAL.add(first[i]);
			else firstHM.put(first[i].getPosition(), first[i]);
		}

		//add BO's from second that exceed the range, check others against the hash
		for (int i=0; i< second.length; i++) {
			if (second[i].getPosition().intValue() >= stop) boAL.add(second[i]);
			//within the overlap
			else {
				//check for it in the first
				BaseObservation left = firstHM.get(second[i].getPosition());
				//found in first!
				if (left != null){
					//check if both the same and add 
					if (left.isConverted() == second[i].isConverted()) boAL.add(left);
					else {						
						//compare scores, if close then skip
						double fraction = left.getScore()/ second[i].getScore();
						if (fraction >= 1.5) {
							boAL.add(left);
						}
						else if (fraction <= 0.667) boAL.add(second[i]);
					}
					//remove the BO
					firstHM.remove(second[i].getPosition());
				}
				//not found in first, check for high quality score
				else if (second[i].getScore() >= minimumStandAloneBaseScore) boAL.add(second[i]);
			}
		}

		//scan remaining BOs in the first to see if they should be added
		if (firstHM.size() !=0){
			Iterator<BaseObservation> c = firstHM.values().iterator();
			while (c.hasNext()){
				BaseObservation bo = c.next();
				if (bo.getScore() >= minimumStandAloneBaseScore) boAL.add(bo);
			}
		}

		//return, no need to sort, this is done later
		first = new BaseObservation[boAL.size()];
		boAL.toArray(first);
		return first;
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new NovoalignBisulfiteParser(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File[] fastas = null;
		File forExtraction = null;
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': forExtraction = new File(args[i+1]); i++; break;
					case 'f': fastas = IO.extractFiles(new File(args[++i])); break;
					case 's': saveDirectory = new File(args[++i]); saveDirectory.mkdir(); break;
					case 'v': versionedGenome = args[i+1]; i++; break;
					case 'p': printBed = true; break;
					case 'z': useParsedFiles = true; break;
					case 'd': removeDuplicateReads = true; break;
					case 'b': minimumBaseScore = Integer.parseInt(args[++i]);break;
					case 'c': minimumStandAloneBaseScore = Integer.parseInt(args[++i]);break;
					case 'x': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'q': minimumPosteriorProbability = Float.parseFloat(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//fetch sam files
		File[][] tot = new File[4][];
		tot[0] = IO.extractFiles(forExtraction,".sam");
		tot[1] = IO.extractFiles(forExtraction,".sam.gz");
		tot[2] = IO.extractFiles(forExtraction,".sam.zip");
		tot[3] = IO.extractFiles(forExtraction,".bam");
		samFiles = IO.collapseFileArray(tot);
		if (samFiles == null || samFiles.length ==0 || samFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.sam(.zip/.gz) file(s)!\n");

		if (versionedGenome == null) Misc.printErrAndExit("\nPlease provide a versioned genome (e.g. H_sapiens_Mar_2006).\n");
		if (fastas == null || fastas.length ==0) Misc.printErrAndExit("\nError: cannot find any fasta sequence files?\n");
		if (samFiles == null || samFiles.length ==0 || samFiles[0].canRead() == false) Misc.printErrAndExit("\nError: cannot find your alignment file(s)!\n");
		if (saveDirectory == null) Misc.printErrAndExit("\nPlease enter a directory to use in saving your results.\n");
		saveDirectory.mkdirs();
		nonConvertedCPointDataDirectory = new File (saveDirectory, "NonConvertedC");
		nonConvertedCPointDataDirectory.mkdir();
		convertedCPointDataDirectory = new File (saveDirectory, "ConvertedC");
		convertedCPointDataDirectory.mkdir();

		//load fastaFiles into hash
		chromosomeFastaFiles = new HashMap<String,File>();
		Pattern chrom = Pattern.compile("(.+)\\.fa.*");
		for (int i=0; i< fastas.length; i++){
			Matcher mat = chrom.matcher(fastas[i].getName());
			if (mat.matches()) chromosomeFastaFiles.put(mat.group(1), fastas[i]);
		}
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Novoalign Bisulfite Parser: May 2016                      **\n" +
				"**************************************************************************************\n" +
				"Parses Novoalign -b2 and -b4 single and paired bisulfite sequence alignment files into\n" +
				"PointData file formats. Generates several summary statistics on converted and non-\n" +
				"converted C contexts. Flattens overlapping reads in a pair to call consensus bps.\n" +
				"Note: for paired read RNA-Seq data run through the SamTranscriptomeParser first.\n" +

				"\nOptions:\n"+
				"-a Alignment file or directory containing non merged novoalignments in SAM/BAM\n"+
				"      (xxx.sam(.zip/.gz OK) or xxx.bam) format. Multiple files are combine. \n" +
				"-f Fasta file directory, chromosome specific xxx.fa/.fasta(.zip/.gz OK) files.\n" +
				"-s Save directory.\n"+
				"-v Versioned Genome (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +

				"\nDefault Options:\n"+
				"-p Print bed file parsed data.\n"+
				"-x Maximum alignment score. Defaults to 300, smaller numbers are more stringent.\n"+
				"-q Minimum mapping quality score. Defaults to 13, bigger numbers are more stringent.\n" +
				"      This is a phred-scaled posterior probability that the mapping position of read\n" +
				"      is incorrect. For RNASeq data, set this to 0.\n" +
				"-b Minimum base quality score for reporting a non/converted C, defaults to 13.\n"+
				"-c Minimum base quality score for reporting a overlapping non/converted C not found\n" +
				"      in the other pair, defaults to 13.\n"+
				"-d Remove duplicate reads prior to generating PointData. Defaults to not removing\n" +
				"      duplicates.\n" +

				"\nExample: java -Xmx25G -jar pathToUSeq/Apps/NovoalignBisulfiteParser -x 240 -a\n" +
				"      /Novo/Run7/ -f /Genomes/Hg19/Fastas/ -v H_sapiens_Feb_2009 -s /Novo/Run7/NBP \n\n" +


		"**************************************************************************************\n");

	}

	public String getGenomicSequence() {
		return genomicSequence;
	}

	public int getGenomicSequenceLengthMinus3() {
		return genomicSequenceLengthMinus3;
	}

	public int getMinimumBaseScore() {
		return minimumBaseScore;
	}

	public String getChromosome() {
		return chromosome;
	}

	public String getStrand() {
		return strand;
	}

	public void incrementBPPassQual() {
		bpPassQual++;
	}

	public void incrementBPFailQual() {
		bpFailQual++;
	}	

}
