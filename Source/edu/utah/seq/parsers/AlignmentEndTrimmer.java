package edu.utah.seq.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import util.gen.Misc;



public class AlignmentEndTrimmer {
	private File refFile = null;
	private File inputFile = null;
	private File outputFile = null;
	private HashMap<String,String> refHash = null;
	
	private int mScore = 1;
	private int mmScore = 2;
	private int minLength = 10;
	private boolean verbose = false;
	private boolean rnaEditing = false;
	
	//mismatch scoring
	private int maxMM = 0;
	private boolean mmMode = false;
	
	public static HashMap<Character,Character> revCompHash = new HashMap<Character,Character>() {
		private static final long serialVersionUID = 1L;

	{
		put('A','T');
		put('T','A');
		put('C','G');
		put('G','C');
		put('N','N');
	}};
	
	public static void main(String[] args) {
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		new AlignmentEndTrimmer(args);
	}
	
	
	public AlignmentEndTrimmer(String[] args) {
		System.out.println("Parsing command line arguments");
		this.processArgs(args);
		
		System.out.println("Loading reference fasta");
		this.parseFasta();
		
		System.out.println("Processing data");
		this.processData();
		
		System.out.println("Finished");
		
	}
	
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
						case 'i': this.inputFile = new File (args[++i]); break;
						case 'r': this.refFile = new File(args[++i]); break;
						case 'o': this.outputFile = new File(args[++i]); break;
						case 'v': this.verbose = true; break;
						case 'm': this.mScore = Integer.parseInt(args[++i]); break;
						case 'n': this.mmScore = Integer.parseInt(args[++i]); break;
						case 'l': this.minLength = Integer.parseInt(args[++i]); break;
						case 'e': this.rnaEditing = true; break;
						case 's': this.mmMode = true; break;
						case 'x': this.maxMM = Integer.parseInt(args[++i]); break;
						case 'h': printDocs(); System.exit(0);
						default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}	
		}
		
		if (this.mScore < 0) {
			System.out.println("Match score must be zero or greater\n");
			System.exit(1);
		}
		
		if (this.mmScore < 0) {
			System.out.println("Mismatch score must be zero or greater\n");
			System.exit(1);
		}
		
		if (this.minLength < 0) {
			System.out.println("Min read length must be zero or greater\n");
			System.exit(1);
		}
		
		//Make sure refFile exists
		if (this.refFile == null) {
			System.out.println("Reference file has not been set, exiting");
			System.exit(1);
		}
		
		if (!this.refFile.canRead()) {
			System.out.println("The specified reference file cannot be read, exiting");
			System.exit(1);
		}
		
		//Make sure input File exists
		if (this.inputFile == null) {
			System.out.println("Input file has not been set, exiting");
			System.exit(1);
		}
		
		if (!this.inputFile.canRead()) {
			System.out.println("The specified file cannot be read, exiting");
			System.exit(1);
		}
		
	}
	
	public void processData() {
		
			//Set up file reader
		    SamReader sfr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(this.inputFile);
		    
			
			//Set up file writer
		    SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
		    sfwf.setCreateIndex(true);
		    SAMFileWriter sfw = sfwf.makeSAMOrBAMWriter(sfr.getFileHeader(), false, this.outputFile);
		    
		    
			
			
			//counters
			int totalReads = 0;
			int trimmedReads = 0;
			int droppedReads = 0;
			int dropTrimReads = 0;
			int dropMmReads = 0;
			
			//Containers
			HashSet<String> notFound = new HashSet<String>();
			HashMap<String, SAMRecord> mateList = new HashMap<String,SAMRecord>();
			HashSet<String> removedList = new HashSet<String>();
			HashMap<String,SAMRecord> editedList = new HashMap<String,SAMRecord>();
			
			for (SAMRecord sr: sfr) {
				//Messaging
				if (totalReads % 1000000 == 0 && totalReads != 0) {
					System.out.println(String.format("Finished processing %d reads. %d were trimmed, %d were set as unmapped. Currently storing mates for %d reads.",totalReads,trimmedReads,droppedReads,mateList.size()));
					
				}
				totalReads += 1;
				
				String keyToCheck = sr.getReadName() + ":" + String.valueOf(sr.getIntegerAttribute("HI"));
			
				//Make sure chromsome is available
				String chrom  = sr.getReferenceName();
				if (!this.refHash.containsKey(chrom)) {
					if (!notFound.contains(chrom)) {
						notFound.add(chrom);
						System.out.println(String.format("Chromosome %s not found in reference file, skipping trimming step", chrom));
					}
				} else if (!sr.getReadUnmappedFlag()) {
					String refSeq = null;
					String obsSeq = null;
					List<CigarElement> cigar = null;
					
					//Get necessary sequence information depending on orientation
					if (sr.getReadNegativeStrandFlag()) {
						refSeq = this.revComp(this.refHash.get(chrom).substring(sr.getAlignmentStart()-1,sr.getAlignmentEnd()));
						obsSeq = this.revComp(sr.getReadString());
						cigar = this.reverseCigar(sr.getCigar().getCigarElements());
						
					} else {
						refSeq = this.refHash.get(chrom).substring(sr.getAlignmentStart()-1,sr.getAlignmentEnd());
						obsSeq = sr.getReadString();
						cigar = sr.getCigar().getCigarElements();
					}
				
					
					//Get alignments
					String[] alns = this.createAlignmentStrings(cigar, refSeq, obsSeq, totalReads);
					
					//Identify Trim Point
					int idx = this.identifyTrimPoint(alns,sr.getReadNegativeStrandFlag());
					
					//Check error rate
					boolean mmPassed = false;
					if (mmMode) {
						mmPassed = this.isPoorQuality(alns, sr.getReadNegativeStrandFlag(), idx);
					}
					
					
					
					//Create new cigar string
					if (idx < minLength || mmPassed) {
						
						
						sr.setAlignmentStart(0);
						sr.setReadUnmappedFlag(true);
						sr.setProperPairFlag(false);
						sr.setReferenceIndex(-1);
						sr.setMappingQuality(0);
						sr.setNotPrimaryAlignmentFlag(false);
						
						
						
						if (sr.getReadPairedFlag() && !sr.getMateUnmappedFlag()) {
							if (mateList.containsKey(keyToCheck)) {
								mateList.put(keyToCheck, this.changeMateUnmapped(mateList.get(keyToCheck)));
							} else {
								removedList.add(keyToCheck);
							}
						} 
						droppedReads += 1;
						if (idx < minLength) {
							dropTrimReads += 1;
						} else {
							dropMmReads += 1;
						}
					} else if (idx+1 != alns[0].length()) {
						trimmedReads++;
						Cigar oldCig = sr.getCigar();
						Cigar newCig = this.createNewCigar(alns, cigar, idx, sr.getReadNegativeStrandFlag());
						sr.setCigar(newCig);
						
						if (sr.getReadNegativeStrandFlag()) {
							int newStart = this.determineStart(oldCig, newCig, sr.getAlignmentStart());
							sr.setAlignmentStart(newStart);
						}
						
						
						if (this.verbose) {
							this.printAlignments(sr, oldCig, alns, idx);
						}
						
						if (sr.getReadPairedFlag() && !sr.getMateUnmappedFlag()) {
							if (mateList.containsKey(keyToCheck)) {
								mateList.put(keyToCheck, this.changeMatePos(mateList.get(keyToCheck),sr.getAlignmentStart()));
							} else {
								editedList.put(keyToCheck,sr);
							}
						}
					}
				}
				
				//System.out.println(sr.getReadName());
				if (sr.getReadPairedFlag() && !sr.getMateUnmappedFlag()) {
					//String rn = sr.getReadName();
					if (mateList.containsKey(keyToCheck)) {
						if (editedList.containsKey(keyToCheck)) {
							sr = this.changeMatePos(sr,editedList.get(keyToCheck).getAlignmentStart());
							editedList.remove(keyToCheck);
						} else if (removedList.contains(keyToCheck)) {
							sr = this.changeMateUnmapped(sr);
							removedList.remove(keyToCheck);
						}
						sfw.addAlignment(sr);
						sfw.addAlignment(mateList.get(keyToCheck));
						mateList.remove(keyToCheck);
					} else {
						mateList.put(keyToCheck, sr);
					}
				} else {
					sfw.addAlignment(sr);
					if (mateList.containsKey(keyToCheck)) {
						mateList.remove(keyToCheck);
					}
				}
			}
			
			System.out.println(String.format("Finished processing %d reads. %d were trimmed, %d were set as unmapped. Of the unmapped, %d were too short and %d had too many mismatches. Currently storing mates for %d reads.",
					totalReads,trimmedReads,droppedReads,dropTrimReads, dropMmReads, mateList.size()));
			System.out.println(String.format("Reads left in hash: %d.  Writing to disk.",mateList.size()));
			for (SAMRecord sr2: mateList.values()) {
				sfw.addAlignment(sr2);
				
			}
			
			sfw.close();
			try {
				sfr.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}		
	}
	
	public SAMRecord changeMateUnmapped(SAMRecord sr) {
		sr.setMateAlignmentStart(0);
		sr.setMateUnmappedFlag(true);
		sr.setProperPairFlag(false);
		sr.setMateReferenceIndex(-1);

		return sr;	
	}
	
	public SAMRecord changeMatePos(SAMRecord sr, int pos) {
		//sr.setMateAlignmentStart(pos);
		return sr;
	}
	
	public int determineStart(Cigar origCigar, Cigar newCigar, int oldPos) {
		int newPos = oldPos;
		if (newCigar.getCigarElement(0).getOperator() == CigarOperator.S) {
			
			boolean diff = false;
			List<CigarElement> newEL = this.reverseCigar(newCigar.getCigarElements());
			List<CigarElement> oldEL = this.reverseCigar(origCigar.getCigarElements());
			
			
			for (int i=0;i<oldEL.size();i++) {
				CigarOperator con = oldEL.get(i).getOperator();
				int col = oldEL.get(i).getLength();
				if (!diff) {
					CigarOperator cnn = newEL.get(i).getOperator();
					int cnl = newEL.get(i).getLength();
					if (con != cnn) {
						diff = true;
					} else if (col != cnl) {
						diff = true;
						col = col - cnl;
					}
				}
				
				if (diff) {
					if (con == CigarOperator.M || con == CigarOperator.N || con == CigarOperator.D) {
						newPos += col;
					}
				}
			}
		}
		return newPos;
		
	}
	
//	public int determineStart(Cigar origCigar, Cigar newCigar, int oldPos) {
//		int newPos = oldPos;
//		int remaining = 0;
//		CigarElement ce = newCigar.getCigarElement(0);
//		if (ce.getOperator() == CigarOperator.S) {
//			remaining = ce.getLength();
//			
//			boolean finished = false;
//			int index = 0;
//			outerloop:
//			for (CigarElement o: origCigar.getCigarElements()) {
//				CigarOperator co = o.getOperator();
//		
//				for (int i=0;i<o.getLength();i++) {
//					if (remaining == 0) {
//						break outerloop;
//					}
//					finished = false;
//					if (co != CigarOperator.D) {
//						remaining -= 1;
//					}
//					if (co == CigarOperator.M || co == CigarOperator.N || co == CigarOperator.D) {
//						newPos += 1;
//					} 
//				}
//				finished = true;
//				index += 1;
//			}
//			if (finished) {
//				for (int i = index; i < origCigar.getCigarElements().size();i++) {
//					CigarElement last = origCigar.getCigarElement(i);
//					if (last.getOperator() == CigarOperator.N || last.getOperator() == CigarOperator.D) {
//						newPos += last.getLength();
//					} else {
//						break;
//					}
//				}
//			}
//			
//		} 
//		return newPos;
//	}
	
	
	public void printAlignments(SAMRecord sr, Cigar cig, String[] alns, int idx) {
		System.out.println(String.format("Read name: %s",sr.getReadName()));
		System.out.println(String.format("Read reversed: %b",sr.getReadNegativeStrandFlag()));
		System.out.println(String.format("Original Cigar: %s",cig.toString()));
		System.out.println(String.format("Updated Cigar : %s",sr.getCigarString()));
		System.out.println(String.format("Orig Rf Seq: %s", alns[0]));
		System.out.println(String.format("Orig Rd Seq: %s", alns[1]));
		System.out.println(String.format("Trim Rf Seq: %s", alns[0].substring(0,idx+1)));
		System.out.println(String.format("Trim Rd Seq: %s", alns[1].substring(0,idx+1)));
		System.out.println(String.format("Trim Idx: %d",idx));
		System.out.println();
	}
	
	public Cigar createNewCigar(String[] alignments, List<CigarElement> cigar, int maxIdx, boolean isReverse) {
		List<CigarElement> ce = new ArrayList<CigarElement>();
		
		char[] aln1 = alignments[0].toCharArray();
		char[] aln2 = alignments[1].toCharArray();
		
		//Hard trims at the beginning aren't included at the beginning of the alignment string, so add them back.
		if (cigar.get(0).getOperator() == CigarOperator.H) {
			ce.add(cigar.get(0));
		}
		
		int cigLen = 0;
		
		CigarOperator cigOp = null;
		
		for (int i = 0;i<aln1.length;i++) {
			int extLen = 0;
			CigarOperator currOp = null;
			if (i>maxIdx) {
				currOp = CigarOperator.S;
				if (aln2[i] != '-') {
					extLen = 1;
				} 
			} else if (aln1[i] == 'S' && aln2[i] == 'S') {
				currOp = CigarOperator.S;
				extLen = 1;
			} else if (aln1[i] != '-') {
				if (aln2[i] != '-') {
					currOp = CigarOperator.M;
				} else {
					currOp = CigarOperator.D;
				}
				extLen = 1;
			} else if (aln2[i] != '-') {
				currOp = CigarOperator.I;
				extLen = 1;
			} else {
				currOp = CigarOperator.N;
				extLen =1;
			}
			
			if (currOp != cigOp) {
				if (cigOp != null) {
					ce.add(new CigarElement(cigLen,cigOp));
				}
				cigLen = 0;
				cigOp = currOp;
			}
			cigLen += extLen;
		}
		
		ce.add(new CigarElement(cigLen,cigOp));
		
		//Hard trims at the end are included so add them back. 
		if (cigar.get(cigar.size()-1).getOperator() == CigarOperator.H) {
			ce.add(cigar.get(cigar.size()-1));
		}
		
		
		if (isReverse) {
			ce = this.reverseCigar(ce);
		}
		
		return new Cigar(ce);
		
	}
	
	
	public boolean isPoorQuality(String[] alignments, boolean isReversed, int trimPoint) {
		char[] aln1 = alignments[0].toCharArray();
		char[] aln2 = alignments[1].toCharArray();
		
		int totalMM = 0;
		
		for (int i=0; i<trimPoint; i++) {
			if (aln1[i] == '-' || aln2[i] == '-' || aln1[i] == 'S' || aln2[i] == 'S') {
				
			} else if (aln1[i] != aln2[i]) {
				if (this.rnaEditing) {
					if (aln1[i] == 'A' && aln2[i] == 'G') {
//					if ((isReversed && aln1[i] == 'T' && aln2[i] == 'C') ||
//						 !isReversed && aln1[i] == 'A' && aln2[i] == 'G') {
//						//OK
					} else {
						totalMM++;
					}
				} else {
					totalMM++;
				}
			} 
		}
		
		if (totalMM > this.maxMM) {
			return true;
		} else {
			return false;
		}
		
		
	}
	
	public int identifyTrimPoint(String[] alignments, boolean isReversed) {
		char[] aln1 = alignments[0].toCharArray();
		char[] aln2 = alignments[1].toCharArray();
		
		int[] scores = new int[aln1.length];
		int lastScore = 0;
		
		for (int i=0; i<scores.length; i++) {
			if (aln1[i] == '-' || aln2[i] == '-' || aln1[i] == 'S' || aln2[i] == 'S') {
				
			} else if (aln1[i] != aln2[i]) {
				if (this.rnaEditing) {
					if (aln1[i] == 'A' && aln2[i] == 'G') {
//					if ((isReversed && aln1[i] == 'T' && aln2[i] == 'C') ||
//						 !isReversed && aln1[i] == 'A' && aln2[i] == 'G') {
						lastScore = lastScore + this.mScore;
					} else {
						lastScore = lastScore - this.mmScore;
					}
				} else {
					lastScore = lastScore - this.mmScore;
				}
			} else {
				lastScore = lastScore + this.mScore;
			}
			scores[i] = lastScore;
		}
		
		return findMaxScore(scores);
	}
	
	public int findMaxScore(int[] scores) {
		int max = 0;
		int pos = 0;
		
		for (int i=0; i<scores.length; i++) {
			if (scores[i] > max) {
				max = scores[i];
				pos = i;
			}
		}
		
		return pos;
	}
	
	
	public String[] createAlignmentStrings(List<CigarElement> cigar, String refSeq, String obsSeq, int totalReads) {
		//Build an alignment from the cigar string
		StringBuilder aln1 = new StringBuilder("");
		StringBuilder aln2 = new StringBuilder("");
		int pos1 = 0;
		int pos2 = 0;
		
		for (int i=0;i<cigar.size();i++) {
		//for (CigarElement ce: cigar) {
			CigarElement ce = cigar.get(i);
			int cel = ce.getLength();

			switch(ce.getOperator()) {
				case M:
					aln1.append(refSeq.substring(pos1, pos1+cel));
					aln2.append(obsSeq.substring(pos2, pos2+cel));
					pos1 += cel;
					pos2 += cel;
					break;
				case N:
					aln1.append(this.createString('-', cel));
					aln2.append(this.createString('-', cel));
					pos1 += cel;
					break;
				case S:
					aln1.append(this.createString('S', cel));
					aln2.append(this.createString('S', cel));
					pos2 += cel;
					break;
				case I:
					aln1.append(this.createString('-', cel));
					aln2.append(obsSeq.substring(pos2, pos2+cel));
					pos2 += cel;
					break;
				case D:
					if (i < cigar.size()-1) { 
						aln1.append(refSeq.substring(pos1, pos1+cel));
						aln2.append(this.createString('-', cel));
						pos1 += cel;
					}
					break;
				case H:
					break;
				default:
					System.out.println(String.format("Don't know how to handle this CIGAR tag: %s. Record %d",ce.getOperator().toString(),totalReads));
					break;
			}
		}
		
		return new String[]{aln1.toString(),aln2.toString()};
		
	}
	
	public String revComp(String original) {
		char[] origChar = original.toCharArray();
		StringBuilder converted = new StringBuilder("");
		for (int i=0; i<original.length();i++) {
			converted.append(AlignmentEndTrimmer.revCompHash.get(origChar[i]));
		}
		
		return converted.reverse().toString();
		
	}
	
	public String createString(char character, int count) {
		StringBuilder sb = new StringBuilder("");
		for (int i=0; i<count; i++) {
			sb.append(character);
		}
		
		return sb.toString();
	}
	

	public void parseFasta() {
		HashMap<String,String> refHash = new HashMap<String,String>();
		try {
			
			BufferedReader br = null;
			if (this.refFile.getName().endsWith(".gz")) {
				br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(this.refFile))));
			} else {
				br = new BufferedReader(new FileReader(this.refFile));
			}
			
			
			Pattern p = Pattern.compile("^>(\\w+)");
			String line = null;
			String chrom = null;
			StringBuilder sequence = new StringBuilder("");
			while (( line = br.readLine()) != null ) {
				Matcher m = p.matcher(line);
				if (m.matches()) {
					if (chrom != null) {
						refHash.put(chrom, sequence.toString());
					}
					chrom = m.group(1);
					sequence = new StringBuilder("");
				} else {
					sequence.append(line.toUpperCase());
				}
			}
			br.close();
			this.refHash = refHash;
		} catch (IOException ioex) {
			System.out.println("Error processing reference file: " + ioex.getMessage());
			ioex.printStackTrace();
			System.exit(1);
		}
	}
	
	public List<CigarElement> reverseCigar(List<CigarElement> original) {
		List<CigarElement> newList = new ArrayList<CigarElement>();
		for (int i=original.size()-1;i>=0;i--) {
			newList.add(original.get(i));
		}
		return newList;
	}
	
	
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Alignment End Trimmer: April 2014                     **\n" +
				"**************************************************************************************\n" +
				"This application can be used to trim alignments according to the density of mismatches.\n"+
				"Each base of the alignment is compared to the reference sequence from the start of the\n" +
				"alignment to the end.  If the bases match, the score is increased by -m. If the bases\n" +
				"don't match, the score is decreased by -n.  The alignment position with the highest \n" +
				"score is used as the new alignment end point. The cigar string, alignment position,\n" +
				"mpos and flags are all updated to reflect trimming. \n"+
				
				"\nNotes:\n"+
				"1) Insertions, deletions and skips are currently not counted as matches or mismatches\n" +
				

				"\nRequired:\n"+
				"-i Path to the orignal alignment, sam/bam/sam.gz OK.\n"+
				"-r Path to the reference sequence, gzipped OK.\n" +
				"-o Name of the trimmed alignment output.  Output is bam and bai.\n" +
				
				
				"\nOptional:\n" +
				"-m Score of match. Default 1\n" +
				"-n Score of mismatch. Default 2\n" +
				"-v Verbose output.  This will write out detailed information for every trimmed read.\n" +
				"    It is suggested to use this option only on small test files.\n" +
				"-l Min length.  If the trimmed length is less than this value, the read is switched\n" +
				"    to unaligned. Default 10bp\n" +
				"-e Turn on RNA Editing mode.  A>G (forward reads) and T>C (reverse reads) are considered\n" +
				"   matches.\n" +
				"-s Turn on mismatch scoring mode. Reads with more than -x mismatches are dropped. If \n" +
				"   RNA Editing mode is on, A>G (forward reads) and T>C (reverse reads) are considered \n" +
				"   matches.\n" +
				"-x Max number of mismatches allowed in max scoring mode. Default 0\n" +

				"\nExamples: \n" +
				"1) java -Xmx4G -jar /path/to/AlignmentEndTrimmer -i 1000X1.bam -o 100X1.trim.bam\n" +
				"           -r /path/to/hg19.fasta\n" +
				"2) java -Xmx4G -jar /path/to/AlignmentEndTrimmer -i 1000X1.bam -o 100X1.trim.bam\n" +
				"           -r /path/to/hg19.fasta -m 0.5 -n 3\n" +
				"3) java -Xmx4G -jar /path/to/AlignmentEndTrimmer -i 1000X1.test.bam \n" +
				"           -o 100X1.test.trim.bam -r /path/to/hg19.fasta -v\n" +

				"**************************************************************************************\n");

	}

}
