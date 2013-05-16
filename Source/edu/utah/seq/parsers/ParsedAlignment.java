package edu.utah.seq.parsers;

import java.io.DataInputStream;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.data.Point;
import edu.utah.seq.data.sam.SamAlignment;
import net.sf.samtools.SAMRecord;
import util.bio.seq.Seq;
import util.gen.Misc;

public class ParsedAlignment implements Comparable<ParsedAlignment>{
	private String readID;
	private int position;
	private String sequence;
	private char[] bases;
	private String baseQualities;
	private int[] baseScores;
	private String[] baseObservations;
	private boolean passesParsing = true;
	private boolean containsIndels = false;
	private BaseObservation[] bo;
	private SAMRecord sam;
	private NovoalignBisulfiteParser parser;
	private String cigar;
	private int minScore;
	
	char[] refSeqChar = null;
	char[] sequenceChar = null;

	public static final Pattern CIGAR_LETTERS = Pattern.compile("[IDSHM]");
	public static final Pattern SPACE = Pattern.compile(" ");
	public static final Pattern CIGAR_SUB = Pattern.compile("(\\d+)([IDNMS])");
	public static final Pattern CIGAR_BAD_CHARACTERS = Pattern.compile(".*[^\\dMNID].*");
	public static final Pattern INSERTION = Pattern.compile("(\\d+)\\+(\\w+)");
	public static final Pattern DELETION = Pattern.compile("(\\d+)-(\\w)");
	public static final Pattern INDEL = Pattern.compile("[-\\+]");
	public static final Pattern CONVERTED_BASE = Pattern.compile("(\\d+)[CG]#[TA]");

	/**Used by the AllelicMethylationDetector*/
	public ParsedAlignment (SAMRecord sam, NovoalignBisulfiteParser parser) {

		try {
			this.sam = sam;
			readID = sam.getReadName();
			position = sam.getAlignmentStart()-1;
			sequence = sam.getReadString();
			baseQualities = sam.getBaseQualityString();
			cigar = sam.getCigarString();
			this.parser = parser;
			minScore = parser.getMinimumBaseScore();
			
			//trim sam info for hard and soft masking
			trimMaskingOfReadToFitAlignment();

			//set isGA based on strand
			boolean isGA; 
			//is it part of a paired alignment? the need to set it to the ori of the 1st read
			if (sam.getReadPairedFlag()){
				if (sam.getFirstOfPairFlag() && sam.getReadNegativeStrandFlag() == false) isGA = false;
				else if (sam.getSecondOfPairFlag() && sam.getReadNegativeStrandFlag()) isGA = false;
				else isGA = true;
			}
			else if (sam.getReadNegativeStrandFlag()) isGA = true;
			else isGA = false;
			String refSeq = new String (parser.getGenomicSequence().substring(position, position + SamAlignment.countLengthOfAlignment(cigar)));

			String baseCalls = makeBaseCalls (sequence, refSeq.toUpperCase(), cigar, isGA);
			containsIndels = INDEL.matcher(baseCalls).find();

			baseObservations = SPACE.split(baseCalls);

			scan();
		} catch (NumberFormatException e){
			System.err.println("Found malformed SAM CIGAR field, skipping ->\n"+sam.getSAMString());
			passesParsing = false;
		}
	}



	/**Used by the NBP*/
	public ParsedAlignment(DataInputStream dis, boolean samFormat, NovoalignBisulfiteParser parser) throws Exception{
		readID = dis.readUTF();
		position = dis.readInt();
		sequence = dis.readUTF().toUpperCase();
		baseQualities = dis.readUTF();
		this.parser = parser;
		minScore = parser.getMinimumBaseScore();
		String baseCalls ;
		if (samFormat){
			String x = dis.readUTF();
			boolean isGA = true;
			if (x.startsWith("f")) isGA = false;
			cigar = x.substring(1);
			//watch for reads that overlap end
			String seq = parser.getGenomicSequence(); 
			int end = position + SamAlignment.countLengthOfAlignment(cigar);
			if (end >= seq.length() || position < 0) {
				passesParsing = false;
				return;
			}
			String refSeq = new String (seq.substring(position, end));		
			baseCalls = makeBaseCalls (sequence, refSeq.toUpperCase(), cigar, isGA);
			
		}
		else baseCalls = dis.readUTF();
		
		//any indels
		containsIndels = INDEL.matcher(baseCalls).find();

		baseObservations = SPACE.split(baseCalls);

		if (samFormat) scanSAM();
		else scan();
	}
	
	private void scanSAM(){
			ArrayList<BaseObservation> boALX = null;
			if (baseObservations[0].equals("CT")) boALX = processLineCTObject(position, sequenceChar, baseScores, baseObservations);
			else if (baseObservations[0].equals("GA")) boALX = processLineGAObject(position, sequenceChar, baseScores, baseObservations);
			bo = new BaseObservation[boALX.size()];
			boALX.toArray(bo);
		
	}

	private void scan(){
		//scan for insertions and deletions
		boolean foundInsertion = false;
		boolean foundDoubleInsertion = false;
		boolean foundDeletions = false;
		for (int x=0; x< baseObservations.length; x++){
			//insertion?
			Matcher mat = INSERTION.matcher(baseObservations[x]);
			if (mat.matches()){
				if (foundInsertion){
					foundDoubleInsertion = true;
					continue;
				}
				int val = Integer.parseInt(mat.group(1)) -1;
				int n = mat.group(2).length();
				String l = sequence.substring(0, val);
				String r = sequence.substring(val+n);
				sequence = l + r;
				l = baseQualities.substring(0, val);
				r = baseQualities.substring(val+n);
				baseQualities = l + r;
				foundInsertion = true;
			}
			//deletion?
			Matcher matD = DELETION.matcher(baseObservations[x]);
			if (matD.matches()){
				foundDeletions = true;
				int val = Integer.parseInt(matD.group(1)) -1;
				String base = matD.group(2);
				sequence = sequence.substring(0, val) + base + sequence.substring(val);
				baseQualities = baseQualities.substring(0, val) + "C" + baseQualities.substring(val);
			}
		}
		//continue?
		if ((foundDeletions && foundInsertion) || foundDoubleInsertion) passesParsing = false;

		else {
			bases = sequence.toCharArray();
			baseScores = Seq.convertScores(baseQualities);
			ArrayList<BaseObservation> boALX = null;
			if (baseObservations[0].equals("CT")) boALX = processLineCTObject(position, bases, baseScores, baseObservations);
			else if (baseObservations[0].equals("GA")) boALX = processLineGAObject(position, bases, baseScores, baseObservations);
			bo = new BaseObservation[boALX.size()];
			boALX.toArray(bo);
		}
	}
	
	/**For reads with H or S masking, strips off these notations from the CIGAR and for S trims the read sequence and base qualities. H's are already trimmed.*/
	public void trimMaskingOfReadToFitAlignment(){	
		//remove hard masking references in cigar
		if (cigar.contains("H")){
			StringBuilder sb = new StringBuilder();
			//(\\d+)([IDH])
			Matcher mat = CIGAR_SUB.matcher(cigar);
			while (mat.find()){
				if (mat.group(2).equals("H") == false) sb.append(mat.group(0));
			}
			cigar = sb.toString();
		}

		//look for soft masking
		if (cigar.contains("S")){
			//at beginning? ^(\\d+)[SH].+
			Matcher mat = SamAlignment.CIGAR_STARTING_MASK.matcher(cigar);
			if (mat.matches()){
				int basesToTrim = Integer.parseInt(mat.group(1));
				sequence = sequence.substring(basesToTrim);
				baseQualities = baseQualities.substring(basesToTrim);
				cigar = cigar.substring(mat.group(1).length()+1);
			}

			//at end?  .+(\\d+)[SH]$
			mat = SamAlignment.CIGAR_STOP_MASKED.matcher(cigar);
			if (mat.matches()){
				int basesToTrim = Integer.parseInt(mat.group(1));
				int stopIndex = sequence.length() - basesToTrim;
				sequence = sequence.substring(0, stopIndex);
				baseQualities = baseQualities.substring(0, stopIndex);
				stopIndex = cigar.length() - mat.group(1).length() - 1;
				cigar = cigar.substring(0, stopIndex);
			}
		}
	}

	/**Generates a native novoalign base call string for bisulfite alignments from SAM junk.*/
	public  String makeBaseCalls (String samSequence, String refSeq, String cigar, boolean isGA){
		StringBuilder sb = new StringBuilder();
		try {
			refSeqChar = refSeq.toCharArray();
			sequenceChar = samSequence.toCharArray();
			baseScores = Seq.convertScores(baseQualities);


			if (isGA) sb.append("GA");
			else sb.append("CT");

			//check to see if cigar contains any unsupported characters, note call SamAlignment.trimMaskingOfReadToFitAlignment to strip off S and H masking
			Matcher mat = CIGAR_BAD_CHARACTERS.matcher(cigar);
			if (mat.matches()) {
				System.err.println("\nUnsupported cigar string "+cigar);
				System.exit(0);
			}

			//any deletions (add to seq) or insertions (subtract seq)?
			int pos =0;
			mat = CIGAR_SUB.matcher(cigar);
			while (mat.find()){
				String call = mat.group(2);
				int numberBases = Integer.parseInt(mat.group(1));
				//System.out.println("\t"+numberBases+" "+call);

				//an insertion, add blanks to refSeq
				if (call.equals("I")) {
					//System.out.println("Insertion of "+numberBases+" bases in read at pos "+pos);
					//add blanks to refSeq char array
					char[] modSeqChar = new char[refSeqChar.length+ numberBases];
					//copy first bit
					System.arraycopy(refSeqChar, 0, modSeqChar, 0, pos);
					//add deletion
					for (int i=0; i< numberBases; i++) modSeqChar[pos+i] = 'I';
					//copy last bit
					System.arraycopy(refSeqChar, pos, modSeqChar, pos+numberBases, refSeqChar.length - pos);
					refSeqChar = modSeqChar;
				}
				//a deletion or gap due to splicing
				else if (call.equals("D") || call.equals("N")) {
					//add blanks to sequence char array
					char[] modSeqChar = new char[sequenceChar.length+ numberBases];
					int[] modBaseScores = new int[baseScores.length+ numberBases];
					
					//copy first bit
					System.arraycopy(sequenceChar, 0, modSeqChar, 0, pos);
					System.arraycopy(baseScores, 0, modBaseScores, 0, pos);

					//add deletion
					char callChar = call.charAt(0);
					for (int i=0; i< numberBases; i++) {
						modSeqChar[pos+i] = callChar;
						modBaseScores[pos+i] = 35;
					}
					
					//copy last bit
					System.arraycopy(sequenceChar, pos, modSeqChar, pos+numberBases, sequenceChar.length - pos);
					System.arraycopy(baseScores, pos, modBaseScores, pos+numberBases, baseScores.length - pos);
					sequenceChar = modSeqChar;
					baseScores = modBaseScores;
				}
				//increment counter
				pos += numberBases;

			}


			//no need to look for hard and soft masking, these are already stripped
			int start = 0; 
			int stop = sequenceChar.length;
			
			pos = start;
			for (int i=start; i< stop; i++){
				pos++;

				//a mis match?
				if (refSeqChar[i] != sequenceChar[i]){
					//N in reference
					if (refSeqChar[i] == 'N'){
						if ( (sequenceChar[i] == 'D' || sequenceChar[i] == 'I') == false) continue;
					}
					
					//N in sequence due to splicing
					if (sequenceChar[i] == 'N') continue;
					
					sb.append(" ");
					
					//is it a G in ref -> A in sequence?
					if (isGA && refSeqChar[i] == 'G' && sequenceChar[i] == 'A'){
						sb.append(pos);
						sb.append("G#A");
					}
					//is it a C in ref -> T in sequence
					else if (isGA == false && refSeqChar[i] == 'C' && sequenceChar[i] == 'T'){
						sb.append(pos);
						sb.append("C#T");
					}
					//is it a deletion?
					else if (sequenceChar[i] == 'D'){
						sb.append(pos);
						sb.append("-");
						sb.append(refSeqChar[i]);
					}
					//is it an insertion?
					else if (refSeqChar[i] == 'I'){
						sb.append(pos);
						sb.append("+");
						sb.append(sequenceChar[i]);
						pos--;
					}
					//snp
					else {
						sb.append(pos);
						sb.append(refSeqChar[i]);
						sb.append(">");
						sb.append(sequenceChar[i]);
					}

				}

			}

		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("trimmed cigar "+cigar);
			System.out.println("samSequence "+samSequence);
			System.out.println("refSeq "+refSeq);
			System.out.println("isGA "+isGA);

			System.exit(0);
		}



		return sb.toString();
	}

	public ArrayList<BaseObservation> processLineCTObject(int position, char[] bases, int[] scores, String[] baseObs){
		ArrayList<BaseObservation> boAL = new ArrayList<BaseObservation>();
		int start = -1;
		String genSeq = "";
		String bedLine = "";
		try {
		//look for non converted C's
		for (int i=0; i< bases.length; i++){
			if (bases[i] == 'C'){
				//check score?
				if (containsIndels == false) {
					if (scores[i] < minScore) {
						parser.incrementBPFailQual();
						continue;
					}
				}
				parser.incrementBPPassQual();
				start = i + position;
				//watch out for out of bounds sequence due to partial matches to sequence termini
				if (start < 2 || start > parser.getGenomicSequenceLengthMinus3()) continue;
				//fetch genomic sequence, use new String to avoid reference
				genSeq = new String(parser.getGenomicSequence().substring(start-2,start+3).toUpperCase());
				//check for center C, often these are either snps or mis calls
				if (genSeq.charAt(2) != 'C') continue;
				bedLine = parser.getChromosome()+ "\t"+ start + "\t" + (start+1)+"\t"+genSeq+"\t"+scores[i]+"\t"+parser.getStrand();
				//save Point
				Point point = new Point(start,1f);
				BaseObservation bo = new BaseObservation(start, bedLine, point, genSeq, scores[i], false); 
				boAL.add(bo);
			}
		}
		//look for converted C's: 'CT 9C#T 35C#T 57C#T', skip first token
		for (int i=1; i< baseObs.length; i++){
			Matcher mat = CONVERTED_BASE.matcher(baseObs[i]);
			if (mat.matches()){
				int relPos = Integer.parseInt(mat.group(1)) - 1;
				//check score?
				if (containsIndels == false) {
					if (scores[relPos] < minScore) {
						parser.incrementBPFailQual();
						continue;
					}
				}
				parser.incrementBPPassQual();
				
				start = relPos + position;
				//watch out for out of bounds sequence due to partial matches to sequence termini
				if (start < 2 || start > parser.getGenomicSequenceLengthMinus3()) continue;
				genSeq = new String(parser.getGenomicSequence().substring(start-2,start+3).toUpperCase());
				//check for center C, often these are either snps or mis calls
				if (genSeq.charAt(2) != 'C') continue;
				bedLine = parser.getChromosome()+ "\t"+ start + "\t" + (start+1)+"\t"+genSeq+"\t"+scores[relPos]+"\t"+parser.getStrand();
				Point point = new Point(start,1f);
				BaseObservation bo = new BaseObservation(start, bedLine, point, genSeq, scores[relPos], true); 
				boAL.add(bo);
			}
		}
		} catch (Exception e){
			e.printStackTrace();
			System.err.println(parser.getChromosome()+"\t"+parser.getStrand()+ "\t"+ position);
			for (char b : bases) System.err.print(" " +b);
			for (int b : scores) System.err.print(" " +b);
			for (String b : baseObs) System.err.print(" " +b);
			System.err.println("\tStart "+start);
			System.err.println("\tGenSeq "+genSeq);
			System.err.println("\tBedLine "+bedLine);
			System.exit(1);
		}
		return boAL;
	}

	public ArrayList<BaseObservation> processLineGAObject(int position, char[] bases, int[] scores, String[] baseObs){
		ArrayList<BaseObservation> boAL = new ArrayList<BaseObservation>();
		//look for non converted C's
		for (int i=0; i< bases.length; i++){
			if (bases[i] == 'G'){
				//check score?
				if (containsIndels == false) {
					if (scores[i] < minScore) {
						parser.incrementBPFailQual();
						continue;
					}
				}
				parser.incrementBPPassQual();

				int start = i + position;
				//watch out for out of bounds sequence due to partial matches to sequence termini
				if (start < 2 || start > parser.getGenomicSequenceLengthMinus3()) continue;
				//fetch genomic sequence, use new String to avoid reference
				String genSeq = new String(parser.getGenomicSequence().substring(start-2,start+3).toUpperCase());
				genSeq = Seq.reverseComplementDNA(genSeq);
				//check for center C, often these are either snps or mis calls
				if (genSeq.charAt(2) != 'C') continue;
				String bedLine = parser.getChromosome()+ "\t"+ start + "\t" + (start+1)+"\t"+genSeq+"\t"+scores[i]+"\t"+parser.getStrand();
				//save Point
				Point point = new Point(start,1f);
				BaseObservation bo = new BaseObservation(start, bedLine, point, genSeq, scores[i], false); 
				boAL.add(bo);
			}
		}
		//look for non converted G's: 'GA 4G#A 5G#A 29G#A 30G#A 31G#A 40G#A', skip first token
		for (int i=1; i< baseObs.length; i++){
			Matcher mat = CONVERTED_BASE.matcher(baseObs[i]);
			if (mat.matches()){
				int relPos = Integer.parseInt(mat.group(1)) - 1;
				//check score?
				if (containsIndels == false) {
					if (scores[relPos] < minScore) {
						parser.incrementBPFailQual();
						continue;
					}
				}
				parser.incrementBPPassQual();

				int start = relPos + position;
				//watch out for out of bounds sequence due to partial matches to sequence termini
				if (start < 2 || start > parser.getGenomicSequenceLengthMinus3()) continue;
				String genSeq = new String(parser.getGenomicSequence().substring(start-2,start+3).toUpperCase());
				genSeq = Seq.reverseComplementDNA(genSeq);
				//check for center C, often these are either snps or mis calls
				if (genSeq.charAt(2) != 'C') continue;
				String bedLine = parser.getChromosome()+ "\t"+ start + "\t" + (start+1)+"\t"+genSeq+"\t"+scores[relPos]+"\t"+parser.getStrand();
				Point point = new Point(start,1f);
				BaseObservation bo = new BaseObservation(start, bedLine, point, genSeq, scores[relPos], true); 
				boAL.add(bo);
			}
		}

		return boAL;
	}



	public int compareTo(ParsedAlignment other) {
		return readID.compareTo(other.readID);
	}

	public String getReadID() {
		return readID;
	}

	public int getPosition() {
		return position;
	}

	public boolean isPassesParsing() {
		return passesParsing;
	}

	public BaseObservation[] getBo() {
		return bo;
	}

	public char[] getBases() {
		return bases;
	}
}