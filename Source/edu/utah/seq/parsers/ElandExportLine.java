package edu.utah.seq.parsers;
import java.util.HashMap;
import java.util.regex.*;
import edu.utah.seq.data.*;
import util.gen.*;
import util.bio.seq.*;

public class ElandExportLine {

	//fields
	private String chromosome;
	private int position;
	private float score;
	private String strand;
	private static final Pattern splitPat = Pattern.compile("\\s");
	private static final Pattern underScore = Pattern.compile("_");
	private int readLength;
	private String sequence;
	private String seqQuality;
	public static final HashMap<String, Integer> asciiValue = ElandExportLine.asciScoreSolexa();


	/**Parses an Eland Extended output line, will shift the stranded position by the offSet.
	 * HWI-EAS240	FC2087UAAXX	8	1	237	408			TACACATGAATTCAACTTAAATTCCTTGTTAAAATT	ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZVUUVVU	chr6.fasta		147077758	F	36	73						Y
	 * HWI-EAS240	FC2087UAAXX	8	1	311	379			TACTTGGTTAAGTAGGATGTTATTCTGCTTCTACAC	ZZZZZZZZZZZZZZZZZZZZZZZZVZRZZZUVLLLL	chr3.fasta		154412997	R	36	65						Y
	 * HWI-EAS240	FC2087UAAXX	8	1	219	459			TACAACATGTACAAGCCTAAATCCTTTTAGCCAGAG	ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZJUOSOUQ	NM											Y
	 * HWI-EAS240	FC2087UAAXX	8	1	582	773			TCATGTTAAAGTGATTGTGATGTTTGAAAACCAATG	ZZZZZZZZZZZZZZZZZZZZZVZZZZJZZZSUSVVG	chr6.fasta		123151134	R	36	54						Y
	 * HWI-EAS240	FC2087UAAXX	8	1	820	387			TAAATTGTAACATAAAATTCTTATGAAATTACCTCA	ZKZZZZZZZZZZZZZZZZZZZZZZIZOZZZOODOOH	chr16.fasta		58060756	R	36	39						N  
	 *
	 *Note, Reverse match annotation are relative to the reverse (seq is reverse, quality reversed, base changes are base changes needed to match consensus, 
	 *     the sequence given is the oligo sequence not the reference, position is in forward space though.)
	 * Only lines that have passed the qc flag and have an alignment will be parsed so check the parsed() before using.
	 */
	public ElandExportLine(String line, int offSet, float minScore){
		String[] tokens = splitPat.split(line);
		if (tokens.length < 15) System.err.println("Error: line does not contain enough columns -> "+line);
		else {
			//trim tokens
			Misc.trim(tokens);
			//worth parsing? pased QC and was aligned
			if (tokens[12].equals("") || (tokens.length == 22 && tokens[21].equals("N"))) return;
			//parse score
			score = Float.parseFloat(tokens[15]);
			if (score < minScore) return;
			//parse the chromosome
			int index = tokens[10].indexOf('.');
			//anything found?
			try {
				if (index != -1) {
					//parse readlength
					readLength = tokens[8].length();
					//parse chromosome
					//anything in the secondary column for a multi fasta?
					if (tokens[11].length() !=0){
						//process splice chrX_71333479_71333942
						String[] splice = underScore.split(tokens[11]);
						chromosome = splice[0];
						//parse position and subtract one to put into interbase coordinates
						position = Integer.parseInt(splice[1]) -1; 
					}
					else {
						chromosome = tokens[10].substring(0, index);
						//parse position and subtract one to put into interbase coordinates
						position = Integer.parseInt(tokens[12]) -1;
						//center
						int halfLength = (int)(Math.round(((double)readLength)/2.0));
						position +=halfLength;
					}
					//strand and offset
					if (tokens[13].equals("F")) {
						strand = "+";
						position += offSet;
					}
					else if (tokens[13].equals("R")) {
						strand = "-";
						position -= offSet;
						if (position < 0) position = 0;
					}
					else {
						System.err.println("Error: cannot parse strand orientation! See -> "+line);
						chromosome = null;
					}
				}
				else System.err.println("Error: missing '.' (e.g. chrX.fasta) from chromosome text cannot parse! See -> "+line);
			} catch (Exception e){
				e.printStackTrace();
				System.out.println(line);
				System.exit(0);
			}
		}
	}

	/**Parses an Eland Extended output line, will shift the stranded position by the offSet.
	 * HWI-EAS240	FC2087UAAXX	8	1	237	408			TACACATGAATTCAACTTAAATTCCTTGTTAAAATT	ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZVUUVVU	chr6.fasta		147077758	F	36	73						Y
	 * HWI-EAS240	FC2087UAAXX	8	1	311	379			TACTTGGTTAAGTAGGATGTTATTCTGCTTCTACAC	ZZZZZZZZZZZZZZZZZZZZZZZZVZRZZZUVLLLL	chr3.fasta		154412997	R	36	65						Y
	 * HWI-EAS240	FC2087UAAXX	8	1	219	459			TACAACATGTACAAGCCTAAATCCTTTTAGCCAGAG	ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZJUOSOUQ	NM											Y
	 * HWI-EAS240	FC2087UAAXX	8	1	582	773			TCATGTTAAAGTGATTGTGATGTTTGAAAACCAATG	ZZZZZZZZZZZZZZZZZZZZZVZZZZJZZZSUSVVG	chr6.fasta		123151134	R	36	54						Y
	 * HWI-EAS240	FC2087UAAXX	8	1	820	387			TAAATTGTAACATAAAATTCTTATGAAATTACCTCA	ZKZZZZZZZZZZZZZZZZZZZZZZIZOZZZOODOOH	chr16.fasta		58060756	R	36	39						N  
	 *
	 *Note, Reverse match annotation are relative to the reverse (seq is reverse, quality reversed, base changes are base changes needed to match consensus, 
	 *     the sequence given is the oligo sequence not the reference, position is in forward space though.)
	 * Only lines that have passed the qc flag and have an alignment will be parsed so check the parsed() before using.
	 */
	public ElandExportLine(String line, int offSet, float minScore, int chromColumn){
		String[] tokens = splitPat.split(line);
		if (tokens.length < 15) System.err.println("Error: line does not contain enough columns -> "+line);
		else {
			//trim tokens
			Misc.trim(tokens);
			//worth parsing? pased QC and was aligned
			if (tokens[12].equals("") || (tokens.length == 22 && tokens[21].equals("N"))) return;
			//parse score
			score = Float.parseFloat(tokens[15]);
			if (score < minScore) return;

			//aligned?
			if (tokens[chromColumn].length() == 0) return;

			//parse readlength
			readLength = tokens[8].length();

			//parse chromosome
			chromosome = tokens[11];

			//parse position and subtract one to put into interbase coordinates
			position = Integer.parseInt(tokens[12]) -1;

			//center
			int halfLength = (int)(Math.round(((double)readLength)/2.0));
			position +=halfLength;

			//strand and offset
			if (tokens[13].equals("F")) {
				strand = "+";
				position += offSet;
			}
			else if (tokens[13].equals("R")) {
				strand = "-";
				position -= offSet;
				if (position < 0) position = 0;
			}
			else {
				System.err.println("Error: cannot parse strand orientation! See -> "+line);
				chromosome = null;
			}
		}
	}

	/**Stripped constructor/ method for ElandSequenceParser*/
	public ElandExportLine(String line, int[][] gatc){
		String[] tokens = splitPat.split(line);
		//for (int i=0; i< tokens.length; i++) System.out.println("\t"+i+" "+tokens[i]);
		position = Integer.parseInt(tokens[12]) -1;
		if (position < 0) {
			System.out.println("\tNegative alignment position, circular genome? additional non reference sequence on ends? Skipping \n\t"+line);
			return;
		}
		//parse read length
		readLength = tokens[8].length();
		//strand
		if (tokens[13].equals("F")) {
			strand = "+";
			sequence = tokens[8];
			seqQuality = tokens[9];
		}
		else if (tokens[13].equals("R")) {
			strand = "-";
			//reverse the sequence, and quality values
			sequence = Seq.reverseComplementDNA(tokens[8]);
			seqQuality = Misc.reverse(tokens[9]);
		}
		else {
			System.err.println("Error: cannot parse strand orientation! See -> "+line);
			chromosome = null;
			return;
		}
		//parse scores
		int[] qscores = convertScores(seqQuality);

		//add to base tracks
		for (int i=0; i< readLength; i++){
			char base = sequence.charAt(i);
			if (base == 'G') gatc[0][position+i] += qscores[i];
			else if (base == 'A') gatc[1][position+i] += qscores[i];
			else if (base == 'T') gatc[2][position+i] += qscores[i];
			else if (base == 'C') gatc[3][position+i] += qscores[i];
		}
	}

	/**Converts the ascii quality scores to numeric scores.*/
	public static int[] convertScores(String seqQual){
		int[] scores = new int[seqQual.length()];
		for (int i=0; i< seqQual.length(); i++){
			String sub = seqQual.substring(i, i+1);
			Integer val = asciiValue.get(sub);
			if (val != null) scores[i] = val.intValue();
			//else System.out.println("\tError converting seq quality character -> "+sub+" from "+seqQual);
		}
		return scores;
	}

	/**Returns a map of asci text character to it's associated Solexa base quality score.*/
	public static HashMap<String,Integer> asciScoreSolexa(){
		String[] asci = {"@","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","[","\\","]","^","_"};
		HashMap<String,Integer> map = new HashMap<String,Integer>();
		for (int i=0; i< asci.length; i++) map.put(asci[i], new Integer(i));
		return map;
	}

	/**
	 * For parsing a stand alone Eland run.
	 * Score set to 1
	 * /home/solexa/Hg17SpikeIns/reads.txt-1	ttccactgagcctgctccagaactt	U2	0	0	1	chr1	53226799	F	..	23A	25C

	 */
	public ElandExportLine(String line, int offSet, boolean parseIt){
		String[] tokens = splitPat.split(line);
		if (tokens.length >= 9){
			chromosome = tokens[6];
			//parse position and subtract one to put into interbase coordinates
			position = Integer.parseInt(tokens[7]) -1;
			//center
			int halfLength = (int)(Math.round(((double)tokens[1].length())/2.0));
			position +=halfLength;
			//strand
			if (tokens[8].equals("F")) {
				strand = "+";
				position += offSet;
			}
			else if (tokens[8].equals("R")) {
				strand = "-";
				position -= offSet;
				if (position < 0) position = 0;
			}
			else {
				System.err.println("Error: cannot parse strand orientation! See -> "+line);
				chromosome = null;
			}
			//parse score
			score = 1;
			//parse readlength
			readLength = tokens[1].length();
		}
	}


	public boolean parsed(){ 
		if (chromosome != null) return true;
		return false;
	}

	public String fetchChromStrandName(){
		return chromosome +"_"+ strand;
	}

	public int fetchReadLength(){
		return readLength;
	}

	public Point fetchPoint(){
		if (chromosome == null) return null;
		return new Point(position,score);
	}

	public String getChromosome() {
		return chromosome;
	}
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public int getPosition() {
		return position;
	}
	public void setPosition(int position) {
		this.position = position;
	}
	public float getScore() {
		return score;
	}
	public void setScore(float score) {
		this.score = score;
	}
}
