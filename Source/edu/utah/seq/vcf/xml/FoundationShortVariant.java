package edu.utah.seq.vcf.xml;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import htsjdk.samtools.reference.ReferenceSequence;
import util.bio.seq.Seq;
import util.gen.Misc;
import util.gen.Num;

public class FoundationShortVariant {
	
	//fields
	private LinkedHashMap<String, String> parsedAttributes;
	private String chromosome;
	private int position;  //one based, subtract 1 to match IGB
	private String reference;
	private String alternate;
	private String type; //snv, ins, del, multi
	private String effectedGene;
	private String proteinEffect;
	private String functionalEffect;
	private String status;
	private double af; //allele frequency
	private int dp; //depth
	private FoundationXml2Vcf parser;
	private boolean failedParsing = false;
	
	//for the vcf header
	static String infoEG = "##INFO=<ID=EG,Number=1,Type=String,Description=\"Effected gene\">";
	private static String infoFE = "##INFO=<ID=FE,Number=1,Type=String,Description=\"Functional effect on gene\">";
	static String infoST = "##INFO=<ID=ST,Number=1,Type=String,Description=\"Status of effect\">";
	private static String infoPE = "##INFO=<ID=PE,Number=1,Type=String,Description=\"Protein effect on gene\">";
	private static String infoDP = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">";
	private static String infoAF = "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">";
	
	//for extracting 
	private static Pattern snvPat = Pattern.compile("[-_\\d\\+]+([GATC]+)>([GATC]+)");
	private static Pattern delPat = Pattern.compile("[-_\\d\\+]+DEL([GATC\\d]+)");
	private static Pattern insPat = Pattern.compile("[-_\\d\\+]+INS([GATC\\d]+)");
	private static Pattern multiPat = Pattern.compile("([-_\\d\\+]+)>([GATC\\d]+)");    
	private static Pattern intPat = Pattern.compile("\\d+");
	private static Pattern minusPat = Pattern.compile("(\\d+)-(\\d+)");
	private static Pattern plusPat = Pattern.compile("(\\d+)\\+(\\d+)");
	
	public FoundationShortVariant (LinkedHashMap<String, String> parsedAttributes, FoundationXml2Vcf parser){
		//examples:
		//cds-effect=81C>G, depth=1079, functional-effect=missense, gene=DICER1, percent-reads=5.0, position=chr14:95599715, protein-effect=F27L, status=unknown, strand=-, subclonal=true, transcript=NM_030621
		//cds-effect=4835_4835delT, depth=462, functional-effect=frameshift, gene=ATRX, percent-reads=78.0, position=chrX:76889174, protein-effect=L1613fs*11, status=likely, strand=-, transcript=NM_000489
		this.parsedAttributes = parsedAttributes;
		this.parser = parser;
		parseChromPos();
		parseRefAlt();
		parseReadSupport();
		parseEffect();
		
		//check the reference
		checkReference();
	}

	/** Returns CHROM POS ID REF ALT QUAL FILTER INFO (EG FE ST PE DP AF)*/
	public String toVcf(){
		StringBuilder sb = new StringBuilder();
		//CHROM
		sb.append(chromosome); sb.append("\t");
		//POS
		sb.append(position); sb.append("\t");
		//ID
		sb.append(".\t");
		//REF
		sb.append(reference); sb.append("\t");
		//ALT
		sb.append(alternate); sb.append("\t");
		//QUAL
		sb.append(".\t");
		//FILTER
		if (failedParsing) sb.append("ci\t");
		else sb.append(".\t");
		//INFO (EG FE ST PE DP AF)
		sb.append("EG=");
		sb.append(effectedGene);
		sb.append(";FE=");
		sb.append(functionalEffect);
		sb.append(";ST=");
		sb.append(status);
		sb.append(";PE=");
		sb.append(proteinEffect);
		sb.append(";DP=");
		sb.append(dp);
		sb.append(";AF=");
		sb.append(Num.formatNumber(af, 2));
		
		return sb.toString();
		
	}
	
	public static void appendInfoLines(TreeSet<String> sb){
		sb.add(infoEG); 
		sb.add(infoFE); 
		sb.add(infoST); 
		sb.add(infoPE); 
		sb.add(infoDP); 
		sb.add(infoAF); 
	}
	
	private void parseEffect() {
		effectedGene = parseString("gene");
		functionalEffect = parseString("functional-effect");
		status = parseString("status").toLowerCase();
		proteinEffect = parseString("protein-effect");
	}
	
	private void parseReadSupport() {
		dp = (int)Math.round(parseDouble("depth"));
		double percentReads = parseDouble("percent-reads");
		af = percentReads/100.0;
	}

	private void parseRefAlt() {
		//parse the strand info
		String strand = parseString("strand");
		if (strand.equals("+") == false && strand.equals("-") == false) Misc.printErrAndExit("\nError: failed to parse + or - strand info from "+parsedAttributes);
		//parse effect
		parseEffect(strand);
	}

	private void parseEffect(String strand) {
		String eff = parseString("cds-effect").toUpperCase();
		
		//snv?   394C>T
		Matcher snv = snvPat.matcher(eff);
		if (snv.matches()){
			type = "snv";
			reference = snv.group(1);
			alternate = snv.group(2);
			if (strand.equals("-")) {
				reference = Seq.reverseComplementDNA(reference);
				alternate = Seq.reverseComplementDNA(alternate);
			}
			//check length
			if (reference.length() != alternate.length()) type = "multi";
//if (reference.length()!=1 || alternate.length()!= 1) System.out.println("\tSNV "+chromosome+":"+position+" "+reference+" "+alternate+"\n"+this.parsedAttributes);

			return;
		}
		
		//deletion?   4835_4835delT   660-107_664del112
		Matcher del = delPat.matcher(eff);
		if (del.matches()){
			type = "del";
			String basesDeleted = del.group(1);
			int lengthOfDeletion = basesDeleted.length();
			//is it a number
			Matcher mat = intPat.matcher(basesDeleted);
			boolean number = false;
			if (mat.matches()) {
				number = true;
				lengthOfDeletion = Integer.parseInt(mat.group(0));
			}
			else {
				if (strand.equals("-")) basesDeleted = Seq.reverseComplementDNA(basesDeleted);
				ReferenceSequence p = parser.getFasta().getSubsequenceAt(chromosome, position+1, position+lengthOfDeletion);
				String fasta = new String(p.getBases());
				if (fasta.equals(basesDeleted) == false){
					failedParsing = true;
					System.err.println("\tWARNING: the deleted bases do not match the fasta ('"+fasta+"')?\n\t"+parsedAttributes+"\n\t"+toVcf());
				}
			}

			//get the bases to delete 
			long start = position;
			long end = position + lengthOfDeletion;
			ReferenceSequence rs = parser.getFasta().getSubsequenceAt(chromosome, start, end);
			reference = new String(rs.getBases());
			alternate = reference.substring(0, 1);
//if (number==false) System.out.println("\tDEL "+chromosome+":"+position+" "+reference+" "+alternate+" "+this.parsedAttributes);
//if (reference.contains("GATGATGATGAAGAG")) Misc.printErrAndExit("\t\tCHECK "+chromosome+":"+position+" "+reference+" "+alternate+" "+this.parsedAttributes);
			return;
		}
		
		//insertion?   386_387insGCAGCAGCA  3558_3559ins148
		Matcher ins = insPat.matcher(eff);
		if (ins.matches()){
			type="ins";
			String basesInserted = ins.group(1);
			
			//is it a number?
			Matcher mat = intPat.matcher(basesInserted);
			if (mat.matches()) {
				int len = Integer.parseInt(mat.group(0));
				ReferenceSequence inser = parser.getFasta().getSubsequenceAt(chromosome, position, position+len);
				basesInserted = new String(inser.getBases());			
			}
			else if (strand.equals("-")) basesInserted = Seq.reverseComplementDNA(basesInserted);
			
			//get the base at position 
			ReferenceSequence rs = parser.getFasta().getSubsequenceAt(chromosome, position, position);
			reference = new String(rs.getBases());
			alternate = reference+basesInserted;
//System.out.println("\tINS "+chromosome+":"+position+" "+reference+" "+alternate+" "+this.parsedAttributes);
			return;
		}
		
		//multiple changes  1322_1336>AT, hmm not sure this can be resolved to genomic coordinates with the given information
		//another odd duck  325-43_518>GC, 594_599+11>GGC, 14344-42_14388>31 
		Matcher multi = multiPat.matcher(eff);
		if (multi.matches()){
			type = "multi";
			//best we can do is pull both start and stop, get a length, use this to fetch the seq from the genome and set it as the Ref, won't be correct in all cases!
			//split on _
			String[] loc = Misc.UNDERSCORE.split(multi.group(1));
			if (loc.length !=2) Misc.printErrAndExit("\nError: failed to parse the coordinates for the multi variant from "+parsedAttributes);
			
			int cDotFirst = parseCoordinate(loc[0]);
			int cDotSecond = parseCoordinate(loc[1]);
			int len = cDotSecond - cDotFirst;
			ReferenceSequence rs = parser.getFasta().getSubsequenceAt(chromosome, position, position+len);
			reference = new String(rs.getBases());
			
			//is alternate a number or bases
			Matcher mat = intPat.matcher(multi.group(2));
			if (mat.matches()) {
				int size =  Integer.parseInt(multi.group(2));
				rs = parser.getFasta().getSubsequenceAt(chromosome, position, position+size);
				alternate = new String(rs.getBases());
			}
			else alternate = multi.group(2);
			
//System.out.println("\tMulti "+chromosome+":"+position+" "+reference+" "+alternate+" "+this.parsedAttributes);
			return;
		}
		
		//should never hit this
		Misc.printErrAndExit("\nError: could not match the variant type from "+parsedAttributes);
	}
	
	/**Compares the given ref to the seq in the fasta.*/
	private void checkReference(){
		ReferenceSequence rs = parser.getFasta().getSubsequenceAt(chromosome, position, position+reference.length()-1);
		String test = new String(rs.getBases());
		if (test.equals(reference) == false){
			failedParsing = true;
			System.err.println("\tWARNING: reference does not match seq from fasta ('"+test+"')?\n\t"+parsedAttributes+"\n\t"+toVcf());
		}
	}
	
	private int parseCoordinate(String val){
		//is it just a number?  e.t. 1322
		Matcher mat = intPat.matcher(val);
		if (mat.matches()) return Integer.parseInt(val);
		
		//is it with a minus sign
		mat = minusPat.matcher(val);
		if (mat.matches()) return Integer.parseInt(mat.group(1)) - Integer.parseInt(mat.group(2));
		
		//is it with a plus sign
		mat = plusPat.matcher(val);
		if (mat.matches()) return Integer.parseInt(mat.group(1)) + Integer.parseInt(mat.group(2));
		
		//should never hit this
		Misc.printErrAndExit("\nError: could not extract the coordinates for "+parsedAttributes);
		return -1;
	}

	private void parseChromPos() {
		String pos = parseString("position");
		String[] chromPos = Misc.COLON.split(pos);
		chromosome = chromPos[0].replace("chr", "");
		position = Integer.parseInt(chromPos[1]);
	}
	
	private double parseDouble(String key){
		String value = parsedAttributes.get(key);
		if (value == null) Misc.printErrAndExit("\nError: failed Variant parsing couldn't find '"+key+"' from "+parsedAttributes);
		return Double.parseDouble(value);
	}
	
	private String parseString(String key){
		String value = parsedAttributes.get(key);
		if (value == null) Misc.printErrAndExit("\nError: failed Variant parsing couldn't find '"+key+"' from "+parsedAttributes);
		//replace whitespace, not allowed in INFO
		value = Misc.WHITESPACE.matcher(value).replaceAll("_");
		return value;
	}

	public boolean isFailedParsing() {
		return failedParsing;
	}
	
	
}
