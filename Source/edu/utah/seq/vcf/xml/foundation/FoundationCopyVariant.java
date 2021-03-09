package edu.utah.seq.vcf.xml.foundation;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import htsjdk.samtools.reference.ReferenceSequence;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class FoundationCopyVariant {
	
	//fields
	private LinkedHashMap<String, String> parsedAttributes;
	private String chromosome;
	private int position;  //start
	private String reference;
	private int end;
	private String type; //DEL or DUP
	private String exonsEffected;
	private double log2Ratio;
	private double absCopyNumber = Double.MIN_NORMAL;
	private String highConfidence; 
	private String effectedGene;
	private String status;
	private FoundationXml2Vcf parser;
	private boolean failedParsing = false;
	
	//for the vcf header
	static String altDel = "##ALT=<ID=DEL,Description=\"Deletion/ decrease in copies relative to control\">";
	static String altDup = "##ALT=<ID=DUP,Description=\"Duplication/ increase in copies relative to control\">";
	static String infoEE = "##INFO=<ID=EE,Number=1,Type=String,Description=\"Number of exons effected\">";
	static String infoEnd = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">";
	static String infoLogRatio = "##INFO=<ID=LR,Number=1,Type=Float,Description=\"Log2 ratio of CNV relative to control\">";
	static String infoACN = "##INFO=<ID=ACN,Number=1,Type=Integer,Description=\"Absolute copy number after correcting for % normal in tumor\">";
	static String infoCon = "##INFO=<ID=CON,Number=1,Type=String,Description=\"High confidence in call, true or false\">";
	static String infoType = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant, DUP, DEL, or BND\">";
	static String infoImp = "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">";
	
	//for extracting 
	private static Pattern coorPat = Pattern.compile("chr([\\dXYM]+):(\\d+)-(\\d+)");
	
	public FoundationCopyVariant (LinkedHashMap<String, String> parsedAttributes, FoundationXml2Vcf parser){
		//copy-number="7" equivocal="true" gene="KRAS" number-of-exons="5 of 5" position="chr12:25362722-25398327" ratio="1.47" status="known" type="amplification"/>
		this.parsedAttributes = parsedAttributes;
		this.parser = parser;
		
		parseEffect();
		if (failedParsing == false) parseChromPos();
		if (failedParsing == false) parseRefAlt();
		if (failedParsing == false) parseStats();
	}

	/** Returns CHROM POS ID REF ALT QUAL FILTER INFO */
	public String toVcf(int id){
		StringBuilder sb = new StringBuilder();
		//CHROM
		sb.append(chromosome); sb.append("\t");
		//POS
		sb.append(position); sb.append("\tFoundation_");
		//ID
		sb.append(id);
		sb.append("\t");		
		//REF
		sb.append(reference); sb.append("\t<");
		//ALT
		sb.append(type); sb.append(">\t");
		//QUAL
		sb.append(".\t");
		//FILTER
		sb.append(".\t");
		//INFO (EG EE ST END LR ACN CON TYPE IMPRECISE)
		sb.append("EG=");
		sb.append(effectedGene);
		sb.append(";EE=");
		sb.append(exonsEffected);
		sb.append(";ST=");
		sb.append(status);
		sb.append(";END=");
		sb.append(end);
		sb.append(";LR=");
		sb.append(Num.formatNumber(log2Ratio, 2));
		//any copy-number info?
		if (absCopyNumber != Double.MIN_NORMAL){
			sb.append(";ACN=");
			sb.append(Num.formatNumber(absCopyNumber, 1));
		}
		sb.append(";CON=");
		sb.append(highConfidence);
		sb.append(";SVTYPE=");
		sb.append(type);
		sb.append(";IMPRECISE");
		
		return sb.toString();
		
	}
	
	public static void appendInfoAltLines(TreeSet<String> alt, TreeSet<String> info){
		alt.add(altDel); 
		alt.add(altDup); ;
		info.add(FoundationShortVariant.infoEG);
		info.add(infoEE); 
		info.add(FoundationShortVariant.infoST); 
		info.add(infoEnd); 
		info.add(infoLogRatio); 
		info.add(infoACN); 
		info.add(infoCon); 
		info.add(infoType); 
		info.add(infoImp); 
	}
	
	private void parseEffect() {
		effectedGene = parseString("gene");
		exonsEffected = parseString("number-of-exons").toLowerCase();
		status = parseString("status").toLowerCase();
		highConfidence = parseString("equivocal").toLowerCase();
	}
	
	private void parseStats() {
		log2Ratio = parseDouble("ratio");
		//watch for cases were this is missing
		String cn = parsedAttributes.get("copy-number");
		if (cn == null) return;
		absCopyNumber = parseDouble("copy-number");
	}

	private void parseRefAlt() {
		//fetch the first base
		ReferenceSequence p = parser.getFasta().getSubsequenceAt(chromosome, position, position);
		reference = new String(p.getBases());
		
		//set the type, a little risky
		String t = parseString("type").toLowerCase();
		if (t.startsWith("amp")) type = "DUP";
		else if (t.startsWith("los")) type = "DEL";
		else Misc.printErrAndExit("\nError: unrecognized CNV type '"+t+"' from "+parsedAttributes);
	}

	private void parseChromPos() {
		String pos = parseString("position");
		//chr([\\dXYM]+):(\\d+)-(\\d+)
		Matcher mat = coorPat.matcher(pos);
		if (mat.matches()){
			chromosome = mat.group(1);
			position = Integer.parseInt(mat.group(2));
			end = Integer.parseInt(mat.group(3));
		}
		else {
			//is this one of be bad CDKN2A/B combo genes? 
			if (pos.equals("NA") && effectedGene.equals("CDKN2A/B") && parser.getGenomeVersion().equals("hg19")) {
				chromosome = "9";
				position = 21967751;
				end = 22009312;
			}
			else {
				System.err.println("\tError: failed to match the chr:pos-end from '"+pos+"' for "+parsedAttributes);
				failedParsing = true;
			}
		}
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
