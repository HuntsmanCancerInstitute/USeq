package edu.utah.seq.vcf.xml;

import java.util.*;
import java.util.regex.Pattern;
import util.gen.Misc;

public class FoundationRearrangeVariant {
	
	//fields
	private LinkedHashMap<String, String> parsedAttributes;
	private MateInfo target;
	private MateInfo other;
	private int id;
	private String status;
	private String inFrame;
	private int supportingReadPairs = 0;
	private String desciptionOfRearrangement;
	private FoundationXml2Vcf parser;
	private boolean failedParsing = false;
	
	//for the vcf header
	private static String altBnd = "##ALT=<ID=BND,Description=\"Structural rearrangement : fusion, truncation, rearrangement, deletion, duplication\">";

	private static String infoSRP = "##INFO=<ID=SRP,Number=1,Type=Integer,Description=\"Supporting read pairs\">";
	private static String infoMG = "##INFO=<ID=MG,Number=1,Type=String,Description=\"Other gene involved in rearrangement\">";
	private static String infoDESC = "##INFO=<ID=DESC,Number=1,Type=String,Description=\"Description of the type of rearrangement observed\">";
	private static String infoMATEID = "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"ID of the mate vcf record\">";
	private static String infoMC = "##INFO=<ID=MC,Number=1,Type=String,Description=\"Cooridnates of the mate\">";
	private static String infoIF = "##INFO=<ID=IF,Number=1,Type=String,Description=\"In-frame status\">";
	
	//for extracting 
	static Pattern coorPat = Pattern.compile("chr([\\dXYM]+):(\\d+)-(\\d+)");
	
	public FoundationRearrangeVariant (LinkedHashMap<String, String> parsedAttributes, FoundationXml2Vcf parser, int id){
		//<rearrangement 
		//pos1="chr6:152080638-152080754" pos2="chr6:151934404-151934799" 
		//status="unknown" supporting-read-pairs="11" targeted-gene="ESR1">
		//description="fusion" in-frame="Yes" other-gene="C6orf97" 
		
		this.parsedAttributes = parsedAttributes;
		this.parser = parser;
		this.id = id;
		
		//parse general data that applies to both
		parseSharedInfo();
		
		//create two MateInfo objects to represent the target gene and other
		createMates();
	}

	private void createMates() {
		//parse genes, the other is often "N/A"
		String targetName = parseString("targeted-gene");
		String otherName = parseString("other-gene");
		
		//parse coordinates
		String targetCoor = parseString("pos1");
		String otherCoor = parseString("pos2");
		
		//create ids
		String targetId = "Foundation_bnd_T"+id;
		String otherId = "Foundation_bnd_M"+id;
		
		//make em, String coordinates, String id, String gene, FoundationXml2Vcf parser
		target = new MateInfo(targetCoor, targetId, targetName, parser);
		other = new MateInfo(otherCoor, otherId, otherName, parser);
	}

	/** Returns two vcf lines, target and other with cross referencing mate info*/
	public String[] toVcf(int id){
		//create common info stuff ST, DESC, IF, SVTYPE,IMPRECISE
		StringBuilder sb = new StringBuilder();
		sb.append(";ST="); sb.append(status);
		sb.append(";DESC="); sb.append(desciptionOfRearrangement);
		sb.append(";IF="); sb.append(inFrame);
		sb.append(";SRP="); sb.append(supportingReadPairs);
		sb.append(";SVTYPE=BND;IMPRECISE");
		String commonInfo = sb.toString();
		
		//get target partial vcf
		String targetVcf = target.getPartialVcf(other);
		String otherVcf = other.getPartialVcf(target);
		
		return new String[]{targetVcf+commonInfo, otherVcf+commonInfo};
	}
	
	public static void appendInfoAltLines(TreeSet<String> alt, TreeSet<String> info){
		//alt
		alt.add(altBnd); 
		//info
		info.add(FoundationShortVariant.infoEG); 
		info.add(infoMG); 
		info.add(infoDESC); 
		info.add(infoIF); 
		info.add(infoSRP); 
		info.add(FoundationCopyVariant.infoType); 
		info.add(infoMATEID); 
		info.add(infoMC); 
		info.add(FoundationShortVariant.infoST); 
		info.add(infoSRP); 
		info.add(infoSRP); 
		info.add(FoundationCopyVariant.infoEnd); 
		info.add(FoundationCopyVariant.infoImp); 

	}
	
	private void parseSharedInfo() {
		desciptionOfRearrangement = parseString("description");
		inFrame = parseString("in-frame").toLowerCase();
		status = parseString("status").toLowerCase();
		supportingReadPairs = (int)parseDouble("supporting-read-pairs");
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
