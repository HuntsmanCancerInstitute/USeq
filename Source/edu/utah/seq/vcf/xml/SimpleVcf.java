package edu.utah.seq.vcf.xml;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Misc;

public class SimpleVcf  implements Comparable<SimpleVcf>{

	//fields
	private String originalRecord;
	private String chr;
	private int pos;  //zero based
	private String id;
	private String ref;
	private String alt;
	private String qual;
	private String filter;
	private String info;
	private int end;
	private int length;
	private int padPos;
	private int padEnd;
	private boolean insertion = false;
	private boolean deletion = false;
	private boolean shortVariant = true;
	private boolean print = true;
	private SimpleVcf match = null;
	private ArrayList<SimpleVcf> overlap = new ArrayList<SimpleVcf>();
	private static final Pattern endPat = Pattern.compile(".+END=(\\d+);.+");
	private static final Pattern afPat = Pattern.compile(".+AF=([\\d+\\.]+).*");
	public static final String ncFilter = "##FILTER=<ID=NC,Description=\"This Foundation variant was not confirmed in subsequent recalling.\">";
	public static final String nrFilter = "##FILTER=<ID=NR,Description=\"This variant was not reported by Foundation.\">";
	public static String infoRAF = "##INFO=<ID=RAF,Number=A,Type=Float,Description=\"Recalled variant allele frequency\">";
	
	//#CHROM POS ID REF ALT QUAL FILTER INFO
	public SimpleVcf (String vcfLine, int bpPadding){
		originalRecord = vcfLine;
		//parse required fields
		String[] t = Misc.TAB.split(vcfLine);
		chr = t[0];
		pos = Integer.parseInt(t[1])-1;
		id = t[2];
		ref = t[3];
		alt = t[4];
		qual = t[5];
		filter = t[6];
		info = t[7];
	
		//short variant? basically anything without the SVTYPE info field
		if (info.contains("SVTYPE")){
			shortVariant = false;
			Matcher mat = endPat.matcher(info);
			if (mat.matches()) end = Integer.parseInt(mat.group(1));
			else Misc.printErrAndExit("\nERROR: failed to match END= from "+vcfLine);
		}
		else {
			//define the end as the max effected genomic region
			int refLen = ref.length();
			int altLen = ref.length();
			if (refLen > altLen) {
				end = pos+refLen;
				deletion = true;
			}
			else if (refLen < altLen){
				end = pos+altLen;
				insertion = true;
			}
			else end = pos+altLen;
		}
		length = end - pos;
		//set bp padded positions
		padPos = pos - bpPadding;
		padEnd = end + bpPadding;
	}
	
	/**This swaps out the pos ref alt end and appends the RAF
	 * @param o */
	public void swapInfoWithOverlap(SimpleVcf o) {
		//Re assign
		pos = o.pos;
		ref = o.ref;
		alt = o.alt;
		end = o.end;
		//append on the RAF
		appendRAF(o);
	}
	
	/**This appends on the AF from the matching SimpleVcf to this vcfLine*/
	public void appendRAF(SimpleVcf o){
		Matcher mat = afPat.matcher(o.info);
		if (mat.matches()) info = info + ";RAF="+mat.group(1);
		else Misc.printErrAndExit("\nERROR: failed to match AF= from "+o.getVcfLine());
	}
	
	/**Appends not confirmed flag to the FILTER field*/
	public void appendFilter(String flag){
		//any existing flags?
		if (filter.equals(".") || filter.toLowerCase().equals("pass")) filter = flag;
		else filter = filter+";"+flag;
		
	}
	
	public String toString(){
		//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLES.....
		return chr+ "\t"+ (1+pos)+ "\t"+ id+ "\t"+ ref+ "\t"+ alt+ "\t"+ qual+ "\t"+ filter+ "\t"+ info;
	}
	
	public boolean compareToExact(SimpleVcf o) {
		if (chr.equals(o.chr) == false) return false;
		if (pos != o.pos) return false;
		if (alt.equals(o.alt) == false) return false;
		return true;
	}
	public boolean compareToOverlap(SimpleVcf o) {
		if (chr.equals(o.chr) == false) return false;
		if (padPos > o.padEnd) return false;
		if (o.padPos > padEnd) return false;
		return true;
	}
	
	public int compareTo(SimpleVcf other) {
		//sort by chromosome
		int x = chr.compareTo(other.chr);
		if (x !=0) return x;
		//sort by position
		if (pos < other.pos) return -1;
		if (pos > other.pos) return 1;
		//sort by length of effected region
		if (length < other.length) return -1;
		if (length > other.length) return 1;
		return 0;
	}
	
	//getters and setters
	public String getVcfLine() {
		return toString();
	}
	public String getChr() {
		return chr;
	}
	public int getPos() {
		return pos;
	}
	public String getRef() {
		return ref;
	}
	public String getAlt() {
		return alt;
	}
	public String getFilter() {
		return filter;
	}
	public String getInfo() {
		return info;
	}
	public boolean isInsertion() {
		return insertion;
	}
	public boolean isDeletion() {
		return deletion;
	}
	public SimpleVcf getMatch() {
		return match;
	}
	public void setMatch(SimpleVcf match) {
		this.match = match;
	}
	public ArrayList<SimpleVcf> getOverlap() {
		return overlap;
	}
	public void setOverlap(ArrayList<SimpleVcf> overlap) {
		this.overlap = overlap;
	}
	public boolean isShortVariant() {
		return shortVariant;
	}
	public boolean isPrint() {
		return print;
	}
	public void setPrint(boolean print) {
		this.print = print;
	}

	public String getOriginalRecord() {
		return originalRecord;
	}
}
