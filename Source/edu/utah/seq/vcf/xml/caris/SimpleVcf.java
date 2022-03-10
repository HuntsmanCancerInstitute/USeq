package edu.utah.seq.vcf.xml.caris;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.Gzipper;
import util.gen.IO;
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
	private String format;
	private String sample;
	private int end;
	private int padPos;
	private int padEnd;
	public static final String xrv = "##INFO=<ID=XRV,Number=.,Type=String,Description=\"Caris xml reported genomic alteration\">\n";
	public static final String xrp = "##INFO=<ID=XRP,Number=.,Type=String,Description=\"Caris xml reported pathogenicity\">\n";
	private boolean print = true;
	private SimpleVcf match = null;
	private ArrayList<SimpleVcf> overlap = new ArrayList<SimpleVcf>();
	private static final Pattern afPat = Pattern.compile(".+AF=([\\d+\\.]+).*");
	public static final String ncFilter = "##FILTER=<ID=NC,Description=\"This clinical test variant was not confirmed in subsequent recalling.\">";
	public static final String nrFilter = "##FILTER=<ID=NR,Description=\"This variant was not reported by Caris.\">";
	public static final String mdFilter = "##FILTER=<ID=MD,Description=\"This clinical variant overlapped a recall, had the same type, and was modified using info from the recall.\">";
	public static String infoRAF = "##INFO=<ID=RAF,Number=A,Type=Float,Description=\"Recalled variant allele frequency\">";


	public SimpleVcf (String vcfLine, int bpPadding){
		originalRecord = vcfLine;
		//parse fields
		String[] t = Misc.TAB.split(vcfLine);
		chr = t[0];
		pos = Integer.parseInt(t[1])-1;
		id = t[2];
		ref = t[3];
		alt = t[4];
		qual = t[5];
		filter = t[6];
		info = t[7];
		format = t[8];
		//two samples or one? 
		// Caris just has one sample, the tumor
		if (t.length==10) sample = t[9];
		// Strelka has two, the normal and tumor, just take the tumor
		else if (t.length==11) sample = t[10];
		else Misc.printErrAndExit("\nERROR: sample mismatch for caris or strelka vcf record:\n"+ vcfLine);

		//define the end as the max effected genomic region
		int refLen = ref.length();
		//find max alt length, there many be more than one, not good!
		String[] splitAlts = Misc.COMMA.split(alt);
		int altLen = splitAlts[0].length();
		for (int i=1; i< splitAlts.length; i++){
			if (splitAlts[i].length() > altLen) altLen = splitAlts[i].length();
		}
		if (refLen > altLen) {
			end = pos+refLen;
		}
		else if (refLen < altLen){
			end = pos+altLen;
		}
		else end = pos+altLen;

		//set bp padded positions
		padPos = pos - bpPadding;
		padEnd = end + bpPadding;
	}

	public void appendInfo(String toAppend) {
		info = toAppend+info;
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

	public void appendID(SimpleVcf o) {
		if (o.getId().equals(".") == false) id = id+";"+o.getId();
		else id = o.getId();
	}

	public void appendINFO(SimpleVcf o) {
		info = info+";"+o.info;
	}

	public String fetchCarisKey() throws IOException {
		//chromosome_ref_alt_readDepth_geneName_hgvsCodingChange
		StringBuilder sb = new StringBuilder();
		sb.append(chr); sb.append("_");
		sb.append(ref); sb.append("_");
		sb.append(alt); sb.append("_");


		//parse DP from INFO
		//DP=1134;TI=NM_003482.3;GI=KMT2D;FC=Nonsense;PC=Q3812*;DC=c.11434C>T;LO=EXON;EXON=39;CI="Pathogenic Variant"
		//parse DP, transcript and hgvs
		//cant use transcript id, these differ between the vcf file and the xml file

		String[] splitInfo = Misc.SEMI_COLON.split(info);
		String dp = null;
		String hgvs = null;
		String gi = null;
		for (String i: splitInfo) {
			if (i.startsWith("DC=")) hgvs = i.substring(3);
			else if (i.startsWith("GI=")) gi = i.substring(3);
			else if (i.startsWith("DP=")) dp = i.substring(3);
		}
		if (hgvs == null || gi == null || dp == null) {
			// just a few Caris vcf records with no annotations, just skip
			//IO.el("\t\tERROR: failed to parse one or more of the DC=, TI=, DP= INFO values from "+originalRecord);
			return null;
		}

		sb.append(dp); sb.append("_");
		sb.append(gi); sb.append("_");
		sb.append(hgvs); 

		return sb.toString();
	}


	public String toString(){
		//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
		StringBuilder sb = new StringBuilder();
		sb.append(chr); sb.append("\t");
		sb.append((1+pos)); sb.append("\t");
		sb.append(id); sb.append("\t");
		sb.append(ref); sb.append("\t");
		sb.append(alt); sb.append("\t");
		sb.append(qual); sb.append("\t");
		sb.append(filter); sb.append("\t");
		sb.append(info); 
		if (format!=null && sample!=null) {
			sb.append("\t");
			sb.append(format); 
			sb.append("\t");
			sb.append(sample);
		}
		 
		return sb.toString();
	}

	public int compareTo(SimpleVcf other) {
		//sort by chromosome
		int x = chr.compareTo(other.chr);
		if (x !=0) return x;
		//sort by position
		if (pos < other.pos) return -1;
		if (pos > other.pos) return 1;
		return 0;
	}

	public boolean compareToExact(SimpleVcf o) {
		if (chr.equals(o.chr) == false) return false;
		if (pos != o.pos) return false;
		if (ref.equals(o.ref) == false) return false;
		if (alt.equals(o.alt) == false) return false;
		return true;
	}
	public boolean compareToOverlap(SimpleVcf o) {
		if (chr.equals(o.chr) == false) return false;
		if (padPos > o.padEnd) return false;
		if (o.padPos > padEnd) return false;
		return true;
	}

	/**Sorts a VCF and writes it out compressed but not indexed.  This loads everything into memory so it can blow up. Careful.
	 * @return number of vcf records*/
	public static int sortVcf(File unsortedVcfFile, File sortedOutVcfFile) throws IOException { 
		Gzipper out = new Gzipper(sortedOutVcfFile);
		BufferedReader in = IO.fetchBufferedReader(unsortedVcfFile);
		String line;
		ArrayList<SimpleVcf> records = new ArrayList<SimpleVcf>();
		//load em watching for blanks
		while ((line = in.readLine())!=null){
			if (line.trim().length() == 0) continue;
			if (line.startsWith("#")) out.println(line);
			else records.add(new SimpleVcf(line, 0));
		}
		//sort em
		SimpleVcf[] toSort = new SimpleVcf[records.size()];
		records.toArray(toSort);
		Arrays.sort(toSort);
		//print em
		for (SimpleVcf v: toSort) out.println(v.getOriginalRecord());
		//cleanup
		in.close();
		out.close();
		return toSort.length;
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
	public String getOriginalRecord() {
		return originalRecord;
	}
	public void setAlt(String alt) {
		this.alt = alt;
	}
	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public boolean isPrint() {
		return print;
	}

	public SimpleVcf getMatch() {
		return match;
	}

	public void setPrint(boolean print) {
		this.print = print;
	}

	public ArrayList<SimpleVcf> getOverlap() {
		return overlap;
	}

	public void setOverlap(ArrayList<SimpleVcf> overlap) {
		this.overlap = overlap;
	}

	public void setMatch(SimpleVcf match) {
		this.match = match;
	}
}
