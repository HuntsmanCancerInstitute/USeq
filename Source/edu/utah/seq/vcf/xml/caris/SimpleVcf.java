package edu.utah.seq.vcf.xml.caris;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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
	private String pathogenicity = null;
	public static final String xrv = "##INFO=<ID=XRV,Number=.,Type=String,Description=\"Caris xml reported genomic alteration\">\n";
	public static final String xrp = "##INFO=<ID=XRP,Number=.,Type=String,Description=\"Caris xml reported pathogenicity\">\n";
	

	public SimpleVcf (String vcfLine){
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
		sample = t[9];
	}
	
	public void appendInfo(String toAppend) {
		info = toAppend+info;
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
		sb.append(info); sb.append("\t");
		sb.append(format); sb.append("\t");
		sb.append(sample); 
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
			else records.add(new SimpleVcf(line));
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
}
