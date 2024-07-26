package edu.utah.seq.data.cbio;

import java.io.File;
import java.io.IOException;

import util.gen.Misc;

public class CBioVar implements Comparable<CBioVar> {

	private String sampleId = "";
	private String chrom;
	private int startMinOne; 
	private int stop = -1; 
	private String ref; 
	private String alt;
	private String proteinChange;
	
	//extras for Vcf
	private String gene;
	private String variantType;
	private String hgvsC;
	private String clinvar;
	private int count=0;
	
	
	public CBioVar (String sampleId, String chrom, int startMinOne, int stop, String ref, String alt, String proteinChange) {
		this.sampleId = sampleId;
		this.chrom = chrom;
		if (chrom.equals("23")) chrom = "X";
		else if (chrom.equals("24")) chrom = "Y";
		this.startMinOne = startMinOne;
		this.stop = stop;
		this.ref = ref;
		this.alt = alt;
		this.proteinChange = proteinChange;
	}
	
	//"Gene", "Chromosome", "Start Pos", "End Pos", "Ref", "Var", "Variant Type", "Protein Change", "HGVSc", "Sample Type"
	public CBioVar (String gene, String chrom, String startPos, String endPos, String ref, String alt, String proteinChange, String variantType, String hgvsC, String clinvar) {
		this.gene = gene;
		this.chrom = chrom;
		//fix and add chr
		if (this.chrom.equals("23")) this.chrom = "chrX";
		else if (this.chrom.equals("24")) this.chrom = "chrY";
		else this.chrom = "chr"+ this.chrom;
		int pos = Integer.parseInt(startPos);
		startMinOne = pos-1;
		if (endPos.length()>0) stop = Integer.parseInt(endPos);
		else stop = pos;
		this.ref = ref;
		this.alt = alt;
		this.proteinChange = proteinChange;
		this.variantType = variantType;
		this.hgvsC = hgvsC;
		this.clinvar = clinvar;
		if (this.clinvar.length()==0) this.clinvar = "NA";
	}
	
	public String fetchKey() {
		StringBuilder sb = new StringBuilder(sampleId); sb.append("\t");
		sb.append(chrom); sb.append("\t");
		sb.append(startMinOne); sb.append("\t");
		sb.append(stop); sb.append("\t");
		sb.append(ref); sb.append("\t");
		sb.append(alt);
		return sb.toString();
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(chrom); sb.append("\t");
		sb.append(startMinOne); sb.append("\t");
		sb.append(stop); sb.append("\t");
		sb.append(ref); sb.append("\t");
		sb.append(alt); sb.append("\t");
		sb.append(variantType); sb.append("\t");
		sb.append(proteinChange); sb.append("\t");
		sb.append(hgvsC); sb.append("\t");
		sb.append(clinvar); sb.append("\t");
		return sb.toString();
	}
	
	public String toVcf(String id) {
		//CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
		StringBuilder sb = new StringBuilder();
		sb.append(chrom); sb.append("\t");
		//pos is 1 base coordinates
		sb.append(startMinOne+1); sb.append("\t");
		sb.append(id); sb.append("\t");
		sb.append(ref); sb.append("\t");
		sb.append(alt); sb.append("\t");
		//qual and filter
		sb.append(".\t.\t");
		//info
		sb.append("gn=");sb.append(gene);
		sb.append(";vt=");sb.append(variantType); 
		sb.append(";pc=");sb.append(proteinChange); 
		sb.append(";hg=");sb.append(hgvsC); 
		sb.append(";cv=");sb.append(clinvar); 
		sb.append(";ct=");sb.append(count);
		return sb.toString();
	}
	public static String getVcfHeader(File fasta) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("##fileformat=VCFv4.1\n");
		sb.append("##reference=file:"+fasta.getCanonicalPath()+"\n");
		sb.append("##INFO=<ID=gn,Number=1,Type=String,Description=\"Associated gene\">\n");
		sb.append("##INFO=<ID=vt,Number=1,Type=String,Description=\"Variant type\">\n");
		sb.append("##INFO=<ID=pc,Number=1,Type=String,Description=\"Protein change\">\n");
		sb.append("##INFO=<ID=hg,Number=1,Type=String,Description=\"HGVSc\">\n");
		sb.append("##INFO=<ID=cv,Number=1,Type=String,Description=\"ClinVar\">\n");
		sb.append("##INFO=<ID=ct,Number=1,Type=Integer,Description=\"Number of occurrences\">\n");
		sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		return sb.toString();
	}

	public String getSvg() {
		float sub = stop - startMinOne;
		int pos = 0;
		if (sub == 1.0f) pos = startMinOne;
		else {
			pos = Math.round(sub/2) + startMinOne;
		}
		return chrom+"\t"+pos+"\t1";
	}
	
	public String getBed(boolean skipSampleName) {
		String si = "";
		if (skipSampleName==false) si = sampleId+"_";
		return chrom+"\t"+startMinOne+"\t"+stop+"\t"+si+ref+"_"+alt+"_"+proteinChange+ "\t0\t.";
	}
	
	/**Sorts by chromsome, start position, length (smallest to largest).*/
	public int compareTo(CBioVar otherCoor){
		//sort by chromosome
		int compare = otherCoor.chrom.compareTo(chrom);
		if (compare !=0) return compare * -1;;
		//sort by start position
		if (startMinOne<otherCoor.startMinOne) return -1;
		if (startMinOne>otherCoor.startMinOne) return 1;
		// if same start, sort by length, smaller to larger
		int len = stop-startMinOne;
		int otherLen = otherCoor.stop-otherCoor.startMinOne;
		if (len<otherLen) return -1;
		if (len>otherLen) return 1;
		return 0;
	}

	public boolean isGoodForVcf() {
		if (chrom.equals("chrNA") || startMinOne<1 || ref.length()==0 || ref.equals("0") || alt.length()==0) return false;
		return true;
	}
	
	public void checkFixRefAlt(String workingSequence) {
		//check the ref
		if (ref.equals("-")==false) {
			String seq = workingSequence.substring(startMinOne, startMinOne+ref.length());
			if (seq.equals(ref)==false) Misc.printErrAndExit("\nRefSeq '"+seq+"' doesn't match for "+ toString());
		}
		
		//insertion
		if (ref.equals("-")) {
			//IO.pl(toVcf("pre"));
			//get the seq before
			String seq = workingSequence.substring(startMinOne-1, startMinOne).toUpperCase();
			//assign it to ref
			ref = seq;
			//add it to the alt
			alt = ref+alt;
			//reset start 
			startMinOne--;
			//should be 1bp in size for insertion point
			stop = startMinOne+1;
			//IO.pl(toVcf("pos")+"\n"+getBed(true)+"\n");
		}

		//deletion
		else if (alt.equals("-")) {
			//add a base to the ref and alt
			//IO.pl(toVcf("pre"));
			//get the seq before
			String seq = workingSequence.substring(startMinOne-1, startMinOne).toUpperCase();
			//assign it to ref
			ref = seq + ref;
			alt = seq;
			//should be the effected bps in the reference
			stop = startMinOne+ ref.length();
			//reset start 
			startMinOne--;
			//reset stop
			stop = startMinOne+ ref.length();
			//IO.pl(toVcf("pos")+"\n"+getBed(true)+"\n");
		}
		//OK so some seq is in the ref and alt, but it's not a snv, append ref base to both
		else if (ref.length()>1) {
			//IO.pl(toVcf("pre"));
			//get the seq before
			String seq = workingSequence.substring(startMinOne-1, startMinOne).toUpperCase();
			//assign it to ref and alt
			ref = seq + ref;
			alt = seq + alt;
			//stop should be the effected bps in the reference
			stop = startMinOne+ ref.length();
			//reset start 
			startMinOne--;
			//reset stop
			stop = startMinOne+ ref.length();
			//IO.pl(toVcf("pos")+"\n"+getBed(true)+"\n");
		}
	}

	public int getStop() {
		return stop;
	}

	public void setStop(int stop) {
		this.stop = stop;
	}

	public String getRef() {
		return ref;
	}

	public void setRef(String ref) {
		this.ref = ref;
	}

	public String getAlt() {
		return alt;
	}

	public void setAlt(String alt) {
		this.alt = alt;
	}

	public String getSampleId() {
		return sampleId;
	}

	public String getChrom() {
		return chrom;
	}

	public int getStartMinOne() {
		return startMinOne;
	}

	public String getProteinChange() {
		return proteinChange;
	}

	public String getGene() {
		return gene;
	}

	public String getVariantType() {
		return variantType;
	}

	public String getHgvsC() {
		return hgvsC;
	}

	public String getClinvar() {
		return clinvar;
	}

	public int getCount() {
		return count;
	}

	public void setCount(int count) {
		this.count = count;
	}
}
