package edu.utah.seq.data.cbio;

import util.bio.annotation.Coordinate;

public class CBioVar implements Comparable<CBioVar> {

	private String sampleId;
	private String chrom;
	private int startMinOne; 
	private int stop; 
	private String ref; 
	private String alt;
	private String proteinChange;
	
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
	
	public String fetchKey() {
		StringBuilder sb = new StringBuilder(sampleId); sb.append("\t");
		sb.append(chrom); sb.append("\t");
		sb.append(startMinOne); sb.append("\t");
		sb.append(stop); sb.append("\t");
		sb.append(ref); sb.append("\t");
		sb.append(alt);
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

}
