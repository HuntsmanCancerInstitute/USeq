package edu.utah.seq.analysis.multi;

import org.apache.poi.ss.usermodel.Row;

import util.gen.Num;

public class GeneResult implements Comparable<GeneResult>{
	private String name;
	private float fdr;
	private float log2Ratio;
	private float sortBy;
	
	public GeneResult (String name, float fdr, float log2Ratio){
		this.name = name;
		this.fdr = fdr;
		this.log2Ratio = log2Ratio;
		sortBy = Math.abs(log2Ratio);
	}
	
	public void addInfo(Row r, int columnIndex, boolean addStraightFDR){
		//name
		r.createCell(columnIndex).setCellValue(name);
		//lg2Rto deseq
		r.createCell(columnIndex +1).setCellValue(log2Ratio);
		//fdr deseq
		if (addStraightFDR){
			r.createCell(columnIndex +2).setCellValue(Num.antiNeg10log10(fdr));
			r.createCell(columnIndex +3).setCellValue(fdr);
		}
		//fdr deseq
		else r.createCell(columnIndex +2).setCellValue(fdr);
	}

	public int compareTo(GeneResult other) {
		if (this.sortBy < other.sortBy) return 1;
		if (this.sortBy > other.sortBy) return -1;
		return 0;
	}


}
