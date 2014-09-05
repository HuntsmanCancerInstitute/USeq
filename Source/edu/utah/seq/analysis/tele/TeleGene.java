package edu.utah.seq.analysis.tele;

import util.bio.parsers.UCSCGeneLine;
import util.gen.Num;

public class TeleGene {
	
	private UCSCGeneLine gene;
	private TeleStats treatment;
	private TeleStats control;
	
	//unspliced stats
	private double transMisSplicePVal = 0;
	private double misSpliceLog2Rto = 0;
	
	//skew stats
	private double pAdjSkewedReadCount;
	
	public TeleGene (UCSCGeneLine gene, TeleStats treatment, TeleStats control){
		this.gene = gene;
		this.treatment = treatment;
		this.control = control;
	}
	
	/**Bunch of info*/
	public String toString(){
		StringBuilder sb = new StringBuilder();
		addIGBExcelHyperLink(sb); sb.append("\t");
		addPNGExcelHyperLink(sb); sb.append("\t");
		sb.append(gene.getTotalExonicBasePairs()); sb.append("\t");
		sb.append(treatment.getNumberExonicAlignments()); sb.append("\t");
		sb.append(treatment.getNumberUnsplicedAlignments()); sb.append("\t");
		sb.append(control.getNumberExonicAlignments()); sb.append("\t");
		sb.append(control.getNumberUnsplicedAlignments()); sb.append("\t");
		
		sb.append(misSpliceLog2Rto); sb.append("\t");
		sb.append(transMisSplicePVal); sb.append("\t");
		sb.append(treatment.getMedianWindowIndex()); sb.append("\t");
		
		sb.append(treatment.getMedianWindow()); sb.append("\t");
		sb.append(treatment.getMedianBackground()); sb.append("\t");
		sb.append(treatment.getMedianSkewLog2Rto()); sb.append("\t");
		sb.append(control.getMedianWindow()); sb.append("\t");
		sb.append(control.getMedianBackground()); sb.append("\t");
		sb.append(control.getMedianSkewLog2Rto()); sb.append("\t");
		sb.append(getMedianSkewLog2Rto()); sb.append("\t");
		
		sb.append(treatment.getCountWindow()); sb.append("\t");
		sb.append(treatment.getCountBackground()); sb.append("\t");
		sb.append(treatment.getCountSkewLog2Rto()); sb.append("\t");
		sb.append(treatment.getBackgroundCoeffVar()); sb.append("\t");
		sb.append(control.getCountWindow()); sb.append("\t");
		sb.append(control.getCountBackground()); sb.append("\t");
		sb.append(control.getCountSkewLog2Rto()); sb.append("\t");
		sb.append(control.getBackgroundCoeffVar()); sb.append("\t");
		sb.append(getCountSkewLog2Rto()); sb.append("\t");
		sb.append(pAdjSkewedReadCount);
		return sb.toString();
	}
	
	private void addPNGExcelHyperLink(StringBuilder sb){
		sb.append("=HYPERLINK(\"Genes/");
		sb.append(gene.getDisplayName());
		sb.append("/");
		sb.append(gene.getDisplayName());
		sb.append("_Exonic.png\",\"");
		sb.append(gene.getName());
		sb.append("\")");
	}
	
	private void addIGBExcelHyperLink(StringBuilder text){
		int start = gene.getTxStart();
		int end = gene.getTxEnd();
		start = start - 10000;
		if (start < 0) start = 0;
		end = end + 10000;
		text.append("=hyperlink(\"http://localhost:7085/UnibrowControl?seqid=");
		text.append(gene.getChrom());
		text.append("&start=");
		text.append(start);
		text.append("&end=");
		text.append(end);
		text.append("\",\"");
		text.append(gene.getDisplayName());
		text.append("\")");
	}

	public double getMedianSkewLog2Rto(){
		return treatment.getMedianSkewLog2Rto() - control.getMedianSkewLog2Rto();
	}
	
	public double getCountSkewLog2Rto(){
		return treatment.getCountSkewLog2Rto() - control.getCountSkewLog2Rto();
	}
	
	public UCSCGeneLine getGene() {
		return gene;
	}

	public TeleStats getTreatment() {
		return treatment;
	}

	public TeleStats getControl() {
		return control;
	}

	public void setMisSpliceLog2Rto(double lrto) {
		this.misSpliceLog2Rto = lrto;
		
	}

	public double getTransMisSplicePVal() {
		return transMisSplicePVal;
	}

	public void setTransMisSplicePVal(double transMisSplicePVal) {
		this.transMisSplicePVal = transMisSplicePVal;
	}

	public double getMisSpliceLog2Rto() {
		return misSpliceLog2Rto;
	}

	public double getpAdjSkewedReadCount() {
		return pAdjSkewedReadCount;
	}

	public void setpAdjSkewedReadCount(double pAdjSkewedReadCount) {
		this.pAdjSkewedReadCount = pAdjSkewedReadCount;
	}

	
	
}
