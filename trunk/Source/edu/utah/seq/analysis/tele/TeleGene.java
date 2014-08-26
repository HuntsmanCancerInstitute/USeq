package edu.utah.seq.analysis.tele;

import java.util.ArrayList;

import util.bio.parsers.UCSCGeneLine;
import util.gen.Misc;

public class TeleGene {
	
	//fields
	private String name;
	private ArrayList<TeleTranscript> scoredTranscripts;
	
	//0 All Exon, 1 Splice Exon, 2 All Intron, 3 MisSpliced
	private int[] masterTypesTreatment;
	private int[] masterTypesControl;
	
	//chiSquareTest
	private double transMisSplicePVal = 0;
	private double misSpliceLog2Rto = 0;
	
	//constructor
	public TeleGene( ArrayList<TeleTranscript> scoredTranscripts, int[] masterTypesTreatment, int[] masterTypesControl) {
		this.scoredTranscripts = scoredTranscripts;
		this.masterTypesTreatment = masterTypesTreatment;
		this.masterTypesControl = masterTypesControl;
		name = scoredTranscripts.get(0).getTranscript().getDisplayName();
	}
	
	// methods
	
	/**GeneName tTypes cTypes transMisSplicePVal misSpliceLog2Rto*/
	public String toStringTabLine(){
		StringBuilder sb = new StringBuilder();
		sb.append(fetchIGBExcelHyperLink()); sb.append("\t");
		sb.append(Misc.intArrayToString(masterTypesTreatment, ",")); sb.append("\t");
		sb.append(Misc.intArrayToString(masterTypesControl, ",")); sb.append("\t");
		sb.append(transMisSplicePVal); sb.append("\t");
		sb.append(misSpliceLog2Rto); 
		return sb.toString();
	}
	
	private String fetchIGBExcelHyperLink(){
		//find min start and max end
		int start = Integer.MAX_VALUE;
		int end = -1;
		for (TeleTranscript tt : scoredTranscripts){
			if (start > tt.getTranscript().getTxStart()) start = tt.getTranscript().getTxStart();
			if (end < tt.getTranscript().getTxEnd()) end = tt.getTranscript().getTxEnd();
		}
		StringBuilder text = new StringBuilder();
		start = start - 10000;
		if (start < 0) start = 0;
		end = end + 10000;

		text.append("=hyperlink(\"http://localhost:7085/UnibrowControl?seqid=");
		text.append(scoredTranscripts.get(0).getTranscript().getChrom());
		text.append("&start=");
		text.append(start);
		text.append("&end=");
		text.append(end);
		text.append("\",\"");
		text.append(scoredTranscripts.get(0).getTranscript().getDisplayName());
		text.append("\")");
		return text.toString();
	}

	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append(name); sb.append("\tGene name\n");
		sb.append(scoredTranscripts.size()); sb.append("\tNumber scored transcripts\n");
		sb.append(Misc.intArrayToString(masterTypesTreatment, "\t")); sb.append("\tTreatment type count (Exon, SplicedExon, Intron, MisSpliced)\n");
		sb.append(Misc.intArrayToString(masterTypesControl, "\t")); sb.append("\tControl type count (Exon, SplicedExon, Intron, MisSpliced)\n");
		sb.append(transMisSplicePVal); sb.append("\t-10Log10(bonn corr pval) chi-square test\n");
		sb.append(misSpliceLog2Rto); sb.append("\tMis splice log2(((tMisSpliced+1)/(tMisSpliced+2+tSpliced)) / ((cMisSpliced+1)/(cMisSpliced+2+cSpliced)))");
		return sb.toString();
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public ArrayList<TeleTranscript> getScoredTranscripts() {
		return scoredTranscripts;
	}

	public void setScoredTranscripts(ArrayList<TeleTranscript> scoredTranscripts) {
		this.scoredTranscripts = scoredTranscripts;
	}

	public int[] getMasterTypesTreatment() {
		return masterTypesTreatment;
	}

	public void setMasterTypesTreatment(int[] masterTypesTreatment) {
		this.masterTypesTreatment = masterTypesTreatment;
	}

	public int[] getMasterTypesControl() {
		return masterTypesControl;
	}

	public void setMasterTypesControl(int[] masterTypesControl) {
		this.masterTypesControl = masterTypesControl;
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

	public void setMisSpliceLog2Rto(double misSpliceLog2Rto) {
		this.misSpliceLog2Rto = misSpliceLog2Rto;
	}
	
}
