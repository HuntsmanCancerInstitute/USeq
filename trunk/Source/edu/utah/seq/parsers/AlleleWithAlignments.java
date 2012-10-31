package edu.utah.seq.parsers;

import java.io.Serializable;
import java.util.*;
import edu.utah.seq.useq.data.RegionScoreText;


import util.bio.annotation.Bed;
import util.bio.seq.Seq;
import util.gen.Misc;

public class AlleleWithAlignments implements Serializable{
	private static final long serialVersionUID = 1L;
	//fields
	String chromosome;
	String strand;
	RegionScoreText allele;
	double[] baseScores;
	String[] parsedBPs;
	Bed[] alignments;
	
	//constructor
	public AlleleWithAlignments (String chromosome, String strand, RegionScoreText allele, double[] baseScores, String[] parsedBPs, Bed[] alignments){
		this.chromosome = chromosome;
		this.strand = strand;
		this.allele = allele;
		this.baseScores = baseScores;
		this.parsedBPs = parsedBPs;
		this.alignments = alignments;
	}
	
	public String toString(){
			String x = chromosome+"\t"+allele.getStart()+"\t"+allele.getStop()+"\t"+allele.getText()+"\t"+allele.getScore()+"\t"+strand+"\t"+parsedBPs.length+"\t"+countNumberUniqueAlignments();
			if (parsedBPs.length!=0) x = x+ "\t" +fetchFractionGATCNBases()+"\t"+Misc.stringArrayToString(parsedBPs, ",")+"\t"+Misc.doubleArrayToString(baseScores, ",");
			else x = x+"\t\t";
			return x;
	}
	
	public String fetchFractionGATCNBases(){
		String seq = Misc.stringArrayToString(parsedBPs, "");
		double[] fractions = Seq.calculateFractionBases(seq);
		return fractions[0] +"\t"+fractions[1] +"\t"+fractions[2] +"\t"+fractions[3] +"\t"+fractions[4];
	}
	
	public String toStringBed(){
		return chromosome+"\t"+allele.getStart()+"\t"+allele.getStop()+"\t"+allele.getText()+"\t"+allele.getScore()+"\t"+strand;
	}
	
	public int countNumberUniqueAlignments(){
		if (alignments.length <2) return alignments.length;
		HashSet<String> uni = new HashSet<String>();
		for (int i=0; i< alignments.length; i++){
			char strand = alignments[i].getStrand();
			if (strand == '+') uni.add(alignments[i].getStart()+"_"+strand);
			else uni.add(alignments[i].getStop()+"_"+strand);
		}
		return uni.size();
	}
	
	public int[] getGATCCounts(){
		int[] b = Seq.countBases(Misc.stringArrayToString(parsedBPs, ""));
		return new int[]{b[0], b[1], b[2], b[3]};
	}
}
