package edu.utah.seq.analysis.ase;

import java.io.Serializable;
import java.util.ArrayList;

import util.gen.Misc;
import util.gen.Num;

public class AllelicExpressionSnp implements Serializable{
	private static final long serialVersionUID = 1L;
	private String chromosome;
	private int position;
	private String name;
	private String reference;
	private String alternate;
	private ArrayList<String> sampleName;
	private ArrayList<int[]> gatcnCounts;
	private ArrayList<Float> sampleGCScore;
	private float fdr = 0;
	private float log2Ratio = 0;
	private ArrayList<String> geneNames;

	public AllelicExpressionSnp(String chromosome, int position, String name, String reference, String alternate){
		this.chromosome = chromosome;
		this.position = position;
		this.name = name;
		this.reference = reference;
		this.alternate = alternate;
	}

	public int[] getRefAltCounts(int[] gatcn){
		int ref = 0; 
		int alt = 0;
		if (reference.equals("G")) ref = gatcn[0];
		else if (reference.equals("A")) ref = gatcn[1];
		else if (reference.equals("T")) ref = gatcn[2];
		else if (reference.equals("C")) ref = gatcn[3];

		if (alternate.equals("G")) alt = gatcn[0];
		else if (alternate.equals("A")) alt = gatcn[1];
		else if (alternate.equals("T")) alt = gatcn[2];
		else if (alternate.equals("C")) alt = gatcn[3];
		return new int[]{ref,alt};
	}

	public void initializeALs(){
		sampleName = new ArrayList<String>();
		gatcnCounts = new ArrayList<int[]>();
		sampleGCScore = new ArrayList<Float>();
	}

	public String toStringCoordinates(){
		StringBuilder sb = new StringBuilder();
		sb.append(chromosome); sb.append("_");
		sb.append(position); sb.append("_");
		sb.append(name); sb.append("_");
		sb.append(reference); sb.append("_");
		sb.append(alternate); 
		
		return sb.toString();
	}
	
	public String getChrPosRefAlt(){
		StringBuilder sb = new StringBuilder();
		sb.append(chromosome); sb.append("_");
		sb.append(position); sb.append("_");
		sb.append(reference); sb.append("_");
		sb.append(alternate); 
		
		return sb.toString();
	}

	/**Chr\tPos\tRef\tAlt\tName\tSampleNames\tRefAltCounts\tTotalRef\tTotalAlt\tSampleGCScores\tFDR\tLog2Rto\tBadHetCall?*/
	public String toString(boolean skipGeneral){
		StringBuilder sb = new StringBuilder();
		if (skipGeneral == false){
			sb.append(chromosome); sb.append("\t");
			sb.append(position); sb.append("\t");
			sb.append(reference); sb.append("\t");
			sb.append(alternate); sb.append("\t");
			sb.append(name); sb.append("\t");
		}
		int numRef = 0;
		int numAlt = 0;
		if (sampleName != null){
			//sampleNames
			String sampleNames = Misc.stringArrayListToString(sampleName, ",");
			sb.append(sampleNames); sb.append("\t");
			
			//ref/alt counts
			int[] refAltCounts = getRefAltCounts(gatcnCounts.get(0));
			sb.append(refAltCounts[0]); sb.append("_"); sb.append(refAltCounts[1]);
			numRef+= refAltCounts[0];
			numAlt+= refAltCounts[1];
			for (int i=1; i< gatcnCounts.size(); i++){
				sb.append(",");
				refAltCounts = getRefAltCounts(gatcnCounts.get(i));
				sb.append(refAltCounts[0]); sb.append("_"); sb.append(refAltCounts[1]);
				numRef+= refAltCounts[0];
				numAlt+= refAltCounts[1];
			}
			sb.append("\t");
			//totals
			sb.append(numRef);sb.append("\t");
			sb.append(numAlt);sb.append("\t");
			
			//sample scores
			float[] s = Num.arrayListOfFloatToArray(sampleGCScore);
			String scores = Misc.floatArrayToString(s, ",");
			sb.append(scores); sb.append("\t");
			//fdr
			sb.append(fdr); sb.append("\t");
			//log2Rto
			sb.append(log2Ratio);
			//bad het call?
			if ((numRef+numAlt !=0) && (numRef ==0 || numAlt== 0)) sb.append("\ttrue");
			else sb.append("\tfalse");
			//gene names?
			if (geneNames != null) {
				sb.append("\t");
				sb.append(Misc.stringArrayListToString(geneNames, ","));
			}
			
			
			
		}
		return sb.toString();
	}
	
	public int[] getTotalRefAltCounts() {
		int numRef = 0;
		int numAlt = 0;
		if (sampleName != null){
			for (int i=0; i< gatcnCounts.size(); i++){
				int[] refAltCounts = getRefAltCounts(gatcnCounts.get(i));
				numRef+= refAltCounts[0];
				numAlt+= refAltCounts[1];
			}
		}
		return new int[]{numRef, numAlt};
	}

	public boolean passMinCoverage(int minimumReadCoverage){
		for (int i=0; i< gatcnCounts.size(); i++){
			if (Num.sumIntArray(gatcnCounts.get(i)) >= minimumReadCoverage) return true;
		}
		return false;
	}

	public String getChromosome() {
		return chromosome;
	}

	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	public int getPosition() {
		return position;
	}

	public void setPosition(int position) {
		this.position = position;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getReference() {
		return reference;
	}

	public void setReference(String reference) {
		this.reference = reference;
	}

	public String getAlternate() {
		return alternate;
	}

	public void setAlternate(String alternate) {
		this.alternate = alternate;
	}

	public ArrayList<String> getSampleName() {
		return sampleName;
	}

	public void setSampleName(ArrayList<String> sampleName) {
		this.sampleName = sampleName;
	}

	public ArrayList<int[]> getGatcnCounts() {
		return gatcnCounts;
	}

	public void setGatcnCounts(ArrayList<int[]> gatcnCounts) {
		this.gatcnCounts = gatcnCounts;
	}

	public ArrayList<Float> getSampleGCScore() {
		return sampleGCScore;
	}

	public void setSampleGCScore(ArrayList<Float> sampleGCScore) {
		this.sampleGCScore = sampleGCScore;
	}

	public float getFdr() {
		return fdr;
	}

	public void setFdr(float fdr) {
		this.fdr = fdr;
	}

	public float getLog2Ratio() {
		return log2Ratio;
	}

	public void setLog2Ratio(float log2Ratio) {
		this.log2Ratio = log2Ratio;
	}

	public static long getSerialversionuid() {
		return serialVersionUID;
	}

	public ArrayList<String> getGeneNames() {
		return geneNames;
	}

	public void setGeneNames(ArrayList<String> geneNames) {
		this.geneNames = geneNames;
	}
}

