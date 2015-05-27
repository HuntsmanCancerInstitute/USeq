package edu.utah.seq.cnv;

import util.bio.parsers.UCSCGeneLine;

public class GeneExonSample {
	private UCSCGeneLine gene;
	//all indexes zero based
	private int globalExonIndex;
	private short exonIndex;
	private short sampleIndex;
	//data
	private int count;
	private float zscore;
	private float residual;
	private float obsExpLgRto;
	private String genotypePosition;
	private boolean passedThresholds = false;

	public GeneExonSample(UCSCGeneLine gene, int globalExonIndex, short exonIndex, short sampleIndex, int count){
		this.gene = gene;
		this.globalExonIndex = globalExonIndex;
		this.exonIndex = exonIndex;
		this.sampleIndex = sampleIndex;
		this.count = count;
	}

	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append(gene.getDisplayNameThenName()); sb.append("\t");
		sb.append(globalExonIndex); sb.append("\t");
		sb.append(exonIndex); sb.append("\t");
		sb.append(sampleIndex); sb.append("\t");
		sb.append(count);
		return sb.toString();
	}
	
	public String fetchRDataLine(){
		return (globalExonIndex+1)+"\t"+(sampleIndex+1)+"\t"+count;
	}

	public float getResidual() {
		return residual;
	}

	public void setResidual(float residual) {
		this.residual = residual;
	}

	public float getObsExpLgRto() {
		return obsExpLgRto;
	}

	public void setObsExpLgRto(float obsExpLgRto) {
		this.obsExpLgRto = obsExpLgRto;
	}

	public UCSCGeneLine getGene() {
		return gene;
	}

	public int getGlobalExonIndex() {
		return globalExonIndex;
	}

	public short getExonIndex() {
		return exonIndex;
	}

	public short getSampleIndex() {
		return sampleIndex;
	}

	public int getCount() {
		return count;
	}

	/**Tab delim: log2(obs/exp), residual, count*/
	public String getDataString() {
		return obsExpLgRto+"\t"+residual+"\t"+zscore+"\t"+count;
	}
	/**_ delim: log2(obs/exp)_residual_count*/
	public String getDataStringUnderscore() {
		return obsExpLgRto+"_"+residual+"_"+zscore+"_"+count;
	}

	public String getGenotypePosition() {
		return genotypePosition;
	}

	public void setGenotypePosition(String genotypePosition) {
		this.genotypePosition = genotypePosition;
	}

	public boolean isPassedThresholds() {
		return passedThresholds;
	}

	public void setPassedThresholds(boolean passedThresholds) {
		this.passedThresholds = passedThresholds;
	}

	public float getZscore() {
		return zscore;
	}

	public void setZscore(float zscore) {
		this.zscore = zscore;
	}
}
