package edu.utah.seq.cnv.wgs;

import edu.utah.seq.useq.data.RegionScoreText;
import util.gen.Num;

public class CopyRegion {

	private String chr;
	private RegionScoreText region;
	private float tumorObs;
	private double normalizedTumorObs;
	private float[] ponObs;
	private float[] normalizedPonObs;
	private double normalizedPonStandardDeviation = Double.MAX_VALUE;
	private double normalizedPonMean = 0;
	
	public CopyRegion(String chr, RegionScoreText region, float tumorObs, float[] ponObs) {
		this.chr = chr;
		this.region = region;
		this.tumorObs = tumorObs;
		this.ponObs = ponObs;
		normalizedPonObs = new float[ponObs.length];
	}

	public double calculateZScore() {
		if (normalizedPonStandardDeviation == Double.MAX_VALUE) {
			normalizedPonStandardDeviation = Num.standardDeviation(normalizedPonObs);
			normalizedPonMean = (double)Num.mean(normalizedPonObs);
		}
		return (normalizedTumorObs-normalizedPonMean)/normalizedPonStandardDeviation;
	}
	
	public boolean isTumorZScoreOutsidePonZScores() {
		double minPonZScore = 0;
		double maxPonZScore = 0;
		for (int i=0; i< normalizedPonObs.length; i++) {
			double mockTumor = normalizedPonObs[i];
			double zs = (mockTumor-normalizedPonMean)/normalizedPonStandardDeviation;
			if (zs > maxPonZScore) maxPonZScore = zs;
			else if (zs < minPonZScore) minPonZScore = zs;
		}
		double tumorZScore = calculateZScore();
		if (tumorZScore < minPonZScore) return true;
		if (tumorZScore > maxPonZScore) return true;
		return false;
	}
	

	public String getChr() {
		return chr;
	}

	public RegionScoreText getRegion() {
		return region;
	}

	public float getTumorObs() {
		return tumorObs;
	}

	public float[] getPonObs() {
		return ponObs;
	}

	public float[] getNormalizedPonObs() {
		return normalizedPonObs;
	}

	public double getNormalizedTumorObs() {
		return normalizedTumorObs;
	}

	public void setNormalizedTumorObs(float normalizedTumorObs) {
		this.normalizedTumorObs = normalizedTumorObs;
	}
}
