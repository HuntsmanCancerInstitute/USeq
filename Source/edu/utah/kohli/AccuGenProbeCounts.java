package edu.utah.kohli;

import util.gen.Num;

public class AccuGenProbeCounts {
	
	private AccuGenProbe probe;
	private int refCounts;
	private int altCounts;
	private double scaledAlleleFraction;
	
	
	public AccuGenProbeCounts(AccuGenProbe probe, int refCounts, int altCounts) {
		this.probe = probe;
		this.refCounts = refCounts;
		this.altCounts = altCounts;
	}
	
	public int getReadDepth(){
		return refCounts+ altCounts;
	}
	
	public double getAlleleFraction() {
		return (double) altCounts/ (double)(refCounts+ altCounts);
	}
	
	public String toString() {
		return refCounts+":"+altCounts+":"+getReadDepth()+":"+Num.formatNumber(getAlleleFraction(), 3);
	}

	public int getRefCounts() {
		return refCounts;
	}

	public int getAltCounts() {
		return altCounts;
	}

	public AccuGenProbe getProbe() {
		return probe;
	}

	public double getScaledAlleleFraction() {
		return scaledAlleleFraction;
	}

	public void setScaledAlleleFraction(double scaledAlleleFraction) {
		this.scaledAlleleFraction = scaledAlleleFraction;
	}

}
