package edu.utah.kohli;

import util.gen.Num;

public class AccuGenProbeCounts {
	
	private int refCounts;
	private int altCounts;
	
	
	public AccuGenProbeCounts(int refCounts, int altCounts) {
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

}
