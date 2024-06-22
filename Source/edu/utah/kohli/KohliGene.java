package edu.utah.kohli;

import java.util.ArrayList;

import util.gen.IO;

public class KohliGene {
	
	private String geneName;
	private ArrayList<CaptureRegion> captureRegions = new ArrayList<CaptureRegion>();

	public String toString() {
		StringBuilder sb = new StringBuilder(geneName);
		sb.append("\n");
		for (CaptureRegion cr: captureRegions) {
			sb.append("\t");
			sb.append(cr.toString());
			sb.append("\n");
		}
		return sb.toString();
	}
	public KohliGene (String geneName) {
		this.geneName = geneName;
	}
	
	public void addCaptureRegion(CaptureRegion cr) {
		captureRegions.add(cr);
	}
	
	public ArrayList<CaptureRegion> getCaptureRegions() {
		return captureRegions;
	}
	public double[] getScaledPonCounts() {
		double[] scaledPonCounts = new double[captureRegions.size()];
		for (int i=0; i< scaledPonCounts.length; i++) {
			scaledPonCounts[i] = captureRegions.get(i).getPonScaledMean();
		}
		return scaledPonCounts;
	}
	public double[] getTestCounts() {
		double[] testCounts = new double[captureRegions.size()];
		for (int i=0; i< testCounts.length; i++) {
			testCounts[i] = captureRegions.get(i).getRawTestCount();
		}
		return testCounts;
	}
	public String getGeneName() {
		return geneName;
	}
	public double compareTestvsPoN() {
		//how about a paired Wilcoxon signed-rank test?
		
		//for each CaptureRegion, calculate z-score of test against PoN
		double totalZScore = 0;
		for (CaptureRegion cr: captureRegions) {
			totalZScore += cr.calculateZScoreFromPoN(cr.getScaledTestCount());
		}
		return totalZScore/(double)captureRegions.size();
	}

}
