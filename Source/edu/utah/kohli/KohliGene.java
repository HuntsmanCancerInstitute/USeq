package edu.utah.kohli;

import java.util.ArrayList;

import trans.main.WilcoxonSignedRankTest;
import util.gen.IO;
import util.gen.Num;

public class KohliGene {
	
	private String geneName;
	private ArrayList<CaptureRegion> captureRegions = new ArrayList<CaptureRegion>();
	private double combinePValueTTest = -1;
	private double fractionPassingAdjPvals;

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
	
	public String getGeneName() {
		return geneName;
	}
	public CopyTestResult compareTestvsPoN() {

		double[] zscores = new double[captureRegions.size()];
		double[] scaledTestCounts = new double[captureRegions.size()];
		double[][] scaledGermlineCounts = new double[captureRegions.size()][];
		
		//for each CaptureRegion in the gene, calculate z-score of test against PoN
		for (int i=0; i< zscores.length; i++) {
			CaptureRegion cr = captureRegions.get(i);
			scaledTestCounts[i] = cr.getScaledTestCount();
			scaledGermlineCounts[i] = cr.getPonScaledCounts();
			zscores[i] = cr.calculateZScoreFromPoN(cr.getScaledTestCount());
		}
		return new CopyTestResult(this, zscores, scaledTestCounts, scaledGermlineCounts);
	}
	
	public double[] fetchWilcoxonSignedRankTestPValues() {
		

		//pull the tumor values from each capture region
		float[] tumor = new float[captureRegions.size()];
		for (int i=0; i< tumor.length; i++) tumor[i] = (float)captureRegions.get(i).getScaledTestCount();
		
		//for each PoN
		int numPoN = captureRegions.get(0).getPonScaledCounts().length;
		double[] pvalues = new double[numPoN];
		for (int i=0; i< numPoN; i++) {
			float[] singlePoN = new float[tumor.length];
			for (int j=0; j< singlePoN.length; j++) singlePoN[j] = (float)(captureRegions.get(j).getPonScaledCounts()[i]);
			//compare the tumor to the normal
			WilcoxonSignedRankTest w = new WilcoxonSignedRankTest(tumor, singlePoN);
			//IO.pl("\nTumor:\t"+Num.floatArrayToString(tumor, "\t"));
			//IO.pl(i+"_PoN:\t"+Num.floatArrayToString(singlePoN, "\t"));
			//IO.pl(w.getPValue());
			pvalues[i] = w.getPValue();
		}
		return pvalues;
		
	}
	
	public double[][] fetchCaptureRegionTestAndPoNScaledScores() {
		
		double[][] crScores = new double[captureRegions.size()][];
		int numScores = captureRegions.get(0).getPonScaledCounts().length + 1;
		for (int i=0; i< crScores.length; i++) crScores[i] = new double[numScores];

		//for each CaptureRegion 
		for (int i=0; i< crScores.length; i++) {
			CaptureRegion cr = captureRegions.get(i);
			crScores[i][0] = cr.getScaledTestCount();
			double[] germline = cr.getPonScaledCounts();
			for (int j=0; j< germline.length; j++) crScores[i][j+1] = germline[j]; 
		}
		return crScores;
	}
	
	public class CopyTestResult {
		
		private KohliGene kohliGene;  //contains CaptureRegion ArrayList
		
		private double[] captureRegionZScores;
		private double[] scaledTestCounts;
		private double[][] scaledGermlineCounts;
		
		public CopyTestResult(KohliGene kohliGene, double[] captureRegionZScores, double[] scaledTestCounts, double[][] scaledGermlineCounts) {
			this.kohliGene = kohliGene;
			this.captureRegionZScores = captureRegionZScores;
			this.scaledTestCounts = scaledTestCounts;
			this.scaledGermlineCounts = scaledGermlineCounts;
		}
		
		public double getMeanZScore() {
			return Num.mean(captureRegionZScores);
		}

		public KohliGene getKohliGene() {
			return kohliGene;
		}

		public double[] getCaptureRegionZScores() {
			return captureRegionZScores;
		}

		public double[] getScaledTestCounts() {
			return scaledTestCounts;
		}

		public double[][] getScaledGermlineCounts() {
			return scaledGermlineCounts;
		}
	}

	public int getNumberTestValuesOutsideCaptureRegions() {
		int numOutside = 0;
		for (CaptureRegion cr: captureRegions) {
			if (cr.isTestOutsidePoN()) numOutside++;
		}
		return numOutside;
	}
	
	public double getMeanCaptureRegionZScore() {
		double meanZScore = 0.0;
		for (CaptureRegion cr: captureRegions) {
			meanZScore+= cr.calculateZScoreFromPoN();
		}
		return meanZScore/(double)captureRegions.size();
	}
	
	public double getPseudoMedianCaptureRegionZScore() {
		double[] zscores = new double[captureRegions.size()];
		for (int i=0; i< zscores.length; i++) zscores[i] = captureRegions.get(i).calculateZScoreFromPoN();
		return Num.pseudoMedian(zscores);
	}
	
	public double getMeanCaptureRegionScaledTestScore() {
		double total = 0.0;
		for (CaptureRegion cr: captureRegions) {
			total+= cr.getScaledTestCount();
		}
		return total/(double)captureRegions.size();
	}
	public double getMeanCaptureRegionPoNScore() {
		double total = 0.0;
		for (CaptureRegion cr: captureRegions) {
			total+= cr.getPonScaledMean();
		}
		return total/(double)captureRegions.size();
	}
	public String getChromosome() {
		return captureRegions.get(0).getChr();
	}
	public double getCombinePValueTTest() {
		return combinePValueTTest;
	}
	public void setCombinePValueTTest(double combinePValueTTest) {
		this.combinePValueTTest = combinePValueTTest;
	}
	public void setFractionPassingAdjPvals(double fractionPassingAdjPvals) {
		this.fractionPassingAdjPvals = fractionPassingAdjPvals;
		
	}
	public double getFractionPassingAdjPvals() {
		return fractionPassingAdjPvals;
	}

}
