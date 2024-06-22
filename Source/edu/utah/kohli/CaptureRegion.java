package edu.utah.kohli;

import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class CaptureRegion {
	
	private String originalInput;
	private String chr;
	private int start;
	private int stop;
	private String gene;
	private boolean ok = true;
	
	//Panel of Normals
	private double[] ponScaledCounts = null;
	private double ponScaledMean = -1;
	private double ponStdDev = -1;
	
	//Test dataset
	private double rawTestCount = -1;
	private double scaledTestCount = -1;
	
	public CaptureRegion (String rep) {
		//Chr_Start_Stop_Info
		originalInput = rep;
		String[] f = Misc.UNDERSCORE.split(rep);
		chr = f[0];
		start = Integer.parseInt(f[1]);
		stop = Integer.parseInt(f[2]);
		gene = f[3];
	}
	
	public double calculateZScoreFromPoN(double testValue) {
		return (testValue-ponScaledMean)/ponStdDev;
	}
	
	public String toString() {	
		StringBuilder sb = new StringBuilder(originalInput);
		sb.append("\t");
		sb.append(ponScaledMean);
		sb.append("\t[");
		sb.append(Num.doubleArrayToString(ponScaledCounts, ","));
		sb.append("]");
		return sb.toString();
	}
	
	/**for each of the pons, save the scaled counts and calculate the mean*/
	public void setPonScaledCounts (double[] ponScaledCounts) {
		this.ponScaledCounts = ponScaledCounts;
		ponScaledMean = Num.mean(ponScaledCounts);
	}
	

	public String getOriginalInput() {
		return originalInput;
	}

	public String getChr() {
		return chr;
	}

	public String getGene() {
		return gene;
	}

	public boolean isOk() {
		return ok;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	public double getPonScaledMean() {
		return ponScaledMean;
	}

	public double getRawTestCount() {
		return rawTestCount;
	}

	public void setRawTestCount(double rawTestCount) {
		this.rawTestCount = rawTestCount;
	}

	public double getScaledTestCount() {
		return scaledTestCount;
	}

	public void setScaledTestCount(double scaledTestCount) {
		this.scaledTestCount = scaledTestCount;
	}

}
