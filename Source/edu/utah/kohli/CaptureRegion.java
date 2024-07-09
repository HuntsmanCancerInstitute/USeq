package edu.utah.kohli;

import java.util.Arrays;
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
	
	//Panel of Normals based on initial 1000 scaling
	private double[] initialPonScaledCounts = null;
	private double initialPonScaledMean = -1;
	private double initialPonStdDev = -1;
	
	//Test dataset
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
	
	public boolean isTestOutsidePoN() {
		double[] minMax = Num.findMinMaxDoubleValues(ponScaledCounts);
		if (scaledTestCount < minMax[0]) return true;
		if (scaledTestCount > minMax[1]) return true;
		return false;
	}
	
	public double calculateZScoreFromPoN(double testValue) {
		return (testValue-ponScaledMean)/ponStdDev;
	}
	
	public double calculateZScoreFromPoN() {
		return (scaledTestCount-ponScaledMean)/ponStdDev;
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
	public void setPonScaledCountsCalculateStats (double[] ponScaledCounts, boolean setInitial) {
		this.ponScaledCounts = ponScaledCounts;
		ponScaledMean = Num.mean(ponScaledCounts);
		ponStdDev = Num.standardDeviation(ponScaledCounts, ponScaledMean);
		if (setInitial) {
			initialPonScaledCounts = Arrays.copyOf(ponScaledCounts, ponScaledCounts.length);
			initialPonScaledMean = ponScaledMean;
			initialPonStdDev = ponStdDev;
		}
	}
	
	/**Replace the working PoN stats with the initially assigned values.*/
	public void copyInitialPoNToCurrentPon() {
		ponScaledCounts = Arrays.copyOf(initialPonScaledCounts, initialPonScaledCounts.length);
		ponScaledMean = initialPonScaledMean;
		ponStdDev = initialPonStdDev;
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

	public double getScaledTestCount() {
		return scaledTestCount;
	}

	public void setScaledTestCount(double scaledTestCount) {
		this.scaledTestCount = scaledTestCount;
	}

	public double[] getPonScaledCounts() {
		return ponScaledCounts;
	}

}
