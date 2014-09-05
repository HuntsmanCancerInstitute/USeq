package edu.utah.seq.analysis.tele;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import util.gen.Num;

public class TeleStats {
	
	private int numberExonicAlignments;
	private int numberUnsplicedAlignments;
	//these are exonic not genomic
	private int[] baseCoverage;
	private ArrayList<String>[] baseCoverageNames;
	private double medianBackground;
	private int backgroundStartIndex;
	private double medianWindow;
	private int medianWindowIndex;
	private int countBackground;
	private double normalizedBackgroundCount;
	private int countWindow;
	private double normalizedWindowCount;
	private double backgroundCoeffVar;

	public TeleStats(int numberExonicAlignments, int numberUnsplicedAlignments, int[] baseCoverage, ArrayList<String>[] baseCoverageNames) {
		this.numberExonicAlignments = numberExonicAlignments;
		this.numberUnsplicedAlignments = numberUnsplicedAlignments;
		this.baseCoverage = baseCoverage;
		this.baseCoverageNames = baseCoverageNames;
	}
	
	public void nullBigArrays(){
		baseCoverage = null;
		baseCoverageNames = null;
	}
	
	public void windowScan(int length5PrimeWindowScan) {
		int[] toScan = new int[backgroundStartIndex];
		System.arraycopy(baseCoverage, 0, toScan, 0, backgroundStartIndex);
		double[] medians = Num.windowMedianScores(toScan, length5PrimeWindowScan);
		int bestIndex = 0;
		double maxMedian = medians[0];
		for (int i=1; i< medians.length; i++){
			if (medians[i]> maxMedian){
				maxMedian = medians[i];
				bestIndex = i;
			}
		}
		medianWindow = maxMedian;
		medianWindowIndex = bestIndex;
	}
	
	public void calculateMedianBackground(double fractionLengthBackground) {
		int num = (int)Math.round( (double)baseCoverage.length * fractionLengthBackground );
		backgroundStartIndex = baseCoverage.length - num;
		int[] slice = new int[num];
		System.arraycopy(baseCoverage, backgroundStartIndex, slice, 0, num);
		Arrays.sort(slice);
		medianBackground = Num.median(slice);
	}
	
	public void windowScanBackground(int length5PrimeWindowScan) {
		int len = baseCoverage.length - backgroundStartIndex;
		int[] toScan = new int[len];
		System.arraycopy(baseCoverage, backgroundStartIndex, toScan, 0, len);
		double[] aves = Num.windowAverageScoresIgnoreZeros(Num.intArrayToDouble(toScan), length5PrimeWindowScan);
		aves = Num.removeNaNAndInfinite(aves);
		double mean = Num.mean(aves);
		double stdDev = Num.standardDeviation(aves, mean);
		backgroundCoeffVar = stdDev/mean;
		
	}
	
	public void calculateMedianWindow(int windowIndex, int windowLength) {
		medianWindowIndex = windowIndex;
		int[] slice = new int[windowLength];
		System.arraycopy(baseCoverage, medianWindowIndex, slice, 0, windowLength);
		Arrays.sort(slice);
		medianWindow = Num.median(slice);
	}
	
	public void countBackgroundReads() {
		HashSet<String> uniqueNames = new HashSet<String>();
		int num = baseCoverageNames.length;
		for (int i=backgroundStartIndex; i< num; i++) if (baseCoverageNames[i]!=null) uniqueNames.addAll(baseCoverageNames[i]);
		countBackground = uniqueNames.size();
		normalizedBackgroundCount = (double)(countBackground+1) / (double)(num-backgroundStartIndex);
	}
	
	public void countWindowReads(int windowSize){
		int stop = windowSize+medianWindowIndex;
		HashSet<String> uniqueNames = new HashSet<String>();
		for (int i=medianWindowIndex; i< stop; i++) {
			if (baseCoverageNames[i]!=null) uniqueNames.addAll(baseCoverageNames[i]);
		}
		countWindow = uniqueNames.size();
		normalizedWindowCount = (double)(countWindow+1)/ (double) windowSize;
	}
	
	public double getCountSkewLog2Rto(){
		return Num.log2(normalizedWindowCount/ normalizedBackgroundCount);
	}
	
	public double getMedianSkewLog2Rto(){
		double t = (medianWindow +1.0) / (medianBackground + 1.0);
		return Num.log2(t);
	}


	public int getNumberExonicAlignments() {
		return numberExonicAlignments;
	}
	public int getNumberUnsplicedAlignments() {
		return numberUnsplicedAlignments;
	}
	public int[] getBaseCoverage() {
		return baseCoverage;
	}
	public ArrayList<String>[] getBaseCoverageNames() {
		return baseCoverageNames;
	}
	public int getBackgroundStartIndex() {
		return backgroundStartIndex;
	}
	public void setBackgroundStartIndex(int backgroundStartIndex) {
		this.backgroundStartIndex = backgroundStartIndex;
	}
	public double getMedianBackground() {
		return medianBackground;
	}
	public void setMedianBackground(double medianBackground) {
		this.medianBackground = medianBackground;
	}
	public double getMedianWindow() {
		return medianWindow;
	}
	public void setMedianWindow(double medianWindow) {
		this.medianWindow = medianWindow;
	}
	public int getMedianWindowIndex() {
		return medianWindowIndex;
	}
	public void setMedianWindowIndex(int medianWindowIndex) {
		this.medianWindowIndex = medianWindowIndex;
	}

	public int getCountBackground() {
		return countBackground;
	}

	public int getCountWindow() {
		return countWindow;
	}

	public double getBackgroundCoeffVar() {
		return backgroundCoeffVar;
	}
}
