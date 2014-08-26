package edu.utah.seq.analysis.tele;

import java.util.ArrayList;
import java.util.HashSet;

import util.gen.Misc;
import util.gen.Num;

public class TeleStats {
	//fields
	private int numberExonicAlignments;
	private float[] baseCoverage;
	private ArrayList<String>[] baseCoverageNames;
	private float[] exonicBaseCoverage;
	private ArrayList<String>[] exonicBaseCoverageNames;
	private int[] types;
	private int[] misSplicePositions;
	private float fivePrimeMedianBaseCoverage;
	private int fivePrimeWindowIndex;
	private float threePrimeMedianBaseCoverage;
	private int fivePrimeAlignmentFragmentCount;
	private int threePrimeAlignmentFragmentCount;
	private float pAdjSkewedReadCount;
	
	//methods
	public double getLog2MedianSkew(){
		double ratio = (1+fivePrimeMedianBaseCoverage) / (1+threePrimeMedianBaseCoverage);
		return Num.log2(ratio);
	}
	public double getLog2CountSkew(double windowSize, double cDNALength, double fractionLengthBackground){
		double ratio = ((double)(1+fivePrimeAlignmentFragmentCount)/windowSize) / ((double)(1+threePrimeAlignmentFragmentCount)/ (cDNALength* fractionLengthBackground));
		return Num.log2(ratio);
	}
	
	/**Returns tab delim, 5'median, 3'median,log2RtoSkew*/
	public void addMedianInfo(StringBuilder sb){
		sb.append(fivePrimeMedianBaseCoverage); sb.append("\t");
		sb.append(threePrimeMedianBaseCoverage); sb.append("\t");
		sb.append(getLog2MedianSkew()); 
	}
	
	/**Returns tab delim, 5'count, 3'count,count log2RtoSkew*/
	public void addCountInfo(StringBuilder sb, double windowSize, double cDNALength, double fractionLengthBackground){
		sb.append(fivePrimeAlignmentFragmentCount); sb.append("\t");
		sb.append(threePrimeAlignmentFragmentCount); sb.append("\t");
		sb.append(getLog2CountSkew(windowSize,cDNALength, fractionLengthBackground)); 
	}
	
	/**This hashes the names and sets the counts for the best 5' window and the background count.*/
	public void setReadCountsForBestAndBackground(int windowSize, int startIndexBkgrnd) {
		HashSet<String> names = new HashSet<String>();
		//5' end
		int stop = fivePrimeWindowIndex + windowSize;
		for (int i=fivePrimeWindowIndex; i< stop; i++){
			ArrayList<String> al = exonicBaseCoverageNames[i];
			if (al != null){
				for (String name: al) names.add(name);
			}
		}
		fivePrimeAlignmentFragmentCount = names.size();
		names.clear();
		
		//3' end
		for (int i=startIndexBkgrnd; i< exonicBaseCoverageNames.length; i++){
			ArrayList<String> al = exonicBaseCoverageNames[i];
			if (al != null){
				for (String name: al) names.add(name);
			}
		}
		threePrimeAlignmentFragmentCount = names.size();
	}
	
	//getters and setters
	public int getNumberExonicAlignments() {
		return numberExonicAlignments;
	}
	public void setNumberExonicAlignments(int numberExonicAlignments) {
		this.numberExonicAlignments = numberExonicAlignments;
	}
	public float[] getBaseCoverage() {
		return baseCoverage;
	}
	public void setBaseCoverage(float[] baseCoverage) {
		this.baseCoverage = baseCoverage;
	}
	public int[] getTypes() {
		return types;
	}
	public void setTypes(int[] types) {
		this.types = types;
	}
	public int[] getMisSplicePositions() {
		return misSplicePositions;
	}
	public void setMisSplicePositions(int[] misSplicePositions) {
		this.misSplicePositions = misSplicePositions;
	}
	public float getFivePrimeMedianBaseCoverage() {
		return fivePrimeMedianBaseCoverage;
	}
	public void setFivePrimeMedianBaseCoverage(float fivePrimeMedianBaseCoverage) {
		this.fivePrimeMedianBaseCoverage = fivePrimeMedianBaseCoverage;
	}
	public float getThreePrimeMedianBaseCoverage() {
		return threePrimeMedianBaseCoverage;
	}
	public void setThreePrimeMedianBaseCoverage(float threePrimeMedianBaseCoverage) {
		this.threePrimeMedianBaseCoverage = threePrimeMedianBaseCoverage;
	}
	public int getFivePrimeWindowIndex() {
		return fivePrimeWindowIndex;
	}
	public void setFivePrimeWindowIndex(int fivePrimeWindowIndex) {
		this.fivePrimeWindowIndex = fivePrimeWindowIndex;
	}

	public float[] getExonicBaseCoverage() {
		return exonicBaseCoverage;
	}

	public void setExonicBaseCoverage(float[] exonicBaseCoverage) {
		this.exonicBaseCoverage = exonicBaseCoverage;
	}

	public void setBaseCoverageNames(ArrayList<String>[] baseCoverageNames) {
		this.baseCoverageNames = baseCoverageNames;
		
	}

	public ArrayList<String>[] getBaseCoverageNames() {
		return baseCoverageNames;
	}

	public void setExonicBaseCoverageNames(ArrayList<String>[] rcNamesExonic) {
		this.exonicBaseCoverageNames = rcNamesExonic;	
	}

	public int getFivePrimeAlignmentFragmentCount() {
		return fivePrimeAlignmentFragmentCount;
	}

	public int getThreePrimeAlignmentFragmentCount() {
		return threePrimeAlignmentFragmentCount;
	}
	public float getpAdjSkewedReadCount() {
		return pAdjSkewedReadCount;
	}
	public void setpAdjSkewedReadCount(float pAdjSkewedReadCount) {
		this.pAdjSkewedReadCount = pAdjSkewedReadCount;
	}

	


	
	
}
