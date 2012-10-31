package util.gen;

import java.io.*;
import java.util.*;

/**
 * Generates histograms.
 */
public class Histogram {
	//fields
	private double minimum;
	private double maximum;
	private double increment;
	private int numberOfBins;
	private Bin[] bins;
	private int numLessThanMin = 0;
	private double numLessThanMinTotalScore = 0;
	private int numMoreThanMax = 0;
	private double numMoreThanMaxTotalScore = 0;
	private double totalBinHits = -1;
	private double totalBinHitsLeftOfZero = -1;
	private double totalBinHitsRightOfZero = -1;
	private boolean meanScoresPresent = false;
	private boolean trimLabelsToSingleInteger = false;
	private boolean skipZeroBins = true;


	public Histogram( double minimum, double maximum, int numberOfBins){
		this.minimum = minimum;
		this.maximum = maximum;
		this.numberOfBins = numberOfBins;
		increment = (maximum-minimum) /numberOfBins;
		makeBins();
	}

	/**Adds the counts of one histogram to this histogram.*/
	public void addCounts(Histogram h) throws Exception{
		//check number of bins is the same
		if (h.getNumberOfBins() != this.numberOfBins) throw new Exception ("Different number of bins found.");
		Bin[] other = h.bins;
		for (int i=0; i< other.length; i++){
			bins[i].addCountToHits(other[i].getHits());
		}
		//less than
		numLessThanMin += h.getNumLessThanMin();
		//more than
		numMoreThanMax += h.getNumMoreThanMax();
	}

	/**Given an array of scores, creates a histogram of set number of bins, trims tails by a set fraction of the data.
	 * Call printScaledHistogram() to print to screen.*/
	public Histogram (double[] scores, int numberOfBins, double tailTrimFraction){
		this.numberOfBins = numberOfBins;
		double[] sorted = new double[scores.length];
		System.arraycopy(scores,0,sorted,0,scores.length);
		Arrays.sort(sorted);
		int index = (int)Math.round(tailTrimFraction * scores.length);
		minimum = sorted[index];
		index = scores.length-index-1;
		maximum = sorted[index];
		increment = (maximum-minimum)/(double)numberOfBins;
		makeBins();
		countAll(scores);
	}

	/**Returns a one tailed, right side, p value given a threshold.
	 * Assumes zero is base.
	 * Returns 1/total bin hits right side if value > last bin. Thus will not return zero.*/
	public double pValueRightSide(double value){
		//has the histogram been totaled?
		if (totalBinHitsRightOfZero == -1) totalBinHitsRightOfZero = numberBinHitsToRightAndIncludingValue(0);
		double totalRightSide = numberBinHitsToRightAndIncludingValue(value);
		if (totalRightSide == 0) return 1.0/totalBinHitsRightOfZero;
		return totalRightSide/totalBinHitsRightOfZero;
	}

	/**Returns a one tailed, left side, p value given a threshold.
	 * Assumes zero is base.
	 * Returns 1/total bin hits left side if value < last bin. Thus will not return zero.*/
	public double pValueLeftSide(double value){
		//has the histogram been totaled?
		if (totalBinHitsLeftOfZero == -1) totalBinHitsLeftOfZero = numberBinHitsToLeftAndIncludingValue(0);
		double totalLeftSide = numberBinHitsToLeftAndIncludingValue(value);
		if (totalLeftSide == 0) return 1.0/totalBinHitsLeftOfZero;
		return totalLeftSide/totalBinHitsLeftOfZero;
	}

	/**Returns a two tailed, p value given a threshold.
	 * Returns 1/total bin hits if value exceeds last bin. Thus will not return zero.*/
	public double pValue(double value){
		//has the histogram been totaled?
		if (totalBinHits == -1) getTotalBinCounts();
		double total;
		if (value >=0) total= numberBinHitsToRightAndIncludingValue(value);
		else total = numberBinHitsToLeftAndIncludingValue(value);
		if (total == 0) return 1.0/totalBinHits;
		return total/totalBinHits;
	}

	/**Returns number of bin hits in the bins to the right of the value's bin
	 * including the numMoreThanMax.*/
	public double numberBinHitsToRightOfValue(double value){
		//find total right of value
		int index = findBinIndex(value);
		if (index == -1) return 0;
		double totalRightSide = 0;
		//start counting in the first bin to the right
		index++;
		for (int i=index; i< numberOfBins; i++){
			totalRightSide += bins[i].getHits();
		}
		totalRightSide += numMoreThanMax;
		return totalRightSide;
	}

	/**Returns number of bin hits in the bins to the right of the value's bin
	 * including the numMoreThanMax.*/
	public double numberBinHitsToRightAndIncludingValue(double value){
		//find total right of value
		int index = findBinIndex(value);
		if (index == -1) return 0;
		double totalRightSide = 0;
		//start counting
		for (int i=index; i< numberOfBins; i++){
			totalRightSide += bins[i].getHits();
		}
		totalRightSide += numMoreThanMax;
		return totalRightSide;
	}

	/**Returns number of bin hits in the bins to the right of the value's bin
	 * including the numMoreThanMax.*/
	public double numberBinHitsToLeftAndIncludingValue(double value){
		//find total right of value
		int index = findBinIndex(value);
		if (index == -1) return 0;
		double totalLeftSide = 0;
		//start counting
		for (int i=index; i>= 0; i--){
			totalLeftSide += bins[i].getHits();
		}
		totalLeftSide += numLessThanMin;
		return totalLeftSide;
	}


	public void makeBins(){
		bins = new Bin[numberOfBins];
		double start = minimum;
		double stop = 0;
		for (int i=0; i< numberOfBins; i++){
			stop = start + increment;
			bins[i] = new Bin(start, stop, this);
			start = stop;
		}
	}

	/**Counts all the values incrementing the proper bin.*/
	public void countAll (float[] values){
		int num = values.length;
		for (int i=0; i<num; i++){
			count(values[i]);
		}
	}

	/**Counts all the values incrementing the proper bin.*/
	public void countAll (double[] values){	
		int num = values.length;
		for (int i=0; i<num; i++){
			count(values[i]);
		}
	}

	/**Counts all the values incrementing the proper bin and scores.*/
	public void countAll (double[] values, double[] scores){	
		int num = values.length;
		meanScoresPresent = true;
		for (int i=0; i<num; i++){
			count(values[i], scores[i]);
		}
	}

	/**Increments the proper bin.*/
	public void count(double value){
		//check to see if value is too small or too tall
		if (value < minimum ) {
			numLessThanMin ++;
			return;
		}
		if (value >= maximum) {
			numMoreThanMax ++;
			return;
		}
		//run through bins
		for (int i=0; i<numberOfBins; i++){
			if (bins[i].count(value)) return;
		}
	}

	/**Increments the proper bin.*/
	public void count(double value, double score){
		//check to see if value is too small or too tall
		if (value < minimum ) {
			numLessThanMin ++;
			numLessThanMinTotalScore += score;
			return;
		}
		if (value >= maximum) {
			numMoreThanMax ++;
			numMoreThanMaxTotalScore += score;
			return;
		}
		//run through bins
		for (int i=0; i<numberOfBins; i++){
			if (bins[i].count(value, score)) return;
		}
	}

	/**Returns the index of the containing bin or -1 if value is too small or big.*/
	public int findBinIndex(double value){
		//check to see if value is too small or too tall
		if (value < minimum || value >= maximum) return -1;

		//run through bins
		for (int i=0; i<numberOfBins; i++){
			if (bins[i].contains(value)) return i;
		}
		return -1;
	}

	/**Given a value, finds the corresponding bin and returns the number of hits to that bin.*/
	public int findBinHits(double value){
		int index = findBinIndex(value);
		if (index == -1) return -1;
		return bins[index].getHits();
	}

	/**Prints scaled histogram.*/
	public void printScaledHistogram(){
		//extract counts and labels from bins
		int[] counts = getBinCounts();
		String[] binLabels = getBinLabels();
		//scale
		double maxCount = Num.findHighestInt(counts);
		double scalar = 100.0/maxCount;

		if (numLessThanMin !=0) System.out.println("<  "+Num.formatNumber(minimum,3)+"\t\t"+numLessThanMin);
		//multiply counts by scalar
		int[] scaledCounts = new int[counts.length];
		for (int i=0; i<scaledCounts.length; i++){
			scaledCounts[i] = (int)Math.round((double)counts[i] * scalar);
		}
		printHistogram(counts, scaledCounts, binLabels, skipZeroBins);
		if (numMoreThanMax !=0) System.out.println(">= "+(Num.formatNumber(maximum,3) )+ "\t\t"+ numMoreThanMax);
	}
	
	/**Prints scaled histogram.*/
	public void printScaledHistogram(PrintWriter out){
		//extract counts and labels from bins
		int[] counts = getBinCounts();
		String[] binLabels = getBinLabels();
		//scale
		double maxCount = Num.findHighestInt(counts);
		double scalar = 100.0/maxCount;

		if (numLessThanMin !=0) out.println("<  "+Num.formatNumber(minimum,3)+"\t\t"+numLessThanMin);
		//multiply counts by scalar
		int[] scaledCounts = new int[counts.length];
		for (int i=0; i<scaledCounts.length; i++){
			scaledCounts[i] = (int)Math.round((double)counts[i] * scalar);
		}
		printHistogram(counts, scaledCounts, binLabels, skipZeroBins, out);
		if (numMoreThanMax !=0) out.println(">= "+(Num.formatNumber(maximum,3) )+ "\t\t"+ numMoreThanMax);
	}
	
	/**Prints histogram.*/
	public void printHistogram(PrintWriter out){
		//extract counts and labels from bins
		int[] counts = getBinCounts();
		String[] binLabels = getBinLabels();
		
		if (numLessThanMin !=0) out.println("<  "+Num.formatNumber(minimum,3)+"\t\t"+numLessThanMin);
		printHistogram(counts, counts, binLabels, skipZeroBins, out);
		if (numMoreThanMax !=0) out.println(">= "+(Num.formatNumber(maximum,3) )+ "\t\t"+ numMoreThanMax);
	}

	/**Prints bin info including mean if scores were provided while counting.*/
	public void printBins(){
		System.out.print("Name\tStart\tMiddle\tEnd\tCounts");
		if (meanScoresPresent){
			String mean = "";
			if (numLessThanMin !=0){
				double ave = numLessThanMinTotalScore/(double)numLessThanMin;
				mean = ave+"";
			}
			System.out.println("\tMean\tMedian");
			System.out.println("<"+bins[0].getStart()+"\t\t\t"+bins[0].getStart()+"\t"+numLessThanMin+"\t"+mean);
			for (int i=0; i< bins.length; i++){
				System.out.println(bins[i]);
			}
			mean = "";
			if (numMoreThanMax !=0){
				double ave = numMoreThanMaxTotalScore/(double)numMoreThanMax;
				mean = ave+"";
			}
			System.out.println(">"+bins[bins.length-1].getStop()+"\t\t\t"+bins[bins.length-1].getStop()+"\t"+numMoreThanMax+"\t"+mean);
		}
		else {
			System.out.println("\n<"+bins[0].getStart()+"\t\t\t"+bins[0].getStart()+"\t"+numLessThanMin);
			for (int i=0; i< bins.length; i++){
				System.out.println(bins[i]);
			}
			System.out.println(">"+bins[bins.length-1].getStop()+"\t\t\t"+bins[bins.length-1].getStop()+"\t"+numMoreThanMax);
		}
	}

	/**Prints a simple histogram scaled using stars but real counts shown.*/
	public static long printHistogram(int[] realCounts, int[] stars, String[] binLabels, boolean skipZeroBins){
		int num = stars.length;
		StringBuffer sb;
		long total =0;

		//find index of first and last realCounts to trim histogram
		int firstIndex = 0;
		int lastIndex = num;
		if (skipZeroBins) {
			for (int i=0; i<num; i++){
				if (realCounts[i] !=0){
					firstIndex = i;
					break;
				}
			}
			for (int i=num-1; i>=0; i--){
				if (realCounts[i] !=0){
					lastIndex = i+1;
					break;
				}
			}
		}

		//make histogram
		for (int i=firstIndex; i<lastIndex; i++){
			//make stars
			sb = new StringBuffer();
			if (stars[i]>150) System.out.println(binLabels[i]+"\t"+realCounts[i]+"\t|***************************************************************" +
			"***************************************************************************************~");
			else {
				for (int j=0; j<stars[i]; j++) sb.append("*");
				System.out.println(binLabels[i]+"\t"+realCounts[i]+"\t|"+sb);
			}
			total+= realCounts[i];
		}
		return total;
	}
	
	/**Prints a simple histogram scaled using stars but real counts shown.*/
	public static long printHistogram(int[] realCounts, int[] stars, String[] binLabels, boolean skipZeroBins, PrintWriter out){
		int num = stars.length;
		StringBuffer sb;
		long total =0;

		//find index of first and last realCounts to trim histogram
		int firstIndex = 0;
		int lastIndex = num;
		if (skipZeroBins) {
			for (int i=0; i<num; i++){
				if (realCounts[i] !=0){
					firstIndex = i;
					break;
				}
			}
			for (int i=num-1; i>=0; i--){
				if (realCounts[i] !=0){
					lastIndex = i+1;
					break;
				}
			}
		}

		//make histogram
		for (int i=firstIndex; i<lastIndex; i++){
			//make stars
			sb = new StringBuffer();
			if (stars[i]>150) out.println(binLabels[i]+"\t"+realCounts[i]+"\t|***************************************************************" +
			"***************************************************************************************~");
			else {
				for (int j=0; j<stars[i]; j++) sb.append("*");
				out.println(binLabels[i]+"\t"+realCounts[i]+"\t|"+sb);
			}
			total+= realCounts[i];
		}
		return total;
	}

	/**Prints a simple histogram to a PrintWriter.*/
	public static void printHistogram(int[] stars, PrintWriter out){
		int num = stars.length;
		StringBuffer sb;
		long total =0;
		int end =1;
		//start backwards looking for a value
		for (int i=num-1; i>=0; i--){
			if (stars[i]!=0){
				end = i+1;
				break;
			}
		}
		//add one to stop to show zero if possible
		if (end!=num) end++;
		//make histogram
		for (int i=0; i<end; i++){
			//make stars
			sb = new StringBuffer();
			if (stars[i]>150) out.println(i+"\t"+stars[i]+"\t|***************************************************************" +
			"***************************************************************************************~");
			else {
				for (int j=0; j<stars[i]; j++) sb.append("*");
				out.println(i+"\t"+stars[i]+"\t|"+sb);
			}
			total+= stars[i];
		}
		out.println("Total: "+total);
	}

	public int[] getReversedBinCounts() {
		int[] rev = new int[numberOfBins];
		int counter = 0;
		for (int i=numberOfBins-1; i>=0; i--){
			rev[counter++] = bins[i].getHits();
		}
		return rev;
	}

	/**Including num less than and num more than.*/
	public double getTotalBinCounts(){
		totalBinHits = 0;
		for (int i=0; i< numberOfBins; i++){
			totalBinHits += bins[i].getHits();
		}
		totalBinHits += numMoreThanMax;
		totalBinHits += numLessThanMin;
		return totalBinHits;
	}

	public int[] getBinCounts(){
		int[] counts = new int[numberOfBins];
		for (int i=0; i< numberOfBins; i++){
			counts[i] = bins[i].getHits();
		}
		return counts;
	}

	public String[] getBinLabels(){
		String[] binLabels = new String[numberOfBins];
		for (int i=0; i< numberOfBins; i++){
			if(trimLabelsToSingleInteger) binLabels[i] = bins[i].getIntLabel();
			else binLabels[i] = bins[i].getLabel();
		}
		return binLabels;
	}


	/**Returns two columns as int[columns][counts] centered on max value with zeros as filler.*/
	public int[][] fetchCenteredFowardAndReverseCounts(){
		int[] revCounts = getReversedBinCounts();
		int[] counts = getBinCounts();
		int maxIndexReverse = Num.findMaxIntIndex(revCounts);
		int maxIndexForward = Num.findMaxIntIndex(counts);
		int diff = Math.abs(maxIndexReverse - maxIndexForward);
		int[] f = new int[counts.length+diff];
		int[] r = new int[f.length];
		int counter = 0;
		if (maxIndexForward< maxIndexReverse){
			for (int i=diff; i<f.length; i++){
				f[i]=counts[counter];
				r[counter] = revCounts[counter];
				counter++;
			}
		}
		else if  (maxIndexForward> maxIndexReverse){
			for (int i=diff; i<f.length; i++){
				r[i]=revCounts[counter];
				f[counter] = counts[counter];
				counter++;
			}
		}
		else {
			f= counts;
			r= revCounts;
		}
		int[][] x = new int[2][];
		x[0] = f;
		x[1] =r;
		return x;
	}
	/**Prints two columns, forward counts and reversed counts centered on max value with zeros as filler.*/
	public void printCenteredForwardAndReverseCounts(){
		int[][] x = fetchCenteredFowardAndReverseCounts();
		System.out.println("Forward\tReverse");
		for (int i=0; i<x[0].length; i++) System.out.println(x[0][i]+"\t"+ x[1][i]);
	}

	/**Prints a simple histogram.*/
	public static void printHistogram(int[] stars){
		int num = stars.length;
		StringBuffer sb;
		long total =0;
		int end =1;
		//start backwards looking for a value
		for (int i=num-1; i>=0; i--){
			if (stars[i]!=0){
				end = i+1;
				break;
			}
		}
		//add one to stop to show zero if possible
		if (end!=num) end++;
		//make histogram
		for (int i=0; i<end; i++){
			//make stars
			sb = new StringBuffer();
			if (stars[i]>150) System.out.println(i+"\t"+stars[i]+"\t|***************************************************************" +
			"***************************************************************************************~");
			else {
				for (int j=0; j<stars[i]; j++) sb.append("*");
				System.out.println(i+"\t"+stars[i]+"\t|"+sb);
			}
			total+= stars[i];
		}
		System.out.println("Total: "+total);
	}

	public double getIncrement() {
		return increment;
	}
	public void setIncrement(double increment) {
		this.increment = increment;
	}
	public double getMaximum() {
		return maximum;
	}
	public void setMaximum(double maximum) {
		this.maximum = maximum;
	}
	public double getMinimum() {
		return minimum;
	}
	public void setMinimum(double minimum) {
		this.minimum = minimum;
	}
	public int getNumberOfBins() {
		return numberOfBins;
	}
	public void setNumberOfBins(int numberOfBins) {
		this.numberOfBins = numberOfBins;
	}

	public int getNumLessThanMin() {
		return numLessThanMin;
	}

	public int getNumMoreThanMax() {
		return numMoreThanMax;
	}

	public boolean isMeanScoresPresent() {
		return meanScoresPresent;
	}

	public boolean isTrimLabelsToSingleInteger() {
		return trimLabelsToSingleInteger;
	}

	public void setTrimLabelsToSingleInteger(boolean trimLabelsToSingleInteger) {
		this.trimLabelsToSingleInteger = trimLabelsToSingleInteger;
	}

	public void setSkipZeroBins(boolean skipZeroBins) {
		this.skipZeroBins = skipZeroBins;
	}
}
