package edu.utah.seq.analysis;

import java.util.ArrayList;
import util.gen.*;

public class BaseContext implements Comparable{

	private String sequence;
	private long numberConvertedReads;
	private long numberNonConvertedReads;
	private double fractionNonConvertedReads = -1;
	private long numberConvertedGenomicContexts = 0;
	private long numberNonConvertedGenomicContexts = 0;
	private double fractionNonConvertedGenomicContexts = -1;
	private Histogram histogram;
	private ArrayList<Float> pValuesAL = new ArrayList<Float>();
	
	public BaseContext(String sequence){
		this.sequence = sequence;
		histogram = new Histogram(0, 1.1, 11);
	}
	
	public String toString(){
		double ratioReads = (double) numberNonConvertedReads/ (double)(numberConvertedReads+ numberNonConvertedReads);
		double ratioGC = (double) numberNonConvertedGenomicContexts/ (double)(numberConvertedGenomicContexts+ numberNonConvertedGenomicContexts);
		return sequence+"\t"+numberNonConvertedReads+"\t"+numberConvertedReads+"\t"+Num.formatNumber(ratioReads, 3)+"\t"+ numberNonConvertedGenomicContexts+"\t"+ numberConvertedGenomicContexts +"\t"+ Num.formatNumber(ratioGC, 3);
	}

	public int compareTo(Object o) {
		BaseContext other = (BaseContext)o;
		if (fractionNonConvertedReads > other.fractionNonConvertedReads) return -1;
		if (fractionNonConvertedReads < other.fractionNonConvertedReads) return 1;
		return 0;
	}
	
	public void setFractionNonConverted(){
		fractionNonConvertedReads = (double) numberNonConvertedReads/ (double)(numberConvertedReads+ numberNonConvertedReads);
		fractionNonConvertedGenomicContexts = (double) numberNonConvertedGenomicContexts/ (double)(numberConvertedGenomicContexts+ numberNonConvertedGenomicContexts);
	}
	
	public void addFractionNonConvertedToHistogram (double fractionNonConverted){
		histogram.count(fractionNonConverted);
	}
	
	public void incrementNumberConvertedReads (long numCon){
		numberConvertedReads += numCon;
	}
	public void incrementNumberNonConvertedReads (long numNonCon){
		numberNonConvertedReads += numNonCon;
	}
	
	public void incrementNumberConvertedGenomicContexts(){
		numberConvertedGenomicContexts++;
	}
	
	public void incrementNumberNonConvertedGenomicContexts(){
		numberNonConvertedGenomicContexts++;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public long getNumberConvertedReads() {
		return numberConvertedReads;
	}

	public void setNumberConvertedReads(long numberConvertedReads) {
		this.numberConvertedReads = numberConvertedReads;
	}

	public long getNumberNonConvertedReads() {
		return numberNonConvertedReads;
	}

	public void setNumberNonConvertedReads(long numberNonConvertedReads) {
		this.numberNonConvertedReads = numberNonConvertedReads;
	}

	public double getFractionNonConvertedReads() {
		if (fractionNonConvertedReads == -1) setFractionNonConverted();
		return fractionNonConvertedReads;
	}

	public long getNumberConvertedGenomicContexts() {
		return numberConvertedGenomicContexts;
	}

	public long getNumberNonConvertedGenomicContexts() {
		return numberNonConvertedGenomicContexts;
	}

	public double getFractionNonConvertedGenomicContexts() {
		if (fractionNonConvertedGenomicContexts == -1) setFractionNonConverted();
		return fractionNonConvertedGenomicContexts;
	}

	public Histogram getHistogram() {
		return histogram;
	}

	public ArrayList<Float> getpValuesAL() {
		return pValuesAL;
	}

}
