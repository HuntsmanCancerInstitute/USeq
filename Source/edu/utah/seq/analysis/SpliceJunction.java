package edu.utah.seq.analysis;

public class SpliceJunction implements Comparable {
	//fields
	private String name;
	private int numberTreatmentObservations;
	private int numberControlObservations;
	private float binomialPValue;
	private float log2Ratio;
	
	public SpliceJunction (String name, int numberTreatmentObservations, int numberControlObservations){
		this.name = name;
		this.numberTreatmentObservations = numberTreatmentObservations;
		this.numberControlObservations = numberControlObservations;
	}

	/**Sorts by binomialPValue, biggest to smallest.*/
	public int compareTo (Object o){
		SpliceJunction other = (SpliceJunction) o;
		if (other.binomialPValue > binomialPValue) return 1;
		if (other.binomialPValue < binomialPValue) return -1;
		return 0;
	}
	
	public String toString(){
		return name+"\t"+binomialPValue+"\t"+log2Ratio+"\t"+numberTreatmentObservations+"\t"+ numberControlObservations;
	}
	
	public float getBinomialPValue() {
		return binomialPValue;
	}

	public void setBinomialPValue(float binomialPValue) {
		this.binomialPValue = binomialPValue;
	}

	public float getLog2Ratio() {
		return log2Ratio;
	}

	public void setLog2Ratio(float log2Ratio) {
		this.log2Ratio = log2Ratio;
	}

	public String getName() {
		return name;
	}

	public int getNumberTreatmentObservations() {
		return numberTreatmentObservations;
	}

	public int getNumberControlObservations() {
		return numberControlObservations;
	}
	
}
