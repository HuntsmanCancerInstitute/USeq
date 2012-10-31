package trans.main;

import util.bio.parsers.UCSCGeneLine;

public class ScoredGene implements Comparable{

	//fields
	private UCSCGeneLine gene = null;
	/*float[replicaNumber][oligoIntensityIndex]*/
	private float[][] treatmentValues = null;
	private float[][] controlValues = null;
	private double median = 0;
	private double misc = 0;
	
	//constructor
	public ScoredGene (UCSCGeneLine geneReference, float[][] treatmentValues, float[][] controlValues){
		gene = geneReference;
		this.treatmentValues = treatmentValues;
		this.controlValues = controlValues;
	}
	
	public ScoredGene (UCSCGeneLine geneReference){
		gene = geneReference;
	}

	public int compareTo (Object obj){
		ScoredGene other = (ScoredGene) obj;
		if (other.median > median) return 1;
		if (other.median < median) return -1;
		return 0;
	}
	
	public double getMedian() {
		return median;
	}

	public void setMedian(double median) {
		this.median = median;
	}

	public float[][] getControlValues() {
		return controlValues;
	}

	public UCSCGeneLine getGene() {
		return gene;
	}

	public float[][] getTreatmentValues() {
		return treatmentValues;
	}

	public double getMisc() {
		return misc;
	}

	public void setMisc(double misc) {
		this.misc = misc;
	}
}
