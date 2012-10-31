package edu.cnv;

/**Represents a particular hit to a CompositCNV from a particular CNV.*/
public class CNVHit{
	
	private CNVGroup cnvGroup;
	private double medianCNV;
	
	public CNVHit (CNVGroup g, double median){
		this.cnvGroup = g;
		this.medianCNV = median;
	}
	
	public String toString(){
		return cnvGroup.getName()+":"+medianCNV;
	}

	public CNVGroup getCnvGroup() {
		return cnvGroup;
	}

	public double getMedianCNV() {
		return medianCNV;
	}
}