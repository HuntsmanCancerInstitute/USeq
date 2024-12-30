package edu.utah.kohli;

public class GeneCallResult {

	//fields
	private boolean isCopyAltered;
	private boolean isAmplified;
	private String ichorCode = null;
	
	//for SmallPanel, tab delimited: Patient, Test Sample, Germline, Used To Norm, Gene, # Capture Regions, # CRs Test Outside PoN, Mean CR Z-Score, Mean CR Scaled Test, Mean CR Scaled PoN, Log2(Test/PoN)
	//for GATK, ; delimited:         numOb=28;lg2Tum=-0.191;lg2Norm=0
	//for iChorCNA; 
	private String statistics;
	
	//constructor
	public GeneCallResult(boolean isCopyAltered, boolean isAmplified, String statistics) {
		this.isCopyAltered = isCopyAltered;
		this.isAmplified = isAmplified;
		this.statistics = statistics;
		//possibly assign ichorCode
		if (statistics.contains("HOMD")) ichorCode = "-2";
		else if (statistics.contains("HETD")) ichorCode = "-1";
		else if (statistics.contains("NEUT")) ichorCode = "+0";
		else if (statistics.contains("GAIN")) ichorCode = "+1";
		else if (statistics.contains("AMP")) ichorCode = "+2";
		else if (statistics.contains("HLAMP")) ichorCode = "+3";
	}

	public boolean isCopyAltered() {
		return isCopyAltered;
	}

	public String getStatistics() {
		return statistics;
	}

	public boolean isAmplified() {
		return isAmplified;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(isCopyAltered); sb.append("\t");
		sb.append(isAmplified); sb.append("\t");
		sb.append(statistics);
		return sb.toString();
	}

	public String getIchorCode() {
		return ichorCode;
	}
}
