package edu.utah.ames.bioinfo;

public class FastqFile {
	
	//fields from fastq file
	private String sequenceIdentifier;
	private String readSequence;
	private String readQualityValues;
	
	//getters and setters
	public String getSequenceIdentifier() {
		return sequenceIdentifier;
	}
	public void setSequenceIdentifier(String sequenceIdentifier) {
		this.sequenceIdentifier = sequenceIdentifier;
	}
	public String getReadSequence() {
		return readSequence;
	}
	public void setReadSequence(String readSequence) {
		this.readSequence = readSequence;
	}
	public String getReadQualityValues() {
		return readQualityValues;
	}
	public void setReadQualityValues(String readQualityValues) {
		this.readQualityValues = readQualityValues;
	}
	
	
}
