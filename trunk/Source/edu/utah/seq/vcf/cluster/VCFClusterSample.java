package edu.utah.seq.vcf.cluster;

import java.util.ArrayList;

public class VCFClusterSample {
	//fields
	private String sampleName;
	private int sampleIndex;
	private ArrayList<Byte> callsAL = new ArrayList<Byte>();
	/**0=homozygous ref, 1=heterozygous, 2=homozygous alt, 3=no call*/
	private byte[] calls = null;
	
	//constructor
	public VCFClusterSample (String sampleName, int sampleIndex){
		this.sampleName = sampleName;
		this.sampleIndex = sampleIndex;
	}
	/**Returns the number of noCalls*/
	public int makeCalls(){
		int numFail = 0;
		calls = new byte[callsAL.size()];
		for (int i=0; i< calls.length; i++) {
			calls[i] = callsAL.get(i).byteValue();
			if (calls[i] == 3) numFail++;
		}
		callsAL = null;
		return numFail;
	}
	
	//getters and setters
	public byte[] getCalls() {
		return calls;
	}
	public String getSampleName() {
		return sampleName;
	}
	public ArrayList<Byte> getCallsAL() {
		return callsAL;
	}
	public int getSampleIndex() {
		return sampleIndex;
	}
}
