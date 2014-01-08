package edu.utah.seq.vcf.splice;

import java.util.ArrayList;

import util.bio.parsers.UCSCGeneLine;
import edu.utah.seq.vcf.VCFRecord;

public class SpliceHit {
	//fields
	private VCFRecord vcf;
	private int vcfAltIndex;
	private UCSCGeneLine transcript;
	private ArrayList<SpliceJunction> affectedSpliceJunctions = null;
	
	//constructors
	public SpliceHit (VCFRecord vcf, int vcfAltIndex, UCSCGeneLine transcript){
		this.vcf = vcf;
		this.vcfAltIndex = vcfAltIndex;
		this.transcript = transcript;
	}
	
	public void saveSpliceJunction(SpliceJunction sj){
		if (affectedSpliceJunctions == null) affectedSpliceJunctions = new ArrayList<SpliceJunction>();
		affectedSpliceJunctions.add(sj);
	}
	
	//getters and setters
	public VCFRecord getVcf() {
		return vcf;
	}
	public void setVcf(VCFRecord vcf) {
		this.vcf = vcf;
	}
	public int getVcfAltIndex() {
		return vcfAltIndex;
	}
	public void setVcfAltIndex(int vcfAltIndex) {
		this.vcfAltIndex = vcfAltIndex;
	}
	public UCSCGeneLine getTranscript() {
		return transcript;
	}
	public void setTranscript(UCSCGeneLine transcript) {
		this.transcript = transcript;
	}

	public ArrayList<SpliceJunction> getAffectedSpliceJunctions() {
		return affectedSpliceJunctions;
	}
	
}
