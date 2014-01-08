package edu.utah.seq.vcf.splice;

import java.util.HashMap;

import util.bio.parsers.UCSCGeneLine;

public class UCSCGene {
	//fields
	private String geneName;
	private HashMap<String, UCSCGeneLine> transcripts = new HashMap<String, UCSCGeneLine>();
	
	public UCSCGene (String geneName){
		this.geneName = geneName;
	}

	public String getGeneName() {
		return geneName;
	}

	public HashMap<String, UCSCGeneLine> getTranscripts() {
		return transcripts;
	}
	
	
}
