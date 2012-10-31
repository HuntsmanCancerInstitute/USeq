package edu.utah.seq.analysis.multi;

import java.util.HashSet;

import util.bio.parsers.UCSCGeneLine;

public class MappedGene {
	private UCSCGeneLine gene;
	private HashSet<String>[] exonReadNames;
	
	public MappedGene(UCSCGeneLine gene){
		this.gene = gene;
		exonReadNames = new HashSet[gene.getExons().length];
	}

	public UCSCGeneLine getGene() {
		return gene;
	}

	public HashSet<String>[] getExonReadNames() {
		return exonReadNames;
	}
	
}
