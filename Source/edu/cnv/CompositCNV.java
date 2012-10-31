package edu.cnv;

import java.util.ArrayList;

/**Represents a merged or intersected grouping of CNVs.*/
public class CompositCNV{
	private String chromosome;
	private int start;
	private int stop;
	private ArrayList<CNVHit> hits = new ArrayList<CNVHit>();
	private String effectedGenes = "";
	private String effectedCytobands = "";

	public CompositCNV(String chromosome, int start, int stop){
		this.chromosome = chromosome;
		this.start = start;
		this.stop = stop;
	}

	public ArrayList<CNVHit> getHits() {
		return hits;
	}

	public void setHits(ArrayList<CNVHit> hits) {
		this.hits = hits;
	}

	public String getEffectedGenes() {
		return effectedGenes;
	}

	public void setEffectedGenes(String effectedGenes) {
		this.effectedGenes = effectedGenes;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	public String getEffectedCytobands() {
		return effectedCytobands;
	}

	public void setEffectedCytobands(String effectedCytobands) {
		this.effectedCytobands = effectedCytobands;
	}
}