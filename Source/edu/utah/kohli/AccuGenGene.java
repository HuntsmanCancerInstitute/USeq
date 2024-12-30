package edu.utah.kohli;

import java.util.ArrayList;

public class AccuGenGene {

	private String geneName;
	private ArrayList<AccuGenProbe> probes = new ArrayList<AccuGenProbe>();
	
	public AccuGenGene (String geneName) {
		this.geneName = geneName;
	}
	
	public double meanZScore(double[] testAFs) {
		return 0;
	}

	public ArrayList<AccuGenProbe> getProbes() {
		return probes;
	}
}
