package edu.utah.seq.vcf.loh;

import java.io.File;

public class BamPileupNormalizationJob {

	//fields
	private File bamPileup;
	private PassingHet[] hets;
	private int minimumDP;
	private double[] stats;
	
	//constructor
	public BamPileupNormalizationJob(File bamPileup, NonLoH nonLoH, double[] stats) {
		this.bamPileup = bamPileup;
		this.stats = stats;
		//must have at least 1/2 the minimum
		minimumDP = nonLoH.getMinimumDP()/2;
		hets = nonLoH.getPassingHets();
	}
	
}
