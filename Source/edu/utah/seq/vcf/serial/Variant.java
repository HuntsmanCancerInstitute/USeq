package edu.utah.seq.vcf.serial;

import util.gen.Misc;

public class Variant {
	
	private String key = null;
	private String consequence = null;
	private boolean foundInSample = false;
	private String gene = null;
	private double vaf = 0.0;
	private int dp = 0;

	public Variant(String[] f) throws Exception {
		//PatientID SampleName SampleNumber VariantID VariantConseq FoundInSample Gene VAF DP
		//    0          1           2           3         4              5         6   7   8
		key = f[3];
		consequence = f[4];
		if (f[5].equals("T")) foundInSample = true;
		else if (f[5].equals("F")) foundInSample = false;
		else throw new Exception("FAILED to find T or F for FoundInSample in line "+Misc.stringArrayToString(f, "\t"));
		gene = f[6];
		vaf = Double.parseDouble(f[7]);
		dp = Integer.parseInt(f[8]);
	}

	public String getKey() {
		return key;
	}

	public String getConsequence() {
		return consequence;
	}

	public boolean isFoundInSample() {
		return foundInSample;
	}

	public String getGene() {
		return gene;
	}

	public double getVaf() {
		return vaf;
	}

	public int getDp() {
		return dp;
	}

}
