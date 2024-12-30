package edu.utah.kohli;

import java.util.ArrayList;

import util.gen.Misc;
import util.gen.Num;

public class AccuGenProbe {
	
	private String originalInput;
	private String chr;
	private int pos; //this is 1 base space so subtract one before using is queries
	private String ref;
	private String alt;
	private String gene;
	private boolean ok = true;
	private ArrayList<Double> ponAfs = new ArrayList<Double>();
	private double ponAfsMean = Double.NaN;
	private double ponStdDev = 0;
	
	public AccuGenProbe (String rep) {
		//Chr_POS_REF_ALT_GENE
		originalInput = rep;
		String[] f = Misc.UNDERSCORE.split(rep);
		chr = f[0];
		pos = Integer.parseInt(f[1]);
		ref = f[2];
		alt = f[3];
		gene = f[4];
	}
	
	public double calculateZScore (double testAf) {
		if (Double.isNaN(ponAfsMean)) {
			double[] afs = Num.arrayListToDoubles(ponAfs);
			ponAfsMean = Num.mean(afs);
			ponStdDev = Num.standardDeviation(afs, ponAfsMean);
		}
		return (testAf-ponAfsMean)/ponStdDev;
	}
	

	public String getOriginalInput() {
		return originalInput;
	}

	public String getChr() {
		return chr;
	}

	public int getPos() {
		return pos;
	}

	public String getRef() {
		return ref;
	}

	public String getAlt() {
		return alt;
	}

	public String getGene() {
		return gene;
	}

	public boolean isOk() {
		return ok;
	}

	public ArrayList<Double> getPonAfs() {
		return ponAfs;
	}

}
