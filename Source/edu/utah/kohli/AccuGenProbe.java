package edu.utah.kohli;

import util.gen.Misc;

public class AccuGenProbe {
	
	private String originalInput;
	private String chr;
	private int pos; //this is 1 base space so subtract one before using is queries
	private String ref;
	private String alt;
	private String gene;
	private boolean ok = true;
	
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

}
