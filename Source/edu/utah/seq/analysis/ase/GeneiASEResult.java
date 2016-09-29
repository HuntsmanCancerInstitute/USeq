package edu.utah.seq.analysis.ase;

import java.io.File;
import java.util.ArrayList;

import util.gen.IO;
import util.gen.Misc;

public class GeneiASEResult {
	
	//fields
	private String gene;
	private double fdr;

	/*
	feat	n.vars	mean.s	median.s	sd.s	cv.s	liptak.s	p.nom	fdr
	12_38716579_C_T	1	0	0	NA	NA	0	1	1
	14_70833767_G_T	1	0	0	NA	NA	0	1	1
	ENSG00000003756_RBM5	1	4.51678588235554	4.51678588235554	NA	NA	4.51678588235554	0.10542	0.860264798206278
	*/
	public GeneiASEResult (String line){
		String[] t = Misc.TAB.split(line);
		if (t.length != 9) Misc.printErrAndExit("\nToo few columns in the results line "+line);
		gene = t[0];
		fdr = Double.parseDouble(t[8]);
	}

	public String getGene() {
		return gene;
	}
	public double getFdr() {
		return fdr;
	}

	public static GeneiASEResult[] loadResults(File resultsFile) {
		String[] lines = IO.loadFile(resultsFile);
		ArrayList<GeneiASEResult> al = new ArrayList<GeneiASEResult>();
		for (String line : lines){
			if (line.startsWith("feat") == false) al.add(new GeneiASEResult(line));
		}
		GeneiASEResult[] r = new GeneiASEResult[al.size()];
		al.toArray(r);
		return r;
	}

	
}
