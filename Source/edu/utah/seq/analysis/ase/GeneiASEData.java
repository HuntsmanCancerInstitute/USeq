package edu.utah.seq.analysis.ase;

import java.io.File;
import java.util.ArrayList;

import util.gen.IO;
import util.gen.Misc;

public class GeneiASEData {
	
	//fields
	private String gene;
	private String snp;
	private int altCount;
	private int refCount;
	
	/*
	gene	snp.id	alt.dp	ref.dp
	12_38716579_C_T	12_38716579_C_T_1012_0	5	5
	14_70833767_G_T	14_70833767_G_T_1178_1	5	5
	ENSG00000007062_PROM1	4_16020161_C_T_839_22	10	10
	*/
	public GeneiASEData (String line){
		String[] t = Misc.TAB.split(line);
		if (t.length != 4) Misc.printErrAndExit("\nToo few columns in the data line "+line);
		gene = t[0];
		String[] c = Misc.UNDERSCORE.split(t[1]);
		//reverse it to alt ref
		snp = c[0]+ "_"+ c[1]+ "_" +c[3]+ "_"+ c[2];
		altCount = Integer.parseInt(t[2]);
		refCount = Integer.parseInt(t[3]);
	}
	
	public String getSnpAltRef(){
		return snp+"_"+ altCount+"_"+ refCount;
	}
	
	public static GeneiASEData[] loadData(File f){
		String[] lines = IO.loadFile(f);
		ArrayList<GeneiASEData> al = new ArrayList<GeneiASEData>();
		for (String line : lines){
			if (line.startsWith("gene") == false) al.add(new GeneiASEData(line));
		}
		GeneiASEData[] r = new GeneiASEData[al.size()];
		al.toArray(r);
		return r;
	}

	public String getGene() {
		return gene;
	}

	public String getSnp() {
		return snp;
	}

	public int getAlt() {
		return altCount;
	}

	public int getRef() {
		return refCount;
	}
}
