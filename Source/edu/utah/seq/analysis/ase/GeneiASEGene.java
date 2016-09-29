package edu.utah.seq.analysis.ase;

import java.util.ArrayList;

import util.gen.Num;

public class GeneiASEGene {

	private GeneiASEResult result;
	private ArrayList<GeneiASEData> data = new ArrayList<GeneiASEData>();

	/**Gene FDR TotalAlt TotalRef A/R #Feat */
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append(result.getGene()); sb.append("\t");
		double totalAlt = 0;
		double totalRef = 0;
		for (GeneiASEData d: data){
			totalAlt+= d.getAlt();
			totalRef+= d.getRef();
		}
		sb.append((int)totalAlt); sb.append("\t");
		sb.append((int)totalRef); sb.append("\t");
		sb.append(Num.formatNumber(totalAlt/(totalAlt+totalRef), 3)); sb.append("\t");
		sb.append(result.getFdr()); sb.append("\t");
		sb.append(data.size()); 
		for (GeneiASEData d: data){
			sb.append("\t");
			sb.append(d.getSnpAltRef());
		}
		return sb.toString();

	}
	public GeneiASEGene (GeneiASEResult result){
		this.result = result;
	}

	public GeneiASEResult getResult() {
		return result;
	}

	public ArrayList<GeneiASEData> getData() {
		return data;
	}
}
