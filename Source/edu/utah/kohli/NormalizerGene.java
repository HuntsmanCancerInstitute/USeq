package edu.utah.kohli;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

import util.gen.Misc;
import util.gen.Num;

public class NormalizerGene {
	
	private String geneName;
	private LinkedHashMap<String,ArrayList<Double>> geneNameGermlineZScores = new LinkedHashMap<String,ArrayList<Double>>();
	private HashMap<String,double[]> geneNameZScoreMinMax = null;
	private ArrayList<String> germlineSampleNames = null;
	
	public NormalizerGene (String geneName) {
		this.geneName = geneName;
	}
	
	public HashMap<String,double[]> getGeneNameZScoreMinMax(){
		if (geneNameZScoreMinMax != null) return geneNameZScoreMinMax;
		
		geneNameZScoreMinMax = new HashMap<String,double[]>();
		
		for (String geneName: geneNameGermlineZScores.keySet()) {
			double min = Double.MAX_VALUE;
			double max = Double.MIN_VALUE;
			
			for (Double d: geneNameGermlineZScores.get(geneName)) {
				if (d>max) max = d;
				if (d<min) min = d;
			}
			geneNameZScoreMinMax.put(geneName, new double[] {min, max});
		}
		
		return geneNameZScoreMinMax;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder(geneName);
		sb.append("\t");
		sb.append(Misc.stringArrayListToString(germlineSampleNames, "\t"));
		sb.append("\n");
		for (String g: geneNameGermlineZScores.keySet()) {
			sb.append(g);
			sb.append("\t");
			sb.append(Num.doubleArrayToString(Num.arrayListOfDoubleToArray(geneNameGermlineZScores.get(g)), "\t"));
			sb.append("\n");
		}
		return sb.toString();
	}
	public void addZScore(String gene, Double zscore) {
		ArrayList<Double> scores = geneNameGermlineZScores.get(gene);
		if (scores == null) {
			scores = new ArrayList<Double>();
			geneNameGermlineZScores.put(gene, scores);
		}
		scores.add(zscore);
	}

	public String getGeneName() {
		return geneName;
	}

	public LinkedHashMap<String, ArrayList<Double>> getGeneNameGermlineZScores() {
		return geneNameGermlineZScores;
	}

	public ArrayList<String> getGermlineSampleNames() {
		return germlineSampleNames;
	}

	public void setGermlineSampleNames(ArrayList<String> germlineSampleNames) {
		this.germlineSampleNames = germlineSampleNames;
	}

}
