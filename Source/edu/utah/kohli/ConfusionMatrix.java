package edu.utah.kohli;

import java.util.HashSet;
import java.util.TreeSet;

import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;


public class ConfusionMatrix {
	
	//fields
	private String[] interrogatedGenes = null;  	//AR, MYC,TP53,BRCA1, etc all genes tested
	private TreeSet<String> keyGeneCalls = null; 	//AR+, MYC+, etc
	private HashSet<String> testGeneCalls = null;	//AR+, MYC-, TP53-
	private float TP = 0;
	private float FP = 0;
	private float FN = 0;
	private float TN = 0;
	
	public ConfusionMatrix (String[] interrogatedGenes, TreeSet<String> keyGeneCalls, HashSet<String> testGeneCalls) {
		this.interrogatedGenes = interrogatedGenes;
		this.keyGeneCalls = keyGeneCalls;
		this.testGeneCalls = testGeneCalls;
		scoreGenes();
	}
	
	public static String toStringHeader() {
		return "TP\tFP\tFN\tTN\tTPR Recall Sensitivity\tFPR\tSpecificity\tPrecision PPV\tFDR\tAccuracy";
	}
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append((int)TP); sb.append("\t"); 
		sb.append((int)FP); sb.append("\t"); 
		sb.append((int)FN); sb.append("\t"); 
		sb.append((int)TN); sb.append("\t"); 
		sb.append(getSensitivityRecallTpr()); sb.append("\t"); 
		sb.append(getFalsePositiveRate()); sb.append("\t");
		sb.append(getSpecificity()); sb.append("\t");
		sb.append(getPrecisionPpv()); sb.append("\t");
		sb.append(getFalseDiscoveryRate()); sb.append("\t");
		sb.append(getAccuracy());
		return sb.toString();
	}
	
	
	public float getAccuracy() {
		return (TP + TN )/ (TP + TN + FP + FN);
	}
	
	public float getPrecisionPpv() {
		return TP / (TP + FP);
	}
	
	public float getSensitivityRecallTpr() {
		return TP / (TP + FN);
	}
	
	public float getSpecificity() {
		return TN / (TN + FP);
	}
	
	public float getFalseDiscoveryRate() {
		return FP / ( FP + TP );
	}
	
	public float getFalsePositiveRate() {
		return FP / (FP + TN);
	}
	
	public float getFScore() {
		return (float) Num.harmonicMean(new double[] {getPrecisionPpv(), getSensitivityRecallTpr()});
	}
	
	//from https://vitalflux.com/cohen-kappa-score-python-example-machine-learning/
	//public float getKappaScore() {
		//float N = TP + FP + FN + TN; 
		//float Po = (TP + TN) / N
	//}

	private void scoreGenes() {
		//for each gene
		for (String geneName: interrogatedGenes) {
			String geneAmp = geneName+"+";
			String geneDel = geneName+"-";
			
			//is it present in the test?
			if (testGeneCalls.contains(geneAmp) || testGeneCalls.contains(geneDel)) {
				//yes, present in the test so it's either a TP or FP
				String testGeneCall = null;
				if (testGeneCalls.contains(geneAmp)) testGeneCall = geneAmp;
				else testGeneCall = geneDel;
				
				//does the test match the key
				if (keyGeneCalls.contains(testGeneCall)) TP++;
				else FP++;
			}
			// OK it is not present in the the test
			else {
				//is it present in the key? thus either a fn or tn
				if (keyGeneCalls.contains(geneAmp) || keyGeneCalls.contains(geneDel)) FN++;
				else TN++;
			}
		}
	}
	
	public static void main (String[] args) {
		String[] interrogatedGenes = new String[]{"AR","BRCA2","CHEK2","MYC","NKX3-1","OPHN1","PIK3CA","PIK3CB","TP53","ZBTB16"};
		
		String wgsKeyResults = "AR+,MYC+,NKX3-1-,OPHN1+,PIK3CA+,PIK3CB+";
		String testResults = "AR+,MYC+,NKX3-1-,OPHN1+,PIK3CB+";
		
		TreeSet<String> wgs = new TreeSet<String>();
		for (String s: Misc.COMMA.split(wgsKeyResults)) wgs.add(s);
		HashSet<String> test = new HashSet<String>();
		for (String s: Misc.COMMA.split(testResults)) test.add(s);
		
		ConfusionMatrix tbtc = new ConfusionMatrix(interrogatedGenes, wgs, test);
		IO.pl(tbtc);
	}

	public float getTP() {
		return TP;
	}

	public float getFP() {
		return FP;
	}

	public float getFN() {
		return FN;
	}

	public float getTN() {
		return TN;
	}

	

}
