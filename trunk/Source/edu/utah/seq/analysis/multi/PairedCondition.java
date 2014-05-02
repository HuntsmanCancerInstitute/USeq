package edu.utah.seq.analysis.multi;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

import util.gen.Misc;
import util.gen.Num;

public class PairedCondition {
	//fields
	private Condition firstCondition;
	private Condition secondCondition;
	
	//deseq analysis
	private String[] geneNames;
	private File diffExpResults;
	private float[][] parsedDiffExpResults;
	private HashSet<String> diffExpGeneNames;
	
	//splice junction info
	HashMap<String, float[]> geneNameSpliceStats;

	public PairedCondition (Condition firstCondition, Condition secondCondition){
		this.firstCondition = firstCondition;
		this.secondCondition = secondCondition;
	}

	/**For every gene name, creates a float[]{FDR, log2VarCorRto} from the results file*/
	public void parseDESeqStatResults(String[] geneNames, float minFDR, float minLog2Rto){
		this.geneNames = geneNames;
		diffExpGeneNames = new HashSet<String>();
		try {
			BufferedReader inStats = new BufferedReader (new FileReader (diffExpResults));
			String line;
			String[] stats;
			Pattern tab = Pattern.compile("\t");
			float maxAdjPVal = 0;
			parsedDiffExpResults = new float[geneNames.length][3];
			for (int x=0; x< geneNames.length; x++){
				//parse stats line: padj, (meanVarCorT- meanVarCorC)
				line= inStats.readLine();
				stats = tab.split(line);
				if (stats.length!=2) Misc.printErrAndExit("One of the DESeq stats R results rows is malformed -> "+line);
				parsedDiffExpResults[x] = new float[2];
				//adjPval
				if (stats[0].equals("NA")) parsedDiffExpResults[x][0] = 0;
				else if (stats[0].equals("Inf") || stats[0].equals("-Inf")) {
					parsedDiffExpResults[x][0] = Float.MIN_VALUE;
				}
				else parsedDiffExpResults[x][0] = Float.parseFloat(stats[0]);
				//reset max?
				if (parsedDiffExpResults[x][0] > maxAdjPVal) maxAdjPVal = parsedDiffExpResults[x][0];
				//parse log2rto
				parsedDiffExpResults[x][1] = Float.parseFloat(stats[1]);
			}
			inStats.close();

			//convert Inf to max values * 1%
			maxAdjPVal = maxAdjPVal * 1.01f;
			for (int i=0; i< geneNames.length; i++){
				float[] scores = parsedDiffExpResults[i];
				//check adjPVal
				if (scores[0] == Float.MIN_VALUE) scores[0] = maxAdjPVal;
				//pass?
				if (scores[0] >= minFDR && Math.abs(scores[1]) >= minLog2Rto) {
					diffExpGeneNames.add(geneNames[i]);
				}
			}
		} catch (Exception e){
			System.err.println("\nProblem parsing DESeq stats results from R, see "+diffExpResults);
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	/**For every gene name, creates a float[]{FDR, log2VarCorRto} from the results file*/
	public void parseDESeq2Results(String[] geneNames, float minFDR, float minLog2Rto){
		this.geneNames = geneNames;
		diffExpGeneNames = new HashSet<String>();
		try {
			BufferedReader inStats = new BufferedReader (new FileReader (diffExpResults));
			//skip header
			inStats.readLine();
			String line;
			String[] stats;
			Pattern tab = Pattern.compile("\t");
			float maxAdjPVal = 0;
			parsedDiffExpResults = new float[geneNames.length][3];
			for (int x=0; x< geneNames.length; x++){
				//parse name baseMean log2FoldChange lfcSE stat	pvalue padj
				//        0      1          2          3     4    5      6
				line= inStats.readLine();
				stats = tab.split(line);
				if (stats.length!=7) Misc.printErrAndExit("One of the DESeq2 stats R results rows is malformed -> "+line);
				//check names
				if (stats[0].equals(geneNames[x]) == false) Misc.printErrAndExit("The gene name from the DESeq2 stats R results doesn't match? -> "+line+"\n"+geneNames[x]);
				parsedDiffExpResults[x] = new float[2];
				//adjPval, log transform
				if (stats[6].equals("Inf")) parsedDiffExpResults[x][0] = Float.MIN_VALUE;
				else {
					parsedDiffExpResults[x][0] = Float.parseFloat(stats[6]);
					if (parsedDiffExpResults[x][0] > maxAdjPVal) maxAdjPVal = parsedDiffExpResults[x][0];
				}
				//parse log2rto
				parsedDiffExpResults[x][1] = Float.parseFloat(stats[2]);
			}
			inStats.close();
			
			//convert Inf to max values * 1%
			maxAdjPVal = maxAdjPVal * 1.01f;
			for (int i=0; i< geneNames.length; i++){
				float[] scores = parsedDiffExpResults[i];
				//check adjPVal
				if (scores[0] == Float.MIN_VALUE) scores[0] = maxAdjPVal;
				//pass?
				if (scores[0] >= minFDR && Math.abs(scores[1]) >= minLog2Rto) {
					diffExpGeneNames.add(geneNames[i]);
				}
			}

		} catch (Exception e){
			System.err.println("\nProblem parsing DESeq2 stats results from R, see "+diffExpResults);
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	/**For every gene name, creates a float[]{FDR, log2VarRto} from the SAMseq results file*/
	public void parseSamSeqStatResults(String[] geneNames, float minFDR, float minLog2Rto){
		this.geneNames = geneNames;
		diffExpGeneNames = new HashSet<String>();
		try {
			BufferedReader inStats = new BufferedReader (new FileReader (diffExpResults));
			String line;
			String[] stats;
			Pattern tab = Pattern.compile("\t");
		
			//Parse header
			inStats.readLine();
			parsedDiffExpResults = new float[geneNames.length][3];
			
			for (int x=0; x< geneNames.length; x++){
				//parse stats line: padj, (meanVarCorT- meanVarCorC)
				line= inStats.readLine();
				stats = tab.split(line);
				if (stats.length!=4) Misc.printErrAndExit("One of the SAMseq stats R results rows is malformed -> "+line);
				parsedDiffExpResults[x] = new float[2];
				//Parse FDR and convert to -10 * log10 ( FDR )
				float FDR = Float.parseFloat(stats[3]);
				if (FDR > 1) {
					FDR = 1;
				} else if (FDR == 0) {
					FDR = Float.MIN_VALUE;
				}
				
				FDR = Num.minus10log10Float(FDR);
				
				float log2Ratio = Float.parseFloat(stats[2]);
				
				parsedDiffExpResults[x][0] = FDR;
				parsedDiffExpResults[x][1] = -(log2Ratio);
				
				if (FDR >= minFDR && Math.abs(log2Ratio) >= minLog2Rto) {
					diffExpGeneNames.add(geneNames[x]);
				}
			}
			
			
			inStats.close();

		} catch (Exception e){
			System.err.println("\nProblem parsing SAMseq stats results from R, see "+diffExpResults);
			e.printStackTrace();
			System.exit(1);
		}
	}

	public String getName(){
		return firstCondition.getName()+"_"+secondCondition.getName();
	}
	
	public String getVsName(){
		return firstCondition.getName()+" vs "+secondCondition.getName();
	}

	public Condition getFirstCondition() {
		return firstCondition;
	}

	public void setFirstCondition(Condition firstCondition) {
		this.firstCondition = firstCondition;
	}

	public Condition getSecondCondition() {
		return secondCondition;
	}

	public void setSecondCondition(Condition secondCondition) {
		this.secondCondition = secondCondition;
	}

	public File getDiffExpResults() {
		return diffExpResults;
	}

	public void setDiffExpResults(File diffExpResults) {
		this.diffExpResults = diffExpResults;
	}

	public float[][] getParsedDiffExpResults() {
		return parsedDiffExpResults;
	}

	public void setParsedDiffExpResults(float[][] parsedDiffExpResults) {
		this.parsedDiffExpResults = parsedDiffExpResults;
	}

	public HashMap<String, float[]> getGeneNameSpliceStats() {
		return geneNameSpliceStats;
	}

	public void setGeneNameSpliceStats(HashMap<String, float[]> geneNameSpliceStats) {
		this.geneNameSpliceStats = geneNameSpliceStats;
	}

	public HashSet<String> getDiffExpGeneNames() {
		return diffExpGeneNames;
	}
}
