package edu.utah.seq.vcf.loh;

import java.io.IOException;

import util.gen.Num;

public class HetWindow {
	
	//these pvals have been -10Log10(pval) transformed!
	private float pvalue;
	private float adjPvalue;
	private HeterozygousVar[] hetVars;  //pvalues in these have not been transformed
	private float meanAfSomatic = -1f;
	private float meanAfGermline = -1f;
	private boolean passesThresholds = false;
	
	public HetWindow (float pvalue, HeterozygousVar[] hetVars) {
		this.pvalue = pvalue;
		this.hetVars = hetVars;
	}
	
	public String toStringBedFormat(int bpPadding) throws IOException {
		StringBuilder sb = new StringBuilder();
		//chr,start,stop,name,score
		addCoordinatesRegion(bpPadding, sb, true, "\t");
		sb.append("\t");
		//Name #Vars_MeanDiff_-10Log10(AdjPval)
		String p = Num.formatNumberNoComma(adjPvalue, 1);
		String diff = Num.formatNumber(getMeanAfDiff(), 3);
		sb.append("LoH_"+hetVars.length+"_"+diff+"_"+p);
		sb.append("\t");
		//score
		p = Num.formatNumberNoComma(pvalue, 0);
		sb.append(p);
		return sb.toString();
	}
	
	public String toStringVcfFormat(int bpPadding) throws IOException {
		
		//build the Block INFO string
		//LoHBlock=PASS,start-stop,#vars,meanAfDiffSomGerm,pval,adjpval;
		StringBuilder sb = new StringBuilder();
		sb.append("LoHBlock=");
		if (passesThresholds) sb.append("PASS,");
		else sb.append("FAIL,");
		addCoordinatesRegion(bpPadding, sb, false, "-"); sb.append(",");
		sb.append(hetVars.length); sb.append(",");
		String diff = Num.formatNumber(getMeanAfDiff(), 3);
		sb.append(diff); sb.append(",");
		String p = Num.formatNumberNoComma(pvalue, 1);
		sb.append(p); sb.append(",");
		p = Num.formatNumberNoComma(adjPvalue, 1);
		sb.append(p); sb.append(";");
		String blockInfo = sb.toString();
		
		//for each variant
		StringBuilder vcfs = new StringBuilder();
		for (int i=0; i< hetVars.length; i++) {
			//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLES
			//   0    1   2  3   4   5      6    7     8       9+
			String[] vcfFields = hetVars[i].getVcfRecord();
			for (int x=0; x<7; x++) {
				vcfs.append(vcfFields[x]);
				vcfs.append("\t");
			}
			//add the info fields
			vcfs.append(blockInfo);
			vcfs.append(hetVars[i].toStringVcf());
			vcfs.append(vcfFields[7]);
			
			//add the rest
			for (int x=8; x< vcfFields.length; x++) {
				vcfs.append("\t");
				vcfs.append(vcfFields[x]);
			}
			vcfs.append("\n");
		}

		return vcfs.toString();
	}
	
	public String toString() {

		StringBuilder sb = new StringBuilder();
		try {
			float diff = getMeanAfDiff();
			addCoordinates(sb);
			sb.append("\n\tMeanAFDiff "+diff +" ("+meanAfSomatic+"-"+meanAfGermline+")");
			sb.append("\n\t-10Log10(combPVal) "+pvalue+ "\t-10Log10(adjP)  "+adjPvalue);
			sb.append("\nLoH block vars:\n");
			for (int i=0; i< hetVars.length; i++) {
				sb.append(hetVars[i]);
				sb.append("\n");
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
		return sb.toString();
	}

	private void addCoordinates(StringBuilder sb) {
		String firstPos = hetVars[0].getVcfRecord()[1];
		String lastPos = hetVars[hetVars.length-1].getVcfRecord()[1];
		sb.append(hetVars[0].getVcfRecord()[0]); //chrom
		sb.append(":");
		sb.append(firstPos);
		sb.append("-");
		sb.append(lastPos);
	}
	
	private void addCoordinatesRegion(int bpPadding, StringBuilder sb, boolean includeChr, String delimiter) {
		Integer firstPos = hetVars[0].getPosition()- bpPadding;
		Integer lastPos = hetVars[hetVars.length-1].getPosition()+bpPadding;
		if (includeChr) {
			sb.append(hetVars[0].getVcfRecord()[0]); //chrom
			sb.append(delimiter);
		}
		
		sb.append(firstPos.toString());
		sb.append(delimiter);
		sb.append(lastPos.toString());
	}

	public float getTransPvalue() {
		return pvalue;
	}
	
	public float getTransAdjPVal() {
		return adjPvalue;
	}
	
	public float getMeanAfDiff() throws IOException {
		if (meanAfSomatic != -1) return meanAfSomatic - meanAfGermline;
		double somAfSum = 0;
		double germAfSum = 0;
		for (HeterozygousVar h: hetVars) {
			somAfSum += h.getAlleleFractionSomatic();
			germAfSum += h.getAlleleFractionGermline();
		}
		double count = hetVars.length;
		meanAfSomatic = (float)(somAfSum/count);
		meanAfGermline = (float)(germAfSum/count);
		return meanAfSomatic - meanAfGermline;
	}
	
	public void setTransPVal(float transPVal ) {
		pvalue = transPVal;
	}

	public void setAdjTransPvalue(float adjPvalue) {
		this.adjPvalue = adjPvalue;
	}

	public HeterozygousVar[] getHetVars() {
		return hetVars;
	}

	public HetWindow compare(HetWindow other) throws IOException {
		
		//all pass
		boolean thisPass = this.isPassesThresholds();
		boolean otherPass = other.isPassesThresholds();
		
		//one pass and the other doesn't
		if (thisPass==true && otherPass==false) return this;
		if (thisPass==false && otherPass==true) return other;
		
		//both pass or both fail return one with > pval
		if (this.adjPvalue > other.adjPvalue) return this;
		if (this.adjPvalue < other.adjPvalue) return other;
		
		//must have same pvalue so return greatest diff
		if (this.getMeanAfDiff() > other.getMeanAfDiff()) return this;
		if (this.getMeanAfDiff() < other.getMeanAfDiff()) return other;
		
		//ok same pval and diff so return null;
		return null;
	}

	public boolean isPassesThresholds() {
		return passesThresholds;
	}

	public void setPassesThresholds(boolean passesThresholds) {
		this.passesThresholds = passesThresholds;
	}

	public float getMeanAfSomatic() {
		return meanAfSomatic;
	}

	public float getMeanAfGermline() {
		return meanAfGermline;
	}




}
