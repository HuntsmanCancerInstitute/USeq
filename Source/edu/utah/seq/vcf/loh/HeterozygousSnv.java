package edu.utah.seq.vcf.loh;

import java.io.IOException;
import java.util.ArrayList;
import edu.utah.seq.parsers.jpileup.BaseCount;
import edu.utah.seq.parsers.jpileup.BpileupLine;
import util.gen.Num;

public class HeterozygousSnv {

	//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLES
	private String[] vcfRecord = null;
	private BpileupLine germlineBP = null;
	private BpileupLine somaticBP = null;
	private int[] somAltRefGermAltRef = null;
	private int position;
	private double pvalue = -1;
	private boolean passing = true;
	private double somaticAf;
	private double germlineAf;
	private HetWindow bestHetWindow = null;
	
	public HeterozygousSnv (String[] vcfRecord) {
		this.vcfRecord = vcfRecord;
		position = Integer.parseInt(vcfRecord[1]);
	}

	public String[] getVcfRecord() {
		return vcfRecord;
	}

	public void setGermlineBPs(ArrayList<BpileupLine> data, int minimumDepth) throws IOException {
		if (data.size() != 1) passing = false;
		else {
			germlineBP = data.get(0);
			int dp = (int)Math.round(germlineBP.getSamples()[0].getPassingReadCoverageSnv());
			if (dp < minimumDepth) passing = false;
		}
	}
	
	public void setSomaticBPs(ArrayList<BpileupLine> data, int minimumDepth) throws IOException {
		if (data.size() != 1) passing = false;
		else {
			somaticBP = data.get(0);
			int dp = (int)Math.round(somaticBP.getSamples()[0].getPassingReadCoverageSnv());
			if (dp < minimumDepth) passing = false;
		}
	}
	
	public String toStringVcf() {
		StringBuilder sb = new StringBuilder("LoHVar=");
		//LoHSnv=afS,afG,diff,sAlt/Ref,gAlt/Ref,pval;
		sb.append(Num.formatNumber(somaticAf,3)); sb.append(",");
		sb.append(Num.formatNumber(germlineAf,3)); sb.append(",");
		String diff = Num.formatNumber(getAlleleFractionDifference(), 3);
		sb.append(diff); sb.append(",");
		sb.append(somAltRefGermAltRef[0]); sb.append("/"); sb.append(somAltRefGermAltRef[1]); sb.append(",");
		sb.append(somAltRefGermAltRef[2]); sb.append("/"); sb.append(somAltRefGermAltRef[3]); sb.append(",");
		
		String p = Num.formatNumberNoComma(Num.minus10log10(pvalue), 1);
		sb.append(p); sb.append(";");
		return sb.toString();
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder(vcfRecord[0]);
		try {
			for (int i=1; i< vcfRecord.length; i++) {
				sb.append(" ");
				sb.append(vcfRecord[i]);
			}
			sb.append("\n");
			fetchSomAltRefGermAltRefCounts();
			sb.append("\tsom "+somAltRefGermAltRef[0]+"/"+somAltRefGermAltRef[1]+" "+somaticAf+"\n");
			sb.append("\tger "+somAltRefGermAltRef[2]+"/"+somAltRefGermAltRef[3]+" "+germlineAf+"\n");
			sb.append("\t-10Log10(snvPVal) "+Num.minus10log10Float(pvalue)+ "\tafDiff " +getAlleleFractionDifference());
			if(bestHetWindow!=null) {
				sb.append("\n");
				sb.append("\t\tmeanAFs Somatic "+bestHetWindow.getMeanAfSomatic()+"\tGermline "+bestHetWindow.getMeanAfGermline()+"\tDiff "+bestHetWindow.getMeanAfDiff());
				sb.append("\n\t\t-10Log10(combP) "+bestHetWindow.getTransPvalue());
				sb.append("\t-10Log10(adjP)  "+bestHetWindow.getTransAdjPVal());
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		return sb.toString();
	}
	
	public int[] fetchSomAltRefGermAltRefCounts() throws IOException {
		if (somAltRefGermAltRef != null) return somAltRefGermAltRef;
		char alt = vcfRecord[4].charAt(0);
		char ref = vcfRecord[3].charAt(0);
		BaseCount somBC =  somaticBP.getSamples()[0];
		int somAlt = somBC.getSnvCount(alt);
		int somRef = somBC.getSnvCount(ref);
		BaseCount germBC = germlineBP.getSamples()[0];
		int germAlt = germBC.getSnvCount(alt);
		int germRef = germBC.getSnvCount(ref);
		somAltRefGermAltRef = new int[] {somAlt, somRef, germAlt, germRef};
		
		//set AFs
		double total = somAlt+somRef;
		somaticAf = (double)somAlt/ total;
		total = germAlt+germRef;
		germlineAf = (double)germAlt/ total;
		
		return somAltRefGermAltRef;
	}
	
	public int getTotalRefAltCount() throws IOException {
		int[] counts = fetchSomAltRefGermAltRefCounts();
		int sum = 0;
		for (int i: counts) sum+=i;
		return sum;
	}

	public BpileupLine getGermlineBP() {
		return germlineBP;
	}

	public BpileupLine getSomaticBP() {
		return somaticBP;
	}

	public double getAlleleFractionGermline() throws IOException {
		return germlineBP.getSamples()[0].getSnvAlleleFreq(vcfRecord[4].charAt(0));
	}

	public double getAlleleFractionSomatic() throws IOException {
		return somaticBP.getSamples()[0].getSnvAlleleFreq(vcfRecord[4].charAt(0));
	}
	
	public double getAlleleFractionDifference() {
		return somaticAf - germlineAf;
	}

	public boolean isPassing() {
		return passing;
	}

	public double getPvalue() {
		return pvalue;
	}

	public void setPvalue(double pvalue) {
		this.pvalue = pvalue;
	}

	public double getSomaticAf() {
		return somaticAf;
	}

	public double getGermlineAf() {
		return germlineAf;
	}

	public int getPosition() {
		return position;
	}

	public HetWindow getBestHetWindow() {
		return bestHetWindow;
	}

	public void setBestHetWindow(HetWindow bestHetWindow) {
		this.bestHetWindow = bestHetWindow;
	}
}
