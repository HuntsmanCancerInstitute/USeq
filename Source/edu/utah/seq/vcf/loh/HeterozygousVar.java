package edu.utah.seq.vcf.loh;

import java.io.IOException;
import java.util.ArrayList;
import edu.utah.seq.parsers.jpileup.BaseCount;
import edu.utah.seq.parsers.jpileup.BpileupLine;
import edu.utah.seq.useq.data.RegionScoreText;
import edu.utah.seq.parsers.jpileup.BamPileupTabixLoaderSingle;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class HeterozygousVar {

	//#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLES
	private String[] vcfRecord = null;
	private boolean isSnv;
	private boolean isDeletion;
	private BpileupLine germlineBP = null;
	private BpileupLine somaticBP = null;
	private int[] somAltRefGermAltRef = null;
	private int position;
	private double pvalue = -1;
	private boolean passing = true;
	private double somaticAf;
	private double germlineAf;
	private HetWindow bestHetWindow = null;
	private RegionScoreText copyRatioRegion = null;
	
	public HeterozygousVar (String[] vcfRecord) throws IOException {
		this.vcfRecord = vcfRecord;
		position = Integer.parseInt(vcfRecord[1]);
		//figure out the type
		char allele = BamPileupTabixLoaderSingle.fetchAllele(vcfRecord);
		if (allele == 'D') {
			isDeletion = true;
			isSnv = false;
		}
		else if (allele == 'I') {
			isDeletion = false;
			isSnv = false;
		} 
		else isSnv = true;
	}

	public String[] getVcfRecord() {
		return vcfRecord;
	}

	public void setGermlineBPs(ArrayList<BpileupLine> data, int minimumDepth) throws IOException {
		if (data== null || data.size() == 0) passing = false;
		else {
			germlineBP = data.get(0);
			int dp = -1;
			//snv?
			if (isSnv) dp = (int)Math.round(germlineBP.getSamples()[0].getPassingReadCoverageSnv());
			//indel
			else dp = (int)Math.round(germlineBP.getSamples()[0].getPassingReadCoverageIndel(isDeletion));
			if (dp < minimumDepth) passing = false;
		}
	}
	
	public void setSomaticBPs(ArrayList<BpileupLine> data, int minimumDepth) throws IOException {
		if (data== null || data.size() == 0) passing = false;
		else {
			somaticBP = data.get(0);
			int dp = -1;
			//snv?
			if (isSnv) dp = (int)Math.round(germlineBP.getSamples()[0].getPassingReadCoverageSnv());
			//indel
			else dp = (int)Math.round(germlineBP.getSamples()[0].getPassingReadCoverageIndel(isDeletion));
			if (dp < minimumDepth) passing = false;
		}
	}
	
	public String toStringVcf() throws IOException {
		StringBuilder sb = new StringBuilder("LoHVar=");
		//LoHVar=afS,afG,diff,sAlt/Ref,gAlt/Ref,pval;
		sb.append(Num.formatNumber(somaticAf,3)); sb.append(",");
		sb.append(Num.formatNumber(germlineAf,3)); sb.append(",");
		String diff = Num.formatNumber(getAlleleFractionDifference(), 3);
		sb.append(diff); sb.append(",");
		sb.append(somAltRefGermAltRef[0]); sb.append("/"); sb.append(somAltRefGermAltRef[1]); sb.append(",");
		sb.append(somAltRefGermAltRef[2]); sb.append("/"); sb.append(somAltRefGermAltRef[3]); sb.append(",");
		String p = Num.formatNumberNoComma(Num.minus10log10(pvalue), 1);
		sb.append(p); sb.append(";");
		if (copyRatioRegion != null) {
			sb.append(fetchCopyRatioForVcf(copyRatioRegion)); sb.append(";");
		}
		return sb.toString();
	}
	
	private String fetchCopyRatioForVcf(RegionScoreText cr) throws IOException {
		
		String lg2T = null;
		String lg2N = null;
		String obs = null;
		
		//numOb=281;lg2Tum=-0.5847;lg2Norm=0.0088;genes=LINC00558,
		String[] fields = Misc.SEMI_COLON.split(cr.getText());
		for (String kv: fields) {
			if (kv.startsWith("numOb")) obs = kv;
			else if (kv.startsWith("lg2Tum")) lg2T = kv;
			else if (kv.startsWith("lg2Norm")) lg2N = kv;
		}
		if (lg2T==null || lg2N==null || obs==null) throw new IOException("FAILED to find all three lg2T,lg2N,obs from "+cr.getText());
		lg2T = lg2T.substring(lg2T.indexOf("=")+1);
		lg2N = lg2N.substring(lg2N.indexOf("=")+1);
		obs = obs.substring(obs.indexOf("=")+1);
		
		String lg2R = Num.formatNumber(cr.getScore(), 3);
		
		//LoHCR=lg2R,lg2T,lg2N,#Obs,Coor
		StringBuilder sb = new StringBuilder("LoHCR=");
		sb.append(lg2R); sb.append(",");
		sb.append(lg2T); sb.append(",");
		sb.append(lg2N); sb.append(",");
		sb.append(obs); sb.append(",");
		sb.append(cr.getStart()); sb.append("-");
		sb.append(cr.getStop());
		
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
			sb.append("\tSom  "+somAltRefGermAltRef[0]+"/"+somAltRefGermAltRef[1]+" "+somaticAf+"\n");
			sb.append("\tGerm "+somAltRefGermAltRef[2]+"/"+somAltRefGermAltRef[3]+" "+germlineAf+"\n");
			sb.append("\t-10Log10(varPVal) "+Num.minus10log10Float(pvalue)+ "\tafDiff " +getAlleleFractionDifference());
			if (copyRatioRegion!= null) {
				sb.append("\n\t"+getCopyRatioInfoScreen());
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		return sb.toString();
	}
	
	private String getCopyRatioInfoScreen() {
		StringBuilder sb = new StringBuilder();
		String score = Num.formatNumber(copyRatioRegion.getScore(), 3);
		String infoNoGene = copyRatioRegion.getText().substring(0, copyRatioRegion.getText().indexOf(";gene"));
		//numOb=1914;lg2Tum=0.3102;lg2Norm=-0.0115;genes=....
		sb.append("CopyRatio lg2Ratio="+score+";coor="+copyRatioRegion.getStart()+"-"+copyRatioRegion.getStop()+";"+infoNoGene);
		return sb.toString();
	}

	public int[] fetchSomAltRefGermAltRefCounts() throws IOException {
		if (somAltRefGermAltRef != null) return somAltRefGermAltRef;
		int somAlt;
		int somRef;
		int germAlt;
		int germRef;
		BaseCount somBC =  somaticBP.getSamples()[0];
		BaseCount germBC = germlineBP.getSamples()[0];
		
		if (isSnv) {
			char alt = vcfRecord[4].charAt(0);
			char ref = vcfRecord[3].charAt(0);
			somAlt = somBC.getSnvCount(alt);
			somRef = somBC.getSnvCount(ref);
			germAlt = germBC.getSnvCount(alt);
			germRef = germBC.getSnvCount(ref);
			
		}
		else {
			somRef = (int)somBC.getPassingReadCoverageSnv(); //don't want to include any Ds
			germRef = (int)germBC.getPassingReadCoverageSnv();
			if (isDeletion) {
				somAlt = somBC.getIndelCount('D');
				germAlt = germBC.getIndelCount('D');
			}
			else {
				somAlt = somBC.getIndelCount('I');
				germAlt = germBC.getIndelCount('I');
			}
		}
		
		//set AFs
		somAltRefGermAltRef = new int[] {somAlt, somRef, germAlt, germRef};
		
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

	public void setCopyRatioRegion(RegionScoreText copyRatioRegion) {
		this.copyRatioRegion= copyRatioRegion;
		
	}
}
