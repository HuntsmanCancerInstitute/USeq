package edu.utah.seq.vcf;

import java.util.ArrayList;
import java.util.HashMap;

//##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
//##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
//##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
//##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
//##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
//##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
//##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
//##INFO=<ID=Dels,Number=1,Type=Float,Description="Fraction of Reads Containing Spanning Deletions">
//##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
//##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
//##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
//##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
//##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
//##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
//##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
//##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
//##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
//##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
//##INFO=<ID=RPA,Number=.,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
//##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
//##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
//##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
//##INFO=<ID=VQSLOD,Number=1,Type=Float,Description="Log odds ratio of being a true variant versus being false under the trained gaussian mixture model">
//##INFO=<ID=culprit,Number=1,Type=String,Description="The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out">

public class VCFInfo {	
	private HashMap<String, String> hashInfoString = new HashMap<String,String>() {{
		put("RU",null);
		put("culprit",null);
		put("DB",null);
		put("DS",null);
		put("STR",null);
		put("BaseQRankSum",null);
		put("Dels",null);
		put("FS",null);
		put("HaplotypeScore",null);
		put("InbreedingCoeff",null);
		put("MQ",null);
		put("MQRankSum",null);
		put("QD",null);
		put("ReadPosRankSum",null);
		put("VQSLOD",null);
		put("AN",null);
		put("DP",null);
		put("END",null);
		put("MQ0",null);
		put("AF",null);
		put("MLEAF",null);
		put("MLEAC",null);
		put("AC",null);
		put("RPA",null);
	}};
	
	private String infoString = null;
	

	public VCFInfo() {
		// TODO Auto-generated constructor stub
	}
	
	public void parseInfoGatk(String infoString) {
		String[] entries = infoString.split(";");
		for (String e: entries) {
			//create simpleInfoString
			this.infoString = infoString;
					
			//Get key-par, account for flag values
			String[] keyValue = e.split("=");
			String key = keyValue[0];
			String value = null;
			
			if (keyValue.length == 1) {
				value = "true";
			} else {
				value = keyValue[1];
			}
			
			if (hashInfoString.containsKey(key)) {
				hashInfoString.put(key,value);
			} else {
				System.out.println("Did not recognize the info entry, skipping..." + keyValue[0]);
				continue;
			}
		}
	}
	
	public void overwriteInfoString(String infoString) {
		this.infoString = infoString;
	}
	
	public String getInfoString() {
		return this.infoString;
	}
	
	public void addInfo(String identifier,String info) {
		this.hashInfoString.put(identifier, info);
		this.infoString += ";" + identifier + "=" + info;
	}
	
	
	
	

}
