package edu.utah.seq.vcf.anno;

import java.io.IOException;
import java.util.ArrayList;

import util.gen.Gzipper;
import util.gen.Misc;
import util.gen.Num;

public class AnnotatedVcfParserDataLine {
	//Report
	String trimmedFileName = null;
	
	//Variant
	String variantId = null;
	String chr = null;
	String pos = null;
	String ref = null;
	String alt = null;
	String qual = null;
	String filter = null;
	
	String priorCallFreq = null;
	String bkz = null;
	String bkzAF = null;
	double varAlleleFreq = -1;
	int varUniOb = -1;
	int totalUniObDepth = -1;
	
	double varPopAlleleFreq = -1;
	
	//ANN
	boolean passesANN = false;
	ArrayList<AnnotatedGene> annoGenes = new ArrayList<AnnotatedGene>();
	
	//CLINVAR
	boolean passesCLINVAR = false;
	boolean rejectedCLINVAR = false;
	String clinSigCLNHGVS = null;
	String clinSig = null;
	String clinSigConf = null;
	
	//Splice
	boolean passesSplice = false;
	String spliceGeneName = null;
	double spliceScoreDiff = -1;
	
	
	//constructor
	public AnnotatedVcfParserDataLine(String trimmedFileName, String[] vcfFields) {
		this.trimmedFileName = trimmedFileName;
		//#CHROM POS ID REF ALT QUAL FILTER INFO ......
		//   0    1   2  3   4   5     6      7
		chr = vcfFields[0];
		pos = vcfFields[1];
		variantId = vcfFields[2];
		ref = vcfFields[3];
		alt = vcfFields[4];
		qual = vcfFields[5];
		filter = vcfFields[6];
	}

	public static final String headerSpreadSheet = "FileName\tIGVLink\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAlleleFreq\t"
			+ "UniOb\tReadDepth\tPriorCallFreq\tBKZ\tBKAF\tPopFreq\tPassClinvar\tClinHGVS\tClinSig\t"
			+ "ClinSigConf\tPassSpliceScan\tSpliceGene\tSpliceDiff";
	
	public void println(Gzipper sumarySpreadSheet) throws IOException {
		//calc varUniOb
		varUniOb = (int)Math.round(varAlleleFreq*(double)totalUniObDepth);
		ArrayList<String> al = new ArrayList<String>();
		//FileName
		al.add(trimmedFileName);
		//IGVLink
		al.add(getIGVLink());

		//CHROM
		al.add(chr);
		//POS
		al.add(pos);
		//ID
		al.add(variantId);
		//REF
		al.add(ref);
		//ALT
		al.add(alt);
		//QUAL
		al.add(qual);
		//FILTER
		al.add(filter);
		//AlleleFreq
		al.add(Num.formatNumber(varAlleleFreq, 3));
		//UniOb
		al.add(varUniOb+"");
		//ReadDepth
		al.add(totalUniObDepth+"");
		//PriorCallFreq
		if (priorCallFreq == null) al.add("null");
		else {
			String[] cf = Misc.COMMA.split(priorCallFreq);
			al.add(cf[0]+" ("+cf[1]+"/"+cf[2]+")");
		}
		
		//PoNBkz
		al.add(bkz);
		//PoNBkzAF
		al.add(bkzAF);
		//PopFreq
		if (varPopAlleleFreq != -1.0) al.add(Num.formatNumber(varPopAlleleFreq, 5));
		else al.add("null");
				
		//CLINVAR
		al.add(passesCLINVAR+"");
		al.add(clinSigCLNHGVS);
		al.add(clinSig);
		al.add(clinSigConf);
		
		//Splice
		al.add(passesSplice+"");
		al.add(spliceGeneName);
		if (spliceScoreDiff != -1.0)al.add(spliceScoreDiff+"");
		else al.add("null");
		
		String main = Misc.stringArrayListToString(al, "\t");
		
		//each ANN
		ArrayList<AnnotatedGene> pass = new ArrayList<AnnotatedGene>();
		ArrayList<AnnotatedGene> fail = new ArrayList<AnnotatedGene>();
		for (AnnotatedGene ag: annoGenes) {
			if (ag.passImpact) pass.add(ag);
			else fail.add(ag);
		}
		
		//any pass?
		if (pass.size() !=0) for (AnnotatedGene ag: pass) sumarySpreadSheet.println(main + ag.toString());
		//only print fail if no pass
		else if (fail.size() !=0) for (AnnotatedGene ag: fail) sumarySpreadSheet.println(main + ag.toString());
		else sumarySpreadSheet.println(main);
		sumarySpreadSheet.println();
	}
	
	public String getIGVLink() {
		//create IGV link  =HYPERLINK('http://localhost:60151/goto?locus=17:7379035-7395388','1')
		StringBuilder sb = new StringBuilder();
		int position = Integer.parseInt(pos);
		sb.append("=HYPERLINK(\"http://localhost:60151/goto?locus=");
		sb.append(chr); 
		sb.append(":");
		int s = position-150;
		if (s < 0) s = 0;
		sb.append(s); sb.append("-");
		sb.append(position+150); sb.append("\",\"");
		sb.append(chr); sb.append(":");
		sb.append(position);
		sb.append("\")");
		return sb.toString();
	}
	
	public static final String legend = "\nColumn Descriptions:\n"+
			"FileName\tTrimmed name of the parsed vcf file\n"+
			"IGVLink\tClicking moves IGV to the variant location\n"+
			"CHROM\tVcf record chromosome\n"+
			"POS\tVcf record position\n"+
			"ID\tVcf record ID\n"+
			"REF\tVcf record reference base\n"+
			"ALT\tVcf record alternative base\n"+
			"QUAL\tVcf record quality store, larger the better, if from VCFBkz, should be > 3\n"+
			"FILTER\tVcf record filter field, anything present other than PASS or . Indicates a potential problem with the variant call\n"+
			"AlleleFreq\tVariant allele frequency\n"+
			"UniOb\tNumber of  unique observations of the variant, e.g. no duplicate or overlapping paired alignments\n"+
			"ReadDepth\tTotal number of unique observations\n"+
			"PriorCallFreq\tFraction of analyzed datasets with this variant, lower the better, if greater than say 5% this may be an artifact.\n"+
			"BKZ\tSmallest AF z-score calculated from background AFs over effected bases. Values < ~4 are suspicous and indicate non reference observations in the background samples, e.g. a noisy background and potential artifactual call.\n"+
			"BKAF\tSorted list (largest to smallest) of background non-reference AFs used to calculate the BKZ.\n"+
			"PopFreq\tMaximum observed population frequency, fractions > 0.01 are indicative of a common variant\n"+
			"PassClinvar\tBoolean indicating whether this record passes the CLINVAR include and exclude criteria\n"+
			"ClinHGVS\tCLINVARs HGVS notation for this variant\n"+
			"ClinSig\tCLINVAR's significance annotation\n"+
			"ClinSigConf\tCLINVAR's conflicting interpretation tabulation\n"+
			"PassSpliceScan\tBoolean indicating whether this record passes the VCFSpliceScanner criteria\n"+
			"SpliceGene\tSpliceScanner impacted gene(s)\n"+
			"SpliceDiff\tDifference in MaxEntScan scores between the reference and variant splice sequences, larger the more significant.\n"+
			"PassAnn\tBoolean indicating whether this record passes the ANN criteria\n"+
			"Gene\tANN impacted gene name\n"+
			"TranscriptId\tANN impacted transcript\n"+
			"Annotation\tANN annotation\n"+
			"Impact\tANN impact\n"+
			"cDot\tANN HGVS.c\n"+
			"pDot\tANN HGVS.p\n"+
			"pPos\tANN Affected amino acid vs the total length";
	
}
