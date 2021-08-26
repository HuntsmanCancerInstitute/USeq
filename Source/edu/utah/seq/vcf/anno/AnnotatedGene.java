package edu.utah.seq.vcf.anno;

public class AnnotatedGene {
	
	boolean passImpact = false;
	String annotation = null;
	String impact = null;
	String geneName = null;
	String transcriptId = null;
	String pDot = null;
	String cDot = null;
	String pPos = null;
	String transcriptMatch = null;

	public AnnotatedGene(boolean passImpact, String[] ann) {
		this.passImpact = passImpact;

		/*
		 * 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
		 *    0          1                2               3          4            5				6				7			   8	   9	   10			11							12					13
		 * GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20),
		 * GC|frameshift_variant|HIGH|SPEN|ENSG00000065526|transcript|ENST00000375759|protein_coding|11/15|c.9729dupC|p.Thr3244fs|9934/12232|9730/10995|3244/3664||INFO_REALIGN_3_PRIME;LOF=(SPEN|ENSG00000065526|5|0.20)
		 */
		annotation = ann[1];
		impact = ann[2];
		geneName = ann[3];
		if (ann.length >6) transcriptId = ann[6];
		if (ann.length >9) cDot = ann[9];
		if (ann.length >10) pDot = ann[10];
		if (ann.length >13) pPos = ann[13];
	}
	
	public static final String headerSpreadSheet = "\tPassAnn\tGene\tTranscriptId\tAnnotation\tImpact\tcDot\tpDot\tpPos";
	public static final String headerSpreadSheetMatch = "\tPassAnn\tGene\tTranscriptId\tTranscriptMatch\tAnnotation\tImpact\tcDot\tpDot\tpPos";
	
	public String toString() {
		StringBuffer sb = new StringBuffer("\t");
		sb.append(passImpact); sb.append("\t");
		sb.append(geneName); sb.append("\t");
		sb.append(transcriptId); sb.append("\t");
		if (transcriptMatch != null) sb.append(transcriptMatch); sb.append("\t");
		sb.append(annotation); sb.append("\t");
		sb.append(impact); sb.append("\t");
		sb.append(cDot); sb.append("\t");
		sb.append(pDot); sb.append("\t");
		sb.append(pPos);
		return sb.toString();
	}

}
