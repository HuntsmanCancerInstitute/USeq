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
		 * C|missense_variant&splice_region_variant|MODERATE|BRAF|BRAF|transcript|NM_001374258.1|protein_coding|16/20|c.1862A>G|p.Asn621Ser|2088/9807|1862/2424|621/807||,C|missense_variant&splice_region_variant|MODERATE|BRAF|BRAF|transcript|NM_004333.6|protein_coding|15/18|c.1742A>G|p.Asn581Ser|1968/6459|1742/2301|581/766||;ALLELEID=174177;CLNDISDB=Human_Phenotype_Ontology:HP:0030078,MONDO:MONDO:0005061,MeSH:D000077192,MedGen:C0152013|Human_Phenotype_Ontology:HP:0030358,MONDO:MONDO:0005233,MeSH:D002289,MedGen:C0007131;CLNDN=Lung_adenocarcinoma|Non-small_cell_lung_carcinoma;CLNHGVS=NC_000007.14:g.140753393T>C;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_pathogenic;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=ClinGen:CA180747|UniProtKB:P15056#VAR_040393;GENEINFO=BRAF:673;MC=SO:0001583|missense_variant;ORIGIN=2;RS=121913370;VCFSS=BRAF:G5S,140753394,CAGATATAT,CAGGTATAT,-0.3,7.88,7.88;MUTATION_EFFECT=Gain-of-function;ONCOGENIC=Oncogenic	GT:DP:FDP:SDP:SUBDP:AU:CU:GU:TU:AF	./.:75:0:0:0:0,0:0,0:0,0:75,77:0	./.:193:0:0:0:0,0:5,5:0,0:188,193:0.0259

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
		if (transcriptMatch != null) {
			sb.append(transcriptMatch); sb.append("\t");
		}
		sb.append(annotation); sb.append("\t");
		sb.append(impact); sb.append("\t");
		sb.append(cDot); sb.append("\t");
		sb.append(pDot); sb.append("\t");
		sb.append(pPos);
		return sb.toString();
	}

}
