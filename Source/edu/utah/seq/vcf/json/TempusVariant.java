package edu.utah.seq.vcf.json;

import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;
import util.gen.Misc;

public class TempusVariant{
	
	//many of these fields will remain null
	private String variantSource = null;
	private String gene = null;
	private String transcript = null;
	private String mutationEffect = null;
	private String cHGVS = null;
	private String pHGVS = null;
	private String pFullHGVS = null;
	private String referenceGenome = null;
	private String variantType = null;
	private String variantDescription = null;
	private String nucleotideAlteration = null;
	private String chromosome = null;
	private String pos = null;
	private String ref = null;
	private String alt = null;
	private String allelicFraction = null;
	private String coverage = null;
	private String copyNumber = null;
	private String clinicalSignificance = null;
	private String disease = null;
	private String structuralVariant = null;
	private String fusionType = null;
	private String gene3 = null;
	private String gene5 = null;
		
	/**Object to represent a Tempus variant, SNV/ INDEL/ CNV/ Fusion*/
	public TempusVariant(String variantSource, String geneName, JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException {
		this.variantSource = variantSource;
		
		//attempt to parse gene, if null pull from constructor
		gene = Json.forceGetString(object, "gene");
		if (gene == null) gene = geneName;

		mutationEffect = Json.forceGetString(object, "mutationEffect");
		transcript = Json.forceGetString(object, "transcript");
		cHGVS = Json.forceGetString(object, "HGVS.c");
		pHGVS = Json.forceGetString(object, "HGVS.p");
		pFullHGVS = Json.forceGetString(object, "HGVS.pFull");
		referenceGenome = Json.forceGetString(object, "referenceGenome");

		variantType = Json.forceGetString(object, "variantType");
		variantDescription = Json.forceGetString(object, "variantDescription");

		nucleotideAlteration = Json.forceGetString(object, "nucleotideAlteration");
		chromosome = Json.forceGetString(object, "chromosome");
		pos = Json.forceGetString(object, "pos");
		ref = Json.forceGetString(object, "ref");
		alt = Json.forceGetString(object, "alt");

		allelicFraction = Json.forceGetString(object, "allelicFraction");
		coverage = Json.forceGetString(object, "coverage");

		copyNumber = Json.forceGetString(object, "copyNumber");

		clinicalSignificance = Json.forceGetString(object, "clinicalSignificance");
		disease = Json.forceGetString(object, "disease");

		structuralVariant = Json.forceGetString(object, "structuralVariant");
		fusionType = Json.forceGetString(object, "fusionType");
		gene3 = Json.forceGetString(object, "gene3");
		gene5 = Json.forceGetString(object, "gene5");
		
		//set coordinates from nucleotideAlteration?
		if (nucleotideAlteration != null){
			String[] t = Misc.COLON.split(nucleotideAlteration);
			chromosome = t[0];
			pos = t[1];
			ref = t[2];
			alt = t[3];
		}
		
		//increment histograms
		if (allelicFraction !=null) tempusJson2Vcf.tumorAF.count(Double.parseDouble(allelicFraction));
		if (coverage !=null) tempusJson2Vcf.tumorDP.count(Double.parseDouble(coverage));
	}
	
	/**Returns all non null values*/
	public String toString(){
		StringBuilder sb = new StringBuilder();
		if (variantSource != null) sb.append("variantSource\t"+variantSource+ "\n");
		if (gene != null) sb.append("gene\t"+gene+ "\n");
		if (transcript != null) sb.append("transcript\t"+transcript+ "\n"); 
		if (mutationEffect != null) sb.append("mutationEffect\t"+mutationEffect+ "\n"); 
		if (cHGVS != null) sb.append("cHGVS\t"+cHGVS+ "\n"); 
		if (pHGVS != null) sb.append("pHGVS\t"+pHGVS+ "\n"); 
		if (pFullHGVS != null) sb.append("pFullHGVS\t"+pFullHGVS+ "\n"); 
		if (referenceGenome != null) sb.append("referenceGenome\t"+referenceGenome+ "\n"); 
		if (variantType != null) sb.append("variantType\t"+variantType+ "\n"); 
		if (variantDescription != null) sb.append("variantDescription\t"+variantDescription+ "\n"); 
		if (nucleotideAlteration != null) sb.append("nucleotideAlteration\t"+nucleotideAlteration+ "\n"); 
		if (chromosome != null) sb.append("chromosome\t"+chromosome+ "\n"); 
		if (pos != null) sb.append("pos\t"+pos+ "\n"); 
		if (ref != null) sb.append("ref\t"+ref+ "\n"); 
		if (alt != null) sb.append("alt\t"+alt+ "\n"); 
		if (allelicFraction != null) sb.append("allelicFraction\t"+allelicFraction+ "\n"); 
		if (coverage != null) sb.append("coverage\t"+coverage+ "\n"); 
		if (copyNumber != null) sb.append("copyNumber\t"+copyNumber+ "\n"); 
		if (clinicalSignificance != null) sb.append("clinicalSignificance\t"+clinicalSignificance+ "\n"); 
		if (disease != null) sb.append("disease\t"+disease+ "\n"); 
		if (structuralVariant != null) sb.append("structuralVariant\t"+structuralVariant+ "\n"); 
		if (fusionType != null) sb.append("fusionType\t"+fusionType+ "\n"); 
		if (gene3 != null) sb.append("gene3\t"+gene3+ "\n"); 
		if (gene5 != null) sb.append("gene5\t"+gene5+ "\n"); 
		return sb.toString();
	}


	public String getVariantSource() {
		return variantSource;
	}

	public String getGene() {
		return gene;
	}

	public String getMutationEffect() {
		return mutationEffect;
	}

	public String getTranscript() {
		return transcript;
	}

	public String getcHGVS() {
		return cHGVS;
	}

	public String getpHGVS() {
		return pHGVS;
	}

	public String getpFullHGVS() {
		return pFullHGVS;
	}

	public String getReferenceGenome() {
		return referenceGenome;
	}

	public String getVariantType() {
		return variantType;
	}

	public String getVariantDescription() {
		return variantDescription;
	}

	public String getNucleotideAlteration() {
		return nucleotideAlteration;
	}

	public String getChromosome() {
		return chromosome;
	}

	public String getPos() {
		return pos;
	}

	public String getRef() {
		return ref;
	}

	public String getAlt() {
		return alt;
	}

	public String getAllelicFraction() {
		return allelicFraction;
	}

	public String getCoverage() {
		return coverage;
	}

	public String getCopyNumber() {
		return copyNumber;
	}

	public String getClinicalSignificance() {
		return clinicalSignificance;
	}

	public String getDisease() {
		return disease;
	}

	public String getStructuralVariant() {
		return structuralVariant;
	}

	public String getFusionType() {
		return fusionType;
	}

	public String getGene3() {
		return gene3;
	}

	public String getGene5() {
		return gene5;
	}
	
}