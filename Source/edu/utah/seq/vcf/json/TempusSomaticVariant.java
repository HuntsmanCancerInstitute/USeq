package edu.utah.seq.vcf.json;

import org.json.JSONException;
import org.json.JSONObject;

import util.gen.Json;
import util.gen.Misc;

public class TempusSomaticVariant{
	private String mutationEffect = null;
	private String HGVSp = null;
	private String HGVSc = null;
	private String transcript = null;
	private String nucleotideAlteration = null;
	private String referenceGenome = null;
	private String allelicPercent = null;
	private String variantDescription = null;
	private Integer coverage = null;
	private static final String acceptedReferenceGenome = "GRCh37/hg19";
	
	public TempusSomaticVariant(JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException {
		mutationEffect = Json.getStringAttribute(object, "mutationEffect");
		HGVSp = Json.getStringAttribute(object, "HGVS.p");
		HGVSc = Json.getStringAttribute(object, "c.619C>T");
		transcript = Json.getStringAttribute(object, "transcript");
		nucleotideAlteration = Json.getStringAttribute(object, "nucleotideAlteration");
		referenceGenome = Json.getStringAttribute(object, "referenceGenome");
		if (referenceGenome.equals(acceptedReferenceGenome) == false) Misc.printErrAndExit("\nError: found different reference genome than what's accepted "+acceptedReferenceGenome +" vs "+referenceGenome);
		//messy Tempus json, not always a String
		String allelicPercent = Json.forceGetString(object, "allelicFraction");
		variantDescription = Json.getStringAttribute(object, "variantDescription");
		coverage = Json.getIntegerAttribute(object, "coverage");
		
		//increment histograms
		if (allelicPercent !=null) tempusJson2Vcf.tumorAF.count(Double.parseDouble(allelicPercent));
		if (coverage !=null) tempusJson2Vcf.tumorDP.count(coverage);
	}

	public String getMutationEffect() {
		return mutationEffect;
	}

	public String getHGVSp() {
		return HGVSp;
	}

	public String getHGVSc() {
		return HGVSc;
	}

	public String getTranscript() {
		return transcript;
	}

	public String getNucleotideAlteration() {
		return nucleotideAlteration;
	}

	public String getReferenceGenome() {
		return referenceGenome;
	}

	public String getAllelicPercent() {
		return allelicPercent;
	}

	public String getVariantDescription() {
		return variantDescription;
	}

	public Integer getCoverage() {
		return coverage;
	}
	
}