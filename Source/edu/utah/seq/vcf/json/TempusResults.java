package edu.utah.seq.vcf.json;

import java.util.ArrayList;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;
import util.gen.Misc;

public class TempusResults {
	
	private Double tumorMutationalBurden = null;
	private Double tumorMutationBurdenPercentile = null;
	private String msiStatus = null;
	private ArrayList<TempusVariant> variants = new ArrayList<TempusVariant>();
	
	/*
    "results": {
    //general
        "tumorMutationalBurden": 1.67,
        "tumorMutationBurdenPercentile": 30,
        "msiStatus": "stable",
   //somatic short vars
        "somaticPotentiallyActionableMutations": [],
        "somaticBiologicallyRelevantVariants": [],
        "somaticVariantsOfUnknownSignificance": [],
   //structural rearrangements
        "somaticPotentiallyActionableCopyNumberVariants": [], 
        "fusionVariants": [], - not seen so try but might fail
   //germline
        "inheritedRelevantVariants": [],
        "inheritedIncidentalFindings": [],
        "inheritedVariantsOfUnknownSignificance": [], dups with inheritedRelevantVariants
   //misc
        "ihcFindings": {},
        "rnaFindings": [],
        "lowCoverageAmplicons": []
    },
	 */
	public TempusResults(JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException {
		JSONObject results = object.getJSONObject("results");
		//messy tempus json output for tumorMutationBurden stuff
		String tmb = Json.forceGetString(results, "tumorMutationalBurden");
		if (tmb !=null) tumorMutationalBurden = Double.parseDouble(tmb);
		String tmbp = Json.forceGetString(results, "tumorMutationBurdenPercentile");
		if (tmbp != null) tumorMutationBurdenPercentile = Double.parseDouble(tmbp);
		msiStatus = Json.getStringAttribute(results, "msiStatus");
		
		//any somaticPotentiallyActionableMutations? these are special where the gene is separated from the variant
		if (results.has("somaticPotentiallyActionableMutations")) {
			JSONArray ja = results.getJSONArray("somaticPotentiallyActionableMutations");
			for (int i=0; i< ja.length(); i++){
				JSONObject geneVars = ja.getJSONObject(i);
				String gene = geneVars.getString("gene");
				if (geneVars.has("variants")){
					JSONArray vars = geneVars.getJSONArray("variants");
					for (int j=0; j< vars.length(); j++) variants.add( new TempusVariant("somaticPotentiallyActionableMutations", gene, vars.getJSONObject(j), tempusJson2Vcf) );
				}
			}
		}
		//any somaticVariantsOfUnknownSignificance
		if (results.has("somaticVariantsOfUnknownSignificance")) {
			JSONArray ja = results.getJSONArray("somaticVariantsOfUnknownSignificance");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("somaticVariantsOfUnknownSignificance", null, ja.getJSONObject(i), tempusJson2Vcf) );
		}
		//any somaticBiologicallyRelevantVariants, mix of cnv and short indels, maybe fusions?
		if (results.has("somaticBiologicallyRelevantVariants")) {
			JSONArray ja = results.getJSONArray("somaticBiologicallyRelevantVariants");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("somaticBiologicallyRelevantVariants", null, ja.getJSONObject(i), tempusJson2Vcf) );
		}
		//any somaticPotentiallyActionableCopyNumberVariants, mix of cnv and short indels, maybe fusions?
		if (results.has("somaticPotentiallyActionableCopyNumberVariants")) {
			JSONArray ja = results.getJSONArray("somaticPotentiallyActionableCopyNumberVariants");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("somaticPotentiallyActionableCopyNumberVariants", null, ja.getJSONObject(i), tempusJson2Vcf) );
		}
		//any fusionVariants, none seen so this might fail
		if (results.has("fusionVariants")) {
			JSONArray ja = results.getJSONArray("fusionVariants");
			if (ja.length() !=0) Misc.printErrAndExit("\nFound fusion variants!!! Contact Nix to fix parser\n");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("fusionVariants", null, ja.getJSONObject(i), tempusJson2Vcf) );
		}
		//any inheritedRelevantVariants
		if (results.has("inheritedRelevantVariants")) {
			JSONArray ja = results.getJSONArray("inheritedRelevantVariants");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("inheritedRelevantVariants", null, ja.getJSONObject(i), tempusJson2Vcf) );
		}
		//any inheritedIncidentalFindings
		if (results.has("inheritedIncidentalFindings")) {
			JSONArray ja = results.getJSONArray("inheritedIncidentalFindings");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("inheritedIncidentalFindings", null, ja.getJSONObject(i), tempusJson2Vcf) );
		}
		//any inheritedVariantsOfUnknownSignificance
		if (results.has("inheritedVariantsOfUnknownSignificance")) {
			JSONArray ja = results.getJSONArray("inheritedVariantsOfUnknownSignificance");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("inheritedVariantsOfUnknownSignificance", null, ja.getJSONObject(i), tempusJson2Vcf) );
		}
	}

	public Double getTumorMutationalBurden() {
		return tumorMutationalBurden;
	}

	public Double getTumorMutationBurdenPercentile() {
		return tumorMutationBurdenPercentile;
	}

	public String getMsiStatus() {
		return msiStatus;
	}

	public ArrayList<TempusVariant> getVariants() {
		return variants;
	}
}
