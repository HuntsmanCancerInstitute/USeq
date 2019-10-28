package edu.utah.seq.vcf.json;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import util.bio.annotation.Bed;
import util.gen.IO;
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
	public TempusResults(JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException, IOException {
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
			int numVars = 0;
			for (int i=0; i< ja.length(); i++){
				JSONObject geneVars = ja.getJSONObject(i);
				String gene = geneVars.getString("gene");
				if (geneVars.has("variants")){
					JSONArray vars = geneVars.getJSONArray("variants");
					for (int j=0; j< vars.length(); j++) {
						variants.add( new TempusVariant("somaticPotentiallyActionableMutations", gene, vars.getJSONObject(j), tempusJson2Vcf) );
						numVars++;
					}
				}
			}
			tempusJson2Vcf.setWorkingNumSomaticPotentiallyActionableMutations(numVars);
		}
		//any somaticVariantsOfUnknownSignificance
		if (results.has("somaticVariantsOfUnknownSignificance")) {
			JSONArray ja = results.getJSONArray("somaticVariantsOfUnknownSignificance");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("somaticVariantsOfUnknownSignificance", null, ja.getJSONObject(i), tempusJson2Vcf) );
			tempusJson2Vcf.setWorkingNumSomaticVariantsOfUnknownSignificance(ja.length());
		}
		//any somaticBiologicallyRelevantVariants, mix of cnv and short indels, maybe fusions?
		if (results.has("somaticBiologicallyRelevantVariants")) {
			JSONArray ja = results.getJSONArray("somaticBiologicallyRelevantVariants");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("somaticBiologicallyRelevantVariant", null, ja.getJSONObject(i), tempusJson2Vcf) );
			tempusJson2Vcf.setWorkingNumSomaticBiologicallyRelevantVariants(ja.length());
		}
		//any somaticPotentiallyActionableCopyNumberVariants, mix of cnv and short indels, maybe fusions?
		if (results.has("somaticPotentiallyActionableCopyNumberVariants")) {
			JSONArray ja = results.getJSONArray("somaticPotentiallyActionableCopyNumberVariants");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("somaticPotentiallyActionableCopyNumberVariant", null, ja.getJSONObject(i), tempusJson2Vcf) );
			tempusJson2Vcf.setWorkingNumSomaticPotentiallyActionableCopyNumberVariants(ja.length());
		}
		//any fusionVariants, none seen so this might fail
		if (results.has("fusionVariants")) {
			/*  These are just so far large chromosome rearrangements, no gene specific coordinates
			JSONArray ja = results.getJSONArray("fusionVariants");
			if (ja.length() !=0) Misc.printErrAndExit("\nFound fusion variants!!! Contact Nix to fix parser\n");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("fusionVariant", null, ja.getJSONObject(i), tempusJson2Vcf) );
			tempusJson2Vcf.setWorkingNumFusionVariants(ja.length());
			*/
		}
		//any inheritedRelevantVariants
		if (results.has("inheritedRelevantVariants")) {
			JSONArray ja = results.getJSONArray("inheritedRelevantVariants");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("inheritedRelevantVariant", null, ja.getJSONObject(i), tempusJson2Vcf) );
			tempusJson2Vcf.setWorkingNumInheritedRelevantVariants(ja.length());
		}
		//any inheritedIncidentalFindings
		if (results.has("inheritedIncidentalFindings")) {
			//need to silently trap a potential json issue with Tempus doc, their inserting a string instead of an object array when the patient elects to receive no germline info
			try {
				JSONArray ja = results.getJSONArray("inheritedIncidentalFindings");
				for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("inheritedIncidentalFinding", null, ja.getJSONObject(i), tempusJson2Vcf) );
				tempusJson2Vcf.setWorkingNumInheritedIncidentalFindings(ja.length());
			} catch (JSONException e) {}
		}
		//any inheritedVariantsOfUnknownSignificance
		if (results.has("inheritedVariantsOfUnknownSignificance")) {
			JSONArray ja = results.getJSONArray("inheritedVariantsOfUnknownSignificance");
			for (int i=0; i<ja.length(); i++) variants.add( new TempusVariant("inheritedVariantsOfUnknownSignificance", null, ja.getJSONObject(i), tempusJson2Vcf) );
			tempusJson2Vcf.setWorkingNumInheritedVariantsOfUnknownSignificance(ja.length());
		}
		
		//cnvs processing
		addCoordinatesToCNVs(tempusJson2Vcf.getCnvGeneNameBed(), tempusJson2Vcf.getFasta());
	}
	
	private void addCoordinatesToCNVs(HashMap<String, Bed> cnvGeneNameBed, IndexedFastaSequenceFile fasta) throws IOException {
		if (cnvGeneNameBed == null) return;
		for (TempusVariant tv: variants) {
			if (tv.getChromosome() == null) {
				if ((tv.getVariantType()!= null && tv.getVariantType().contains("CNV")) || 
						(tv.getVariantDescription()!=null && tv.getVariantDescription().contains("Copy"))) tv.addCnvInfo(cnvGeneNameBed, fasta);
			}
		}
		
	}

	public void addMetaData(LinkedHashMap<String, String> meta) {
		if (tumorMutationalBurden != null) meta.put("tempusTumorMutationalBurden", tumorMutationalBurden.toString());
		if (tumorMutationBurdenPercentile != null) meta.put("tempusTumorMutationBurdenPercentile", tumorMutationBurdenPercentile.toString());
		if (msiStatus != null) meta.put("tempusMsiStatus", msiStatus);
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
