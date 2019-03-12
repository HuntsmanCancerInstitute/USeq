package edu.utah.seq.vcf.json;

import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;

public class TempusResults {
	
	private Double tumorMutationalBurden = null;
	private Double tumorMutationBurdenPercentile = null;
	private String msiStatus = null;
	private TempusSomaticMutation[] somaticPotentiallyActionableMutations = null;
	
	
	/*
    "results": {
    //general
        "tumorMutationalBurden": 1.67,
        "tumorMutationBurdenPercentile": 30,
        "msiStatus": "stable",
        
       	//somatic
        "somaticPotentiallyActionableMutations": [],
        "somaticPotentiallyActionableCopyNumberVariants": [],
        "somaticBiologicallyRelevantVariants": [],
        "somaticVariantsOfUnknownSignificance": [],
        "fusionVariants": [],
        
        //germline
        "inheritedRelevantVariants": [],
        "inheritedIncidentalFindings": [],
        "inheritedVariantsOfUnknownSignificance": [],
        
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
		//any somaticPotentiallyActionableMutations?
		if (results.has("somaticPotentiallyActionableMutations")) somaticPotentiallyActionableMutations = TempusSomaticMutation.getMutations(results.getJSONArray("somaticPotentiallyActionableMutations"), tempusJson2Vcf);
	}

	
}
