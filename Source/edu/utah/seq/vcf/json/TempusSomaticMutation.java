package edu.utah.seq.vcf.json;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;
import util.gen.Misc;

public class TempusSomaticMutation {
	
	private String gene = null;
	private TempusSomaticVariant[] variants = null;
	

	
	/*
        "somaticPotentiallyActionableMutations": [
            {
                "gene": "RB1",
                "display": "RB1",
                "hgncId": "9884",
                "entrezId": "5925",
                "variants": [
                    {
                        "mutationEffect": "p.Q207*",
                        "HGVS.p": "p.Q207*",
                        "HGVS.pFull": "p.Gln207*",
                        "HGVS.c": "c.619C>T",
                        "transcript": "NM_000321",
                        "nucleotideAlteration": "13:48934164:C:T",
                        "referenceGenome": "GRCh37/hg19",
                        "allelicFraction": "77.01",
                        "variantDescription": "Stop gain - LOF",
                        "coverage": 1044,
                        "therapies": []
                    }
                ]
            },

	 */
	

	public TempusSomaticMutation(JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException {
		gene = Json.getStringAttribute(object, "gene");
		TempusJson2Vcf.add(gene, tempusJson2Vcf.somaticGenes);
		
		JSONArray vars = object.getJSONArray("variants");
		variants = new TempusSomaticVariant[vars.length()];
		for (int i=0; i< variants.length; i++) variants[i] = new TempusSomaticVariant(vars.getJSONObject(i), tempusJson2Vcf);
		

	}
	
	public static TempusSomaticMutation[] getMutations(JSONArray so, TempusJson2Vcf tempusJson2Vcf) throws JSONException{
		TempusSomaticMutation[] muts = new TempusSomaticMutation[so.length()];
		for (int i=0; i< muts.length; i++) muts[i] = new TempusSomaticMutation(so.getJSONObject(i), tempusJson2Vcf);
		return muts;
	}
}
