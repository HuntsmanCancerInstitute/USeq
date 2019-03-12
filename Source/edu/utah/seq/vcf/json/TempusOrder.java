package edu.utah.seq.vcf.json;

import org.json.JSONException;
import org.json.JSONObject;

import util.gen.Json;

public class TempusOrder {

	private String physician = null;
	private String accessionId = null;
	private String testCode = null;
	private String testDescription = null;
	
	/*
    "order": {
        "institution": "Huntsman Cancer Institute <c5d1a86d-e03a-4c91-9af2-be7c034ecb86>",
        "physician": "Ignacio Garrido-Laguna",
        "tempusOrder_id": "451fa6ea-3206-47f6-a921-925d6acb6244",
        "accessionId": "TL-18-03CFD6",
        "test": {
            "code": "XT",
            "name": "Targeted Panel",
            "description": "Panel of 595 genes (germline, tumor) + mRNA (tumor)"
        }
    },
	 */
	public TempusOrder(JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException {
		JSONObject order = object.getJSONObject("order");
		physician = Json.getStringAttribute(order, "physician");
		accessionId = Json.getStringAttribute(order, "accessionId");
		JSONObject test = order.getJSONObject("test");
		testCode = Json.getStringAttribute(test, "code");
		testDescription = Json.getStringAttribute(test, "description");
		
		TempusJson2Vcf.add(physician, tempusJson2Vcf.physicians);
		TempusJson2Vcf.add(testCode, tempusJson2Vcf.testCodes);
		TempusJson2Vcf.add(testDescription, tempusJson2Vcf.testDescriptions);
	}

	public String getPhysician() {
		return physician;
	}

	public String getAccessionId() {
		return accessionId;
	}

	public String getTestCode() {
		return testCode;
	}

	public String getTestDescription() {
		return testDescription;
	}
}
