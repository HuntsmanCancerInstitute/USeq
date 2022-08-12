package edu.utah.seq.run.tempus;

import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;

public class TempusJsonOrder {

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
    
       "order": {
        "institution": "Huntsman Cancer Institute",
        "physician": "Kevin Jones",
        "tempusOrderId": "22cxjvwi",
        "accessionId": "TL-22-JAHQRJV8",
        "test": {
            "code": "XT.V4",
            "name": "Targeted Panel",
            "description": "Panel of 648 genes (germline, tumor) + mRNA (tumor)"
        }
	 */
	public TempusJsonOrder(JSONObject main) throws JSONException {
		JSONObject order = main.getJSONObject("order");
		physician = Json.getStringAttribute(order, "physician");
		accessionId = Json.getStringAttribute(order, "accessionId");
		JSONObject test = order.getJSONObject("test");
		testCode = Json.getStringAttribute(test, "code");
		fixTestCode();
		testDescription = Json.getStringAttribute(test, "description");
	}

	private void fixTestCode() {
		if (testCode.contains(".") == false) testCode = testCode+".V1";
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
