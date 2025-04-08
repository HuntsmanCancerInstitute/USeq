package edu.utah.seq.vcf.json.tempusv3;

import java.util.LinkedHashMap;
import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;

public class TempusV3Order {

	private String physician = null;
	private String accessionId = null;
	private String tempusOrderId = null;
	private String testCode = null;
	private String testDescription = null;
	private String referenceGenome = null;
	private String institution = null;
	
	/*
     "order": {
          "institution": "Huntsman Cancer Institute",
          "test": {
               "code": "XT.V4",
               "name": "xT",
               "description": "648 gene panel",
               "referenceGenome": "GRCh37/hg19"
          },
          "physician": "Conan Kinsey",
          "tempusOrderId": "25cnychi",
          "accessionId": "TL-25-MOCK"
     }
	 */
	public TempusV3Order(JSONObject object, TempusV3Json2Vcf tempusJson2Vcf) throws JSONException {
		JSONObject order = object.getJSONObject("order");
		institution = Json.getStringAttribute(order, "institution");
		physician = Json.getStringAttribute(order, "physician");
		accessionId = Json.getStringAttribute(order, "accessionId");
		tempusOrderId = Json.getStringAttribute(order, "tempusOrderId");
		JSONObject test = order.getJSONObject("test");
		testCode = Json.getStringAttribute(test, "code");
		testDescription = Json.getStringAttribute(test, "description");
		referenceGenome = Json.getStringAttribute(test, "referenceGenome");
		
		TempusV3Json2Vcf.add(physician, tempusJson2Vcf.physicians);
		TempusV3Json2Vcf.add(testCode, tempusJson2Vcf.testCodes);
		TempusV3Json2Vcf.add(testDescription, tempusJson2Vcf.testDescriptions);
	}
	
	/**Placed in VCF Header*/
	public void addMetaData(LinkedHashMap<String, String> meta) {
		meta.put("tempusPhysician", physician);
		meta.put("tempusTestCode", testCode);
		meta.put("tempusTestDescription", testDescription);
		meta.put("tempusAccessionId", accessionId);
		meta.put("tempusOrderId", tempusOrderId);
		meta.put("tempusReportingReferenceGenome", referenceGenome);
	}
	
	/**Put into the spreadsheet output.*/
	public void addAttributes(LinkedHashMap<String, String> meta) {
		meta.put("institution", institution);
		meta.put("physician", physician);
		meta.put("accessionId", accessionId);
		meta.put("testCode", testCode);
		meta.put("testDescription", testDescription);
		meta.put("reportingReferenceGenome", referenceGenome);
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
