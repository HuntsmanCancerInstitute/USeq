package edu.utah.seq.vcf.json;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;

public class TempusSpecimen {
	
	private String tempusSampleId = null;
	private String collectionDate = null;
	private String sampleCategory = null;
	private String sampleSite = null;
	private String sampleType = null;
	private String notes = null;
	private Integer tumorPercentage = null;
	
	/*
    "specimens": [
        {
            "tempusSampleId": "6e3be2d4-7a75-4d00-84f2-1380e9a56e64",
            "collectionDate": "2018-09-24T05:00:00",
            "receiptDate": "2018-09-28T19:26:15.431",
            "sampleCategory": "normal",
            "sampleSite": null,
            "sampleType": "Blood",
            "notes": "",
            "institutionData": {
                "caseId": "",
                "blockId": "",
                "tumorPercentage": null
            }
        },
        {
            "tempusSampleId": "f6013666-6655-4f85-9414-9bc9f3d1237c",
            "collectionDate": "2018-08-15T05:00:00",
            "receiptDate": "2018-10-01T18:55:59.375",
            "sampleCategory": "tumor",
            "sampleSite": "Gallbladder",
            "sampleType": "Biopsy",
            "notes": "",
            "institutionData": {
                "caseId": "Castleview Hospital S18-2014",
                "blockId": "C1",
                "tumorPercentage": 80
            }
        }
    ],
	 */
	public TempusSpecimen(JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException {
		tempusSampleId = Json.getStringAttribute(object, "tempusSampleId");
		collectionDate = Json.getStringAttribute(object, "collectionDate");
		sampleCategory = Json.getStringAttribute(object, "sampleCategory");
		sampleSite = Json.getStringAttribute(object, "sampleSite");
		sampleType = Json.getStringAttribute(object, "sampleType");
		if (object.has("institutionData")){
			JSONObject id = object.getJSONObject("institutionData");
			tumorPercentage = Json.getIntegerAttribute(id, "tumorPercentage");
			if (tumorPercentage != null) tempusJson2Vcf.tumorPercentages.count(tumorPercentage);
		}
		
		TempusJson2Vcf.add(sampleCategory, tempusJson2Vcf.sampleCategories);
		TempusJson2Vcf.add(sampleSite, tempusJson2Vcf.sampleSites);
		TempusJson2Vcf.add(sampleType, tempusJson2Vcf.sampleTypes);
		

	}
	
	public static TempusSpecimen[] getSpecimens(JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException{
		JSONArray so = object.getJSONArray("specimens");
		TempusSpecimen[] specimens = new TempusSpecimen[so.length()];
		for (int i=0; i< specimens.length; i++) specimens[i] = new TempusSpecimen(so.getJSONObject(i), tempusJson2Vcf);
		return specimens;
	}

	public String getTempusSampleId() {
		return tempusSampleId;
	}

	public String getCollectionDate() {
		return collectionDate;
	}

	public String getSampleCategory() {
		return sampleCategory;
	}

	public String getSampleSite() {
		return sampleSite;
	}

	public String getSampleType() {
		return sampleType;
	}

	public String getNotes() {
		return notes;
	}

	public Integer getTumorPercentage() {
		return tumorPercentage;
	}

}
