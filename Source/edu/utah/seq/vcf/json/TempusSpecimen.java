package edu.utah.seq.vcf.json;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;
import util.gen.Misc;

public class TempusSpecimen {
	
	private String tempusSampleId = null;
	private String collectionDate = null;
	private String sampleCategory = null;
	private String sampleSite = null;
	private String sampleType = null;
	private String notes = null;
	private String caseId = null;
	//private String caseSourceId = null;
	private String blockId = null;
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
			//caseId
			caseId = Json.getStringAttribute(id, "caseId");
			//don't do this anymore
			/*
			String unParsed = Json.getStringAttribute(id, "caseId");
			if (unParsed != null && unParsed.equals("")== false) {
				String[] split = Misc.WHITESPACE.split(unParsed);
				caseId = split[split.length-1].trim();
				caseSourceId = unParsed.replace(caseId, "").trim();
			}*/
			//blockId
			blockId = Json.getStringAttribute(id, "blockId");
			//tumorPercentage
			tumorPercentage = Json.getIntegerAttribute(id, "tumorPercentage");
			if (tumorPercentage != null) tempusJson2Vcf.tumorPercentages.count(tumorPercentage);
		}
		
		TempusJson2Vcf.add(sampleCategory, tempusJson2Vcf.sampleCategories);
		TempusJson2Vcf.add(sampleSite, tempusJson2Vcf.sampleSites);
		TempusJson2Vcf.add(sampleType, tempusJson2Vcf.sampleTypes);
		

	}
	
	public void addMetaData(LinkedHashMap<String, String> meta, int id) {
		if (sampleCategory != null) meta.put("tempusSampleCategory_"+id, sampleCategory);
		if (sampleSite != null) meta.put("tempusSampleSite_"+id, sampleSite);
		if (sampleType != null) meta.put("tempusSampleType_"+id, sampleType);
		if (tumorPercentage != null) meta.put("tempusTumorPercentage_"+id, tumorPercentage.toString());	
		//if (caseId != null) meta.put("tempusCaseId_"+id, caseId +" : "+caseSourceId);
		if (caseId != null) meta.put("tempusCaseId_"+id, caseId);
		if (blockId != null) meta.put("tempusBlockId_"+id, blockId);
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

	public static void addAttributes(LinkedHashMap<String, String> reportAttributes, TempusSpecimen[] specimens) {
		ArrayList<String> al = new ArrayList<String>();
		Integer tp = null;
		String caseId = null;
		String caseSourceId = null;
		String blockId = null;
		for (int i=0; i< specimens.length; i++ ) {
			String type = specimens[i].getSampleType();
			String site = specimens[i].getSampleSite();
			String tn = specimens[i].getSampleCategory();
			if (site == null && type == null) continue;
			if (site == null) al.add(tn+" : "+type);
			else if (type == null) al.add(tn+" : "+site);
			else al.add(tn+" : "+type+" - "+site);
			//tumor percentage
			if (specimens[i].getTumorPercentage() != null) tp = specimens[i].getTumorPercentage();
			//caseId
			if (specimens[i].getCaseId() != null) {
				caseId = specimens[i].getCaseId();
				//caseSourceId = specimens[i].getCaseSourceId();
				blockId = specimens[i].getBlockId();
			}
		}
		String sum = Misc.stringArrayListToString(al, "; ");
		reportAttributes.put("specimines", sum);
		if (tp != null) reportAttributes.put("tumorPercentage", tp.toString());
		//should only be one of these
		if (caseId != null) {
			reportAttributes.put("caseId", caseId);
			reportAttributes.put("blockId", blockId);
			reportAttributes.put("caseIdSource", caseSourceId);
		}
	}

	public String getCaseId() {
		return caseId;
	}

	/*public String getCaseSourceId() {
		return caseSourceId;
	}*/

	public String getBlockId() {
		return blockId;
	}

}
