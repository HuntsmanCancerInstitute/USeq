package edu.utah.seq.vcf.json.tempusv3;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;
import util.gen.Misc;

public class TempusV3Specimen {
	
	//institution data
	private String caseId = null;
	private String blockId = null;
	private Integer tumorPercentage = null;
	
	private String sampleCategory = null;
	private String collectionDate = null;
	private String notes = null;
	
	private String primarySampleSite = null; //new
	private String normalSampleSite = null; //new
	
	//diagnosis and specimine site codes and txt, these only populate if the map is created
	private String tempusIcdOCodeMorphology = null;
	private String tempusIcdOCodeTopography = null;
	private String tempusIcd10Code = null;
	private String originPathLabDiagnosis = null;
	
	private LinkedHashSet<String> tempusIcdOTxtMorphology = new LinkedHashSet<String>();
	private LinkedHashSet<String> tempusIcdOTxtTopography = new LinkedHashSet<String>();
	private LinkedHashSet<String> tempusIcd10Txt = new LinkedHashSet<String>();
	


	/*

    NEW:
         "specimens": [
          {
               "institutionData": {
                    "blockId": "A1",
                    "caseId": "Intermountain Medical Center - Pathology CLS-25-010466",
                    "tumorPercentage": 20
               },
               "notes": null,
               "primarySampleSite": "Lung, left lower lobe",
               "receiptDate": "2025-02-22T14:00:00",
               "sampleCategory": "tumor",
               "diagnosis": {
                    "tempusIcdOCodeMorphology": "8140/3",
                    "tempusIcdOCodeTopography": "C25.9",
                    "originPathLabDiagnosis": "Adenocarcinoma with enteric morphology",
                    "tempusIcd10Code": "C25,C25.9"
               },
               "normalSampleSite": null,
               "collectionDate": "2025-01-31T00:00:00"
          },
          {
               "institutionData": {
                    "blockId": null,
                    "caseId": null,
                    "tumorPercentage": null
               },
               "notes": null,
               "primarySampleSite": null,
               "receiptDate": "2025-02-21T14:00:00",
               "sampleCategory": "normal",
               "diagnosis": {
                    "tempusIcdOCodeMorphology": null,
                    "tempusIcdOCodeTopography": null,
                    "originPathLabDiagnosis": null,
                    "tempusIcd10Code": null
               },
               "normalSampleSite": "Blood",
               "collectionDate": "2025-02-19T00:00:00"
          }
     ],

    Example of Heme:
    "specimens": [{
          "institutionData": {
               "blockId": null,
               "caseId": null,
               "tumorPercentage": null
          },
          "notes": null,
          "primarySampleSite": "Peripheral Blood",
          "receiptDate": "2025-02-26T14:00:00",
          "sampleCategory": "heme",
          "diagnosis": {
               "tempusIcdOCodeMorphology": "9823/3",
               "tempusIcdOCodeTopography": "C77.9",
               "originPathLabDiagnosis": "Chronic lymphocytic leukemia ",
               "tempusIcd10Code": "C81-C96,C96.9"
          },
          "normalSampleSite": "",
          "collectionDate": "2025-02-25T00:00:00"
     }],
    
	 */
	public TempusV3Specimen(JSONObject object, TempusV3Json2Vcf tempusJson2Vcf) throws JSONException {
		collectionDate = Json.getStringAttribute(object, "collectionDate");
		sampleCategory = Json.getStringAttribute(object, "sampleCategory");
		notes = Json.getStringAttribute(object, "notes");
		primarySampleSite = Json.getStringAttribute(object, "primarySampleSite");
		normalSampleSite = Json.getStringAttribute(object, "normalSampleSite");
		
		//institution
		JSONObject id = object.getJSONObject("institutionData");
		caseId = Json.getStringAttribute(id, "caseId");
		blockId = Json.getStringAttribute(id, "blockId");
		tumorPercentage = Json.getIntegerAttribute(id, "tumorPercentage");
		if (tumorPercentage != null && tempusJson2Vcf!=null) tempusJson2Vcf.tumorPercentages.count(tumorPercentage);
		
		//diagnosis
		JSONObject d = object.getJSONObject("diagnosis");
		tempusIcdOCodeMorphology = Json.getStringAttribute(d, "tempusIcdOCodeMorphology");
		tempusIcdOCodeTopography = Json.getStringAttribute(d, "tempusIcdOCodeTopography");
		tempusIcd10Code = Json.getStringAttribute(d, "tempusIcd10Code");
		originPathLabDiagnosis = Json.getStringAttribute(d, "originPathLabDiagnosis");
		
		if (tempusJson2Vcf!=null) {
			TempusV3Json2Vcf.add(sampleCategory, tempusJson2Vcf.sampleCategories);
			TempusV3Json2Vcf.add(primarySampleSite, tempusJson2Vcf.sampleSites);
			TempusV3Json2Vcf.add(normalSampleSite, tempusJson2Vcf.sampleSites);
			TempusV3Json2Vcf.add(originPathLabDiagnosis, tempusJson2Vcf.diagnosis);
		}
		
	}
	
	/**Added to the vcf header*/
	public void addMetaData(LinkedHashMap<String, String> meta, int id) {
		if (sampleCategory != null) Misc.addConcatinatedValue (meta, "tempusSampleCategory_"+id, sampleCategory, "; ");
		if (primarySampleSite != null) Misc.addConcatinatedValue (meta, "tempusPrimarySampleSite_"+id, primarySampleSite, "; ");
		if (tumorPercentage != null) Misc.addConcatinatedValue (meta, "tempusTumorPercentage_"+id, tumorPercentage.toString(), "; ");	
		if (caseId != null) Misc.addConcatinatedValue (meta, "tempusCaseId_"+id, caseId, "; ");
		if (blockId != null) Misc.addConcatinatedValue (meta, "tempusBlockId_"+id, blockId, "; ");
		
		if (tempusIcdOCodeMorphology != null) Misc.addConcatinatedValue (meta, "tempusIcdOCodeMorphology_"+id, tempusIcdOCodeMorphology, "; ");
		if (tempusIcdOCodeTopography != null) Misc.addConcatinatedValue (meta, "tempusIcdOCodeTopography_"+id, tempusIcdOCodeTopography, "; ");
		if (tempusIcd10Code != null) Misc.addConcatinatedValue (meta, "tempusIcd10Code_"+id, tempusIcd10Code, "; ");
		if (originPathLabDiagnosis != null) Misc.addConcatinatedValue (meta, "tempusOriginPathLabDiagnosis_"+id, originPathLabDiagnosis, "; ");
	}
	
	public static TempusV3Specimen[] getSpecimens(JSONObject object, TempusV3Json2Vcf tempusJson2Vcf) throws JSONException{
		JSONArray so = object.getJSONArray("specimens");
		TempusV3Specimen[] specimens = new TempusV3Specimen[so.length()];
		for (int i=0; i< specimens.length; i++) specimens[i] = new TempusV3Specimen(so.getJSONObject(i), tempusJson2Vcf);
		return specimens;
	}

	/**Added to the spreadsheet output and used in the PMRSearch*/
	public static void addAttributes(LinkedHashMap<String, String> reportAttributes, TempusV3Specimen[] specimens) {
		ArrayList<String> al = new ArrayList<String>();
		Integer tp = null;
		String caseId = null;
		String blockId = null;
		String normalCollectionDate = null;
		String tumorCollectionDate = null;
		String tempusIcdOCodeMorphology = null;
		String tempusIcdOCodeTopography = null;
		String tempusIcd10Code = null;
		
		String tempusIcdOTxtMorphology = null;
		String tempusIcdOTxtTopography = null;
		String tempusIcd10Txt = null;
		
		String originPathLabDiagnosis = null;
		//for each specimen
		for (int i=0; i< specimens.length; i++ ) {
			String site = null;
			String tn = specimens[i].getSampleCategory();  //tumor, heme, normal
			if (tn.contains("normal")) {
				normalCollectionDate = specimens[i].getCollectionDate();
				site = specimens[i].getNormalSampleSite();
			}
			else {
				tumorCollectionDate = specimens[i].getCollectionDate();
				site = specimens[i].getPrimarySampleSite();
				//tumor percentage
				if (specimens[i].getTumorPercentage() != null) tp = specimens[i].getTumorPercentage();
				//caseId
				if (specimens[i].getCaseId() != null) {
					caseId = specimens[i].getCaseId();
					blockId = specimens[i].getBlockId();
				}
				//ICD10 codes
				tempusIcdOCodeMorphology = specimens[i].getTempusIcdOCodeMorphology();
				tempusIcdOCodeTopography = specimens[i].getTempusIcdOCodeTopography();
				tempusIcd10Code = specimens[i].getTempusIcd10Code();
				originPathLabDiagnosis = specimens[i].getOriginPathLabDiagnosis();
				//ICD txt
				tempusIcdOTxtMorphology = fetchIcdTxt(specimens[i].getTempusIcdOTxtMorphology());
				tempusIcdOTxtTopography = fetchIcdTxt(specimens[i].getTempusIcdOTxtTopography());
				tempusIcd10Txt = fetchIcdTxt(specimens[i].getTempusIcd10Txt());
				
			}
			al.add(tn+" : "+site);
		}
		String sum = Misc.stringArrayListToString(al, "; ");
		reportAttributes.put("specimines", sum);
		
		//collection dates?
		if (tumorCollectionDate!=null) reportAttributes.put("tumorCollectionDate", tumorCollectionDate);
		if (normalCollectionDate!=null) reportAttributes.put("normalCollectionDate", normalCollectionDate);
		
		if (tp != null) reportAttributes.put("tumorPercentage", tp.toString());
		//should only be one of these
		if (caseId != null) {
			reportAttributes.put("caseId", caseId);
			reportAttributes.put("blockId", blockId);
		}
		if (tempusIcdOCodeMorphology != null) reportAttributes.put("icdOCodeMorphology", tempusIcdOCodeMorphology);
		if (tempusIcdOTxtMorphology != null) reportAttributes.put("icdOTxtMorphology", tempusIcdOTxtMorphology);
		
		if (tempusIcdOCodeTopography != null) reportAttributes.put("icdOCodeTopography", tempusIcdOCodeTopography);
		if (tempusIcdOTxtTopography != null) reportAttributes.put("icdOTxtTopography", tempusIcdOTxtTopography);
		
		if (tempusIcd10Code != null) reportAttributes.put("icd10Code", tempusIcd10Code);
		if (tempusIcd10Txt != null) reportAttributes.put("icd10Txt", tempusIcd10Txt);

		if (originPathLabDiagnosis != null) reportAttributes.put("diagnosis", originPathLabDiagnosis);
	}
	
	private static String fetchIcdTxt(LinkedHashSet<String> lhs) {
		if (lhs == null || lhs.size() ==0)return null;
		String merge = Misc.linkedSetToString(lhs, "; ");
		return merge;
		
	}

	public String getKey() {
		StringBuilder sb = new StringBuilder(sampleCategory);
		sb.append("_");
		if (caseId!=null) {
			sb.append(caseId);
			sb.append("_");
		}
		if (blockId!=null) {
			sb.append(blockId);
			sb.append("_");
		}
		//don't use sometimes not in the report
		/*if (tumorPercentage!=null) {
			sb.append(tumorPercentage);
			sb.append("_");
		}*/
		if (primarySampleSite!=null) {
			sb.append(primarySampleSite);
			sb.append("_");
		}
		if (collectionDate!=null) {
			sb.append(collectionDate);
		}
		return sb.toString();
	}

	public String getCollectionDate() {
		return collectionDate;
	}
	public String getSampleCategory() {
		return sampleCategory;
	}
	public String getNotes() {
		return notes;
	}
	public Integer getTumorPercentage() {
		return tumorPercentage;
	}
	public String getCaseId() {
		return caseId;
	}
	public String getBlockId() {
		return blockId;
	}
	public String getPrimarySampleSite() {
		return primarySampleSite;
	}

	public String getTempusIcdOCodeMorphology() {
		return tempusIcdOCodeMorphology;
	}

	public String getTempusIcdOCodeTopography() {
		return tempusIcdOCodeTopography;
	}

	public String getTempusIcd10Code() {
		return tempusIcd10Code;
	}

	public String getOriginPathLabDiagnosis() {
		return originPathLabDiagnosis;
	}

	public String getNormalSampleSite() {
		return normalSampleSite;
	}

	public LinkedHashSet<String> getTempusIcdOTxtMorphology() {
		return tempusIcdOTxtMorphology;
	}

	public LinkedHashSet<String> getTempusIcdOTxtTopography() {
		return tempusIcdOTxtTopography;
	}

	public LinkedHashSet<String> getTempusIcd10Txt() {
		return tempusIcd10Txt;
	}


}
