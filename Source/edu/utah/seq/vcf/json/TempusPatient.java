package edu.utah.seq.vcf.json;

import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;

public class TempusPatient {
	
	private String firstName = null;
	private String lastName = null;
	private String tempusId = null;
	private String emr_id = null;
	private String sex = null;
	private String dateOfBirth = null;
	private String diagnosis = null;
	private String diagnosisDate = null;
	
	/*
	  "patient": {
        "firstName": "xxx",
        "lastName": "xxx",
        "tempusId": "87909418-f877-4387-88df-b7732b8feb43",
        "emr_id": "1111111",
        "sex": "FEMALE",
        "DoB": "19xx-0x-0x",
        "diagnosis": "Pancreatic adenocarcinoma",
        "diagnosisDate": "xx/02/201x"
    }
	 */
	public TempusPatient(JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException {
		JSONObject patient = object.getJSONObject("patient");
		firstName = Json.getStringAttribute(patient, "firstName");
		lastName = Json.getStringAttribute(patient, "lastName");
		tempusId = Json.getStringAttribute(patient, "tempusId");
		emr_id = Json.getStringAttribute(patient, "emr_id");
		sex = Json.getStringAttribute(patient, "sex");
		dateOfBirth = Json.getStringAttribute(patient, "DoB");
		diagnosis = Json.getStringAttribute(patient, "diagnosis");
		diagnosisDate = Json.getStringAttribute(patient, "diagnosisDate");
		
		TempusJson2Vcf.add(diagnosis, tempusJson2Vcf.diagnosis);
	}

	public String getFirstName() {
		return firstName;
	}

	public String getLastName() {
		return lastName;
	}

	public String getTempusId() {
		return tempusId;
	}

	public String getEmr_id() {
		return emr_id;
	}

	public String getSex() {
		return sex;
	}

	public String getDateOfBirth() {
		return dateOfBirth;
	}

	public String getDiagnosis() {
		return diagnosis;
	}

	public String getDiagnosisDate() {
		return diagnosisDate;
	}

}
