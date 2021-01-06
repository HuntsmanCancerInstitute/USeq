package edu.utah.seq.vcf.json;

import java.util.LinkedHashMap;
import org.joda.time.LocalDate;
import org.joda.time.Period;
import org.joda.time.PeriodType;
import org.json.JSONException;
import org.json.JSONObject;

import util.gen.IO;
import util.gen.Json;
import util.gen.Misc;

public class TempusPatient {

	private String firstName = null;
	private String lastName = null;
	private String tempusId = null;
	private String emr_id = null;
	private String sex = null;
	private String dateOfBirth = null;
	private String diagnosis = null;
	private String diagnosisDate = null;
	private Integer ageAtDiagnosis = null;

	/*
	  "patient": {
        "firstName": "xxx",
        "lastName": "xxx",
        "tempusId": "87909418-f877-4387-88df-b7732b8feb43",
        "emr_id": "1111111",
        "sex": "FEMALE",
        "DoB": "19xx-0x-0x",  year, month, day
        "diagnosis": "Pancreatic adenocarcinoma",
        "diagnosisDate": "xx/02/201x" month day year, can also be 1999-08-29 so year month day
    }
	 */
	public TempusPatient(JSONObject object, TempusJson2Vcf tempusJson2Vcf) throws JSONException {
		JSONObject patient = object.getJSONObject("patient");
		firstName = Json.getStringAttribute(patient, "firstName");
		lastName = Json.getStringAttribute(patient, "lastName");
		tempusId = Json.getStringAttribute(patient, "tempusId");
		emr_id = Json.getStringAttribute(patient, "emr_id");
		if (emr_id == null) emr_id = Json.getStringAttribute(patient, "emrId");
		sex = Json.getStringAttribute(patient, "sex");
		if (sex != null) sex = sex.toUpperCase();
		dateOfBirth = Json.getStringAttribute(patient, "DoB");
		if (dateOfBirth == null) dateOfBirth = Json.getStringAttribute(patient, "dateOfBirth");
		diagnosis = Json.getStringAttribute(patient, "diagnosis");
		diagnosisDate = Json.getStringAttribute(patient, "diagnosisDate");

		//compute ageAtDiagnosis?
		if (dateOfBirth != null && diagnosisDate != null) {
			String[] dobSplit = Misc.DASH.split(dateOfBirth);
			if (dobSplit.length != 3) throw new JSONException("\nFailed to parse three fields from DoB "+dateOfBirth);
			LocalDate ldob = new LocalDate(Integer.parseInt(dobSplit[0]), Integer.parseInt(dobSplit[1]), Integer.parseInt(dobSplit[2]));
			String[] dodSplit = Misc.FORWARD_SLASH.split(diagnosisDate);
			LocalDate ldod = null;
			if (dodSplit.length != 3) {
				//might be year-month-day
				if (dodSplit.length != 3) dodSplit = Misc.DASH.split(diagnosisDate);
				
				if (dodSplit.length !=3) {
					tempusJson2Vcf.getWorkingReport().getWarningMessages().add("Failed to parse three fields from diagnosisDate "+diagnosisDate);
				}
				else ldod = new LocalDate(Integer.parseInt(dodSplit[0]), Integer.parseInt(dodSplit[1]), Integer.parseInt(dodSplit[2]));
			}
			else {
				ldod = new LocalDate(Integer.parseInt(dodSplit[2]), Integer.parseInt(dodSplit[0]), Integer.parseInt(dodSplit[1]));
			}
			if (ldod != null) {
				Period p = new Period(ldob, ldod, PeriodType.yearMonthDay());
				ageAtDiagnosis = p.getYears();
				tempusJson2Vcf.ageAtDiagnosis.count(ageAtDiagnosis);
				IO.pl("\nAAD "+ageAtDiagnosis);
			}	
		}
		TempusJson2Vcf.add(diagnosis, tempusJson2Vcf.diagnosis);
	}
	
	public void addMetaData(LinkedHashMap<String, String> meta) {
		meta.put("tempusPatientId", tempusId);
		meta.put("tempusPatientGender", sex);
		if (diagnosis != null) meta.put("tempusDiagnosis", diagnosis);
		if (ageAtDiagnosis != null) meta.put("tempusAgeAtDiagnosis", ageAtDiagnosis.toString());
	}
	
	public void addAttributes(LinkedHashMap<String, String> meta, boolean includePhi) {
		meta.put("tempusId", tempusId);
		meta.put("sex", sex);
		if (diagnosis != null) meta.put("diagnosis", diagnosis);
		if (ageAtDiagnosis != null) meta.put("ageAtDiagnosis", ageAtDiagnosis.toString());
		if (includePhi) {
			meta.put("firstName", firstName);
			meta.put("lastName", lastName);
			meta.put("DoB", dateOfBirth);
			meta.put("emr_id", emr_id);
		}
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
