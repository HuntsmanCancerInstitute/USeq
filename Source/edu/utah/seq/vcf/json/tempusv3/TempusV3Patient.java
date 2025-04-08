package edu.utah.seq.vcf.json.tempusv3;

import java.util.LinkedHashMap;
import org.joda.time.LocalDate;
import org.joda.time.Period;
import org.joda.time.PeriodType;
import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;
import util.gen.Misc;

public class TempusV3Patient {

	private String firstName = null;
	private String lastName = null;
	private String tempusId = null;
	private String emrId = null;
	private String sex = null;
	private String dateOfBirth = null;
	
	//might not be used anymore
	private String diagnosisDate = null;
	private Integer ageAtDiagnosis = null;

	/*
  "patient": {
    "firstName": "Larry",
    "lastName": "Mock",
    "tempusId": "4cf15a4e-9988-310a-a046-9de30140e79c",
    "emrId": "06591694",
    "sex": "Male",
    "dateOfBirth": "1953-09-30",
    "diagnosisDate": null
  },

	 */
	public TempusV3Patient(JSONObject object, TempusV3Json2Vcf tempusJson2Vcf) throws JSONException {
		JSONObject patient = object.getJSONObject("patient");
		firstName = Json.getStringAttribute(patient, "firstName");
		lastName = Json.getStringAttribute(patient, "lastName");
		tempusId = Json.getStringAttribute(patient, "tempusId");
		emrId = Json.getStringAttribute(patient, "emrId");
		sex = Json.getStringAttribute(patient, "sex");
		if (sex != null) sex = sex.toUpperCase();
		dateOfBirth = Json.getStringAttribute(patient, "dateOfBirth");
		
		//not sure if these are used anymore? both are null in all jsons
		diagnosisDate = Json.getStringAttribute(patient, "diagnosisDate");

		//compute ageAtDiagnosis?
		if (dateOfBirth != null && diagnosisDate != null) {
			String[] dobSplit = Misc.DASH.split(dateOfBirth);
			if (dobSplit.length != 3) throw new JSONException("\nFailed to parse three fields from DoB "+dateOfBirth);
			// year, month, day
			LocalDate ldob = new LocalDate(Integer.parseInt(dobSplit[0]), Integer.parseInt(dobSplit[1]), Integer.parseInt(dobSplit[2]));
			
			String[] dodSplit = Misc.DASH.split(diagnosisDate);
			if (dodSplit.length != 3) {
				tempusJson2Vcf.getWorkingReport().getWarningMessages().add("Failed to parse three fields from diagnosisDate "+diagnosisDate);
			}
			LocalDate ldod = new LocalDate(Integer.parseInt(dodSplit[0]), Integer.parseInt(dodSplit[1]), Integer.parseInt(dodSplit[2]));
			
			Period p = new Period(ldob, ldod, PeriodType.yearMonthDay());
			ageAtDiagnosis = p.getYears();
			tempusJson2Vcf.ageAtDiagnosis.count(ageAtDiagnosis);
		}
	}
	
	/**Added to VCF header*/
	public void addMetaData(LinkedHashMap<String, String> meta) {
		meta.put("tempusPatientId", tempusId);
		meta.put("tempusPatientGender", sex);
		if (ageAtDiagnosis != null) meta.put("tempusAgeAtDiagnosis", ageAtDiagnosis.toString());
	}
	
	/**Added to spreadsheet output*/
	public void addAttributes(LinkedHashMap<String, String> meta, boolean includePhi) {		
		meta.put("tempusId", tempusId);
		meta.put("sex", sex);
		if (ageAtDiagnosis != null) meta.put("ageAtDiagnosis", ageAtDiagnosis.toString());
		if (includePhi) {
			meta.put("firstName", firstName);
			meta.put("lastName", lastName);
			meta.put("DoB", dateOfBirth);
			meta.put("emrId", emrId);
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

	public String getEmrId() {
		return emrId;
	}

	public String getSex() {
		return sex;
	}

	public String getDateOfBirth() {
		return dateOfBirth;
	}

	public String getDiagnosisDate() {
		return diagnosisDate;
	}

}
