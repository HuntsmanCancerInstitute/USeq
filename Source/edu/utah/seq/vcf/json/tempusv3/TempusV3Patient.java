package edu.utah.seq.vcf.json.tempusv3;

import java.io.IOException;
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
	//private String tempusId = null;  do not use! not the same between diff tests results
	private String emrId = null;
	private String sex = null;
	private String dateOfBirth = null;
	
	//might not be used anymore
	private String diagnosisDate = null;
	private Integer ageAtDiagnosis = null;
	
	//patient info for Subject Match Maker
	private int dobMonth = -1;
	private int dobDay = -1;
	private int dobYear = -1;
	private String gender = null; //M or F
	private String pmrId = null;

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
	public TempusV3Patient(JSONObject object, TempusV3Json2Vcf tempusJson2Vcf, boolean removePhi) throws Exception {
		JSONObject patient = object.getJSONObject("patient");
		firstName = Json.getStringAttribute(patient, "firstName");
		lastName = Json.getStringAttribute(patient, "lastName");
		//tempusId = Json.getStringAttribute(patient, "tempusId");
		emrId = Json.getStringAttribute(patient, "emrId");
		sex = Json.getStringAttribute(patient, "sex");
		if (sex != null) sex = sex.toUpperCase();
		dateOfBirth = Json.getStringAttribute(patient, "dateOfBirth");
		
		//not sure if these are used anymore? both are null in all jsons
		diagnosisDate = Json.getStringAttribute(patient, "diagnosisDate");
		
		if (removePhi) {
			patient.remove("firstName");
			patient.remove("lastName");
			patient.remove("emrId");
			patient.remove("dateOfBirth");
		}

		//compute ageAtDiagnosis?
		if (dateOfBirth != null && diagnosisDate != null) {
			String[] dobSplit = Misc.DASH.split(dateOfBirth);
			if (dobSplit.length != 3) throw new JSONException("\nFailed to parse three fields from DoB "+dateOfBirth);
			// year, month, day
			LocalDate ldob = new LocalDate(Integer.parseInt(dobSplit[0]), Integer.parseInt(dobSplit[1]), Integer.parseInt(dobSplit[2]));
			
			//not always present
			String[] dodSplit = Misc.DASH.split(diagnosisDate);
			if (dodSplit.length == 3 ) {
				LocalDate ldod = new LocalDate(Integer.parseInt(dodSplit[0]), Integer.parseInt(dodSplit[1]), Integer.parseInt(dodSplit[2]));
				Period p = new Period(ldob, ldod, PeriodType.yearMonthDay());
				ageAtDiagnosis = p.getYears();
				if (tempusJson2Vcf!=null) tempusJson2Vcf.ageAtDiagnosis.count(ageAtDiagnosis);
			}
		}
		
		//for the Subject Match Maker
		parseDoB();
		parseGender();
	}
	
	/**lastName firstName dobMonth(1-12) dobDay(1-31) dobYear(1900-2050) gender(M|F) mrn */
	public String fetchSubjectMatchMakerLine() {
		StringBuilder sb = new StringBuilder();
		if (lastName != null) sb.append(lastName);
		sb.append("\t");
		
		if (firstName != null) sb.append(firstName);
		sb.append("\t");
		
		if (dobMonth != -1) sb.append(new Integer(dobMonth).toString());
		sb.append("\t");
		
		if (dobDay != -1) sb.append(new Integer(dobDay).toString());
		sb.append("\t");
		
		if (dobYear != -1) sb.append(new Integer(dobYear).toString());
		sb.append("\t");
		
		if (gender != null) sb.append(gender);
		sb.append("\t");
		
		if (emrId != null) sb.append(emrId);
		return sb.toString();
	}
	
	private void parseGender() throws Exception {
		if (sex == null) return;
		//Male or Female
		if (sex.toUpperCase().equals("MALE")) gender = "M";
		else if (sex.toUpperCase().equals("FEMALE")) gender = "F";
		else {
			gender = "NA";
			throw new Exception("ERROR: parsing gender, must be male or female, case insensitive.");
		}
	}

	private void parseDoB() throws Exception {
		if (dateOfBirth == null || dateOfBirth.length()==0) return;
		//1965-06-15, 1950-07-26
		String[] t = Misc.DASH.split(dateOfBirth);
		if (t.length!=3) throw new Exception("ERROR: parsing date of birth, not 3 fields after splitting on dash, see <dateOfBirth>");
		if (t[0].length()!=0) {
			dobYear = Integer.parseInt(t[0]);
			if (dobYear< 1900 || dobYear > 2050) throw new IOException("ERROR: dob year field '"+t[0]+"' is malformed, must be 1900-2050, see <dateOfBirth>");
		}
		
		if (t[1].length()!=0) {
			dobMonth = Integer.parseInt(t[1]);
			if (dobMonth< 1 || dobMonth > 12) throw new IOException("ERROR: dob month field '"+t[1]+"' is malformed, must be 1-12, see <dateOfBirth>");
		}
		
		if (t[2].length()!=0) {
			dobDay = Integer.parseInt(t[2]);
			if (dobDay< 1 || dobDay > 31) throw new IOException("ERROR: dob day field '"+t[2]+"' is malformed, must be 1-31, see <dateOfBirth>");
		}
	}
	
	/**Added to VCF header*/
	public void addMetaData(LinkedHashMap<String, String> meta) {
		//Misc.addConcatinatedValue (meta,"tempusPatientId", tempusId, "; ");
		Misc.addConcatinatedValue (meta,"tempusPatientSex", sex, "; ");
		if (ageAtDiagnosis != null) Misc.addConcatinatedValue (meta,"tempusAgeAtDiagnosis", ageAtDiagnosis.toString(), "; ");
	}
	
	/**Added to spreadsheet output*/
	public void addAttributes(LinkedHashMap<String, String> meta, boolean includePhi) {		
		//meta.put("tempusPatientId", tempusId);
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
	public void setPmrId(String coreId) {
		pmrId = coreId;
	}
	public String getPmrId() {
		return pmrId;
	}

	public String getGender() {
		return gender;
	}

}
