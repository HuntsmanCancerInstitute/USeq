package edu.utah.seq.run.tempus;

import org.json.JSONException;
import org.json.JSONObject;
import util.gen.Json;

public class TempusJsonPatient {

		private String firstName = null;
		private String lastName = null;
		private String tempusId = null;
		private String emr_id = null;
		private String sex = null;
		private String dateOfBirth = null;
		private String diagnosis = null;
		private String diagnosisDate = null;

		/* Loads and then deletes the patient phi
		  "patient": {
	        "firstName": "xxx",
	        "lastName": "xxx",
	        "tempusId": "87909418-f877-4387-88df-b7732b8feb43",
	        "emr_id": "1111111",
	        "sex": "FEMALE",
	        "DoB": "19xx-0x-0x",  year, month, day
	        "diagnosis": "Pancreatic adenocarcinoma",
	        "diagnosisDate": "xx/02/201x" month day year, can also be 1999-08-29 so year month day
	        
	    "patient": {
        	"firstName": "John",
        	"lastName": "Smith",
        	"tempusId": "80b0b12-05d5-44bb-8a01-f998380301",
        	"emrId": "2222222",
        	"sex": "Male",
        	"dateOfBirth": "1972-07-12",
        	"diagnosis": "High-grade pleomorphic epithelioid and spindle cell malignancy",
        	"diagnosisDate": ""
	    }
		 */
		public TempusJsonPatient(JSONObject object) throws JSONException {
			JSONObject patient = object.getJSONObject("patient");
			firstName = Json.getStringAttribute(patient, "firstName");
			patient.remove("firstName");
			lastName = Json.getStringAttribute(patient, "lastName");
			patient.remove("lastName");
			tempusId = Json.getStringAttribute(patient, "tempusId");
			emr_id = Json.getStringAttribute(patient, "emr_id");
			patient.remove("emr_id");
			if (emr_id == null) {
				emr_id = Json.getStringAttribute(patient, "emrId");
				patient.remove("emrId");
			}
			sex = Json.getStringAttribute(patient, "sex");
			if (sex != null) sex = sex.toUpperCase();
			dateOfBirth = Json.getStringAttribute(patient, "DoB");
			patient.remove("DoB");
			if (dateOfBirth == null) {
				dateOfBirth = Json.getStringAttribute(patient, "dateOfBirth");
				patient.remove("dateOfBirth");
			}
			diagnosis = Json.getStringAttribute(patient, "diagnosis");
			diagnosisDate = Json.getStringAttribute(patient, "diagnosisDate");
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
