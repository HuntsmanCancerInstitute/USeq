package edu.utah.seq.run.ambry;

import java.util.ArrayList;
import java.util.HashMap;

public class AmbryPatient {
	
	//fields
	private String patientKey = null; //first_last_dob_gender_mrn some are missing mrn
	private ArrayList<AmbryResult> ambryResults = new ArrayList<AmbryResult>();
	
	public AmbryPatient (String patientKey) {
		this.patientKey = patientKey;
	}

	public String getPatientKey() {
		return patientKey;
	}

	public ArrayList<AmbryResult> getAmbryResults() {
		return ambryResults;
	}
	
}
