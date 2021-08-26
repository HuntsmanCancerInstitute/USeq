package edu.utah.seq.run.avproj;

import java.util.ArrayList;

public class Patient {
	
	private String patientId = null;
	private String gender = null;
	private ArrayList<TumorSample> tumorSamples = new ArrayList<TumorSample>();
	private ArrayList<NormalSample> normalSamples = new ArrayList<NormalSample>();
	
	public Patient (String patientId, String gender) {
		this.patientId = patientId;
		this.gender = gender;
	}

	public String getPatientId() {
		return patientId;
	}

	public String getGender() {
		return gender;
	}

	public ArrayList<TumorSample> getTumorSamples() {
		return tumorSamples;
	}

	public ArrayList<NormalSample> getNormalSamples() {
		return normalSamples;
	}
}
