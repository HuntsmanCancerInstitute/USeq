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
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(patientId); 
		sb.append("\n\tGender "+gender);
		sb.append("\n\t#Norm "+normalSamples.size());
		sb.append("\n\t#Tum "+tumorSamples.size());
		return sb.toString();
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

	public void setGender(String gender) {
		this.gender = gender;
	}
}
