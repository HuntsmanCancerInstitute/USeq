package edu.utah.seq.run.avproj.adw;

import java.util.ArrayList;
import org.json.JSONObject;
import edu.utah.hci.bioinfo.smm.Subject;

public class PatientADW {
	
	private String patientId = null;
	private ArrayList<TumorSampleADW> tumorSampleADWs = new ArrayList<TumorSampleADW>();
	private ArrayList<NormalSampleADW> normalSampleADWs = new ArrayList<NormalSampleADW>();
	private ArrayList<AvatarAnalysisJob> analysisJobs = new ArrayList<AvatarAnalysisJob>();
	private Subject subjectMatchMaker = null;
	
	public PatientADW (String patientId) {
		this.patientId = patientId;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(patientId); 
		sb.append("\n\t#Norm "+normalSampleADWs.size());
		sb.append("\n\t#Tum "+tumorSampleADWs.size());
		return sb.toString();
	}
	
	public JSONObject fetchJson() {
		JSONObject jo = new JSONObject();
		//Molecular Data Patient ID
		jo.put("molecularDataPatientId", subjectMatchMaker.getCoreIdNewOrMatch());
		//ORIENAvatarKey
		jo.put("orienAvatarKey", patientId);
		//HCIPatientID
		String hciId = subjectMatchMaker.getOtherSubjectIds()[1];
		jo.put("hciPatientId", hciId);
		//sex is what you were born with, gender is what you choose to be.
		jo.put("sex", subjectMatchMaker.getGender());
		return jo;
	}
	
	
	public void addAnalysisJob(AvatarAnalysisJob aj) {
		analysisJobs.add(aj);
	}

	public String getPatientId() {
		return patientId;
	}

	public ArrayList<TumorSampleADW> getTumorSamples() {
		return tumorSampleADWs;
	}

	public ArrayList<NormalSampleADW> getNormalSamples() {
		return normalSampleADWs;
	}

	public ArrayList<AvatarAnalysisJob> getAnalysisJobs() {
		return analysisJobs;
	}

	public void setAnalysisJobs(ArrayList<AvatarAnalysisJob> analysisJobs) {
		this.analysisJobs = analysisJobs;
	}

	public void setSubject(Subject subject) {
		subjectMatchMaker = subject;
		
	}

	public Subject getSubjectMatchMaker() {
		return subjectMatchMaker;
	}
}
