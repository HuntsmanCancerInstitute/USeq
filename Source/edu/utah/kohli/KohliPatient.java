package edu.utah.kohli;

import java.util.ArrayList;

public class KohliPatient {

	private String hciPatientId = null;
	private ArrayList<KohliSample> samples = new ArrayList<KohliSample>();
	
	public KohliPatient (String[] tokens) {
		// #HCIPersonID 	SampleID	Type	DateDrawn
		this.hciPatientId = tokens[0];
		boolean isGermline = false;
		if (tokens[0].toLowerCase().contains("germline")) isGermline = true;
		samples.add(new KohliSample(tokens[1], isGermline, tokens[3]));
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(hciPatientId); sb.append("\n");
		for (KohliSample ks: samples) {
			sb.append("\t");
			sb.append(ks.toString());
			sb.append("\n");
		}
		return sb.toString();
	}

	public String getHciPatientId() {
		return hciPatientId;
	}

	public ArrayList<KohliSample> getSamples() {
		return samples;
	}
}
