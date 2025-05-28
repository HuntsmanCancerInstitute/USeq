package edu.utah.seq.vcf.json.tempusv3;

import java.util.ArrayList;
import java.util.HashMap;

public class TempusV3JsonCollection {
	
	//fields
	private String tempusOrderId = null;
	private boolean passMinAge = false;
	private boolean alreadyProcessed = true;
	
	//all jsons from a particular 'tempusOrderId'
	private ArrayList<TempusV3JsonSummary> jsonSummaries = new ArrayList<TempusV3JsonSummary>();
	
	public TempusV3JsonCollection (String tempusOrderId) {
		this.tempusOrderId = tempusOrderId;
	}

	public ArrayList<TempusV3JsonSummary> getJsonSummaries() {
		return jsonSummaries;
	}

	public boolean isPassMinAge() {
		return passMinAge;
	}

	public void setPassMinAge(boolean passMinAge) {
		this.passMinAge = passMinAge;
	}

	public String getTempusOrderId() {
		return tempusOrderId;
	}

	public boolean isAlreadyProcessed() {
		return alreadyProcessed;
	}

	public void setAlreadyProcessed(boolean alreadyProcessed) {
		this.alreadyProcessed = alreadyProcessed;
	}

	
	
	
	
	
}
