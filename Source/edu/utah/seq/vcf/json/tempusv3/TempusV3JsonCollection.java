package edu.utah.seq.vcf.json.tempusv3;

import java.util.ArrayList;
import java.util.HashMap;

public class TempusV3JsonCollection {
	
	//fields
	private String tempusPatientId = null;
	
	//all jsons from a particular 'tempusOrderId'
	private ArrayList<TempusV3JsonSummary> jsonSummaries = new ArrayList<TempusV3JsonSummary>();
	
	//above jsons split by tumor specimen, sometimes they submit multiple tumors (primary, met1, met2, cfDNA) under one order?
	private HashMap<String, ArrayList<TempusV3JsonSummary>> tumorSummaries = new HashMap<String, ArrayList<TempusV3JsonSummary>>();
	
	public TempusV3JsonCollection (String tempusPatientId) {
		this.tempusPatientId = tempusPatientId;
	}

	public ArrayList<TempusV3JsonSummary> getJsonSummaries() {
		return jsonSummaries;
	}

	public void groupByTumor() {
		
		//for each json report file
		for (TempusV3JsonSummary s: jsonSummaries) {
			
			//for each specimen in the report
			for (TempusV3Specimen specimen: s.getTempusSpecimens()) {
				
				//if not normal? tumor, cfDNA specimen, heme, normal
				if (specimen.getSampleCategory().contains("normal")==false) {
					String key = specimen.getKey();
					ArrayList<TempusV3JsonSummary> al = tumorSummaries.get(key);
					if (al == null) {
						al = new ArrayList<TempusV3JsonSummary>();
						tumorSummaries.put(key, al);
					}
					al.add(s);
				}
			}
		}
		
	}

	public String getTempusPatientId() {
		return tempusPatientId;
	}

	public HashMap<String, ArrayList<TempusV3JsonSummary>> getTumorSummaries() {
		return tumorSummaries;
	}
	
	
	
	
	
}
