package edu.utah.seq.pmr;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Iterator;

import org.json.JSONArray;
import org.json.JSONObject;
import util.gen.IO;

public class AvatarClinicalInfo {
	
	private LinkedHashMap<String, String> patient = new LinkedHashMap<String, String>();
	private LinkedHashMap<String, String> root = new LinkedHashMap<String, String>();
	private LinkedHashMap<String, String> tumorDna = new LinkedHashMap<String, String>();
	private LinkedHashMap<String, String> tumorRna = new LinkedHashMap<String, String>();
	private ArrayList<LinkedHashMap<String, String>> normalDna = new ArrayList<LinkedHashMap<String, String>>();
	
	public static final String TUMOR_DNA_NAME = "TumorDNA";
	public static final String TUMOR_RNA_NAME = "TumorRNA";
	public static final String NORMAL_DNA_NAME = "NormalDNA";
	public static final String PATIENT_NAME = "Patient";
	
	
	public AvatarClinicalInfo (File json) {
		String jString = IO.loadFile(json, " ", true);
		JSONObject mainJsonObject = new JSONObject(jString);
		Iterator it = mainJsonObject.keys();
		while (it.hasNext()) {
			String key = (String)it.next();
			if (key.equals(PATIENT_NAME)) loadLinkedHashMap(patient,  mainJsonObject.getJSONObject(key));
			else if (key.equals(TUMOR_DNA_NAME)) loadLinkedHashMap(tumorDna,  mainJsonObject.getJSONObject(key));
			else if (key.equals(TUMOR_RNA_NAME)) loadLinkedHashMap(tumorRna,  mainJsonObject.getJSONObject(key));
			else if (key.equals(NORMAL_DNA_NAME)) {
				JSONArray allNorm = mainJsonObject.getJSONArray(key);
				int num = allNorm.length();
				for (int i=0; i< num; i++) {
					JSONObject no = allNorm.getJSONObject(i);
					LinkedHashMap<String, String> normal = new LinkedHashMap<String, String>();
					loadLinkedHashMap(normal, no);
					normalDna.add(normal);
				}
			}
			else {
				root.put(key, mainJsonObject.get(key).toString());
			}
		}
		
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("  Patient:\n");
		addLinkedHashMap(patient, sb);
		sb.append("  Root:\n");
		addLinkedHashMap(root, sb);
		if (normalDna.size()!=0) {
			for (LinkedHashMap<String, String> n: normalDna) {
				sb.append("  NormalDNA:\n");
				addLinkedHashMap(n, sb);
			}
		}
		if (tumorDna.size()!=0) {
			sb.append("  TumorDNA:\n");
			addLinkedHashMap(tumorDna, sb);
		}
		if (tumorRna.size()!=0) {
			sb.append("  TumorRNA:\n");
			addLinkedHashMap(tumorRna, sb);
		}
		return sb.toString();
	}
	
	private void addLinkedHashMap(LinkedHashMap<String,String> lhm, StringBuilder sb) {
		for (String key : lhm.keySet()) {
			sb.append("    ");
			sb.append(key);
			sb.append(" : ");
			sb.append(lhm.get(key));
			sb.append("\n");
		}
	}
	


	private void loadLinkedHashMap(LinkedHashMap<String, String> hm, JSONObject j) {
		Iterator it = j.keys();
		while (it.hasNext()) {
			String key = (String)it.next();
			hm.put(key, j.getString(key));
		}
	}

	public LinkedHashMap<String, String> getPatient() {
		return patient;
	}

	public LinkedHashMap<String, String> getRoot() {
		return root;
	}

	public LinkedHashMap<String, String> getTumorDna() {
		return tumorDna;
	}

	public LinkedHashMap<String, String> getTumorRna() {
		return tumorRna;
	}

	public ArrayList<LinkedHashMap<String, String>> getNormalDna() {
		return normalDna;
	}

}
