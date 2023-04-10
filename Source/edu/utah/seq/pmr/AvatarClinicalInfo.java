package edu.utah.seq.pmr;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import org.json.JSONArray;
import org.json.JSONObject;
import util.gen.IO;

public class AvatarClinicalInfo {
	
	private HashMap<String, String> patient = new HashMap<String, String>();
	private HashMap<String, String> root = new HashMap<String, String>();
	private HashMap<String, String> tumorDna = new HashMap<String, String>();
	private HashMap<String, String> tumorRna = new HashMap<String, String>();
	private ArrayList<HashMap<String, String>> normalDna = new ArrayList<HashMap<String, String>>();
	
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
			if (key.equals(PATIENT_NAME)) loadHashMap(patient,  mainJsonObject.getJSONObject(key));
			else if (key.equals(TUMOR_DNA_NAME)) loadHashMap(tumorDna,  mainJsonObject.getJSONObject(key));
			else if (key.equals(TUMOR_RNA_NAME)) loadHashMap(tumorRna,  mainJsonObject.getJSONObject(key));
			else if (key.equals(NORMAL_DNA_NAME)) {
				JSONArray allNorm = mainJsonObject.getJSONArray(key);
				int num = allNorm.length();
				for (int i=0; i< num; i++) {
					JSONObject no = allNorm.getJSONObject(i);
					HashMap<String, String> normal = new HashMap<String, String>();
					loadHashMap(normal, no);
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
		sb.append("Patient:\t"+ patient+"\n");
		sb.append("Root:\t"+ root+"\n");
		sb.append("TumorDNA:\t"+ tumorDna+"\n");
		sb.append("TumorRNA:\t"+ tumorRna+"\n");
		for (HashMap<String, String> n: normalDna) {
			sb.append("NormalDNA:\t"+ n+"\n");
		}
		return sb.toString();
	}
	


	private void loadHashMap(HashMap<String, String> hm, JSONObject j) {
		Iterator it = j.keys();
		while (it.hasNext()) {
			String key = (String)it.next();
			hm.put(key, j.getString(key));
		}
	}

	public HashMap<String, String> getPatient() {
		return patient;
	}

	public HashMap<String, String> getRoot() {
		return root;
	}

	public HashMap<String, String> getTumorDna() {
		return tumorDna;
	}

	public HashMap<String, String> getTumorRna() {
		return tumorRna;
	}

	public ArrayList<HashMap<String, String>> getNormalDna() {
		return normalDna;
	}

}
