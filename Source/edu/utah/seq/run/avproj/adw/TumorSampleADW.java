package edu.utah.seq.run.avproj.adw;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import org.json.JSONObject;


public class TumorSampleADW {
	
	//fields
	private String genericSpecimineId = null;
	private ArrayList<String[]> tumorLinkageDataLines = new ArrayList<String[]>();

	//Tumor Exome
	private String platformName = null;
	private String tumorDnaName = null;
	private ArrayList<File> tumorDnaFastqCram = new ArrayList<File>();
	private String tumorWesCramFileNameToFetch = null;
	
	//Tumor RNA
	private String tumorRnaName = null;
	private ArrayList<File> tumorRnaFastqCram = new ArrayList<File>();
	private String tumorRnaCramFileNameToFetch = null;
	

	public TumorSampleADW (String tumorDnaName, String tumorRnaName, String platformName, String trimmedSpecimineId, String[] tumorLinkageDataLine) {
		this.tumorDnaName = tumorDnaName;
		this.tumorRnaName = tumorRnaName;
		this.platformName = platformName;
		this.genericSpecimineId = trimmedSpecimineId;
		tumorLinkageDataLines.add(tumorLinkageDataLine);
	}
		
	public JSONObject fetchTumorDnaJson(ClinicalMolLinkage linkage) throws IOException {
		JSONObject jo = new JSONObject();
		jo.put("capturePlatform", platformName);
		jo.put("trimmedSpecimineId", genericSpecimineId);
		jo.put("tumorDNASeqFile", tumorWesCramFileNameToFetch);
		jo.put("tumorDNASampleLibraryId", tumorDnaName);
		addLinkageInfo(tumorDnaName, jo, linkage, tumorLinkageDataLines, true);
		return jo;
	}
	
	public JSONObject fetchTumorRnaJson(ClinicalMolLinkage linkage) throws IOException {
		JSONObject jo = new JSONObject();
		jo.put("trimmedSpecimineId", genericSpecimineId);
		jo.put("tumorRNASeqFile", tumorRnaCramFileNameToFetch);
		jo.put("tumorRNASampleLibraryId", tumorRnaName);
		addLinkageInfo(tumorRnaName, jo, linkage, tumorLinkageDataLines, false);
		return jo;
	}
	
	public static void addLinkageInfo (String sampleName, JSONObject jo, ClinicalMolLinkage linkage, ArrayList<String[]> linkageDataLines, boolean isWes) throws IOException {
		HashMap<String, Integer> headerKeyIndex = linkage.getHeaderKeyIndex();
		String[] headerKeys = linkage.getHeaderKeys();
		
		//for each linkage line see if it belongs to this sample requested
		boolean found = false;
		for (String[] splitLine: linkageDataLines) {
			if (splitLine[headerKeyIndex.get("WES")].equals(sampleName) || splitLine[headerKeyIndex.get("RNASeq")].equals(sampleName)) {
				found = true;
				//for each cell that isn't empty
				for (int i=0; i< splitLine.length; i++) {
					String cell = splitLine[i].trim();
					if (cell.length()!=0) {
						if (isWes && headerKeys[i].startsWith("RNA") == false) jo.put(headerKeys[i], cell);
						else if (isWes==false && headerKeys[i].startsWith("WES") == false) jo.put(headerKeys[i], cell);
					}
				}
			}
		}
		if (found == false) throw new IOException("\nFailed to find a clinical linkage data line for "+sampleName);
	}

	public String getPlatformName() {
		return platformName;
	}

	public String getTumorDnaName() {
		return tumorDnaName;
	}

	public ArrayList<File> getTumorDnaFastqCram() {
		return tumorDnaFastqCram;
	}

	public String getTumorRnaName() {
		return tumorRnaName;
	}

	public ArrayList<File> getTumorRnaFastqCram() {
		return tumorRnaFastqCram;
	}

	public String getTumorWesCramFileNameToFetch() {
		return tumorWesCramFileNameToFetch;
	}

	public void setTumorWesCramFileNameToFetch(String tumorWesCramFileNameToFetch) {
		this.tumorWesCramFileNameToFetch = tumorWesCramFileNameToFetch;
	}

	public String getTumorRnaCramFileNameToFetch() {
		return tumorRnaCramFileNameToFetch;
	}

	public void setTumorRnaCramFileNameToFetch(String tumorRnaCramFileNameToFetch) {
		this.tumorRnaCramFileNameToFetch = tumorRnaCramFileNameToFetch;
	}

	public String getGenericSpecimineId() {
		return genericSpecimineId;
	}

	public void setTumorRnaName(String tumorRnaName) {
		this.tumorRnaName = tumorRnaName;
	}

	public ArrayList<String[]> getTumorLinkageDataLines() {
		return tumorLinkageDataLines;
	}




	
	
}
