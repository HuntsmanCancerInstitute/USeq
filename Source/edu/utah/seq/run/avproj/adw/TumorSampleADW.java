package edu.utah.seq.run.avproj.adw;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.json.JSONArray;
import org.json.JSONObject;

import util.gen.Misc;


public class TumorSampleADW {
	
	//fields
	private String genericSpecimineId = null;
	private ArrayList<String[]> tumorLinkageDataLines = new ArrayList<String[]>();

	//Tumor Exome
	private String platformName = null;
	private String tumorDnaSampleName = null;
	private ArrayList<File> tumorDnaFastqFiles = new ArrayList<File>();
	private HashSet<String>  tumorWesFastqPathsToFetch = null;
	
	//Tumor RNA
	private String tumorRnaSampleName = null;
	private ArrayList<File> tumorRnaFastqFiles = new ArrayList<File>();
	private HashSet<String>  tumorRnaFastqPathsToFetch = null;
	

	public TumorSampleADW (String tumorDnaName, String tumorRnaName, String platformName, String trimmedSpecimineId, String[] tumorLinkageDataLine) {
		this.tumorDnaSampleName = tumorDnaName;
		this.tumorRnaSampleName = tumorRnaName;
		this.platformName = platformName;
		this.genericSpecimineId = trimmedSpecimineId;
		tumorLinkageDataLines.add(tumorLinkageDataLine);
	}
		
	public JSONObject fetchTumorDnaJson(ClinicalMolLinkage linkage) throws IOException {
		JSONObject jo = new JSONObject();
		jo.put("capturePlatform", platformName);
		jo.put("trimmedSpecimineId", genericSpecimineId);
		JSONArray ja = new JSONArray();
		String[] paths = mergePairedFastqPaths(tumorWesFastqPathsToFetch);
		ja.put(paths[0]);
		ja.put(paths[1]);
		jo.put("tumorDNASeqPaths", ja);
		jo.put("tumorDNASampleLibraryId", tumorDnaSampleName);
		addLinkageInfo(tumorDnaSampleName, jo, linkage, tumorLinkageDataLines, true);
		return jo;
	}
	
	public JSONObject fetchTumorRnaJson(ClinicalMolLinkage linkage) throws IOException {
		JSONObject jo = new JSONObject();
		jo.put("trimmedSpecimineId", genericSpecimineId);
		JSONArray ja = new JSONArray();
		String[] paths = mergePairedFastqPaths(tumorRnaFastqPathsToFetch);
		for (String p: paths) ja.put(p);
		jo.put("tumorRNASeqPaths", ja);
		jo.put("tumorRNASampleLibraryId", tumorRnaSampleName);
		addLinkageInfo(tumorRnaSampleName, jo, linkage, tumorLinkageDataLines, false);
		return jo;
	}
	
	public static String[] mergePairedFastqPaths(HashSet<String> al) throws IOException {
		return Misc.hashSetToStringArray(al);
	}
	
	public static String fetchFastqPathDir(ArrayList<String[]> al) {
		//  /Avatar_MolecularData_hg38/2023_06_30/Whole_Exome/FASTq/FT-SA212052_R1.fastq.gz
		//  /           1                   2          3        4        5
		// look at just first, they will be the same; skip the actual fastq.gz
		String merged = Misc.stringArrayToString(al.get(0), "/");
		return merged.substring(0,merged.lastIndexOf('/'));
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
		return tumorDnaSampleName;
	}

	public ArrayList<File> getTumorDnaFastqFiles() {
		return tumorDnaFastqFiles;
	}

	public String getTumorRnaName() {
		return tumorRnaSampleName;
	}

	public ArrayList<File> getTumorRnaFastqFiles() {
		return tumorRnaFastqFiles;
	}

	public HashSet<String>  getTumorWesFastqPathsToFetch() {
		return tumorWesFastqPathsToFetch;
	}

	public void setTumorWesFastqPathsToFetch(HashSet<String>  tumorWesFastqPathsToFetch) {
		this.tumorWesFastqPathsToFetch = tumorWesFastqPathsToFetch;
	}

	public HashSet<String>  getTumorRnaFastqPathsToFetch() {
		return tumorRnaFastqPathsToFetch;
	}

	public void setTumorRnaPathsToFetch(HashSet<String>  tumorRnaFastqPathsToFetch) {
		this.tumorRnaFastqPathsToFetch = tumorRnaFastqPathsToFetch;
	}

	public String getGenericSpecimineId() {
		return genericSpecimineId;
	}

	public void setTumorRnaName(String tumorRnaName) {
		this.tumorRnaSampleName = tumorRnaName;
	}

	public ArrayList<String[]> getTumorLinkageDataLines() {
		return tumorLinkageDataLines;
	}




	
	
}
