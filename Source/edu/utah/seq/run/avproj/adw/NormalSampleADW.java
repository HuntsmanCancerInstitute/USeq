package edu.utah.seq.run.avproj.adw;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import org.json.JSONObject;

public class NormalSampleADW {

	private String platformName = null;
	private String normalDnaSampleName = null;
	private ArrayList<File> normalDnaFastqFiles = new ArrayList<File>();
	private ArrayList<String[]>  normalWesFastqPathsToFetch = null;
	private String[] linkageDataLine = null;

	public NormalSampleADW (String normalDnaName, String platformName, String[] linkageDataLine) {
		this.normalDnaSampleName = normalDnaName;
		this.platformName = platformName;
		this.linkageDataLine = linkageDataLine;
	}

	public JSONObject fetchJson(ClinicalMolLinkage linkage) throws IOException {
		JSONObject jo = new JSONObject();
		jo.put("capturePlatform", platformName);
		jo.put("normalDNASeqPaths", TumorSampleADW.mergePairedFastqPaths(normalWesFastqPathsToFetch));
		jo.put("normalDNASampleLibraryId", normalDnaSampleName);
		ArrayList<String[]> al = new ArrayList<String[]>();
		al.add(linkageDataLine);
		TumorSampleADW.addLinkageInfo(normalDnaSampleName, jo, linkage, al, true);
		return jo;
	}

	public String getPlatformName() {
		return platformName;
	}

	public String getNormalDnaName() {
		return normalDnaSampleName;
	}

	public ArrayList<File> getNormalDnaFastqFiles() {
		return normalDnaFastqFiles;
	}

	public ArrayList<String[]>  getNormalWesFastqPathsToFetch() {
		return normalWesFastqPathsToFetch;
	}

	public void setNormalWesFastqPathsToFetch(ArrayList<String[]>  normalWesFastqPathsToFetch) {
		this.normalWesFastqPathsToFetch = normalWesFastqPathsToFetch;
	}

	public String[] getLinkageDataLine() {
		return linkageDataLine;
	}

}
