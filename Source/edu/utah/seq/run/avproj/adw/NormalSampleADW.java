package edu.utah.seq.run.avproj.adw;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import org.json.JSONObject;

public class NormalSampleADW {

	private String platformName = null;
	private String normalDnaName = null;
	private ArrayList<File> normalDnaFastqCram = new ArrayList<File>();
	private String normalWesCramFileNameToFetch = null;
	private String[] linkageDataLine = null;


	public NormalSampleADW (String normalDnaName, String platformName, String[] linkageDataLine) {
		this.normalDnaName = normalDnaName;
		this.platformName = platformName;
		this.linkageDataLine = linkageDataLine;
	}

	public JSONObject fetchJson(ClinicalMolLinkage linkage) throws IOException {
		JSONObject jo = new JSONObject();
		jo.put("capturePlatform", platformName);
		jo.put("normalDNASeqFile", normalWesCramFileNameToFetch);
		jo.put("normalDNASampleLibraryId", normalDnaName);
		ArrayList<String[]> al = new ArrayList<String[]>();
		al.add(linkageDataLine);
		TumorSampleADW.addLinkageInfo(normalDnaName, jo, linkage, al, true);
		return jo;
	}

	public String getPlatformName() {
		return platformName;
	}

	public String getNormalDnaName() {
		return normalDnaName;
	}

	public ArrayList<File> getNormalDnaFastqCram() {
		return normalDnaFastqCram;
	}

	public String getNormalWesCramFileNameToFetch() {
		return normalWesCramFileNameToFetch;
	}

	public void setNormalWesCramFileNameToFetch(String normalWesCramFileNameToFetch) {
		this.normalWesCramFileNameToFetch = normalWesCramFileNameToFetch;
	}

	public String[] getLinkageDataLine() {
		return linkageDataLine;
	}

}
