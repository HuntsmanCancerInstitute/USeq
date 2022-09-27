package edu.utah.seq.run.avproj.adw;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

public class ClinicalMolLinkage {
	
	//fields
	private boolean parsed = false;
	private String header = null;
	private String[] headerKeys = null;
	private HashMap<String, Integer> headerKeyIndex = new HashMap<String, Integer>();
	private HashMap<String, ArrayList<String[]>> avatarIdDataLines = new HashMap<String, ArrayList<String[]>>();
	private HashMap<String, String> keyDiseaseType = new HashMap<String, String>();
	
	public static final Pattern commaQuote = Pattern.compile("\",\"");
	public static final Pattern quote = Pattern.compile("\"");
	public static final Pattern dash = Pattern.compile("-");

	public ClinicalMolLinkage( File cvsFromM2Gen) {
		
		parseIt(cvsFromM2Gen);
		
		makeHeaderLookup();
		
		parseDiseaseType();
		
	}

	private void parseDiseaseType() {
		//pull indexes
		int tumGermType = headerKeyIndex.get("Tumor/Germline");
		int slIds = headerKeyIndex.get("WES");
		int disType = headerKeyIndex.get("Disease Type");
		//for each patient
		for (String avatarKey: avatarIdDataLines.keySet()) {
			//for each data line
			for (String[] line: avatarIdDataLines.get(avatarKey)) {
				//is it a Tumor? and the SLID isn't empty
				if (line[tumGermType].contains("Tumor") && line[slIds].length()!=0) {
					String key = avatarKey+"_"+line[slIds];
					// GYN - Endometrial Cancer; GYN - Ovarian Cancer; SAR - Sarcoma; HEM - Myeloma Spectrum; HEM - Leukemia; NEU - Brain Cancer; GU - Prostate Cancer; NEU - Other Neurologic Cancer; GI - Gastric Cancer; H&N - Head and Neck Cancer; GU - Bladder Cancer; GI - Pancreatic Cancer; HEM - Lymphoma; GI - Other GI Cancer; GU - Kidney Cancer; THO - Lung Cancer; GYN - Other GYN Cancer; END - Thyroid Cancer; GI - Colorectal Cancer; BRE - Breast Cancer; CUT - Melanoma; GU - Other GU Cancer; HEM - Other HEME Disease; THO - Other Thoracic Cancer; CUT - Other Cutaneous Cancer; GI - Esophageal Cancer; GI - Liver Cancer; OTH - Other Cancer, NOS
					String value = line[disType];
					//watch out for H&N
					if (value.startsWith("H&N")) value = "HN"+ value.substring(3);
					// split on - and trim
					String[] v = dash.split(value);
					value = v[0].trim();
					keyDiseaseType.put(key, value);
				}
			}
		}
	}

	private void makeHeaderLookup() {
		headerKeys = Misc.TAB.split(header);
		for (int i=0; i<headerKeys.length; i++) headerKeyIndex.put(headerKeys[i], i);
	}

	private void parseIt(File cvsFromM2Gen) {
		String line = null;
		try {
			BufferedReader in = IO.fetchBufferedReader(cvsFromM2Gen);
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() == 0) continue;
				//swap out "," for tabs
				line = commaQuote.matcher(line).replaceAll("\t");
				line = quote.matcher(line).replaceAll("");
				if (line.startsWith("ORIENAvatarKey")) header = line;
				else addToHash(line);
				
			}
			in.close();
			parsed = true;
		} catch (IOException e) {
			IO.el("Failed to parse "+ line + " in "+cvsFromM2Gen);
			e.printStackTrace();
		}
	}

	private void addToHash(String line) {
		String[] f = Misc.TAB.split(line);
		String avatarId = f[0];
		ArrayList<String[]> al = avatarIdDataLines.get(avatarId);
		if (al == null) {
			al = new ArrayList<String[]>();
			avatarIdDataLines.put(avatarId, al);
		}
		al.add(f);
		
	}
	
	public static void main (String[] args) {
		File test = new File ("/Users/u0028003/HCI/AvatarORIEN/AutoAvatar/20220816_HCI_ClinicalMolLinkage_V4.csv");
		ClinicalMolLinkage cml = new ClinicalMolLinkage(test);
		IO.pl(cml.getKeyDiseaseType());
	}

	public boolean isParsed() {
		return parsed;
	}

	public String getHeader() {
		return header;
	}

	public HashMap<String, ArrayList<String[]>> getAvatarIdDataLines() {
		return avatarIdDataLines;
	}

	public HashMap<String, Integer> getHeaderKeyIndex() {
		return headerKeyIndex;
	}

	public HashMap<String, String> getKeyDiseaseType() {
		return keyDiseaseType;
	}

	public String[] getHeaderKeys() {
		return headerKeys;
	}
}
