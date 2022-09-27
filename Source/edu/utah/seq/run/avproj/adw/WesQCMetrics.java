package edu.utah.seq.run.avproj.adw;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

public class WesQCMetrics {
	
	//fields
	private boolean parsed = false;
	private String header = null;
	private HashMap<String, Integer> headerKeyIndex = new HashMap<String, Integer>();
	private HashMap<String, ArrayList<String[]>> avatarIdDataLines = new HashMap<String, ArrayList<String[]>>();
	private HashMap<String, String> slidPlatform = new HashMap<String, String>();
	private HashSet<String> baitSets = new HashSet<String>();
	
	public static final Pattern commaQuote = Pattern.compile("\",\"");
	public static final Pattern quote = Pattern.compile("\"");

	public WesQCMetrics( File wesQcMetrics) throws IOException {
		
		parseIt(wesQcMetrics);
		
		makeHeaderLookup();
		
		loadSlidBaitSets();
		
	}

	private void loadSlidBaitSets() throws IOException {
		Integer indexOfSLID = headerKeyIndex.get("SLID");
		Integer indexOfBaitSet = headerKeyIndex.get("BaitSet");
		if (indexOfSLID == null || indexOfBaitSet == null) throw new IOException("Failed to parse the index of either the SLID or BaitSet from "+header);
		//for each data line
		for (String avatarId: avatarIdDataLines.keySet()) {
			ArrayList<String[]> data = avatarIdDataLines.get(avatarId);
			for (String[] d: data) {
				String slid = d[indexOfSLID];
				//does it already exist?
				if (slidPlatform.containsKey(slid)) throw new IOException ("Duplicate SLID "+slid);
				String baitSet = d[indexOfBaitSet];
				baitSets.add(baitSet);
				//hg38_lifted_SeqCap_EZ_Exome_v3_capture_baits
				//hg38_lifted_IDT_merged_custom_single_double_probes
				//M2GEN_CustomePlus_Exv2_probes_hg38
				if (baitSet.contains("SeqCap_EZ_Exome_v3")) slidPlatform.put(slid, "NIM");
				else if (baitSet.contains("_IDT_merged")) slidPlatform.put(slid, "IDTv1");
				else if (baitSet.contains("M2GEN_CustomePlus_Exv2")) slidPlatform.put(slid, "IDTv2");
				else throw new IOException("Failed to match a platform to "+baitSet);
			}
		}
	}

	private void makeHeaderLookup() {
		String[] f = Misc.TAB.split(header);
		for (int i=0; i<f.length; i++) headerKeyIndex.put(f[i], i);
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
	
	public static void main (String[] args) throws IOException {
		File test = new File ("/Users/u0028003/HCI/AvatarORIEN/AutoAvatar/20220824_HCI_WES_QC_Metrics.csv");
		WesQCMetrics cml = new WesQCMetrics(test);
		IO.pl(cml.getBaitSets());
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

	public HashSet<String> getBaitSets() {
		return baitSets;
	}

	public HashMap<String, String> getSlidPlatform() {
		return slidPlatform;
	}
}
