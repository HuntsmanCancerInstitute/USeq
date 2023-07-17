package edu.utah.seq.pmr;

import java.io.IOException;

import util.gen.Misc;

public class AvatarDataset {
	
	//fields, some may be NA
	private String avatarId = null;
	private String normalExomeId = null;
	private String tumorExomeId = null;
	private String tumorTranscriptomeId = null;
	private boolean mixedPlatform = false;
	private Dataset dataset = null;
	
	public AvatarDataset (Dataset d) throws IOException {
		dataset = d;
		
		// A032049_SL419345_SL419548_SL420681  some may be NA
		String[] t = Misc.UNDERSCORE.split(dataset.getDatasetId());
		avatarId = t[0];
		normalExomeId = t[1];
		tumorExomeId = t[2];
		tumorTranscriptomeId = t[3];
		
		//mixed platform?
		AvatarClinicalInfo i = d.getAvatarClinicalInfo();
		String mixed = i.getRoot().get("mixedCapturePlatforms");
		if (mixed == null) mixed = "NA";
		if (mixed.equals("false")) mixedPlatform = false;
		else if (mixed.equals("true")) mixedPlatform = true;
		else throw new IOException("ERROR: failed to parse true or false from 'mixedCapturePlatforms' in the Avatar json file -> "+d.getAvatarClinicalInfo().getJsonFile());
	}

	public String getAvatarId() {
		return avatarId;
	}

	public String getNormalExomeId() {
		return normalExomeId;
	}

	public String getTumorExomeId() {
		return tumorExomeId;
	}

	public String getTumorTranscriptomeId() {
		return tumorTranscriptomeId;
	}

	public Boolean getMixedPlatform() {
		return mixedPlatform;
	}

	public Dataset getDataset() {
		return dataset;
	}

	public boolean isComplete() {
		//check all three aren't NA
		if (normalExomeId.equals("NA") || tumorExomeId.equals("NA") || tumorTranscriptomeId.equals("NA") || mixedPlatform == true) return false;
		return true;
	}

}
