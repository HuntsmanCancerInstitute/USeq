package edu.utah.seq.run.avproj;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import util.gen.IO;
import util.gen.Misc;

public class NormalSample {

	private String platformName = null;
	private String normalDnaName = null;
	private ArrayList<File> normalDnaFastq = new ArrayList<File>();


	public NormalSample (String normalDnaName, String platformName) {
		this.normalDnaName = normalDnaName;
		this.platformName = platformName;
	}

	public boolean fetchFastqs (HashMap<String, File> nameFastq, String patientId) {
		for (String name : nameFastq.keySet()) {
			if (name.startsWith(normalDnaName)) normalDnaFastq.add(nameFastq.get(name));
		}
		//check number
		int num = normalDnaFastq.size();
		if (num == 1 || num > 2) Misc.printErrAndExit("ERROR: normal sample for patient "+patientId+ "Incorrect # number DNA fastq samples "+ num);
		
		if (num == 2) return true;
		return false;
	}

	public String getPlatformName() {
		return platformName;
	}

	public String getNormalDnaName() {
		return normalDnaName;
	}

	public ArrayList<File> getNormalDnaFastq() {
		return normalDnaFastq;
	}
}
