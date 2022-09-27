package edu.utah.seq.run.avproj;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import util.gen.IO;
import util.gen.Misc;

public class NormalSample {

	private String platformName = null;
	private String normalDnaName = null;
	private ArrayList<File> normalDnaFastqCram = new ArrayList<File>();
	
	public NormalSample (String normalDnaName, String platformName) {
		this.normalDnaName = normalDnaName;
		this.platformName = platformName;
	}

	public boolean fetchFastqCrams (HashMap<String, File> nameFastq, String patientId) {
		for (String name : nameFastq.keySet()) {
			if (name.startsWith(normalDnaName)) normalDnaFastqCram.add(nameFastq.get(name));
		}
		//check number
		int num = normalDnaFastqCram.size();

		if (num == 1 && normalDnaFastqCram.get(0).getName().endsWith(".cram")) return true;
		if (num == 2 && normalDnaFastqCram.get(0).getName().endsWith(".gz")) return true;
		else {
			IO.pl("Platform "+platformName+"  NormalName "+normalDnaName);
			Misc.printErrAndExit("ERROR: normal sample for patient "+patientId+ " Incorrect # number DNA fastq or cram files "+ num+" NormalDnaName "+normalDnaName);
		}
		return false;
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
}
