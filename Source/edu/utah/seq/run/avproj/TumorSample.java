package edu.utah.seq.run.avproj;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import util.gen.IO;
import util.gen.Misc;

public class TumorSample {

	private String platformName = null;

	private String tumorDnaName = null;
	private ArrayList<File> tumorDnaFastq = new ArrayList<File>();

	private String tumorRnaName = null;
	private ArrayList<File> tumorRnaFastq = new ArrayList<File>();

	public TumorSample (String tumorDnaName, String tumorRnaName, String platformName) {
		this.tumorDnaName = tumorDnaName;
		this.tumorRnaName = tumorRnaName;
		this.platformName = platformName;
	}
	
	public boolean fetchFastqs (HashMap<String, File> nameFastq, String patientId) {
		if (tumorDnaName != null && tumorRnaName != null) {
			for (String name : nameFastq.keySet()) {
				 if (name.startsWith(tumorDnaName)) tumorDnaFastq.add(nameFastq.get(name));
				 else if (name.startsWith(tumorRnaName)) tumorRnaFastq.add(nameFastq.get(name));
			}
		}
		else if (tumorDnaName != null) {
			for (String name : nameFastq.keySet()) {
				 if (name.startsWith(tumorDnaName)) tumorDnaFastq.add(nameFastq.get(name));
			}
		}
		else {
			for (String name : nameFastq.keySet()) {
				 if (name.startsWith(tumorRnaName)) tumorRnaFastq.add(nameFastq.get(name));
			}
		}
		
		//check numbers
		String failedDna = null;
		String failedRna = null;
		int numD = tumorDnaFastq.size();
		if (tumorDnaName != null) {
			if (numD == 1 || numD > 2) failedDna = "Incorrect # tumor DNA fastq samples "+ numD;
		}
		
		int numR = tumorRnaFastq.size();
		if (tumorRnaName != null) {
			if (numR == 1 || numR > 2) failedRna = "Incorrect # tumor RNA fastq samples "+ numR;
		}
		
		if (failedDna != null || failedRna != null) {
			Misc.printErrAndExit("ERROR: tumor sample for patient "+patientId+"\n\tDNA: "+failedDna+"\n\tRNA: "+failedRna);
		}
		
		if (numD == 2 || numR == 2) return true;
		
		return false;
	}

	public String getPlatformName() {
		return platformName;
	}

	public String getTumorDnaName() {
		return tumorDnaName;
	}

	public ArrayList<File> getTumorDnaFastq() {
		return tumorDnaFastq;
	}

	public String getTumorRnaName() {
		return tumorRnaName;
	}

	public ArrayList<File> getTumorRnaFastq() {
		return tumorRnaFastq;
	}
	
	
}
