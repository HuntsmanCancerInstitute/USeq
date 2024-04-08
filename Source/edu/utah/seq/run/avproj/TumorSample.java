package edu.utah.seq.run.avproj;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import util.gen.Misc;

public class TumorSample {

	private String platformName = null;

	private String tumorDnaName = null;
	private ArrayList<File> tumorDnaFastqCram = new ArrayList<File>();
	private String tumorRnaName = null;
	private ArrayList<File> tumorRnaFastqCram = new ArrayList<File>();

	public TumorSample (String tumorDnaName, String tumorRnaName, String platformName) {
		this.tumorDnaName = tumorDnaName;
		this.tumorRnaName = tumorRnaName;
		this.platformName = platformName;
	}
	
	public boolean fetchFastqCrams (HashMap<String, File> nameFastq, String patientId) {
		if (tumorDnaName != null && tumorRnaName != null) {
			for (String name : nameFastq.keySet()) {
				 if (name.startsWith(tumorDnaName)) tumorDnaFastqCram.add(nameFastq.get(name));
				 else if (name.startsWith(tumorRnaName)) tumorRnaFastqCram.add(nameFastq.get(name));
			}
		}
		else if (tumorDnaName != null) {
			for (String name : nameFastq.keySet()) {
				 if (name.startsWith(tumorDnaName)) tumorDnaFastqCram.add(nameFastq.get(name));
			}
		}
		else {
			for (String name : nameFastq.keySet()) {
				 if (name.startsWith(tumorRnaName)) {
					if (nameFastq.get(name)!=null) tumorRnaFastqCram.add(nameFastq.get(name));
				 }
			}
		}
		
		//check numbers
		String failedDna = null;
		String failedRna = null;
		int numD = tumorDnaFastqCram.size();
		if (tumorDnaName != null) {
			if (numD == 1 && tumorDnaFastqCram.get(0).getName().endsWith(".cram") == false) failedDna = "Incorrect # tumor DNA fastq cram samples "+ numD;
			if (numD == 2 && tumorDnaFastqCram.get(0).getName().endsWith(".gz") == false) failedDna = "Incorrect # tumor DNA fastq samples "+ numD;
			if (numD > 2) failedDna = "Incorrect # tumor DNA fastq cram samples "+ numD;
		}
		
		int numR = tumorRnaFastqCram.size();
		if (tumorRnaName != null) {
			if (numR == 0 || numR > 2) failedRna = "Incorrect # tumor RNA fastq samples "+ numR;
			if (numR == 1 && tumorRnaFastqCram.get(0).getName().endsWith(".cram") == false) failedRna = "Incorrect # tumor RNA fastq cram samples "+ numR;
			if (numR == 2 && tumorRnaFastqCram.get(0).getName().endsWith(".gz") == false) failedRna = "Incorrect # tumor RNA fastq samples "+ numR;
		}
		
		if (failedDna != null || failedRna != null) {
			Misc.printErrAndExit("ERROR: tumor sample for patient "+patientId+"\n\tDNA: "+failedDna+"\n\tRNA: "+failedRna);
		}

		return true;
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
	
	
}
