package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import util.gen.IO;

/**Info related to a directory of fastq files.*/
public class FastqCramDataset{

	private boolean cramFastqDirExists = false;
	private File[] cramFastqs = null;
	private String name = null;
	private boolean incorrectNumberCramFastq = false;

	public FastqCramDataset (File fastqDir, String name, ArrayList<String> info, TNSample2 tnSample) throws IOException{
		if (fastqDir != null) {
			this.name = name;
			File dir = new File(fastqDir, name);
			if (dir.exists()) {
				cramFastqDirExists = true;
				//look for 2 or more paired fastq
				cramFastqs = checkNumberFiles(dir, ".gz", 2);
				//if not found look for cram
				if (cramFastqs == null) cramFastqs = checkNumberFiles(dir, ".cram", 1);
				if (cramFastqs == null) {
					info.add("\tFAILED: Did not find two or more fastq files or one cram file within "+dir);
					incorrectNumberCramFastq = true;
					tnSample.setFailed(true);
				}
				else info.add("\tREADY "+dir);
			}
		}
		else info.add("\tNO "+name); 
	}
	
	public File[] checkNumberFiles(File dir, String extension, int requiredNumberFiles) throws IOException {
		File[] f = IO.extractFiles(dir, extension);
		if (f.length < requiredNumberFiles) {
			return null;
		}
		return f;
	}

	public boolean isCramFastqDirExists() {
		return cramFastqDirExists;
	}
	public File[] getCramFastqs() {
		return cramFastqs;
	}
	public String getName() {
		return name;
	}

	public boolean isIncorrectNumberCramFastq() {
		return incorrectNumberCramFastq;
	}
	
	public boolean isGoodToAlign() {
		if (cramFastqDirExists == true && incorrectNumberCramFastq == false) return true;
		return false;
	}
}
