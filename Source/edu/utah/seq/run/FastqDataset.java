package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import util.gen.IO;

/**Info related to a directory of fastq files.*/
public class FastqDataset{

	private boolean fastqDirExists = false;
	private File[] fastqs = null;
	private String name = null;
	private boolean incorrectNumberFastq = false;

	public FastqDataset (File fastqDir, String name, ArrayList<String> info) throws IOException{
		if (fastqDir != null) {
			this.name = name;
			File dir = new File(fastqDir, name);
			if (dir.exists()) {
				fastqDirExists = true;
				//ok it exits check that there are two files
				fastqs = checkNumberFiles(dir, ".gz", 2);
				if (fastqs == null) {
					info.add("\tDid not find two fastq files within "+dir);
					incorrectNumberFastq = true;
				}
				else info.add("\tREADY "+dir);
			}
		}
		else info.add("\tNO "+name); 
	}
	
	public File[] checkNumberFiles(File dir, String extension, int requiredNumberFiles) throws IOException {
		File[] f = IO.extractFiles(dir, extension);
		if (f.length != requiredNumberFiles) {
			return null;
		}
		return f;
	}

	public boolean isFastqDirExists() {
		return fastqDirExists;
	}
	public File[] getFastqs() {
		return fastqs;
	}
	public String getName() {
		return name;
	}

	public boolean isIncorrectNumberFastq() {
		return incorrectNumberFastq;
	}
	
	public boolean isGoodToAlign() {
		if (fastqDirExists == true && incorrectNumberFastq == false) return true;
		return false;
	}
}
