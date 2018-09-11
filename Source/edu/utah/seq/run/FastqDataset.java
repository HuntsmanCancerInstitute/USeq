package edu.utah.seq.run;

import java.io.File;
import java.io.IOException;

/**Info related to a directory of fastq files.*/
public class FastqDataset{

	private boolean fastqDirExists = false;
	private File[] fastqs = null;
	private String name = null;

	public FastqDataset (File fastqDir, String name) throws IOException{
		if (fastqDir != null) {
			this.name = name;
			File dir = new File(fastqDir, name);
			if (dir.exists()) {
				fastqDirExists = true;
				//ok it exits check that there are two files
				fastqs = TNSample.checkNumberFiles(dir, ".gz", 2);
			}
		}
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
}
