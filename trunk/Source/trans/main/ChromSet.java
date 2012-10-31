package trans.main;

import java.io.File;

import util.gen.IO;

public class ChromSet {
	private File saveDirectory;
	private File[] chipSetDirectories;
	private boolean treatmentReplica;
	private long numberOligos = 0;
	
	public void loadNumberOligos(){
		File num = new File (saveDirectory, "totalNumberOligos.txt");
		if (num.exists()){
			String[] lines = IO.loadFileIntoStringArray(num);
			if (lines != null && lines.length !=0) numberOligos = Integer.parseInt(lines[0]);
		}
	}

	public File[] getChipSetDirectories() {
		return chipSetDirectories;
	}

	public void setChipSetDirectories(File[] chipSetDirectories) {
		this.chipSetDirectories = chipSetDirectories;
	}

	public long getNumberOligos() {
		return numberOligos;
	}

	public void setNumberOligos(long numberOligos) {
		this.numberOligos = numberOligos;
	}

	public File getSaveDirectory() {
		return saveDirectory;
	}

	public void setSaveDirectory(File saveDirectory) {
		this.saveDirectory = saveDirectory;
	}

	public boolean isTreatmentReplica() {
		return treatmentReplica;
	}

	public void setTreatmentReplica(boolean treatmentReplica) {
		this.treatmentReplica = treatmentReplica;
	}
}