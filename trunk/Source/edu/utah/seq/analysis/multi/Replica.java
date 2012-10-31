package edu.utah.seq.analysis.multi;

import java.io.*;
import java.util.HashMap;
import net.sf.samtools.SAMFileReader;

public class Replica implements Serializable{
	//fields
	private String nameNumber;
	private File bamFile;
	private HashMap<String, GeneCount> geneCounts = new HashMap<String, GeneCount>();
	private static final long serialVersionUID = 1L;
	private int totalCounts = 0;
	
	//constructor
	public Replica (String nameNumber, File bamFile){
		this.bamFile = bamFile;
		this.nameNumber = nameNumber;
	}

	public String getNameNumber() {
		return nameNumber;
	}

	public void setNameNumber(String nameNumber) {
		this.nameNumber = nameNumber;
	}

	public File getBamFile() {
		return bamFile;
	}

	public void setBamFile(File bamFile) {
		this.bamFile = bamFile;
	}

	public void setGeneCounts(HashMap<String, GeneCount> geneCounts) {
		this.geneCounts = geneCounts;
	}

	public HashMap<String, GeneCount> getGeneCounts() {
		return geneCounts;
	}

	public void setTotalCounts(int totalCounts) {
		this.totalCounts = totalCounts;
	}

	public int getTotalCounts() {
		return totalCounts;
	}

	public void removeFlaggedGenes(String[] badGeneNames) {
		//look for any bad genes
		for (String badGene: badGeneNames){
			if (geneCounts.containsKey(badGene)){
				GeneCount gc = geneCounts.get(badGene);
				//yup bad guy, strip counts and delete it
				totalCounts = totalCounts - gc.getCount();
				geneCounts.remove(badGene);
			}
		}
		
	}
}
