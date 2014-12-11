package edu.utah.seq.data;

import java.io.File;
import java.io.Serializable;

/**Class for saving minimal info about a sam alignment for generating read coverage tracks*/
public class ChromDataSave implements Serializable{

	//fields
	private int firstBase;
	private int lastBase;
	private File binaryFile;
	private String chromosome;
	private String strand;
	private static final long serialVersionUID = 1L;

	//constructors
	public ChromDataSave (ChromData cd){
		this.firstBase = cd.firstBase;
		this.lastBase = cd.lastBase;
		this.chromosome = cd.chromosome;
		this.strand = cd.strand;
		this.binaryFile = cd.binaryFile;
	}
	
	public ChromData getChromData(){
		ChromData cd = new ChromData();
		cd.firstBase = this.firstBase;
		cd.lastBase = this.lastBase;
		cd.binaryFile = this.binaryFile;
		cd.chromosome = this.chromosome;
		cd.strand = this.strand;
		return cd;
	}
}

