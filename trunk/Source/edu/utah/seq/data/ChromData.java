package edu.utah.seq.data;

import java.io.DataOutputStream;
import java.io.File;
import java.io.Serializable;

/**Class for minimal info about sam alignments for generating read coverage tracks*/
public class ChromData {

	//fields
	public int firstBase;
	public int lastBase;
	public DataOutputStream out;
	File binaryFile;
	String chromosome;
	String strand;

	//constructors
	public ChromData (int firstBase, int lastBase, String chromosome, String strand, File binaryFile, DataOutputStream out){
		this.firstBase = firstBase;
		this.lastBase = lastBase;
		this.chromosome = chromosome;
		this.strand = strand;
		this.binaryFile = binaryFile;
		this.out = out;
	}

	public ChromData() {}
}