package edu.utah.seq.methylation454;
import util.gen.*;
import java.io.*;

public class DataLine implements Comparable{

	//fields
	private String sample;
	private String amplicon;
	private String cpG;
	private String read;
	private String strand;
	private String sequence;
	
	//constructor
	public DataLine (String line){
		String[] cells = line.split("\\t");
		if (cells.length != 6) Misc.printExit("\nProblem parsing line -> "+line+"\n");
		sample = cells[0];
		amplicon = cells[1];
		cpG = cells[2];
		read = cells[3];
		strand = cells[4];
		sequence = cells[5];
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		String tab = "\t";
		sb.append(sample); sb.append(tab);
		sb.append(amplicon); sb.append(tab);
		sb.append(cpG); sb.append(tab);
		sb.append(read); sb.append(tab);
		sb.append(strand); sb.append(tab);
		sb.append(sequence);
		return sb.toString();
	}
	
	public int compareTo (Object obj){
		DataLine other = (DataLine) obj;
		//sort by sample
		int t = other.sample.compareTo(sample) * -1;
		if (t !=0 ) return t;
		//sort by amplicon
		t = other.amplicon.compareTo(amplicon);
		if (t !=0 ) return t;
		//sort by cpG
		t = other.cpG.compareTo(cpG) * -1;
		if (t !=0 ) return t;
		//sort by strand
		t = other.strand.compareTo(strand);
		if (t !=0 ) return t;
		//sort by sequence
		t = other.sequence.compareTo(sequence);
		if (t !=0 ) return t;
		return 0;
	}

	public String getAmplicon() {
		return amplicon;
	}

	public String getCpG() {
		return cpG;
	}

	public String getRead() {
		return read;
	}

	public String getSample() {
		return sample;
	}

	public String getSequence() {
		return sequence;
	}

	public String getStrand() {
		return strand;
	}
}
