package edu.utah.seq.data;
import java.io.*;


public class QseqFileSet {
	
	//qseq files
	private File first = null;
	private File second = null;
	private File barcode = null;
	
	//methods
	public File getFirst() {
		return first;
	}
	public void setFirst(File first) {
		this.first = first;
	}
	public File getSecond() {
		return second;
	}
	public void setSecond(File second) {
		this.second = second;
	}
	public File getBarcode() {
		return barcode;
	}
	public void setBarcode(File barcode) {
		this.barcode = barcode;
	}
	
}
