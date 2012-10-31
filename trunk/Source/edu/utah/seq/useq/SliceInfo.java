package edu.utah.seq.useq;

import java.io.IOException;
import java.util.regex.*;


/**Container to hold chrom, strand, firstStartPosition, lastStartPosition, numberRecords, and binaryType from a useq data slice. 
 * Each slice name should look like chr5+8394834-8394854-10000.isft, no spaces, interbase coordinates, strand is +,-,or .
 * 
 * @author david.nix@hci.utah.edu*/
public class SliceInfo {
	
	//fields
	public static final Pattern SLICE_NAME_SPLITTER = Pattern.compile("^(.+)([+-.])(\\d+)-(\\d+)-(\\d+)\\.(\\w+)$");
	/**Currently not used for anything, defaults to ""*/
	private String notes = "";
	
	//required fields
	private String chromosome = null;
	/**+, -, or .*/
	private String strand = null;
	/**First start bp position, interbase coordinates*/
	private int firstStartPosition;
	/**Last start bp postion, interbase coordinates*/
	private int lastStartPosition;
	private int numberRecords;
	/**An ordered list of lower case letters describing the binary contents of a single record slice. 
	 * Boolean (o), Byte(b), Short(s), Integer(i), Long(l), Float(f), Double(d), Char(c), String/Text encoded in UTF-8(t), 
	 * For example, use 'isft' to describe a file that contains records with int,short,float,text*/
	private String binaryType = null;

	
	//constructors
	public SliceInfo (String sliceName) throws IOException {
		parseSliceName(sliceName);
	}
	public SliceInfo (String chromosome, String strand, int firstStartPosition, int lastStartPosition, int numberRecords, String binaryType){
		this.chromosome = chromosome;
		this.strand = strand;
		this.firstStartPosition = firstStartPosition;
		this.lastStartPosition = lastStartPosition;
		this.numberRecords = numberRecords;
		this.binaryType = binaryType;
	}
	
	//methods
	/**Rips and loads the info from the slice name into this object.  Throws and IOException if the name is malformed*/
	public void parseSliceName(String sliceName) throws IOException {
		Matcher mat = SLICE_NAME_SPLITTER.matcher(sliceName);
		if (mat.matches() == false) throw new IOException ("Malformed slice name! Failed to parse the slice info from -> "+sliceName);
		chromosome = mat.group(1);
		strand = mat.group(2);
		//no need to catch NumberFormatException, this is checked by the NAME_SPLITTER Pattern
		firstStartPosition = Integer.parseInt(mat.group(3));
		lastStartPosition = Integer.parseInt(mat.group(4));
		numberRecords = Integer.parseInt(mat.group(5));
		binaryType = mat.group(6);
	}
	
	/**Returns a properly formated name, chr5+8394834-8394854-10000.isft*/
	public String getSliceName(){
		return chromosome+strand+firstStartPosition+"-"+lastStartPosition+"-"+numberRecords+"."+binaryType;
	}

	//getters setters
	public String getNotes() {
		return notes;
	}

	public void setNotes(String notes) {
		this.notes = notes;
	}

	public static Pattern getSliceNameSplitter() {
		return SLICE_NAME_SPLITTER;
	}

	public String getChromosome() {
		return chromosome;
	}

	public String getStrand() {
		return strand;
	}

	public int getFirstStartPosition() {
		return firstStartPosition;
	}

	public int getLastStartPosition() {
		return lastStartPosition;
	}

	public int getNumberRecords() {
		return numberRecords;
	}

	public String getBinaryType() {
		return binaryType;
	}
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public void setStrand(String strand) {
		this.strand = strand;
	}
	public void setFirstStartPosition(int firstStartPosition) {
		this.firstStartPosition = firstStartPosition;
	}
	public void setLastStartPosition(int lastStartPosition) {
		this.lastStartPosition = lastStartPosition;
	}
	public void setNumberRecords(int numberRecords) {
		this.numberRecords = numberRecords;
	}
	public void setBinaryType(String binaryType) {
		this.binaryType = binaryType;
	}
	public boolean isContainedBy(int beginningBP, int endingBP) {
		if (firstStartPosition >= beginningBP && firstStartPosition < endingBP && lastStartPosition >= beginningBP && lastStartPosition < endingBP) return true;
		return false;
	}
}
