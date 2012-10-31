package edu.utah.seq.data;
import java.util.*;
import java.io.*;

import edu.utah.seq.parsers.BarParser;

/**Data file information.*/
public class Info implements Serializable{
	//fields
	private String name;
	
	/**As defined by UCSC (ie hg18, dm2, ce2, mm8), see http://genome.ucsc.edu/FAQ/FAQreleases*/
	private String versionedGenome;
	
	/**As defined by UCSC see http://genome.ucsc.edu/ and particular species.*/
	private String chromosome;
	
	/**+, -, or .*/
	private String strand;
	
	/**
	 * Length of the individual reads*/
	private int readLength;
	
	/**Will be set without having to actually load the data*/
	private int numberObservations = 0;
	private double scoreTotal = 0;
	
	/**Optional, can be null.*/
	private HashMap<String,String> notes;
	
	//change this value when making major changes that should void any saved Objects
	public static final long serialVersionUID = 1;

	
	//constructor
	/**Notes, can be null, everything else is required.
	 * @param versionedGenome - As defined by UCSC (ie hg18, dm2, ce2, mm8), see http://genome.ucsc.edu/FAQ/FAQreleases
	 * @param chromosome - As defined by UCSC see http://genome.ucsc.edu/ and particular species.
	 * @param strand - +, -, or .
	 * @param readLength - Length of the individual reads.
	 * @param notes - optional can be null.*/
	public Info (String name, String versionedGenome, String chromosome, String strand, int readLength, HashMap<String,String> notes){
		this.name = name;
		this.versionedGenome = versionedGenome;
		this.chromosome = chromosome;
		this.strand = strand;
		this.readLength = readLength;
		this.notes = notes;
	}

	public Info (BarParser barParser){
		//attempt to find length of reads
		notes = barParser.getTagValues();
		if (notes.containsKey(BarParser.READ_LENGTH_TAG)) readLength = Integer.parseInt(notes.get(BarParser.READ_LENGTH_TAG));
		//attempt to find total score
		if (notes.containsKey(BarParser.SCORE_TOTAL)) scoreTotal = Double.parseDouble(notes.get(BarParser.SCORE_TOTAL));
		//set number of positionValues
		numberObservations = barParser.getNumberPositionValues();
		name = barParser.getBarFile().toString();
		versionedGenome = barParser.getVersionedGenome();
		chromosome = barParser.getChromosome();
		strand = barParser.getStrand();
	}
	public Info (){}
	
	//methods
	/**Returns a copy of this Info, the notes is set to null;*/
	public Info copy(){
		Info i = new Info(new String(name), new String(versionedGenome), new String(chromosome), new String(strand), readLength, null);
		return i;
	}
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append("Name:\t"+name);
		sb.append("\nVersionedGenome:\t"+versionedGenome);
		sb.append("\nChromosome:\t"+chromosome);
		sb.append("\nStrand:\t"+strand);
		sb.append("\nReadLength:\t"+readLength);
		sb.append("\nNumObs:\t"+numberObservations);
		if (notes != null) sb.append("\nNotes:\t"+notes);
		sb.append("\n");
		return sb.toString();
	}
	public String getChromosome() {
		return chromosome;
	}

	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public HashMap<String,String> getNotes() {
		return notes;
	}

	public void setNotes(HashMap<String,String> notes) {
		this.notes = notes;
	}

	public String getStrand() {
		return strand;
	}

	public void setStrand(String strand) {
		this.strand = strand;
	}

	public String getVersionedGenome() {
		return versionedGenome;
	}

	public void setVersionedGenome(String versionedGenome) {
		this.versionedGenome = versionedGenome;
	}

	public int getReadLength() {
		return readLength;
	}

	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}

	public int getNumberObservations() {
		return numberObservations;
	}

	public void setNumberObservations(int numberObservations) {
		this.numberObservations = numberObservations;
	}

	public double getScoreTotal() {
		return scoreTotal;
	}

	public void setScoreTotal(double scoreTotal) {
		this.scoreTotal = scoreTotal;
	}

}
