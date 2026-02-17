package edu.utah.seq.maf;

import java.util.ArrayList;
import java.util.TreeSet;

import util.gen.Misc;

public class RasStats {

	//fields
	private String rasName = null;
	private ArrayList<String[]> mafRasLines = null;
	private TreeSet<String> observedRaspDots = null;
	private int numberSamplesRasAny = -1;
	private int numberSamplesRasCanonical = -1;
	private int numberSamplesRasNonCanonical = -1;
	private int numberSamplesMAP2K_RasCanonical = -1;
	private int numberSamplesMAP2K_RasNonCanonical = -1;
	private String header = null;
	
	private double canOvPval = -1;
	private double canUnPval = -1;
	private double nonCanOvPval = -1;
	private double nonCanUnPval = -1;
	
	public void setPValues(double canOvPval, double canUnPval, double nonCanOvPval, double nonCanUnPval) {
		this.canOvPval = canOvPval;
		this.canUnPval = canUnPval;
		this.nonCanOvPval = nonCanOvPval;
		this.nonCanUnPval = nonCanUnPval;
	}
	//constructor
	public RasStats(String rasName) {
		this.rasName = rasName;
		header = "Obs "+rasName+" pDots\tNum "+rasName+" Mut Spls\tNum Can "+rasName+" Mut Spls\tNum Non Can "+rasName+
				" Mut Spls\tNum MAP "+rasName+" Can Mut Spls\tFrac MAP "+rasName+" Can Mut\tMAP "+rasName+" Can Over PVal\tMAP "+
				rasName+" Can Under PVal\tNum MAP "+rasName+" Non Can Mut Spls\tFrac MAP "+rasName+" Non Can Mut\tMAP "+rasName+
				" Non Can Over PVal\tMAP "+rasName+" Non Can Under PVal";
	}
	
	public void appendInfo(StringBuilder sb, double numberSamples) {
		sb.append(Misc.treeSetToString(observedRaspDots, ", ")); sb.append("\t");
		sb.append(numberSamplesRasAny); sb.append("\t");
		sb.append(numberSamplesRasCanonical); sb.append("\t");
		sb.append(numberSamplesRasNonCanonical); sb.append("\t");
		sb.append(numberSamplesMAP2K_RasCanonical); sb.append("\t");
		sb.append((double)numberSamplesMAP2K_RasCanonical/numberSamples); sb.append("\t");
		sb.append(canOvPval); sb.append("\t");
		sb.append(canUnPval); sb.append("\t");
		sb.append(numberSamplesMAP2K_RasNonCanonical); sb.append("\t");
		sb.append((double)numberSamplesMAP2K_RasNonCanonical/numberSamples); sb.append("\t");
		sb.append(nonCanOvPval); sb.append("\t");
		sb.append(nonCanUnPval); 
	}

	public String getRasName() {
		return rasName;
	}

	public ArrayList<String[]> getMafRasLines() {
		return mafRasLines;
	}

	public void setMafRasLines(ArrayList<String[]> mafRasLines) {
		this.mafRasLines = mafRasLines;
	}

	public TreeSet<String> getObservedRaspDots() {
		return observedRaspDots;
	}

	public void setObservedRaspDots(TreeSet<String> observedRaspDots) {
		this.observedRaspDots = observedRaspDots;
	}

	public int getNumberSamplesRasAny() {
		return numberSamplesRasAny;
	}

	public void setNumberSamplesRasAny(int numberSamplesRasAny) {
		this.numberSamplesRasAny = numberSamplesRasAny;
	}

	public int getNumberSamplesRasCanonical() {
		return numberSamplesRasCanonical;
	}

	public void setNumberSamplesRasCanonical(int numberSamplesRasCanonical) {
		this.numberSamplesRasCanonical = numberSamplesRasCanonical;
	}

	public int getNumberSamplesRasNonCanonical() {
		return numberSamplesRasNonCanonical;
	}

	public void setNumberSamplesRasNonCanonical(int numberSamplesRasNonCanonical) {
		this.numberSamplesRasNonCanonical = numberSamplesRasNonCanonical;
	}

	public int getNumberSamplesMAP2K_RasCanonical() {
		return numberSamplesMAP2K_RasCanonical;
	}

	public void setNumberSamplesMAP2K_RasCanonical(int numberSamplesMAP2K_RasCanonical) {
		this.numberSamplesMAP2K_RasCanonical = numberSamplesMAP2K_RasCanonical;
	}

	public int getNumberSamplesMAP2K_RasNonCanonical() {
		return numberSamplesMAP2K_RasNonCanonical;
	}

	public void setNumberSamplesMAP2K_RasNonCanonical(int numberSamplesMAP2K_RasNonCanonical) {
		this.numberSamplesMAP2K_RasNonCanonical = numberSamplesMAP2K_RasNonCanonical;
	}

	public String getHeader() {
		return header;
	}
}
