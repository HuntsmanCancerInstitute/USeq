package trans.main;
import java.util.*;
import java.io.*;

import util.gen.*;

/**
 * Container to hold information about a particular oligo.
 */
public class Oligo implements Serializable{
	//fields
	private int index;
	private int start;				//bp position
	private int matches;				//number of matches to genome
	private ArrayList intensities;		//Float intensity value, one for each chip, 
	//	first treatment values, then control values, can find out how many of each by looking at interval
	private String sequence;
	//for making binding peaks
	private double smoothedScore;		//used to represent this oligo
	private boolean skip = false; 	//used in looking for binding sites
	private boolean inPeak = false;
	private double rightSlope = 0;
	private double leftSlope = 0;
	
	public Oligo(int index, int start, int hits, ArrayList intensities, String sequence){
		this.index = index;
		this.start = start;
		this.matches = hits;
		this.intensities = intensities;
		this.sequence = sequence;
	}
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append("index\t");
		sb.append(index);
		sb.append("\nstart\t");
		sb.append(start);
		sb.append("\nsequence\t");
		sb.append(sequence);
		if (matches !=0 ){
			sb.append("\nMatches To Genome\t");
			sb.append(matches);
		}
		sb.append("\nintensities");
		sb.append(intensities);
		return sb.toString();
	}
	public String toString(int numTreatInt, int numControlInt, String divider){
		StringBuffer sb = new StringBuffer();
		sb.append(sequence);
		sb.append(divider);
		if (matches!=0){
			sb.append(matches);
			sb.append(divider);
		}
		sb.append(start);
		sb.append(divider);
		sb.append(Misc.floatArrayToString(getTreatmentIntensities(numTreatInt), ","));
		sb.append(divider);
		sb.append(Misc.floatArrayToString(getControlIntensities(numControlInt), ","));
		sb.append(divider);
		return sb.toString();
	}
	
	public float[] getTreatmentIntensities(int numTreatInt){
		float[] treatments = new float[numTreatInt];
		for (int i=0; i<numTreatInt; i++) treatments[i] = ((Float)intensities.get(i)).floatValue();
		return treatments;
	}
	public float[] getControlIntensities(int numControlInt){
		float[] controls = new float[numControlInt];
		int max = intensities.size();
		int counter =0;
		for (int i=max-numControlInt; i<max; i++) controls[counter++] = ((Float)intensities.get(i)).floatValue();
		return controls;
	}
	
	public int getStart() {
		return start;
	}
	public ArrayList getIntensities() {
		return intensities;
	}
	public int getIndex() {
		return index;
	}
	public String getSequence() {
		return sequence;
	}
	public double getSmoothedScore() {
		return smoothedScore;
	}
	public void setSmoothedScore(double smoothedScore) {
		this.smoothedScore = smoothedScore;
	}
	public boolean isInPeak() {
		return inPeak;
	}
	public void setInPeak(boolean inPeak) {
		this.inPeak = inPeak;
	}
	public boolean skip() {
		return skip;
	}
	public void setSkip(boolean skip) {
		this.skip = skip;
	}
	public double getLeftSlope() {
		return leftSlope;
	}
	public void setLeftSlope(double leftSlope) {
		this.leftSlope = leftSlope;
	}
	public double getRightSlope() {
		return rightSlope;
	}
	public void setRightSlope(double rightSlope) {
		this.rightSlope = rightSlope;
	}
	public int getMatches() {
		return matches;
	}
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
}
