package trans.roc;
import java.util.*;

import util.gen.*;
/**
 * Container for a positive BAC used by the Sgr scripts.
 *
 */
public class Positive {
	private String chromosome;
	private int start;
	private int stop;
	private ArrayList scores = new ArrayList();
	private int totalNumberWindows = 0;
	private boolean senseStrand = true;
	private int index = 0;
	
	public Positive(String chromosome, int start, int stop){
		this.chromosome =  chromosome;
		this.start = start;
		this.stop = stop;
	}
	
	public Positive(int index, String chromosome, int start, int stop){
		this.index = index;
		this.chromosome =  chromosome;
		this.start = start;
		this.stop = stop;
	}
	
	public boolean matches(Sgr sgr){
		if (sgr.getChromosome().equals(chromosome) && sgr.getPosition()>=start && sgr.getPosition()<=stop) {
			return true;
		}
		return false;
	}
	
	public boolean contains(int position){
		if (position>= start && position <= stop) return true;
		return false;
	}

	public ArrayList getScores() {
		return scores;
	}
	public String toStringSimple(){
		StringBuffer sb = new StringBuffer();
		sb.append(chromosome);
		sb.append(":");
		sb.append(start);
		sb.append("-");
		sb.append(stop);
		return sb.toString();
	}
	
	public String toStringSimpleFloat(){
		StringBuffer sb = new StringBuffer();
		sb.append(chromosome);
		sb.append("\t");
		sb.append(start);
		sb.append("\t");
		sb.append(stop);
		float[] f = (float[])scores.get(0);
		for (int i=0; i<f.length; i++){
			sb.append("\t");
			sb.append(f[i]);
		}
		return sb.toString();
	}
	
	/**Adds the bpOffSetEnds to the stop of each Positive.
	 * If the start < modified stop, the Pos is saved.*/
	public static Positive[] filter(int bpOffSetEnds, Positive[] pos){
		ArrayList al = new ArrayList();
		for (int i=0; i< pos.length; i++){
			int start = pos[i].getStart();
			int end = pos[i].getStop() + bpOffSetEnds;
			if (start< end) {
				pos[i].stop = end;
				al.add(pos[i]);
			}
		}
		Positive[] filtered = new Positive[al.size()];
		al.toArray(filtered);
		return filtered;
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		int numScores = scores.size();
		double[] s = new double[numScores];
		for (int i=0; i<numScores; i++){
			s[i]= ((Double)scores.get(i)).doubleValue();
		}
		Arrays.sort(s);	//needed for median calc
		sb.append(chromosome);
		sb.append(":");
		sb.append(start);
		sb.append("-");
		sb.append(stop);
		sb.append("\nNumber of Scores: ");
		sb.append(numScores);
		sb.append("\nMean: ");
		sb.append(Num.mean(s));
		sb.append("\nMedian: ");
		sb.append(Num.median(s));
		sb.append("\nStandard Deviation: ");
		sb.append(Num.standardDeviation(s));
		sb.append("\nScores: ");
		sb.append(scores);
		sb.append("\n");
		return sb.toString();
	}
	public String getChromosome() {
		return chromosome;
	}
	public int getStart() {
		return start;
	}
	public int getStop() {
		return stop;
	}
	public int getTotalNumberWindows() {
		return totalNumberWindows;
	}
	public void setTotalNumberWindows(int totalNumberWindows) {
		this.totalNumberWindows = totalNumberWindows;
	}
	public int getLength(){
		return stop-start;
	}
	public void setScores(ArrayList scores) {
		this.scores = scores;
	}

	public boolean isSenseStrand() {
		return senseStrand;
	}

	public void setSenseStrand(boolean senseStrand) {
		this.senseStrand = senseStrand;
	}

	public int getIndex() {
		return index;
	}

	public void setIndex(int index) {
		this.index = index;
	}
}
