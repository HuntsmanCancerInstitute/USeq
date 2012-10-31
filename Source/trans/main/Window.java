package trans.main;
import java.io.*;

import util.gen.*;

/**
 * Class to hold information about a window, a defined portion of the genome typically < 1kb.
 */
public class Window implements Comparable, Serializable{

	//fields
	private String chromosome;
	private int start1stOligo;			//bp start of first oligo
	private int startLastOligo;			//bp start of last oligo
	private int numberOligos;			//number of oligos covered by the window, this is not the total (ie all treatments + all controls)
	private double[] scores;			//variable length see ScanChipNoPerm
	private double sortBy;				//used in sorting

	public Window (String chromosome, int start, int end, int numberOligos, double[] scores){
		this.chromosome = chromosome;
		this.start1stOligo = start;
		this.startLastOligo = end;
		this.numberOligos = numberOligos;
		this.scores = scores;
	}
	
	public int compareTo(Object obj){
		Window other = (Window)obj;
		//by score
		if (other.sortBy>sortBy) return 1;
		if (other.sortBy<sortBy) return -1;
		return 0;
	}
	
	public boolean overlap (Window other){
		//check chrom
		if (other.getChromosome().equals(chromosome) == false ) return false;
		//to left
		if (other.startLastOligo< start1stOligo) return false;
		//to right
		if (other.start1stOligo> startLastOligo) return false;
		//by default
		return true;
	}
	/**Returns chromosome, start first oligo, start last oligo*/
	public String stringRep(int sizeOfOligo){
		StringBuffer sb = new StringBuffer();
		sb.append(chromosome);
		sb.append("\t");
		sb.append(start1stOligo);
		sb.append("\t");
		sb.append(startLastOligo+sizeOfOligo);
		return sb.toString();
	}

	public String getChromosome() {
		return chromosome;
	}
	public int getStartLastOligo() {
		return startLastOligo;
	}
	public int getStart1stOligo() {
		return start1stOligo;
	}
	public void setSortBy(double sortBy) {
		this.sortBy = sortBy;
	}
	public double getSortBy() {
		return sortBy;
	}
	public void setStartLastOligo(int end) {
		this.startLastOligo = end;
	}
	public void setStart1stOligo(int start) {
		this.start1stOligo = start;
	}
	public double[] getScores() {
		return scores;
	}
	public void setScores(double[] scores) {
		this.scores = scores;
	}
	/**
	 * @return Returns the numberOligo positions covered by a window, not all the treatments and controls.
	 */
	public int getNumberOligos() {
		return numberOligos;
	}

}
