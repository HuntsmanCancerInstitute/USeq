package trans.roc;
import java.io.*;

/**
 * For holding information about a region for ROC analysis.
 */
public class RocWindow implements Comparable, Serializable{
	//fields
	private String chromosome;
	private int start;
	private int end;
	private int middle = -1;
	private double score;
	private int type = 0; //-1 for neg, +1 pos, 0 undetermined
	
	//constructors
	public RocWindow(String chromosome, int start, int end,double score){
		this.score = score;
		this.start = start;
		this.end = end;
		this.chromosome = chromosome;
	}
	
	//methods
	public int compareTo(Object obj){
		RocWindow other = (RocWindow)obj;
		//sort by chromosome
		int compare = other.chromosome.compareTo(chromosome);
		if (compare !=0) return compare * -1;
		//sort by start position
		if (other.start>start) return -1;
		if (other.start<start) return 1;
		return 0;
	}
	
	public boolean intersects(RocWindow other){
		if (other.getChromosome().equals(chromosome)== false) return false;
		// is other left of this
		if (other.getEnd() < start) return false;
		// is other right of this
		if (other.getStart() > end ) return false;
		// must overlap
		return true;
	}
	/**Returns a binding region line, tab delimited: chrom, start, stop, score*/
	public String simpleSummaryLine(){
		StringBuffer sb = new StringBuffer();
		sb.append(chromosome);
		sb.append("\t");
		sb.append(start);
		sb.append("\t");
		if (end != start) {
			sb.append(end);
			sb.append("\t");
		}
		sb.append(score);
		return sb.toString();
	}
	
	public static void printSGRWindows(RocWindow[] w){
		int num = w.length;
		for (int i=0; i<num; i++){
			System.out.println(w[i].getChromosome()+"\t"+w[i].getMiddle()+"\t"+w[i].getScore());
		}
	}
	public String getChromosome() {
		return chromosome;
	}
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public int getMiddle() {
		if (middle == -1) {
			middle = (int)Math.round( ((double)end- (double)start)/2.0 ) + start;
		}
		return middle;
	}
	public void setMiddle(int middle) {
		this.middle = middle;
	}
	public double getScore() {
		return score;
	}
	public void setScore(double score) {
		this.score = score;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}

	public int getType() {
		return type;
	}

	public void setType(int type) {
		this.type = type;
	}
}
