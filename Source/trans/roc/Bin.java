package trans.roc;

/**
 * Score collector for Sgr apps.
*/
public class Bin {
	private int numPositives = 0;
	private int numNegatives = 0;
	private double start;
	private double stop;
	private double fdr;
	
	public Bin(){};
	public Bin(double start, double stop){
		this.start = start;
		this.stop = stop;
	}
	
	public boolean contains(double score){
		if (score>=start && score<stop) return true;
		return false;
	}
	
	public void incrementPositives() {
		numPositives++;
	}
	public void incrementNegatives(){
		numNegatives++;
	}
	public String toString(){
		return numPositives+"\t#Pos Scores \t"+numNegatives+"\t#Neg Scores\tthat were >=\t"+start+"\t<\t"+stop+"\t%FDR: "+fdr;
	}
	public int getNumNegatives() {
		return numNegatives;
	}
	public void setNumNegatives(int numNegatives) {
		this.numNegatives = numNegatives;
	}
	public int getNumPositives() {
		return numPositives;
	}
	public void setNumPositives(int numPositives) {
		this.numPositives = numPositives;
	}
	public double getStart() {
		return start;
	}
	public void setStart(double start) {
		this.start = start;
	}
	public double getStop() {
		return stop;
	}
	public void setStop(double stop) {
		this.stop = stop;
	}
	public double getFdr() {
		return fdr;
	}
	public void setFdr(double fdr) {
		this.fdr = fdr;
	}
}
