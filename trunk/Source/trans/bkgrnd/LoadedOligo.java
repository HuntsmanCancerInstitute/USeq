package trans.bkgrnd;
import util.gen.*;

/**
 * BPMap line contain gc and tm calculations.
 */
public class LoadedOligo implements Comparable{
	//fields
	private String chromosome;
	private int startPosition;
	private String sequence;
	private int numberMisMatches;
	private double intensityPM;
	private double intensityMM;
	private double tm;
	private double gc;
	
	//constructors
	public LoadedOligo (String chromosome, int startPosition, String sequence, double intensityPM, 
			double intensityMM, double tm, double gc){
		this.chromosome = chromosome;
		this.startPosition = startPosition;
		this.sequence = sequence;
		this.intensityPM = intensityPM;
		this.intensityMM = intensityMM;
		this.tm = tm;
		this.gc = gc;
	}

	public LoadedOligo(String line){
		String[] tokens = line.split("\\s");
		//oligo sequence, number of 1bp mismatches, chromosome, start position, PM intensity, MM intensity, tm, gc
		sequence = tokens[0];
		numberMisMatches = Num.parseInt(tokens[1], -1);
		chromosome = tokens[2];
		startPosition = Integer.parseInt(tokens[3]);
		intensityPM = Double.parseDouble(tokens[4]);
		intensityMM = Double.parseDouble(tokens[5]);
		tm = Double.parseDouble(tokens[6]);
		gc = Double.parseDouble(tokens[7]);
	}
	
	public int compareTo(Object obj){
		LoadedOligo other = (LoadedOligo)obj;
		//sort by chromosome
		int comp = other.chromosome.compareTo(chromosome);
		if (comp != 0) return comp*-1;
		//sort by base position
		if (other.startPosition>startPosition) return -1;
		if (other.startPosition<startPosition) return 1;
		return 0;
	}
	public String getChromosome() {
		return chromosome;
	}
	public double getGc() {
		return gc;
	}
	public double getIntensityMM() {
		return intensityMM;
	}
	public double getIntensityPM() {
		return intensityPM;
	}
	public String getSequence() {
		return sequence;
	}
	public int getStartPosition() {
		return startPosition;
	}
	public double getTm() {
		return tm;
	}
	public int getNumberMisMatches() {
		return numberMisMatches;
	}
	public void setIntensityPM(double intensityPM) {
		this.intensityPM = intensityPM;
	}
}
