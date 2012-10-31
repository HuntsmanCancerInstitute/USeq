package trans.main;

import java.util.ArrayList;

/**
 * For holding information about a genomic interval, what genes are up, down, overlap, or contain this genomic inteval.
 */
public class PosRegion implements Comparable{
	//fields
	private String chromosome;
	private int start;
	private int end;
	private String name;
	private PosRegion region;
	private ArrayList overlapping = new ArrayList();
	private int closestNonOverlappingInterval = -1;
	private Interval closestInterval;
	
	//constructors
	public PosRegion(String name, String chromosome, int start, int end){
		this.name = name;
		this.start = start;
		this.end = end;
		this.chromosome = chromosome;
	}
	
	//methods
	public int compareTo(Object obj){
		PosRegion other = (PosRegion)obj;
		//sort by chromosome
		int compare = other.chromosome.compareTo(chromosome);
		if (compare !=0) return compare * -1;;
		//sort by start position
		if (other.start>start) return -1;
		if (other.start<start) return 1;
		return 0;
	}
	/**Returns tab delimited line: rank, chrom, start, stop, dist to closes gene, gene text(s), # neighbors*/
	public String summaryLine(){
		StringBuffer sb = new StringBuffer();
		sb.append(name);
		sb.append("\t");
		sb.append(chromosome);
		sb.append("\t");
		sb.append(start);
		sb.append("\t");
		sb.append(end);
		return sb.toString();
	}
	public String getChromosome() {
		return chromosome;
	}
	public int getStart() {
		return start;
	}
	public int getEnd() {
		return end;
	}
	public int getClosestNonOverlappingInterval() {
		return closestNonOverlappingInterval;
	}
	public void setClosestNonOverlappingInterval(int closestOverlap) {
		this.closestNonOverlappingInterval = closestOverlap;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public ArrayList getOverlapping() {
		return overlapping;
	}
	public void setOverlapping(ArrayList overlapping) {
		this.overlapping = overlapping;
	}
	public PosRegion getRegion() {
		return region;
	}
	public void setRegion(PosRegion region) {
		this.region = region;
	}
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public Interval getClosestInterval() {
		return closestInterval;
	}
	public void setClosestInterval(Interval closestInterval) {
		this.closestInterval = closestInterval;
	}
}
