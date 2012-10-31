package trans.graphics;

import java.util.ArrayList;
import util.gen.Misc;

public class GenomicRegion {
	//fields
	private String chromosome;
	private int start;
	private int stop;
	private double score;
	private String notes;
	private ArrayList intersectingRegions = new ArrayList();
	private GenomicRegionGlyph glyph;
	private int rank;
	
	//constructor
	public GenomicRegion (String line, int rank){
		this.rank = rank;
		String[] t = line.split("\\t");
		try{
			chromosome = t[0];
			start = Integer.parseInt(t[1]);
			stop = Integer.parseInt(t[2]);
			if (t.length>3) score = Double.parseDouble(t[3]);
			if (t.length> 4) {
				StringBuffer sb = new StringBuffer(t[4]);
				for (int i=5; i< t.length; i++){
					sb.append(" ");
					sb.append(t[i]);
				}
				notes = sb.toString();
			}
		}catch(Exception e){
			Misc.printExit("\nParsing error: something is wrong with the " +
					"following region, it should be tab delimited: chrom, start, stop, " +
					"score, and optionally notes ->"+line+"\n");
		}
	}
	
	/**Returns tab delimited rank+1, chrom, st, stp, scr, notes.*/
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append(rank+1);
		sb.append("\t");
		sb.append(chromosome);
		sb.append("\t");
		sb.append(start);
		sb.append("\t");
		sb.append(stop);
		sb.append("\t");
		sb.append(score);
		sb.append("\t");
		String stripped ="";
		if (notes !=null) stripped = notes.replaceAll("\\t"," ");
		sb.append(stripped);
		return sb.toString();
	}
	
	/**checks to see if regions overlap by the minimum maxGap, can set maxGap negative to require an overlap.*/
	public boolean overlap (GenomicRegion other, int maxGap){
		//check chromosome
		if (chromosome.equals(other.chromosome) == false) return false;
		//check distances
		int oneLast = stop;
		int overlap;
		//to left
		if (oneLast< other.start){
			//check distance
			overlap = other.start-oneLast;
			if ( overlap <=maxGap )return true;
			else return false;
		}
		int twoLast = other.stop;
		//to right
		if (start> twoLast) {
			//check distance
			overlap = start-twoLast;
			if ( overlap <=maxGap ) return true;
			else return false;
		}
		//by default they overlap
		//entirely contained within?
		if ((start<=other.start && oneLast>=twoLast) || (other.start<=start && twoLast>= oneLast) ) {
			return true;
		}
		//partial overlap
		//two right of one
		if (start< other.start) overlap = other.start- oneLast; //want it negative
		//two left of one
		else overlap = start-twoLast; //want it negative
		if (overlap<= maxGap) return true;
		return false;
	}
	
	/**Returns the length of the biggest GenomicRegion*/
	public static int findBiggestGenomicRegion(GenomicRegion[] gr){
		int max = gr[0].stop - gr[0].start;
		for (int i=1; i< gr.length; i++){
			int len = gr[i].stop - gr[i].start;
			if (len > max) max = len;
		}
		return max;
	}
	
	
	public String getChromosome() {
		return chromosome;
	}
	
	public ArrayList getIntersectingRegions() {
		return intersectingRegions;
	}
	
	public String getNotes() {
		return notes;
	}
	
	public double getScore() {
		return score;
	}
	
	public int getStart() {
		return start;
	}
	
	public int getStop() {
		return stop;
	}
	
	public GenomicRegionGlyph getGlyph() {
		return glyph;
	}
	
	public void setGlyph(GenomicRegionGlyph glyph) {
		this.glyph = glyph;
	}

	public int getRank() {
		return rank;
	}

	public void setRank(int rank) {
		this.rank = rank;
	}
}