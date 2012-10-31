package trans.anno;

import java.util.*;

import trans.main.*;
import util.bio.annotation.*;

/**
 * For holding information about a genomic interval, what genes are up, down, overlap, or contain this genomic inteval.
 */
public class BindingRegion implements Comparable{
	//fields
	private String chromosome;
	private int start;
	private int end;
	private int rank;
	private int neighborhood; //bp
	private double score;
	private BindingRegion closestBindingRegion;	//used by IntersectionBindingPeaks
	private double gcContent = -1;	//-1 = not set, don't change
	private String name;
	
	//5' gene, to the left
	private GeneGroup geneGrp5Prime;
	private int distanceTo5PrimeGeneGrp;
	//3' gene, to the right
	private GeneGroup geneGrp3Prime;
	private int distanceTo3PrimeGeneGrp;
	//overlap 5', left
	private ArrayList overlap5PrimeGeneGrps = new ArrayList();
	//overlap 3', right
	private ArrayList overlap3PrimeGeneGrps = new ArrayList();
	//contained by or contains
	private ArrayList containingGeneGrps = new ArrayList();
	
	//all genes within a given neighborhood
	private ArrayList neighboringGeneGrps = new ArrayList();
	
	//constructors
	public BindingRegion(int rank, double score, String chromosome, int start, int end, int sizeNeighborhood){
		this.rank = rank;
		this.score = score;
		this.start = start;
		this.end = end;
		this.chromosome = chromosome.substring(3);
		this.neighborhood = sizeNeighborhood;
	}
	
	public BindingRegion(int rank, String chromosome, int start, int end){
		this.rank = rank;
		this.start = start;
		this.end = end;
		this.chromosome = chromosome;
	}
	
	public BindingRegion(Interval interval, int rank){
		this.rank = rank;
		chromosome = interval.getChromosome();
		start = interval.getStart1stOligo();
		end = interval.getSizeOfOligoMinusOne()+ interval.getStartLastOligo();
	}
	
	//methods
	public int compareTo(Object obj){
		BindingRegion other = (BindingRegion)obj;
		//sort by chromosome		
		int compare = other.chromosome.compareTo(chromosome);
		if (compare !=0) return compare * -1;
		//sort by start position
		if (other.start>start) return -1;
		if (other.start<start) return 1;
		return 0;
	}
	
	public boolean intersects(BindingRegion other){
		if (other.getChromosome().equals(chromosome)== false) return false;
		// is other left of this
		if (other.getEnd() < start) return false;
		// is other right of this
		if (other.getStart() > end ) return false;
		// must overlap
		return true;
	}
	
	/**Returns -1 if on diff chromosomes, 0 if they overlap, or the # of bases between the two.*/
	public int intersectReturnGap(BindingRegion other){
		if (other.getChromosome().equals(chromosome)== false) return -1;
		// is other left of this
		if (other.getEnd() < start) return start - other.getEnd();
		// is other right of this
		if (other.getStart() > end ) return other.getStart() - end;
		// must overlap
		return 0;
	}
	
	/**Assumes regions are on the same chromosome, and the stop base is included, not interbase numbering.
	 * Returns 0 if the regions abut, positive ints for each base in the gap (ie 1= 1bp between an stop of 12 and a start of 14),
	 * negative ints for each base of intersection (ie -1= 1bp overlap),
	 * max negative int is the length of the smaller region (ie stop - start +1).*/
	public int bpIntersectionSameChromosome(BindingRegion other){
		// is other left of this?
		if (other.getEnd() < start) return start - other.getEnd()-1;
		//is other right of this?
		if (other.getStart() > end ) return other.getStart() - end-1;
		// must overlap
		//left side
		if (start>=other.getStart() && start<=other.getEnd() && end > other.getEnd()) return start-other.getEnd()-1;
		//right side
		if (end >= other.getStart() && end <= other.getEnd() && start < other.getStart()) return other.getStart()-end -1;
		//contained within
		int lengthThis = this.getLength();
		int lengthOther = other.getLength();
		if (lengthThis > lengthOther) return -1*lengthOther;
		return -1*lengthThis;
		
	}
	
	/**Returns -1 for no overlap or the # bases of intersection.*/
	public int intersectReturnBP(BindingRegion other){
		if (other.getChromosome().equals(chromosome)== false) return -1;
		// is other left of this
		if (other.getEnd() < start) return -1;
		// is other right of this
		if (other.getStart() > end ) return -1;
		// must overlap
		//left side
		if (start>=other.getStart() && start<=other.getEnd() && end > other.getEnd()) return other.getEnd() - start +1;
		//right side
		if (end >= other.getStart() && end <= other.getEnd() && start < other.getStart()) return end-other.getStart() +1;
		//contained within
		int lengthThis = this.getLength();
		int lengthOther = other.getLength();
		if (lengthThis > lengthOther) return lengthOther;
		return lengthThis;
	}
	
	/**Returns -1 for no overlap, 0 for complete overlap, or a positive int for the # bases of overlap.*/
	public int overlap(BindingRegion other){
		if (other.getChromosome().equals(chromosome)== false) return -1;
		// is other left of this
		if (other.getEnd() < start) return -1;
		// is other right of this
		if (other.getStart() > end ) return -1;
		// must overlap
		//left side
		if (start>=other.getStart() && start<=other.getEnd() && end > other.getEnd()) return other.getEnd() - start +1;
		//right side
		if (end >= other.getStart() && end <= other.getEnd() && start < other.getStart()) return end-other.getStart() +1;
		//contained within
		return 0;
	}
	/**Returns a binding region line, tab delimited: chrom, start, stop, score*/
	public String simpleSummaryLine(){
		StringBuffer sb = new StringBuffer();
		sb.append(chromosome);
		sb.append("\t");
		sb.append(start);
		sb.append("\t");
		if (end != start) sb.append(end);
		return sb.toString();
	}
	
	/**Returns a binding region line, tab delimited: chrom, start, stop, score*/
	public String simpleSummaryLineWithRank(){
		StringBuffer sb = new StringBuffer();
		sb.append(rank);
		sb.append("\t");
		if (name != null) {
			sb.append(name);
			sb.append("\t");
		}
		sb.append(chromosome);
		sb.append("\t");
		sb.append(start);
		sb.append("\t");
		if (end != start) sb.append(end);
		return sb.toString();
	}
	
	/**Returns tab delimited line: rank, chrom, start, stop, dist to closes gene, gene text(s), # neighbors*/
	public String summaryLine(){
		StringBuffer sb = new StringBuffer();
		sb.append(rank);
		sb.append("\t");
		sb.append(chromosome);
		sb.append("\t");
		sb.append(start);
		sb.append("\t");
		sb.append(end);
		sb.append("\t");
		//find dist to closest gene
		//any overlaps or contained withins
		StringBuffer closestGeneNames = new StringBuffer();
		boolean overlaps = false;
		if (overlap5PrimeGeneGrps.size()!=0) {
			closestGeneNames.append(extractGeneGrps(overlap5PrimeGeneGrps, false)+" ");
			overlaps = true;
		}
		if (overlap3PrimeGeneGrps.size()!=0) {
			closestGeneNames.append(extractGeneGrps(overlap3PrimeGeneGrps, false)+" ");
			overlaps = true;
		}
		if (containingGeneGrps.size()!=0) {
			closestGeneNames.append(extractGeneGrps(containingGeneGrps, false)+" ");
			overlaps = true;
		}
		//any overlaps? if not then find closest 5' or 3'
		if (overlaps){
			sb.append("0");
			sb.append("\t");
			sb.append(closestGeneNames);
		}
		else {
			//no 5'
			if (geneGrp5Prime==null){
				sb.append(distanceTo3PrimeGeneGrp);
				sb.append("\t");
				sb.append(geneGrp3Prime.getName()+" ");
			}
			//no 3'
			else if (geneGrp3Prime==null){
				sb.append(distanceTo5PrimeGeneGrp);
				sb.append("\t");
				sb.append(geneGrp5Prime.getName()+" ");
			}
			else if (distanceTo5PrimeGeneGrp < distanceTo3PrimeGeneGrp) {
				sb.append(distanceTo5PrimeGeneGrp);
				sb.append("\t");
				sb.append(geneGrp5Prime.getName()+" ");
			}
			else {
				sb.append(distanceTo3PrimeGeneGrp);
				sb.append("\t");
				sb.append(geneGrp3Prime.getName()+" ");
			}
		}
		sb.append("\t");
		//number neighbors
		sb.append(neighboringGeneGrps.size());
		
		return sb.toString();
	}
	
	public String toString(){
		//region line
		StringBuffer sb = new StringBuffer("GenomicRegion: ");
		sb.append(rank);
		sb.append(", chr");
		sb.append(chromosome);
		sb.append(":");
		sb.append(start);
		sb.append("-");
		sb.append(end);
		sb.append(", (+/-5kb) chr");
		sb.append(chromosome);
		sb.append(":");
		sb.append(start-5000);
		sb.append("-");
		sb.append(end+5000);
		sb.append("\n");
		//5' 
		if (geneGrp5Prime!=null){
			sb.append("5' Gene: ");
			sb.append(distanceTo5PrimeGeneGrp);
			sb.append("bp\t");
			sb.append(geneGrp5Prime.fetchSumaryLine());
			sb.append("\n");
		}
		//5' Overlap
		if (overlap5PrimeGeneGrps.size() !=0) {
			sb.append("5' Overlap:\n");
			sb.append (extractGeneGrps(overlap5PrimeGeneGrps, true));
			sb.append("\n");
		}
		//contain
		if (containingGeneGrps.size()!=0) {
			sb.append("Contained By:\n");
			sb.append (extractGeneGrps(containingGeneGrps, true));
			sb.append("\n");
		}
		//3' Overlap
		if (overlap3PrimeGeneGrps.size()!=0) {
			sb.append("3' Overlap:\n");
			sb.append (extractGeneGrps(overlap3PrimeGeneGrps, true));
			sb.append("\n");
		}
		//3'	
		if (geneGrp3Prime != null){
			sb.append("3' Gene: ");
			sb.append(distanceTo3PrimeGeneGrp);
			sb.append("bp\t");
			sb.append(geneGrp3Prime.fetchSumaryLine());
		}
		//neighbors
		int numNeighbors = neighboringGeneGrps.size();
		if (numNeighbors !=0){
			sb.append("\nNeighbors: ");
			sb.append( ((GeneGroup)neighboringGeneGrps.get(0)).getName() );
			for (int i=1; i<numNeighbors; i++){
				sb.append(", ");
				sb.append( ((GeneGroup)neighboringGeneGrps.get(i)).getName() );
			}
		}
		return sb.toString();
	}
	/**Returns summary line or text for each GeneGrp in ArrayList*/
	public static String extractGeneGrps(ArrayList geneGrps, boolean returnSummaryLine){
		StringBuffer sb = new StringBuffer();
		int num = geneGrps.size();
		GeneGroup gene;
		gene = (GeneGroup) geneGrps.get(0);
		if (returnSummaryLine){
			sb.append("\t"+gene.fetchSumaryLine());
			for (int i=1; i< num; i++){
				sb.append("\n\t");
				gene = (GeneGroup) geneGrps.get(i);
				sb.append(gene.fetchSumaryLine());
			}
		}
		else{
			sb.append(gene.getName());
			for (int i=1; i< num; i++){
				sb.append(",");
				gene = (GeneGroup) geneGrps.get(i);
				sb.append(gene.getName());
			}
		}
		return sb.toString();
	}
	
	
	
	/**Compares a binding region to a gene grp looking for where it is in relation to the binding 
	 * region, 5', 3', overlap, contained.
	 * Saves all that overlap or are contained as well as the closest 5' and the closest 3'.*/
	public void compare(GeneGroup geneGrp){
		if (chromosome.equalsIgnoreCase(geneGrp.getChromosome())== false) return;
		//start stop
		int[] se = new int[]{geneGrp.getStart(),geneGrp.getEnd()} ;
		
		//where is geneGrp in relation to this binding region?
		//left
		if (start>se[1]){
			//neighbors?
			if ((start-se[1]) <= neighborhood) neighboringGeneGrps.add(geneGrp);
			//best?
			if (geneGrp5Prime == null){
				geneGrp5Prime = geneGrp;
				distanceTo5PrimeGeneGrp = start-se[1];
			}
			else{
				int newDistance = start-se[1];
				if (newDistance < distanceTo5PrimeGeneGrp){
					geneGrp5Prime = geneGrp;
					distanceTo5PrimeGeneGrp = newDistance;
				}
			}
		}
		//right
		else if (end<se[0]){
			//neighbors?
			if ((se[0]-end) <= neighborhood) neighboringGeneGrps.add(geneGrp);
			//best?
			if (geneGrp3Prime == null){
				geneGrp3Prime = geneGrp;
				distanceTo3PrimeGeneGrp = se[0]-end;
			}
			else{
				int newDistance = se[0]-end;
				if (newDistance < distanceTo3PrimeGeneGrp){
					geneGrp3Prime = geneGrp;
					distanceTo3PrimeGeneGrp = newDistance;
				}
			}
		}
		//left overlap
		else if (se[0]<start && start<=se[1] && end>=se[1]){
			neighboringGeneGrps.add(geneGrp);
			overlap5PrimeGeneGrps.add(geneGrp);
		}
		
		//right overlap
		else if (start <= se[0] && end >= se[0] && end < se[1]){
			neighboringGeneGrps.add(geneGrp);
			overlap3PrimeGeneGrps.add(geneGrp);
		}
		//contain by or contains
		else if ( (se[0]<=start && end <= se[1]) || (start<=se[0] && end>= se[1]) ) {
			neighboringGeneGrps.add(geneGrp);
			containingGeneGrps.add(geneGrp);
		}
		else {
			System.out.println("Problem, no comparison appropriate!");
			System.exit(0);
		}
	}
	
	public String getChromosome() {
		return chromosome;
	}
	public ArrayList getContainingGeneGrps() {
		return containingGeneGrps;
	}
	public int getDistanceTo3PrimeGeneGrp() {
		return distanceTo3PrimeGeneGrp;
	}
	public int getDistanceTo5PrimeGeneGrp() {
		return distanceTo5PrimeGeneGrp;
	}
	public GeneGroup getGeneGrp3Prime() {
		return geneGrp3Prime;
	}
	public GeneGroup getGeneGrp5Prime() {
		return geneGrp5Prime;
	}
	public ArrayList getOverlap3PrimeGeneGrps() {
		return overlap3PrimeGeneGrps;
	}
	public ArrayList getOverlap5PrimeGeneGrps() {
		return overlap5PrimeGeneGrps;
	}
	public int getRank() {
		return rank;
	}
	public int getStart() {
		return start;
	}
	public int getEnd() {
		return end;
	}
	public int getNeighborhood() {
		return neighborhood;
	}
	public void setNeighborhood(int neighborhood) {
		this.neighborhood = neighborhood;
	}
	public ArrayList getNeighboringGeneGrps() {
		return neighboringGeneGrps;
	}
	public double getScore() {
		return score;
	}
	public BindingRegion getClosestBindingRegion() {
		return closestBindingRegion;
	}
	public void setClosestBindingRegion(BindingRegion closestBindingRegion) {
		this.closestBindingRegion = closestBindingRegion;
	}
	/**End-Start*/
	public int getLength(){
		return end-start;
	}

	public double getGcContent() {
		return gcContent;
	}

	public void setGcContent(double gcContent) {
		this.gcContent = gcContent;
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
}
