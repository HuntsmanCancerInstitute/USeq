package util.bio.annotation;

import java.util.ArrayList;
import java.util.HashMap;

public class ScoredCoordinate extends Coordinate{
	//fields
	private double score;  //or rank
	
	//constructor
	public ScoredCoordinate (String chromosome, int start, int stop, double score){
		super(chromosome, start, stop);
		this.score = score;
	}

	//methods
	/**Returns an array of ScoredCoordinate where the score is the index position of the Coordinate[].*/
	public static ScoredCoordinate[] makeRanked (Coordinate[] c){
		ScoredCoordinate[] sc = new ScoredCoordinate[c.length];
		for (int i=0; i< sc.length; i++){
			sc[i] = new ScoredCoordinate(c[i].getChromosome(), c[i].getStart(), c[i].getStop(), i);
		}
		return sc;
	}
	
	public String toString(){
		return chromosome+"\t"+start+"\t"+stop+"\t"+score;
	}
	
	/**Split by chromosome into a HashMap of chromosome:ScoredCoordinate[].
	 * Don't forget to sort!*/
	public static HashMap splitByChromosome(ScoredCoordinate[] sortedCoordinates){
		HashMap chrSpec = new HashMap();
		ArrayList al = new ArrayList();
		String currChrom = sortedCoordinates[0].getChromosome();
		for (int i=0; i< sortedCoordinates.length; i++){
			if (sortedCoordinates[i].getChromosome().equals(currChrom) == false){
				ScoredCoordinate[] sub = new ScoredCoordinate[al.size()];
				al.toArray(sub);
				chrSpec.put(currChrom, sub);
				al.clear();
				currChrom = sortedCoordinates[i].getChromosome();
			}
			al.add(sortedCoordinates[i]);
		}
		//add last to hash
		ScoredCoordinate[] sub = new ScoredCoordinate[al.size()];
		al.toArray(sub);
		chrSpec.put(currChrom, sub);
		return chrSpec;
	}
	
	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}
}
