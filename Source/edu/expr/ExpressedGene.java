package edu.expr;
import java.util.*;
import util.bio.annotation.*;
import util.gen.*;

public class ExpressedGene implements Comparable{
	//fields
	private String name;
	private Coordinate coordinates;
	private float[] values;
	private RankedFloatArray rfa;
	
	public ExpressedGene (String name, Coordinate coordinates, float[] values){
		this.name = name;
		this.coordinates = coordinates;
		this.values = values;
	}

	public String toString(){
		return name+"\t"+coordinates;
	}
	/**Sorts by Coordinate */
	public int compareTo(Object other){
		ExpressedGene exOther = (ExpressedGene) other; 
		return coordinates.compareTo(exOther.coordinates);
	}
	
	/**Splits a ExpressedGene[] by chromosome to ExpressedGene[chrom][lines].
	 * Assumes the ExpressedGene[] is sorted by chromosome!*/
	public static ExpressedGene[][] splitByChromosome(ExpressedGene[] exprGenes){
		//scan for number of chromosomes
		HashSet chroms = new HashSet();
		for (int x=0; x< exprGenes.length; x++){
			Coordinate coor = exprGenes[x].getCoordinates();
			if (chroms.contains(coor.getChromosome()) == false) chroms.add(coor.getChromosome());
		}
		//make array
		ExpressedGene[][] exprGenChrom = new ExpressedGene[chroms.size()][];
		//assume sorted so make sub chrom arrays and add
		String currChrom = exprGenes[0].getCoordinates().getChromosome();
		ArrayList sub = new ArrayList(5000);
		int counter =0;
		for (int x=0; x<exprGenes.length; x++){
			//new chrom?
			String chrom = exprGenes[x].getCoordinates().getChromosome();
			if (chrom.equals(currChrom)) sub.add(exprGenes[x]);
			else{
				//make array
				ExpressedGene[] s = new ExpressedGene[sub.size()];
				sub.toArray(s);
				sub.clear();
				exprGenChrom[counter++] = s;
				currChrom = chrom;
				sub.add(exprGenes[x]);
			}
		}
		//make last
		ExpressedGene[] s = new ExpressedGene[sub.size()];
		sub.toArray(s);
		exprGenChrom[counter] = s;
		return exprGenChrom;
	}

	
	public Coordinate getCoordinates() {
		return coordinates;
	}

	public String getName() {
		return name;
	}

	public float[] getValues() {
		return values;
	}

	public RankedFloatArray getRfa() {
		return rfa;
	}

	public void setRfa(RankedFloatArray rfa) {
		this.rfa = rfa;
	}

	public void setValues(float[] values) {
		this.values = values;
	}
}
