package trans.roc;

import java.io.*;
import trans.misc.*;
import java.util.*;
import util.gen.*;

/**
 * Sgr line container.
 */
public class Sgr implements Serializable{
	//fields
	private String chromosome = null;
	private int position = 0;
	private double score = 0;

	//constructors
	public Sgr (String line){
		String[] tokens = line.split("\\s");
		try {
			if (tokens.length > 2){
				chromosome = tokens[0];
				position = Integer.parseInt(tokens[1]);
				score = Double.parseDouble(tokens[2]);
			}
		} catch (Exception e){
			System.err.println("Problem parsing sgr line for "+line+"\n");
			chromosome = null;
			position = -1;
			score = -1;
		}
	}
	public Sgr (String line, String chromosome){
		String[] tokens = line.split("\\s");
		this.chromosome = chromosome;
		position = Integer.parseInt(tokens[0]);
		score = Double.parseDouble(tokens[1]);
	}

	//main methods
	/**Splits a Sgr[] containing many chromosome lines to a Sgr[chrom][lines].
	 * Assumes the Sgr[] is sorted by chromosome!*/
	public static Sgr[][] splitByChromosome(Sgr[] sortedSgrs){
		//scan for number of chromosomes
		HashSet chroms = new HashSet();
		for (int x=0; x< sortedSgrs.length; x++){
			if (chroms.contains(sortedSgrs[x].getChromosome()) == false) chroms.add(sortedSgrs[x].getChromosome());
		}
		//make array
		Sgr[][] sgrChrom = new Sgr[chroms.size()][];
		//assume sorted so make sub chrom arrays and add
		String currChrom = sortedSgrs[0].getChromosome();
		ArrayList sub = new ArrayList(5000);
		int counter =0;
		for (int x=0; x<sortedSgrs.length; x++){
			//new chrom?
			if (sortedSgrs[x].getChromosome().equals(currChrom)) sub.add(sortedSgrs[x]);
			else{
				//make array
				Sgr[] s = new Sgr[sub.size()];
				sub.toArray(s);
				sub.clear();
				sgrChrom[counter++] = s;
				currChrom = sortedSgrs[x].getChromosome();
				sub.add(sortedSgrs[x]);
			}
		}
		//make last
		Sgr[] s = new Sgr[sub.size()];
		sub.toArray(s);
		sgrChrom[counter] = s;
		return sgrChrom;
	}

	/**Splits a Sgr[] by chromosome. Assumes the Sgr[] is sorted!*/
	public static GrGraph[] splitByChrom(Sgr[] sortedSgrs){
		ArrayList<GrGraph> graphs = new ArrayList<GrGraph>();
		//assume sorted so make sub chrom arrays and add
		String currChrom = sortedSgrs[0].getChromosome();
		ArrayList<Integer> position = new ArrayList<Integer>(10000);
		ArrayList<Float> score = new ArrayList<Float>(10000);
		for (int x=0; x<sortedSgrs.length; x++){
			//new chrom?
			if (sortedSgrs[x].getChromosome().equals(currChrom) == false) {
				//make array
				graphs.add( new GrGraph(currChrom, null, Num.arrayListOfIntegerToInts(position), Num.arrayListOfFloatToArray(score)));
				position.clear();
				score.clear();
				currChrom = sortedSgrs[x].getChromosome();
			}
			position.add(new Integer(sortedSgrs[x].position));
			score.add(new Float (sortedSgrs[x].score));
		}
		//make last
		graphs.add( new GrGraph(currChrom, null, Num.arrayListOfIntegerToInts(position), Num.arrayListOfFloatToArray(score)));
		//convert
		GrGraph[] grs = new GrGraph[graphs.size()];
		graphs.toArray(grs);
		return grs;
	}


	/**Loads an xxx.sgr(.zip) file into an ArrayList of Sgr*/
	public static ArrayList<Sgr> loadSgrFile(File sgrFile){
		ArrayList<Sgr> sgrAL = new ArrayList<Sgr>(100000);
		//load file
		try{
			String line;
			BufferedReader in = IO.fetchBufferedReader(sgrFile);
			while ((line=in.readLine())!=null){
				line = line.trim();
				if (line.length()==0 || line.startsWith("#")) continue;
				Sgr test = new Sgr(line);
				//skip score?
				//if (test.getScore() == -1) continue;
				//transform?
				//test.setScore(Num.log2(test.getScore()/2));
				if (test.getChromosome() == null) {
					System.out.println("\nError: this is not an sgr line! -> "+line+"\n\tIn file "+sgrFile);
					return null;
				}
				sgrAL.add(test);
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}
		return sgrAL;
	}



	//getters setters
	public String getChromosome() {
		return chromosome;
	}
	public int getPosition() {
		return position;
	}
	public double getScore() {
		return score;
	}
	public String toString(){
		return chromosome+"\t"+position+"\t"+score;
	}
	public void setScore(double score) {
		this.score = score;
	}




}
