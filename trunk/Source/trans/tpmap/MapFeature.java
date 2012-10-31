package trans.tpmap;
import java.io.*;

/**
 * Tpmap line broken into some components, chromosome, start, pmX, pmY, and optionally, mmX, mmY. 
 * No getters or setters to minimize size and increase spead recall.
 * 
 */
public class MapFeature implements Comparable, Serializable{
	//fields
	public static final long serialVersionUID = 3;
	public String chromosome;
	public int matches = -1;
	public int start;
	public int pmX = -1;
	public int pmY = -1;
	public int mmX = -1;
	public int mmY = -1;
	
	/**Constructor - bpmapLine - seq ori chrom start pmX pmY mmX mmY, (mm coords are optional).
	 */
	public MapFeature(String tpmapLine){
		try{
			String[] tokens = tpmapLine.split("\\s+");
			//attempt to part number of hits to genome
			if (tokens[1].equals("t") == false && tokens[1].equals("f") == false) matches = Integer.parseInt(tokens[1]);
			chromosome = tokens[2];
			start = Integer.parseInt(tokens[3]);
			pmX = Integer.parseInt(tokens[4]);
			pmY = Integer.parseInt(tokens[5]);
			//mismatches?
			if (tokens.length>=8){
				mmX = Integer.parseInt(tokens[6]);
				mmY = Integer.parseInt(tokens[7]);
			}
		} catch (Exception e){
			System.out.println("\nProblem with parsing tpmap feature line -> "+tpmapLine+"\n");
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	//main methods
	/**Sorts by chromosome and then position.*/
	public int compareTo(Object obj){
		MapFeature other = (MapFeature)obj;
		//sort by chromosome
		int comp = other.chromosome.compareTo(chromosome);
		if (comp != 0) return comp*-1;
		//sort by base position
		if (other.start>start) return -1;
		if (other.start<start) return 1;
		return 0;
	}
}
