package util.bio.annotation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class NamedCoordinate extends Coordinate{
	//fields
	private String name;  //or rank
	
	//constructor
	public NamedCoordinate (String chromosome, int start, int stop, String name){
		super(chromosome, start, stop);
		this.name = name;
	}

	//methods
	public String toString(){
		return chromosome+"\t"+start+"\t"+stop+"\t"+name;
	}
	
	/**Split by chromosome into a HashMap of chromosome:NamedCoordinate[].
	 * Don't forget to sort!*/
	public static HashMap<String,NamedCoordinate[]> splitByChromosome(NamedCoordinate[] sortedCoordinates){
		HashMap<String,NamedCoordinate[]> chrSpec = new HashMap<String,NamedCoordinate[]>();
		ArrayList<NamedCoordinate> al = new ArrayList<NamedCoordinate>();
		String currChrom = sortedCoordinates[0].getChromosome();
		for (int i=0; i< sortedCoordinates.length; i++){
			if (sortedCoordinates[i].getChromosome().equals(currChrom) == false){
				NamedCoordinate[] sub = new NamedCoordinate[al.size()];
				al.toArray(sub);
				chrSpec.put(currChrom, sub);
				al.clear();
				currChrom = sortedCoordinates[i].getChromosome();
			}
			al.add(sortedCoordinates[i]);
		}
		//add last to hash
		NamedCoordinate[] sub = new NamedCoordinate[al.size()];
		al.toArray(sub);
		chrSpec.put(currChrom, sub);
		return chrSpec;
	}
	
	/**Parses a tab delimited chr,start,stop,...,text file.
	 * @param subStart - bases to be subtracted from region starts
	 * @param subEnd - bases to be subtracted from region ends
	 * */
	public static NamedCoordinate[] parseFile(File picksFile, int subStart, int subEnd, int indexOfNameColumn){
		NamedCoordinate[] coor =null;
		try{
			BufferedReader in = new BufferedReader(new FileReader(picksFile));
			String line;
			String[] tokens;
			ArrayList<NamedCoordinate> al = new ArrayList<NamedCoordinate>();
			//chrom, start, stop,.. text
			//0 1 2
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				tokens = line.split("\\s+");
				if (tokens.length < 4) continue;
				al.add(new NamedCoordinate(tokens[0], Integer.parseInt(tokens[1])-subStart, Integer.parseInt(tokens[2])- subEnd, tokens[indexOfNameColumn]));
			}
			coor = new NamedCoordinate[al.size()];
			al.toArray(coor);
		}catch (IOException e){
			e.printStackTrace();
		}
		return coor;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}
	

}
