package edu.utah.seq.useq.data;
import java.util.*;
import java.io.*;
import util.gen.*;

/**
 * Simple start stop object. Assumes interbase coordinates. 
 * @author david.nix@hci.utah.edu*/
public class Region implements Comparable<Region>, Serializable {
	//fields
	private static final long serialVersionUID = 1L;
	protected int start;
	protected int stop;

	//constructors
	public Region (int start, int stop){
		this.start = start;
		this.stop = stop;
	}

	public String toString(){
		return start+"\t"+stop;
	}
	/**Assumes coordinates are interbase.*/
	public boolean intersects (int start, int stop){
		if (stop <= this.start || start >= this.stop) return false;
		return true;
	}

	/**Checks to see if each start is <= the stop*/
	public static boolean checkStartStops(Region[] ss){
		for (int i=0; i< ss.length; i++){
			if (ss[i].start> ss[i].stop) return false;
		}
		return true;
	}
	/**Returns a Region[] for regions defined with baseHitCounts !=0.
	 * Coordinates are interbase.*/
	public static Region[] makeStartStops(short[] baseHitCount){
		ArrayList<Region> ss = new ArrayList<Region>();
		//find first non zero base
		int i=0;
		for (; i< baseHitCount.length; i++) if (baseHitCount[i]!=0) break;
		if (i == baseHitCount.length) return null;
		int start = i;
		int val = baseHitCount[start];
		//find different base
		for (; i< baseHitCount.length; i++){
			if (baseHitCount[i]!=val) {
				//make SS
				if (val!=0) ss.add(new Region(start,i));
				start = i;
				val = baseHitCount[start];
			}
		}
		//add last
		if (val!=0) ss.add(new Region(start,i));
		Region[] ssA = new Region[ss.size()];
		ss.toArray(ssA);
		return ssA;
	}

	/**Sorts by start base, then by length, smaller to larger for both.*/
	public int compareTo(Region se){
		if (start<se.start) return -1;
		if (start>se.start) return 1;
		// if same start, sort by length, smaller to larger
		int len = stop-start;
		int otherLen = se.stop-se.start;
		if (len<otherLen) return -1;
		if (len>otherLen) return 1;
		return 0;
	}

	/**Returns the total number of bases, assumes interbase coordinates.*/
	public static int totalBP (Region[] ss){
		int total = 0;
		for (int i=0; i< ss.length; i++) total += ss[i].getLength();
		return total;
	}

	/**Returns the starts in number of bases, not coordinates for the array.*/
	public static int[] startsInBases (Region[] ss){
		int[] indexes = new int[ss.length];
		int index = 0;
		for (int i=0; i< ss.length; i++) {
			index += ss[i].getLength();
			indexes[i] = index;
		}
		return indexes;
	}

	/**Assumes interbase coordinates.*/
	public int getLength(){
		return stop-start;
	}

	public int getStop() {
		return stop;
	}
	public int getStart() {
		return start;
	}
	public int[] getStartStop(){
		return new int[]{start,stop};
	}

	public boolean isContainedBy(int beginningBP, int endingBP) {
		if (start >= beginningBP && stop < endingBP) return true;
		return false;
	}
	/**Parses a tab delimited file (chr, start, stop, ...), zip/ gz OK. 
	 * @param bed file, skips empty lines and those starting with '#'
	 * @param subStart and subEnd are the number to subtract from the ends of each region
	 * @return a HashMap<Chr,sorted Region[]> or null in none are found
	 * */
	public static HashMap<String,Region[]> parseStartStops (File bedFile, int subStart, int subEnd, int minSize){
		HashMap<String,ArrayList<Region>> ss = new HashMap<String,ArrayList<Region>>();
		try{
			BufferedReader in = IO.fetchBufferedReader(bedFile);
			String line;
			String[] tokens;
			ArrayList<Region> al = new ArrayList<Region>();
			//chrom, start, stop
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				tokens = line.split("\\s+");
				if (tokens.length < 3) continue;
				//does chrom already exist?
				if (ss.containsKey(tokens[0])) al = ss.get(tokens[0]);
				else {
					al = new ArrayList<Region>();
					ss.put(tokens[0], al);
				}
				int start = Integer.parseInt(tokens[1]);
				int stop = Integer.parseInt(tokens[2]);
				if (start > stop) throw new Exception("\nFound a start that is greater than stop!  Cannot parse file "+bedFile+", bad line->\n\t"+line);
				if (start < 0) throw new Exception("\nFound a start with a negative value!  Cannot parse file "+bedFile+", bad line->\n\t"+line);
				int length = stop-start;
				if (length < minSize) continue;
				al.add(new Region(start-subStart, stop- subEnd));
			}
		}catch (Exception e){
			e.printStackTrace();
			return null;
		}
		if (ss.size() == 0) return null;
		//make hashmap
		HashMap<String,Region[]> ssReal = new HashMap<String,Region[]>();
		Iterator<String> it = ss.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			ArrayList<Region> al = ss.get(chrom);
			Region[] array = new Region[al.size()];
			al.toArray(array);
			Arrays.sort(array);
			ssReal.put(chrom, array);
		}
		return ssReal;
	}
	
	
	public int getMiddle(){
		if (start == stop) return start;
		double length = stop - start;
		double halfLength = length/2.0;
		return (int)Math.round(halfLength) + start;
	}
	
	public static int findLastBase(Region[] r){
		int lastBase = -1;
		for (int i=0; i< r.length; i++){
			if (r[i].stop> lastBase) lastBase = r[i].stop;
		}
		return lastBase;
	}
}
