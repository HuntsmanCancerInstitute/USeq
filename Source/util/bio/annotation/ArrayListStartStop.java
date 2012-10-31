package util.bio.annotation;
import java.io.*;
import java.util.*;

import util.gen.IO;

import edu.utah.seq.useq.data.Region;



/**
 * Simple start stop object for sorting. 
 */
public class ArrayListStartStop extends Region {
	private ArrayList arrayList;
	
	//constructors
	public ArrayListStartStop (int start, int end, ArrayList arrayList){
		super (start, end);
		this.arrayList = arrayList;
	}
	
	/**Assumes coordinates are interbase.*/
	public boolean intersects (int begin, int stop){
		if (stop <= this.start || begin >= this.stop) return false;
		return true;
	}
	
	/**Parses a tab delimited file (chr, start, stop, ...).
	 * @param bed file, skips empty lines and those starting with '#'
	 * @param subStart and subEnd are the number to subtract from the ends of each region
	 * @return a HashMap<Chr,sorted Region[]> or null in none are found
	 * */
	public static HashMap<String,ArrayListStartStop[]> parseArrayListStartStops (File bedFile, int subStart, int subEnd, boolean addStrandToChromName){
		HashMap<String,ArrayList<ArrayListStartStop>> ss = new HashMap<String,ArrayList<ArrayListStartStop>>();
		try{
			BufferedReader in = IO.fetchBufferedReader(bedFile); 
			String line;
			String[] tokens;
			ArrayList<ArrayListStartStop> al = new ArrayList<ArrayListStartStop>();
			//chrom, start, stop
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				tokens = line.split("\\s+");
				if (tokens.length < 3) continue;
				//does chrom already exist?
				String chrom;
				if (addStrandToChromName){
					if (tokens.length == 6) chrom = tokens[0]+tokens[5];
					else chrom = tokens[0]+".";
				}
				else chrom = tokens[0];
				if (ss.containsKey(chrom)) al = ss.get(chrom);
				else {
					al = new ArrayList<ArrayListStartStop>();
					ss.put(chrom, al);
				}
				int start = Integer.parseInt(tokens[1])-subStart;
				int stop = Integer.parseInt(tokens[2])- subEnd;
				if (start > stop) throw new Exception ("\nStart is less than stop, aborting bed file parsing. "+line);
				al.add(new ArrayListStartStop(start, stop, new ArrayList()));
			}
		}catch (Exception e){
			e.printStackTrace();
		}
		if (ss.size() == 0) return null;
		//make hashmap
		HashMap<String,ArrayListStartStop[]> ssReal = new HashMap<String,ArrayListStartStop[]>();
		Iterator<String> it = ss.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			ArrayList<ArrayListStartStop> al = ss.get(chrom);
			ArrayListStartStop[] array = new ArrayListStartStop[al.size()];
			al.toArray(array);
			Arrays.sort(array);
			ssReal.put(chrom, array);
		}
		return ssReal;
	}

	public ArrayList getArrayList() {
		return arrayList;
	}

	public void setArrayList(ArrayList arrayList) {
		this.arrayList = arrayList;
	}

}
