package util.bio.annotation;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

import trans.anno.BindingRegion;
import util.gen.IO;
import util.gen.Misc;

/**Base class for holding info about a genomic coordinate.*/
public class Coordinate implements Comparable<Coordinate>, Serializable{
	//fields
	String chromosome;
	int start;
	int stop;
	private static final String className = "Coordinate";
	
	//constructor
	public Coordinate (String chromosome, int start, int stop){
		this.chromosome = chromosome;
		this.start = start;
		this.stop = stop;
	}
	/**Sorts by chromsome, start position, length (smallest to largest).*/
	public int compareTo(Coordinate otherCoor) {
		//sort by chromosome
		int compare = otherCoor.chromosome.compareTo(chromosome);
		if (compare !=0) return compare * -1;
		//sort by start position
		if (start<otherCoor.start) return -1;
		if (start>otherCoor.start) return 1;
		// if same start, sort by length, smaller to larger
		int len = stop-start;
		int otherLen = otherCoor.stop-otherCoor.start;
		if (len<otherLen) return -1;
		if (len>otherLen) return 1;
		return 0;
	}
	public String toString(){
		return new String(chromosome +"\t"+ start+ "\t"+ stop);
	}
	
	/**Returns chr+":"+(start+1)+"-"+stop; */
	public String getTabixSearchCoordinates() {
		return new String(chromosome +":"+ (start+1)+ "-"+ stop);
	}
	
	/**Assumes coordinates are inclusive.*/
	public boolean intersects (int beginning, int end){
		if (end < start || beginning > stop) return false;
		return true;
	}
	
	public static boolean writeToFile(Coordinate[] c, File file){
		try{
			PrintWriter out = new PrintWriter( new FileWriter(file));
			for (int i=0; i< c.length; i++) out.println(c[i]);
			out.close();		
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	/**Assumes coordinates are inclusive.*/
	public boolean intersects (Coordinate other){
		if (other.chromosome.equals(this.chromosome) == false) return false;
		if (other.stop < this.start || other.start > this.stop) return false;
		return true;
	}
	
	
	
	/**Returns -1 for no overlap, 0 for complete overlap, or a positive int for the # bases of overlap.
	 * Assumes coordinates are inclusive.*/
	public int bpsIntersection(Coordinate other){
		if (other.getChromosome().equals(chromosome)== false) return -1;
		// is other left of this
		if (other.getStop() < start) return -1;
		// is other right of this
		if (other.getStart() > stop ) return -1;
		// must overlap
		//left side
		if (start>=other.getStart() && start<=other.getStop() && stop > other.getStop()) return other.getStop() - start +1;
		//right side
		if (stop >= other.getStart() && stop <= other.getStop() && start < other.getStart()) return stop-other.getStart() +1;
		//contained within
		return 0;
	}
	
	/**Clips the forReduction based on the mask.  Will return null if the mask completely covers the forReduction.
	 * Assumes interbase coordinates.*/
	public static Coordinate[] subtract (Coordinate forReduction, Coordinate mask){
		//do they intersect?
		if (mask.getStop() <= forReduction.getStart() || mask.getStart() >= forReduction.getStop()) {
			//System.out.println("no int");
			return new Coordinate[] {forReduction};
		}
		//contained?
		if (mask.getStart() <= forReduction.getStart() &&  mask.getStop()>= forReduction.getStop()) {
			//System.out.println("contained");
			return null;
		}
		//inside center punch
		if (mask.getStart() > forReduction.getStart() && mask.getStop() < forReduction.getStop()) {
			Coordinate c1 = new Coordinate (forReduction.getChromosome(), forReduction.getStart(), mask.getStart());
			Coordinate c2 = new Coordinate (forReduction.getChromosome(), mask.getStop(), forReduction.getStop());
			//System.out.println(c1+" center "+c2);
			return new Coordinate[] {c1, c2};
		}
		//left
		if (mask.getStop() > forReduction.getStart() && mask.getStop() < forReduction.getStop()) {
			forReduction.setStart(mask.getStop());
			//System.out.println("left");
		}
		//right
		else if (mask.getStart() >= forReduction.getStart() && mask.getStart() <= forReduction.getStop()) {
			forReduction.setStop(mask.getStart());
			//System.out.println("right");
		}
		return new Coordinate[] {forReduction};
	}
	
	/**Finds the last base.*/
	public static int findLastBase(Coordinate[] chromSpecific){
		int lastBase = -1;
		for (int i=0; i< chromSpecific.length; i++){
			if (chromSpecific[i].getStop() > lastBase) lastBase = chromSpecific[i].getStop();
		}
		return lastBase;
	}
	
	/**Parses a tab delimited chr,start,stop file (zip or gz OK)
	 * @param subStart - bases to be subtracted from region starts
	 * @param subEnd - bases to be subtracted from region ends
	 * */
	public static Coordinate[] parseFile(File picksFile, int subStart, int subEnd){
		Coordinate[] coor =null;
		try{
			BufferedReader in = IO.fetchBufferedReader(picksFile);
			String line;
			String[] tokens;
			ArrayList<Coordinate> al = new ArrayList<Coordinate>();
			//chrom, start, stop
			//0 1 2
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() ==0 || line.startsWith("#")) continue;
				tokens = line.split("\\s+");
				if (tokens.length < 3) continue;
				
				al.add(new Coordinate(tokens[0], Integer.parseInt(tokens[1])-subStart, Integer.parseInt(tokens[2])- subEnd));
			}
			coor = new Coordinate[al.size()];
			al.toArray(coor);
		}catch (IOException e){
			e.printStackTrace();
		}
		return coor;
	}
	
	/**Split by chromosome into a HashMap of chromosome:Coordinate[].
	 * Don't forget to sort!*/
	public static HashMap<String,Coordinate[]> splitByChromosome(Coordinate[] sortedCoordinates){
		HashMap<String,Coordinate[]> chrSpec = new HashMap();
		ArrayList<Coordinate> al = new ArrayList<Coordinate>();
		String currChrom = sortedCoordinates[0].getChromosome();
		for (int i=0; i< sortedCoordinates.length; i++){
			if (sortedCoordinates[i].getChromosome().equals(currChrom) == false){
				Coordinate[] sub = new Coordinate[al.size()];
				al.toArray(sub);
				chrSpec.put(currChrom, sub);
				al.clear();
				currChrom = sortedCoordinates[i].getChromosome();
			}
			al.add(sortedCoordinates[i]);
		}
		//add last to hash
		Coordinate[] sub = new Coordinate[al.size()];
		al.toArray(sub);
		chrSpec.put(currChrom, sub);
		return chrSpec;
	}
	
	/**Writes a binary Coordinate[].
	 * @return true if sucessful, false if something bad happened.*/
	public static boolean writeBinary(Coordinate[] c, File file){
		try {
			int numCorr = c.length;
			FileOutputStream fos = new FileOutputStream(file);
			DataOutputStream dos = new DataOutputStream( new BufferedOutputStream (fos));
			
			//write text of array
			dos.writeUTF(className);
			
			//write number of Coordinate
			dos.writeInt(numCorr);
			
			//write Coordinate
			for (int i=0; i<numCorr; i++) { 
				dos.writeUTF(c[i].chromosome);
				dos.writeInt(c[i].start);
				dos.writeInt(c[i].stop); 
			}
			dos.close();
			fos.close();
			return true;
			
		} catch (IOException ioe) {
			ioe.printStackTrace();
			return false; 
		}
	}
	
	
	/**Reads a binary Coordinate[] file.
	 * @return null if something bad happened.*/
	public static Coordinate[] readBinary(File file){
		try {
			FileInputStream fis = new FileInputStream(file);
			DataInputStream dis = new DataInputStream( new BufferedInputStream(fis ));
			//read array text
			if (dis.readUTF().equals(className) == false) return null;
			//read number of Coordinate
			int numCoor = dis.readInt();
			//make array
			Coordinate[] c = new Coordinate[numCoor];
			for (int i=0; i< numCoor; i++){
				String chrom = dis.readUTF();
				int start = dis.readInt();
				int stop = dis.readInt();
				c[i] = new Coordinate(chrom, start, stop);
			}
			
			dis.close();
			fis.close();
			return c;
		}
		catch (Exception ioe){
			ioe.printStackTrace();
			return null;   
		} 
	}
	
	public static Coordinate[] insertGaps(Coordinate[] sortedNonOverlappingRegions, int bpPad) {
		ArrayList<Coordinate> toReturn = new ArrayList<Coordinate>();
		Coordinate left = sortedNonOverlappingRegions[0];
		
		for (int i=1; i< sortedNonOverlappingRegions.length; i++) {
			Coordinate right = sortedNonOverlappingRegions[i];
			
			Coordinate[] tLR = insertGapInPair(left, right, bpPad);
			
			//is left not null? save it
			if (tLR[0] != null) toReturn.add(tLR[0]);
			
			//is right not null? make it the new left
			if (tLR[1] != null) left = tLR[1];
	
			else {
				//right was deleted so pull next right
				i++;
				if (i< sortedNonOverlappingRegions.length) left = sortedNonOverlappingRegions[i];
				else left = null;
			}
		}
		
		//add last?
		if (left!=null) toReturn.add(left);
		
		Coordinate[] cor = new Coordinate[toReturn.size()];
		toReturn.toArray(cor);
		return cor;
	}
	
	public static void main(String[] args) {
		
		Coordinate a = new Coordinate("chrX", 100, 200);
		Coordinate b = new Coordinate("chrX", 400, 500);
		Coordinate c = new Coordinate("chrX", 600, 700);
		Coordinate d = new Coordinate("chrX", 900, 1000);
		
		Coordinate[] cor = new Coordinate[] {a,b,c,d};
		Coordinate[] trm = insertGaps(cor, 150);
		
		for (Coordinate cs: trm) IO.pl(cs);
		
	}
	
	/**Inserts the bpPad gap, may lead to complete deletion of the left or right region. 
	 * Assumes they don't overlap and are the same chr.*/
	private static Coordinate[] insertGapInPair(Coordinate left, Coordinate right, int bpPad) {

		int halfBpPad = (int)Math.round((double)bpPad/ 2.0);
		Coordinate trimmedLeft = null;
		Coordinate trimmedRight = null;
		
		//calc gap
		int leftStop = left.getStop();
		int rightStart = right.getStart();
		int diff = rightStart-leftStop;
		//far apart? if so the just return unaltered
		if (diff > bpPad) return new Coordinate[] {left, right};
		
		int midPos = (int)Math.round((double)diff/2.0)+ leftStop;
		
		//trim left
		int testLeftStop = midPos-halfBpPad;
		//OK?
		if (left.getStart() < testLeftStop) trimmedLeft = new Coordinate(left.getChromosome(), left.getStart(), testLeftStop);
		else return new Coordinate[] {null, right};
		
		//trim right
		int testRightStart = midPos+halfBpPad;
		if (testRightStart< right.getStop()) trimmedRight = new Coordinate(right.getChromosome(), testRightStart, right.getStop()); 
		else return new Coordinate[] {left, null};
		
		return new Coordinate[] {trimmedLeft, trimmedRight};
		
		
	}
	
	public String getChromosome() {
		return chromosome;
	}
	
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	
	public int getStart() {
		return start;
	}
	
	public void setStart(int start) {
		this.start = start;
	}
	
	public int getStop() {
		return stop;
	}
	
	public void setStop(int stop) {
		this.stop = stop;
	}
	
	/**rounds down*/
	public int getMiddle(){
		float d = stop - start;
		d = d/2.0f;
		return start + (int)d;
	}
	
	/**Assumes interbase coordinates.*/
	public int getLength(){
		return stop - start;
	}

	
}

