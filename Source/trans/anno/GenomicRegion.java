package trans.anno;
import trans.tpmap.*;
import java.io.*;
import java.util.*;

import trans.main.*;
import util.bio.annotation.Coordinate;
import util.gen.IO;
import util.gen.Misc;
/**
 * Represents a genomic region. See Coordinate for a more simple class.
 *
 */
public class GenomicRegion implements Serializable{

	//fields
	private String chromosome;
	private int start;
	private int end;
	private String notes;
	private HashSet hits;
	private boolean[] gcContent;	//one per base, true = g or c, false = a or t
	
	
	//constructors
	public GenomicRegion(String chromosome, int start, int end, String notes){
		this.start = start;
		this.end = end;
		this.chromosome = chromosome;
		this.notes = notes;
	}
	
	//methods
	
	/**Returns stop - start +1, thus last base is included, not interbase numbering!*/
	public int getLength(){
		return end - start +1;
	}
	
	/**Looks for a binary version of the file xxx.corr, if found loads binary, otherwise it parses the txt file 
	 * and then writes the binary for future use.*/
	public static GenomicRegion[] loadWriteBinaryRegions(File regionsFile){
		if (regionsFile.getName().endsWith(".corr")){
			return loadBinaryCoordinates(regionsFile);
		}
		File binary = new File (IO.getFullPathName(regionsFile)+".corr");
		if (binary.exists()) return loadBinaryCoordinates(binary);
		else {
			GenomicRegion[] rs = GenomicRegion.parseRegions(regionsFile);
			System.out.println("\tSaving binary version of regions file -> "+binary);
			GenomicRegion.writeBinaryCoordinates(binary, rs);
			return rs;
		}
	}
	
	/**Reads binary chrom start stop file into GenomicRegion[]*/
	public static GenomicRegion[] loadBinaryCoordinates(File file){
		Coordinate[] c = Coordinate.readBinary(file);
		GenomicRegion[] rs = new GenomicRegion[c.length];
		for (int i=0; i< c.length; i++){
			rs[i] = new GenomicRegion(c[i].getChromosome(), c[i].getStart(), c[i].getStop(), null);
		}
		c = null;
		return rs;
	}
	
	/**Writes the chrom start stop of a GenomicRegion[] as a binary file.*/
	public static void writeBinaryCoordinates(File file, GenomicRegion[] r){
		Coordinate[]  startStop = new Coordinate[r.length];
		for (int i=0; i< r.length; i++){
			startStop[i] = new Coordinate(r[i].getChromosome(), r[i].getStart(), r[i].getEnd());
		}
		Coordinate.writeBinary(startStop, file);
		startStop = null;
	}
	
	/**Given a chromosome SORTED array of Regions returns a HashMap 
	 * containing chromosome: sorted chromo specific GenomicRegion[]*/
	public static HashMap splitByChromosome(GenomicRegion[] sorted){
		ArrayList al = new ArrayList();
		HashMap map = new HashMap();
		String chromosome = sorted[0].getChromosome();
		//for each GenomicRegion
		for (int i=0; i< sorted.length; i++){
			//same chromosome?
			if (sorted[i].getChromosome().equals(chromosome)){
				al.add(sorted[i]);
			}
			//different chromosome!
			else {
				//set new chrom in map
				GenomicRegion[] regions = new GenomicRegion[al.size()];
				al.toArray(regions);
				map.put(chromosome, regions);
				//reset params
				al = new ArrayList();
				al.add(sorted[i]);
				chromosome = sorted[i].getChromosome();
			}
		}
		//add last
		GenomicRegion[] regions = new GenomicRegion[al.size()];
		al.toArray(regions);
		map.put(chromosome, regions);
		return map;
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append(chromosome);
		sb.append("\t");
		sb.append(start);
		sb.append("\t");
		sb.append(end);
		sb.append("\t");
		if (notes!=null)sb.append(notes);
		return sb.toString();
	}
	
	/**Returns the fraction of GC provided the GenomicRegion.gcContent[] is loaded.
	 * Otherwise returns zero.*/
	public double calculateGCContent(){
		double total = 0;
		for (int i=0; i<gcContent.length; i++){
			if (gcContent[i]) total++;
		}
		return total/(gcContent.length);
	}
	
	/**Writes bed file.*/
	public static void writeRegions(HashMap<String, GenomicRegion[]> regions, File bedFile){
		try{
			PrintWriter out = new PrintWriter (new FileWriter( bedFile));
			//for each chromosome
			Iterator<String> it = regions.keySet().iterator();
			while (it.hasNext()){
				//fetch CNVs
				String chrom = it.next();
				GenomicRegion[] r = regions.get(chrom);
				for (int i=0; i< r.length; i++) out.println(r[i]);
			}
			out.close();
		}catch (Exception e){
			e.printStackTrace();
		}
	}
	
	//methods
	public boolean intersect(Interval interval){
		if (interval.getChromosome().equals(chromosome)== false) return false;
		// is other left of this
		if ((interval.getStartLastOligo()+interval.getSizeOfOligoMinusOne()) < start) return false;
		// is other right of this
		if (interval.getStart1stOligo() > end ) return false;
		// must overlap
		return true;
	}
	public boolean intersect(TPMapLine line){
		if (line.getChromosome().equals(chromosome)== false) return false;
		// is other left of this
		if (line.getStart() < start) return false;
		// is other right of this
		if (line.getStart() > end ) return false;
		// must overlap
		return true;
	}
	/**Returns true only if region contains or is contained or is covered by >= fractionAcceptibleCoverage by a window.
	 * ie 0.5 for 50% or region covered by a window.*/
	public boolean intersect(Window window, int sizeOligo, double fractionAcceptibleCoverage){
		if (window.getChromosome().equals(chromosome)== false) return false;
		// is other left of this
		int winEnd = window.getStartLastOligo()+sizeOligo;
		if (winEnd < start) return false;
		// is other right of this
		if (window.getStart1stOligo() > end ) return false;
		// must overlap
			//region contained within window
			if (window.getStart1stOligo()<= start && winEnd >= end) {
				return true;
			}
			//window contained by region
			if (window.getStart1stOligo()>= start && winEnd <= end) {
				return true;
			}
			// window covers >= 1/2 of region
				double coveredBases = 0;
				//window left of region
				if (window.getStart1stOligo()< start) coveredBases = winEnd - start;
				//window right of region
				else coveredBases = end - window.getStart1stOligo();
				if (coveredBases/(double)(end-start) >= fractionAcceptibleCoverage) return true;
			
		return false;
	}
	
	/**Assumes regions are on the same chromosome, and the stop base is included, not interbase numbering.
	 * Returns 0 if the regions abut, positive ints for each base in the gap (ie 1= 1bp between an stop of 12 and a start of 14),
	 * negative ints for each base of intersection (ie -1= 1bp overlap),
	 * max negative int is the length of the smaller region (ie stop - start +1).*/
	public int bpIntersectionSameChromosome(GenomicRegion other){
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
	
	/**Returns -1 for no overlap, or a positive int for the # bases of overlap.
	 * Doesn't check chromosome!*/
	public int overlap(Interval interval){
		// is other left of this
		int endInterval = interval.getStartLastOligo()+interval.getSizeOfOligoMinusOne();
		if (endInterval < start) return -1;
		// is other right of this
		int startInterval = interval.getStart1stOligo();
		if (startInterval > end ) return -1;
		// must overlap
		//left side
		if (start>=startInterval && start<=endInterval && end > endInterval) return endInterval - start +1;
		//right side
		if (end >= startInterval && end <= endInterval && start < startInterval) return end-startInterval +1;
		//contained within, return length of shortest
		int lengthInterval = endInterval-startInterval +1;
		int lengthThis = end - start +1;
		if (lengthInterval < lengthThis ) return lengthInterval;
		return lengthThis;
	}
	
	/**Parses a file for a tab delimited list of at minimum, chrom, start, stop, (optional) notes into a
	 * GenomicRegion[].*/
	public static GenomicRegion[] parseRegions(File file){
		GenomicRegion[] regions =null;
		String line = null;
		try{
			BufferedReader in = IO.fetchBufferedReader(file);
			String[] tokens;
			ArrayList regionsAL = new ArrayList();
			
			String chromosome;
			int start;
			int stop;
			String notes = null;
			
			//chrom, start, stop, notes
			//  0      1      2     3      
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.length() == 0 || line.startsWith("#")) continue;
				tokens = line.split("\\s+");
				if (tokens.length<3) continue;
				chromosome = tokens[0];
				start = Integer.parseInt(tokens[1]);
				stop = Integer.parseInt(tokens[2]);
				if (stop < start) throw new Exception("Stop position is less than start!");
				if (tokens.length>3){
					StringBuffer sb = new StringBuffer();
					for (int i=3; i<tokens.length; i++) {
						sb.append(tokens[i]);
						sb.append(" ");
					}
					notes = sb.toString();
				}
				regionsAL.add(new GenomicRegion(chromosome, start, stop, notes));
			}
			regions = new GenomicRegion[regionsAL.size()];
			regionsAL.toArray(regions);
			
		}catch (Exception e){
			System.out.println ("Problem parsing this line-> "+line);
			e.printStackTrace();
		}
		return regions;
	}
	
	public String getChromosome() {
		return chromosome;
	}
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public HashSet getHits() {
		return hits;
	}
	public void setHits(HashSet hits) {
		this.hits = hits;
	}
	public String getNotes() {
		return notes;
	}
	public void setNotes(String notes) {
		this.notes = notes;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}

	public boolean[] getGcContent() {
		return gcContent;
	}

	public void setGcContent(boolean[] gcContent) {
		this.gcContent = gcContent;
	}
}
