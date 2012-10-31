package util.bio.annotation;
import java.io.*;

import util.bio.parsers.gff.*;

/**
 * Transcript relevant information, mRNA.
 */
public class Transcript implements Serializable{
	//fields
	private int start;
	private int end;
	private String name;
	private int orientation;
	
	//constructors
	public Transcript(int start, int end, String name, int orientation) {
		this.start = start;
		this.end = end;
		this.name = name;
		this.orientation = orientation;
	}
	
	public Transcript (Gff3Feature f){
		start = f.getStart();
		end = f.getEnd();
		name = f.getId();
		orientation = GeneGroup.convertPlusToNumOrientation(f.getStrand());
	}
	
	public String toString(){
			return "Transcript: "+start+"-"+ end;
		}
	//getters
	public String getName() {
		return name;
	}


	public int getOrientation() {
		return orientation;
	}
	public int getEnd() {
		return end;
	}
	public int getStart() {
		return start;
	}
}
