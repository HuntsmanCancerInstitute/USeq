package util.bio.annotation;

import java.util.*;
import java.io.*;

import util.bio.parsers.gff.Gff3Feature;

/**
Translation specific information, protein but referenced with genomic coordinates thus some
bits code for protein, some don't physically exist.
 */
public class Translation implements Serializable{
	//fields
	private int start;
	private int end;
	TransGroup transGrpRef;
	int orientation;	//1 is forward (start stop 5..25), -1 is reverse (start stop 25..5)
	ArrayList codingSegments = new ArrayList();	//of int[2]'s; parts of the "translation" that code for protein
	ArrayList introns = new ArrayList(); 	//of int[2];

	//constructors
	public Translation(int start, int end, int orientation){
		this.start = start;
		this.end = end;
		this.orientation = orientation;
	}
	public Translation (Gff3Feature f){
		this.start = f.getStart();
		this.end = f.getEnd();
		orientation = GeneGroup.convertPlusToNumOrientation(f.getStrand());
	}
	
	/**Sets TransGroup ref but also fires the makeSegments method*/
	public void setTransGrp(TransGroup ref) {
		transGrpRef = ref;
		makeSegments();
	}

	//primary methods
	public void makeSegments() {
		if (transGrpRef.getExons()==null){
			System.out.println("ERROR: "+transGrpRef);
			codingSegments = new ArrayList();
			codingSegments.add(new int[]{0,0});
			introns = new ArrayList();
			return;
		}
		codingSegments = GeneGroup.extractStartEnds(transGrpRef.getExons());
		//modify flanking exons
		boolean flag = true;
		while (flag) {
			int[] exon = (int[]) codingSegments.get(0);
			if (exon[1] < start) codingSegments.remove(0);
			else {
				exon[0] = start;
				flag = false;
			}
		}
		flag = true;
		while (flag) {
			int[] exon = (int[]) codingSegments.get(codingSegments.size() - 1);
			if (exon[0] > end)
				codingSegments.remove(codingSegments.size() - 1);
			else {
				exon[1] = end;
				flag = false;
			}
		}
		introns = GeneGroup.fetchIntrons(codingSegments);
	}

	public String toString() {
		return "Translation: " + start + "-" + end;
	}

	public static String toStringALInts(ArrayList al) {
		StringBuffer sb = new StringBuffer();
		int len = al.size();
		for (int i = 0; i < len; i++) {
			int[] ints = (int[]) al.get(i);
			sb.append(ints[0] + "-" + ints[1] + " ");
		}
		return sb.toString();
	}
	//getters
	public ArrayList getCodingSegments() {
		return codingSegments;
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
	public ArrayList getIntrons() {
		return introns;
	}
}
