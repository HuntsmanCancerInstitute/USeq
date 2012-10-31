package util.bio.annotation;
import java.io.*;
import java.util.*;

import edu.utah.seq.useq.data.Region;



/**
A GeneGroup is the top organizing level of a "Gene". It may contain multiple transGrps
(one transcript, possibly a translation, each with one or more exons) and a composite/ summary GeneRep (a gene representation broken into
coding and non-coding regions) .
Thus GeneGroup -> one GeneRep + one or more TransGroups -> transcript, exon, introns....
*/
public class GeneGroup implements Comparable, Serializable {
	//fields
	private String name;
	private String chromosome;
	private int start;  
	private int end;
	private String type;
	private TransGroup[] transGrps;
	private GeneRep geneRep;
	private int orientation; //1 is forward (start stop 5..25), -1 is reverse (start stop 25..5)
	private String attributes; //String rep of attributes hash
	
	//Constructor
	public GeneGroup(String name, String chromosome, int start, int end, String type, TransGroup[] transGrps, int orientation, String attributes){
		this.name=name;
		this.chromosome=chromosome;
		this.start = start;
		this.end = end;
		this.type=type;
		this.transGrps=transGrps;
		this.orientation = orientation;
		this.attributes = attributes;
		geneRep = new GeneRep (transGrps, this);
	}
	//primary methods
	public int compareTo(Object other){
		GeneGroup otherGG = (GeneGroup)other;
		//check non overlaping
		if (start<otherGG.start && end<otherGG.end) return -1;
		if (start>otherGG.start && end>otherGG.end) return 1;
		//overlapping then set shorter = 1
			// this should put larger on bottom
			int len = end-start;
			int otherLen = otherGG.end-otherGG.start;
			if (len<otherLen) return -1;
			if (len>otherLen) return 1;
			return 0;
		}
	
	public String fetchSumaryLine(){
		StringBuffer sb = new StringBuffer();
		sb.append(name);
		sb.append("\tchr");
		sb.append(chromosome);
		sb.append(":");
		sb.append(start);
		sb.append("-");
		sb.append(end);
		sb.append("\tori: ");
		if (orientation==1) sb.append("+");
		else sb.append("-");
		return sb.toString();
	}
	public String toString(){
		StringBuffer x = new StringBuffer(
			"\nGeneGrp: "+name+", MapFeature Type: "+type+", Location: "
			+chromosome+ ":"+start+"-"+end+", Orientation: "+convertNumOrientationToStrand(orientation)
			+"\nAttributes: "+attributes);
		x.append(geneRep.toString());
		for (int i=0; i<transGrps.length; i++){
			x.append(transGrps[i].toString());
		}
		x.append("\n");
		return x.toString();
	}
	//getters
	public int getOrientation() {
		return orientation;
	}
	public GeneRep getGeneRep(){
		return geneRep;
	}
	public TransGroup[] getTransGrps(){
		return transGrps;
	}
	public String getLabel(){
		return name;
	}
	public String getChromosome() {
		return chromosome;
	}
	public String getAttributes() {
		return attributes;
	}
	public String getName() {
		return name;
	}
	public int getEnd() {
		return end;
	}
	public int getStart() {
		return start;
	}
	
	
	//static methods for general use
	public static int convertPlusToNumOrientation(String orientation){
		if (orientation.equals("+")) return 1;
		if (orientation.equals("-")) return -1;
		return 0;
	}
	
	public static String convertNumOrientationToStrand(int orientation){
		if (orientation==1) return "+";
		if (orientation==-1) return "-";
		return ".";
	}
	public static int[][] extractStartEndsInts(ExonIntron[] exons){
		int len = exons.length;
		int[][] extra = new int[len][2];
		for (int i=0; i<len; i++){
			extra[i]= new int[]{exons[i].getStart(),exons[i].getEnd()};
		}
		return extra;
	}
	public static ArrayList extractStartEnds(ExonIntron[] exons){
		int len = exons.length;
		ArrayList al = new ArrayList(len);
		for (int i=0; i<len; i++){
			al.add(new int[]{exons[i].getStart(),exons[i].getEnd()});
		}
		return al;
	}
	/**Takes an ArrayList of exon startStop int[2]'s, returns an ArrayList of int[2]'s representing
	 the start,stop for each intron*/
	public static ArrayList fetchIntrons(ArrayList exons) {
		int numIntrons = exons.size() - 1;
		ArrayList introns = new ArrayList();
		for (int i = 0; i < numIntrons; i++) {
			int[] ex1 = (int[]) exons.get(i);
			int[] ex2 = (int[]) exons.get(i + 1);
			introns.add(new int[] { ex1[1] + 1, ex2[0] - 1 });
		}
		return introns;
	}
	
	public static String ALToString(ArrayList al) {
		StringBuffer sb = new StringBuffer();
		int len = al.size();
		for (int i = 0; i < len; i++) {
			int[] ss = (int[]) al.get(i);
			sb.append("\n  " + ss[0] + "-" + ss[1]);
		}
		return sb.toString();
	}
	
	/**Sorts an ArrayList of int[2]'s representing start, stop*/
	public static ArrayList sortSegments(ArrayList segs){
		int num = segs.size();
		ArrayList sorted = new ArrayList (num);
		Region[] se = new Region[num];
		for (int i=0; i<num; i++){
			int[] stst = (int[])segs.get(i);
			se[i] = new Region(stst[0], stst[1]);
		}
		Arrays.sort(se);
		for (int i=0; i<num; i++){
			sorted.add(se[i].getStartStop());
		}
		return sorted;
	}
	
	/**Takes an ArrayList of int[2]'s representing start stop coords
	 * kills the dups, and those cotained within another seg, will
	 * merge overlapping segments.
	 */
	public static ArrayList joinSegments(ArrayList segs) {
		ArrayList mergedSegs = new ArrayList();
		int len = segs.size();
		int[] seg1 = new int[2];
		while (len > 0) {
			boolean test = true;
			seg1 = (int[]) segs.remove(0);
			len--;
			for (int i = 0; i < len; i++) {
				int[] seg2 = (int[]) segs.get(i);
				//if it entirely overlaps, contained within
				if (seg1[0] >= seg2[0] && seg1[1] <= seg2[1]) {
					test = false;
					break;
				} //don't save
				else if (seg1[0] <= seg2[0] && seg1[1] >= seg2[1]) {
					segs.remove(i);
					len--;
					i--;
					continue;
				}
				//if immediately adjacent fuse
				if (seg1[1] == (seg2[0] - 1)) {
					seg1 = new int[] { seg1[0], seg2[1] };
					segs.remove(i);
					len--;
					i = -1; //start over
					continue;
				} else if (seg1[0] == (seg2[1] + 1)) {
					seg1 = new int[] { seg2[0], seg1[1] };
					segs.remove(i);
					len--;
					i = -1; //start over
					continue;
				}
				//if it doesn't overlap then continue
				if (seg1[1] < seg2[0] || seg1[0] > seg2[1]) {
					continue;
				}
				//if it partially overlaps fuse together
				if (seg1[1] >= seg2[0] && seg1[0] <= seg2[0]) {
					seg1 = new int[] { seg1[0], seg2[1] };
					segs.remove(i);
					len--;
					i = -1; //start over
				} else if (seg1[0] <= seg2[1] && seg1[1] >= seg2[1]) {
					seg1 = new int[] { seg2[0], seg1[1] };
					segs.remove(i);
					len--;
					i = -1; //start over 
				}
			}
			//reached stop add exon
			if (test) {
				mergedSegs.add(seg1);
			}
		}
		return mergedSegs;
	}
	public String getType() {
		return type;
	}
	public void setOrientation(int orientation) {
		this.orientation = orientation;
	}
}
