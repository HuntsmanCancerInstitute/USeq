package util.bio.annotation;
import java.io.*;
import java.util.*;
import util.bio.parsers.*;
import util.bio.parsers.gff.Gff3Feature;
import util.bio.seq.MakeSpliceJunctionFasta;
import util.gen.Misc;

/**
 * ExonIntron info. 
 */
public class ExonIntron implements Comparable, Serializable {
	//fields
	private int start;
	private int end;
	private String name;
	private String parent;
	private String sequence; 
	private int index;
	private ExonIntron[] exonsIntrons;
	private int maxNumberSplices;

	//constructors
	public ExonIntron (int start, int end){
		this.start = start;
		this.end = end;
	}
	public ExonIntron (int start, int end, String name){
		this.start = start;
		this.end = end;
		this.name = name;
	}
	public ExonIntron (Gff3Feature f){
		start = f.getStart();
		end = f.getEnd();
		name = f.getId();
		parent = f.getParent();
	}

	//primary methods

	public ExonIntron clone(){
		ExonIntron clone = new ExonIntron(this.start, this.end);
		clone.index = this.index;
		if (this.name != null) clone.name = new String(this.name);
		if (this.parent != null) clone.parent = new String(this.parent);
		if (this.sequence != null) clone.sequence = new String(this.sequence);
		clone.exonsIntrons = this.exonsIntrons;
		clone.maxNumberSplices = this.maxNumberSplices;
		return clone;
	}

	public static ExonIntron[] clone(ExonIntron[] ei){
		ExonIntron[] clones = new ExonIntron[ei.length];
		for (int i=0; i< clones.length; i++){
			clones[i] = ei[i].clone();
		}
		return clones;
	}
	
	/**This just fetches adjacent sequence to the request length.  No expansion into all permutations of downstream exons.
	 * Use this to build known splice junctions from a transcript.*/
	public ArrayList<SubSequence> fetch3PrimeSubSeqNoExpand (int numberBases){
		ArrayList<SubSequence> concatsAL;
		//entirely contained so no need to fetch additional adjacent exons just return trimmed sequence
		if (sequence.length() >= numberBases) {
			concatsAL = new ArrayList<SubSequence>();
			String s = new String(sequence.substring(0, numberBases));
			concatsAL.add(new SubSequence(s, start, start+numberBases));
		}
		//is this the last exon? then likewise nothing to fetch just return seq
		else if (index == (exonsIntrons.length -1)) {
			concatsAL = new ArrayList<SubSequence>();
			concatsAL.add(new SubSequence(sequence, start, end));
		}
		//need to fetch more sequence
		else {
			int diff = numberBases - sequence.length();
			//fetch adjacent exon
			ExonIntron ei = fetch3PrimeExonIntron();
			//fetch the subseqs
			concatsAL = ei.fetch3PrimeSubSeqNoExpand(diff);
			//prepend this exon 
			concatsAL.add(0, new SubSequence(sequence, start, end));
		}
		return concatsAL;
	}
	
	/**This just fetches adjacent sequence to the request length.  No expansion into all permutations of downstream exons.
	 * Use this to build known splice junctions from a transcript.*/
	public ArrayList<SubSequence> fetch5PrimeSubSeqNoExpand (int numberBases){
		ArrayList<SubSequence> concatsAL;
		//entirely contained so no need to fetch additional adjacent exons just return trimmed sequence
		if (sequence.length() >= numberBases) {
			concatsAL = new ArrayList<SubSequence>();
			int diff = sequence.length() - numberBases;
			String s = new String(sequence.substring(diff, sequence.length()));
			concatsAL.add(new SubSequence(s, end-numberBases, end));
		}
		//is this the 1st exon? then likewise nothing to fetch just return seq
		else if (index == 0) {
			concatsAL = new ArrayList<SubSequence>();
			concatsAL.add(new SubSequence(sequence, start, end));
		}
		//need to make some concatinates
		else {
			int diff = numberBases - sequence.length();
			//fetch adjacent exon
			ExonIntron ei = fetch5PrimeExonIntron();
			//fetch the subseqs
			concatsAL = ei.fetch5PrimeSubSeqNoExpand(diff);
			concatsAL.add(new SubSequence(sequence, start, end));
		}
		return concatsAL;
	}

	/**This builds all possible downstream junctions.  Should be behind a thread so you can kill it if it goes too long.
	 * Returns null if the maxNumberSplices is exceeded*/
	public ArrayList<SubSequence>[] fetch3PrimeSubSeq (int numberBases){
		ArrayList<ArrayList<SubSequence>> concatsAL = new ArrayList<ArrayList<SubSequence>>();
		//entirely contained so no need to fetch additional adjacent exons just return trimmed sequence
		if (sequence.length() >= numberBases) {
			ArrayList<SubSequence> subSeq = new ArrayList<SubSequence>();
			String s = new String(sequence.substring(0, numberBases));
			subSeq.add(new SubSequence(s, start, start+numberBases));
			concatsAL.add(subSeq);
		}
		//is this the last exon? then likewise nothing to fetch just return seq
		else if (index == (exonsIntrons.length -1)) {
			ArrayList<SubSequence> subSeq = new ArrayList<SubSequence>();
			subSeq.add(new SubSequence(sequence, start, end));
			concatsAL.add(subSeq);
		}
		//need to make some concatinates
		else {
			int diff = numberBases - sequence.length();

			//for each adjacent exon
			for (ExonIntron ei : fetch3PrimeExonsIntrons()){
				//for each sequence
				ArrayList<SubSequence>[] seqs = ei.fetch3PrimeSubSeq(diff);
				//check if too many
				if (seqs == null || seqs.length > maxNumberSplices) {
					if (seqs == null) System.out.println ("seqs null");
					else System.out.println("seqs length "+seqs.length+" max num "+maxNumberSplices);
					return null;
				}
				for (ArrayList<SubSequence> subSeq : seqs) {
					//append on seq
					subSeq.add(0, new SubSequence(sequence, start, end));
					//add to concatsAL
					concatsAL.add(subSeq);
				}
			}
		}
		ArrayList<SubSequence>[] concats = new ArrayList[concatsAL.size()];
		concatsAL.toArray(concats);
		return concats;
	}

	public ArrayList<SubSequence>[] fetch5PrimeSubSeq (int numberBases){
		ArrayList<ArrayList<SubSequence>> concatsAL = new ArrayList<ArrayList<SubSequence>>();
		//entirely contained so no need to fetch additional adjacent exons just return trimmed sequence
		if (sequence.length() >= numberBases) {
			int diff = sequence.length() - numberBases;
			ArrayList<SubSequence> subSeq = new ArrayList<SubSequence>();
			String s = new String(sequence.substring(diff, sequence.length()));
			subSeq.add(new SubSequence(s, end-numberBases, end));
			concatsAL.add(subSeq);
		}
		//is this the 1st exon? then likewise nothing to fetch just return seq
		else if (index == 0) {
			ArrayList<SubSequence> subSeq = new ArrayList<SubSequence>();
			subSeq.add(new SubSequence(sequence, start, end));
			concatsAL.add(subSeq);
		}
		//need to make some concatinates
		else {
			int diff = numberBases - sequence.length();
			//for each adjacent exon
			for (ExonIntron ei : fetch5PrimeExonsIntrons()){
				//for each sequence
				ArrayList<SubSequence>[] seqs = ei.fetch5PrimeSubSeq(diff);
				//check if too many
				if (seqs == null || seqs.length > maxNumberSplices) return null;
				for (ArrayList<SubSequence> subSeq : seqs) {
					//append on seq to end
					subSeq.add(new SubSequence(sequence, start, end));
					//add to concatsAL
					concatsAL.add(subSeq);
				}
			}
		}
		ArrayList<SubSequence>[] concats = new ArrayList[concatsAL.size()];
		concatsAL.toArray(concats);
		return concats;
	}

	/**Returns null if no exons/introns to the 3'*/
	public ExonIntron[] fetch3PrimeExonsIntrons(){
		int num = exonsIntrons.length - index - 1;
		if (num == 0) return null;
		ExonIntron[] ei = new ExonIntron[num];
		int i = 0;
		for (int x = index +1; x < exonsIntrons.length; x++){
			ei[i++] = exonsIntrons[x];
		}
		return ei;
	}

	/**Returns null if no exons/introns to the 5'*/
	public ExonIntron[] fetch5PrimeExonsIntrons(){
		if (index == 0) return null;
		ExonIntron[] ei = new ExonIntron[index];
		for (int x = 0; x < index; x++){
			ei[x] = exonsIntrons[x];
		}
		return ei;
	}
	
	/**Returns null if no exons/introns to the 3'*/
	public ExonIntron fetch3PrimeExonIntron(){
		int num = exonsIntrons.length - index - 1;
		if (num == 0) return null;
		return exonsIntrons[index+1];
	}
	
	/**Returns null if no exons/introns to the 5'*/
	public ExonIntron fetch5PrimeExonIntron(){
		if (index == 0) return null;
		return exonsIntrons[index -1];
	}

	/**Assumes interbase coordinates.*/
	public String getSequence(String chromSeq){
		return new String( chromSeq.substring(start, end));
	}

	/**Assumes interbase coordinates.*/
	public boolean contains(int position){
		if (position< start || position >= end) return false;
		return true;
	}

	/**Assumes interbase coordinates*/
	public boolean intersects(int start, int stop){
		if (stop < this.start || start >= this.end) return false;
		return true;
	}
	public static String toString(ExonIntron[] e){
		StringBuilder sb = new StringBuilder();
		sb.append(e[0].getStart());
		sb.append(":");
		sb.append(e[0].getEnd());
		for (int i=1; i<e.length; i++){
			sb.append("_");
			sb.append(e[0].getStart());
			sb.append(":");
			sb.append(e[0].getEnd());
		}
		return sb.toString();
	}
	public static ExonIntron[] arrayList2Array(ArrayList<ExonIntron> al){
		ExonIntron[] e = new ExonIntron[al.size()];
		al.toArray(e);
		return e;
	}
	public String toString(){
		return "ExonIntron: "+start+"-"+end;
	}
	public int compareTo(Object other){
		ExonIntron otherExon = (ExonIntron)other;
		//sort by start position
		if (start<otherExon.start) return -1;
		if (start>otherExon.start) return 1;
		// if same start, sort by length, smaller to larger
		int len = end-start;
		int otherLen = otherExon.end-otherExon.start;
		if (len<otherLen) return -1;
		if (len>otherLen) return 1;
		return 0;
	}

	/**Assumes interbase coordinates.*/
	public static ExonIntron[] merge(ExonIntron[] a, ExonIntron[] b){
		//find min and max
		int[] minMaxOne = minMax(a);
		int[] minMaxTwo = minMax(b);
		int min = minMaxOne[0];
		int max = minMaxOne[1];
		if (minMaxTwo[0]< min) min = minMaxTwo[0];
		if (minMaxTwo[1]> max) max = minMaxTwo[1];
		//load boolean array
		int length = max-min+1;
		boolean[] fetchFalse = new boolean[length];
		Arrays.fill(fetchFalse, true);
		for (int i=0; i< a.length; i++){
			int start = a[i].start -min;
			int stop = a[i].end-min;		
			for (int j=start; j< stop; j++) fetchFalse[j] = false;
		}
		for (int i=0; i< b.length; i++){
			int start = b[i].start -min;
			int stop = b[i].end-min;		
			for (int j=start; j< stop; j++) fetchFalse[j] = false;
		}
		//retrieve blocks, ends included

		int[][] blocks = ExportIntergenicRegions.fetchFalseBlocks(fetchFalse, 0, 0);
		//convert to interbase coordinate exon introls
		ExonIntron[] ei = new ExonIntron[blocks.length];
		int minPlusOne = min +1;
		for (int i=0; i< ei.length; i++){
			ei[i] = new ExonIntron(blocks[i][0]+min, blocks[i][1]+minPlusOne);
		}
		return ei;
	}

	public static boolean intersect(ExonIntron[] a, ExonIntron[] b){
		for (int i=0; i< a.length; i++){
			int start = a[i].start;
			int end = a[i].end;
			for (int j=0; j< b.length; j++){
				if (b[j].intersects(start, end)) return true;
			}
		}
		return false;
	}

	public static int[] minMax(ExonIntron[] ei){
		int min = ei[0].getStart();
		int max = ei[0].getEnd();
		for (int i=1; i< ei.length; i++){
			if (ei[i].start < min) min = ei[i].start;
			if (ei[i].end > max) max = ei[i].end;
		}
		return new int[]{min,max};
	}

	//getters
	/**Assumes interbase coordinants.*/
	public int getLength(){
		return (end-start);
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
	public String getParent() {
		return parent;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public void setStart(int start) {
		this.start = start;
	}
	public String getSequence() {
		return sequence;
	}
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	public int getIndex() {
		return index;
	}
	public void setIndex(int index) {
		this.index = index;
	}
	public ExonIntron[] getExonsIntrons() {
		return exonsIntrons;
	}
	public void setExonsIntrons(ExonIntron[] exonsIntrons) {
		this.exonsIntrons = exonsIntrons;
	}
	public void setName(String name) {
		this.name = name;
	}
	public void setParent(String parent) {
		this.parent = parent;
	}
	public int getMaxNumberSplices() {
		return maxNumberSplices;
	}
	public void setMaxNumberSplices(int maxNumberSplices) {
		this.maxNumberSplices = maxNumberSplices;
	}
}
