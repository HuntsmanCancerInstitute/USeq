package edu.utah.seq.data;
import java.util.*;

public class QseqLane implements Comparable<QseqLane>{
	
	//fields
	private int laneNumber = 0;
	private ArrayList<QseqFileSet> qseqFiles =  new ArrayList<QseqFileSet>();
	
	//constructors
	public QseqLane (String laneNumberString){
		laneNumber = Integer.parseInt(laneNumberString);
	}
	
	//methods
	public int compareTo(QseqLane other) {
		if (other.laneNumber > this.laneNumber) return 1;
		if (other.laneNumber < this.laneNumber) return -1;
		return 0;
	}
	public int getLaneNumber() {
		return laneNumber;
	}
	public ArrayList<QseqFileSet> getQseqFiles() {
		return qseqFiles;
	}

	

}
