package edu.utah.seq.parsers;

import edu.utah.seq.data.Point;

public class BaseObservation implements Comparable<BaseObservation>{
	//fields
	private Integer position;
	private String bedLine;
	private Point point;
	private String genSeq;
	private double score;
	private boolean converted;

	public BaseObservation (int position, String bedLine, Point point, String genSeq, double score, boolean converted){
		this.position = new Integer(position);
		this.bedLine = bedLine;
		this.point = point;
		this.genSeq = genSeq;
		this.score = score;
		this.converted = converted;
	}
	
	public String toString(){
		return position +"\t"+ genSeq +"\t"+ score +"\t"+ converted;
	}

	public int compareTo(BaseObservation other){
		return position.compareTo(other.position);
	}

	public Integer getPosition() {
		return position;
	}

	public String getBedLine() {
		return bedLine;
	}

	public Point getPoint() {
		return point;
	}

	public String getGenSeq() {
		return genSeq;
	}

	public double getScore() {
		return score;
	}

	public boolean isConverted() {
		return converted;
	}
}