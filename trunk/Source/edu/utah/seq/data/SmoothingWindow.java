package edu.utah.seq.data;
import java.io.*;

/**Coordinates are interbase!*/
public class SmoothingWindow implements Serializable{
	//fields
	int start;
	int stop;
	float[] scores;
	//change this value when making major changes that should void any saved Objects
	public static final long serialVersionUID = 1;
	
	//constructors
	public SmoothingWindow(){};
	
	public SmoothingWindow(int start, int stop, float[] scores){
		this.start = start;
		this.stop = stop;
		this.scores = scores;
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append(start);
		sb.append("\t");
		sb.append(stop);
		if (scores !=null){
			for (int i=0; i< scores.length; i++){
				sb.append("\t");
				sb.append(scores[i]);
			}
		}
		return sb.toString();
	}

	public float[] getScores() {
		return scores;
	}

	public void setScores(float[] scores) {
		this.scores = scores;
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
}
