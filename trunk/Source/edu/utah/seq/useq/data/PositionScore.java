package edu.utah.seq.useq.data;

/** @author david.nix@hci.utah.edu*/
public class PositionScore extends Position {
	//fields
	protected float score;

	//constructors
	public PositionScore (int position, float score){
		super(position);
		this.score = score;
	}

	public String toString(){
		return position+"\t"+score;
	}

	public float getScore() {
		return score;
	}

	public void setScore(float score) {
		this.score = score;
	}
}
