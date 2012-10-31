package edu.utah.seq.useq.data;

/** @author david.nix@hci.utah.edu*/
public class RegionScore extends Region {

	//fields
	private static final long serialVersionUID = 1L;
	protected float score;

	//constructor
	public RegionScore (int start, int stop, float score){
		super(start, stop);
		this.score = score;
	}

	//methods
	public String toString(){
		return start+"\t"+stop+"\t"+score;
	}

	public float getScore() {
		return score;
	}

	public void setScore(float score) {
		this.score = score;
	}

}
