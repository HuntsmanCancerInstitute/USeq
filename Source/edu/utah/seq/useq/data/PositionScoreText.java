package edu.utah.seq.useq.data;

/**
 * Simple position, score, text object 
 * @author david.nix@hci.utah.edu*/
public class PositionScoreText extends PositionScore {
	//fields
	protected String text;


	//constructors
	public PositionScoreText (int position, float score, String text){
		super(position,score);
		this.text = text;
	}

	public String toString(){
		return position+"\t"+score+"\t"+text;
	}

	public String getText() {
		return text;
	}

	public void setText(String text) {
		this.text = text;
	}
}
