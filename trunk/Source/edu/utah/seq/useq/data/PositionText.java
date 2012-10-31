package edu.utah.seq.useq.data;

/**
 * Simple position text object. 
 * @author david.nix@hci.utah.edu*/
public class PositionText extends Position {
	//fields
	protected String text;


	//constructors
	public PositionText (int position, String text){
		super(position);
		this.text = text;
	}

	public String toString(){
		return position+"\t"+text;
	}

	public String getText() {
		return text;
	}

	public void setText(String text) {
		this.text = text;
	}
}
