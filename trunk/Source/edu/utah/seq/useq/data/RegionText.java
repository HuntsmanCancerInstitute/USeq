package edu.utah.seq.useq.data;

/** @author david.nix@hci.utah.edu*/
public class RegionText extends Region{

	//fields
	protected String text;  //or rank
	private static final long serialVersionUID = 1L;

	//constructor
	public RegionText (int start, int stop, String text){
		super(start, stop);
		this.text = text;
	}

	public String toString(){
		return start+"\t"+stop+"\t"+text;
	}
	public String getText() {
		return text;
	}
	public void setText(String text) {
		this.text = text;
	}	
}
