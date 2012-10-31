package edu.utah.seq.data;

public class InfoPoint extends Point {
	
	//fields
	private Info info;
	
	public InfoPoint (int position, float score, Info info){
		super(position, score);
		this.info = info;
	}

	public Info getInfo() {
		return info;
	}

	public void setInfo(Info info) {
		this.info = info;
	}

}
