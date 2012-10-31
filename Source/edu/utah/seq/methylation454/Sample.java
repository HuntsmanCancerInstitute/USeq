package edu.utah.seq.methylation454;

public class Sample {
	//fields
	private String id;
	private Amplicon[] amplicons;
	
	//constructor
	public Sample (String id){
		this.id = id;
	}
	
	//methods
	public String toString(){
		StringBuffer sb = new StringBuffer();
		String tab = "\t";
		sb.append(id); sb.append(tab);
		for (int i=0; i< amplicons.length; i++){
			sb.append(amplicons[i]);
			sb.append("\n");
			sb.append(tab);
		}
		return sb.toString();
	}
	
	public Amplicon[] getAmplicons() {
		return amplicons;
	}
	public void setAmplicons(Amplicon[] amplicons) {
		this.amplicons = amplicons;
	}
	public String getId() {
		return id;
	}
}
