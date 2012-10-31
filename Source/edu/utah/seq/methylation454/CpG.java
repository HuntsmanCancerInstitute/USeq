package edu.utah.seq.methylation454;


public class CpG {
	
	//fields
	private String id;
	private BaseRead[] reads;
	private double fractionC;
	private int numberCs = -1;
	
	//constructor
	public CpG (String id){
		this.id = id;
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		String tab = "\t\t\t";
		sb.append(id); sb.append("\t");
		for (int i=0; i< reads.length; i++){
			sb.append(reads[i]);
			sb.append("\n");
			sb.append(tab);
		}
		return sb.toString();
	}
	public String getId() {
		return id;
	}
	public BaseRead[] getReads() {
		return reads;
	}
	public void setReads(BaseRead[] reads) {
		this.reads = reads;
	}
	public int getNumberReads(){
		return reads.length;
	}
	public double fractionC(){
		getNumberCs();
		return (double) numberCs/(double)reads.length;
	}
	public int getNumberCs(){
		if (numberCs == -1){
			numberCs = 0;
			for (int i=0; i< reads.length; i++) if (reads[i].getSequence().equals("C")) numberCs++;
		}
		return numberCs;
	}
}
