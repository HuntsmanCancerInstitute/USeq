package edu.utah.seq.methylation454;

public class BaseRead {
	//fields
	private String id;
	private String strand;		//+ or -
	private String sequence;	
	
	//constructor
	public BaseRead (DataLine line){
		this.id = line.getRead();
		this.strand = line.getStrand();
		this.sequence = line.getSequence();
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		String tab = "\t";
		sb.append(id);  sb.append(tab);
		sb.append(strand); sb.append(tab);
		sb.append(sequence);
		return sb.toString();
	}

	public String getId() {
		return id;
	}

	public String getSequence() {
		return sequence;
	}

	public String getStrand() {
		return strand;
	}
}
