package edu.utah.seq.methylation454;

public class Amplicon {
	//fields
	private String id;
	private CpG[] cpGs;
	
	public Amplicon (String id){
		this.id = id;
	}
	
	//methods
	public String toString(){
		StringBuffer sb = new StringBuffer();
		String tab = "\t\t";
		sb.append(id); sb.append ("\t");
		for (int i=0; i< cpGs.length; i++){
			sb.append(cpGs[i]);
			sb.append("\n");
			sb.append(tab);
		}
		return sb.toString();
	}

	
	public CpG[] getCpGs() {
		return cpGs;
	}
	public String getId() {
		return id;
	}

	public void setCpGs(CpG[] cpGs) {
		this.cpGs = cpGs;
	}
}
