package edu.expr;

public class MiRNA {
	
	//fields
	private String name;
	private MiRNATarget[] targets;
	
	//constructor
	public MiRNA (String name, MiRNATarget[] targets) {
		this.name = name;
		this.targets = targets;
	}

	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append (name);
		sb.append ("\n");
		for (int i=0; i< targets.length; i++){
			sb.append("\t");
			sb.append(targets[i].toString());
			sb.append("\n");
		}
		return sb.toString();
	}
	
	public String getName() {
		return name;
	}

	public MiRNATarget[] getTargets() {
		return targets;
	}
	
}
