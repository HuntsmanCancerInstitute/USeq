package edu.utah.kegg.api;

public class KeggApiPathway {
	
	private String id;
	private String name;
	
	public KeggApiPathway(String id, String name) {
		this.name = name;
		this.id = id;
	}
	
	public String toString() {
		return id+":"+name;
	}
	
	public static String toString(KeggApiPathway[] classes, String separator) {
		StringBuilder sb = new StringBuilder(classes[0].toString());
		for (int i=1; i<classes.length; i++) {
			sb.append(separator);
			sb.append(classes[i].toString());
		}
		return sb.toString();
	}

	public String getName() {
		return name;
	}

	public String getId() {
		return id;
	}

}
