package edu.utah.kegg.api;

public class KeggApiGene {
	
	private String name;
	private String id;
	
	public KeggApiGene(String name, String id) {
		this.name = name;
		this.id = id;
	}
	
	public String toString() {
		return name+":"+id;
	}
	
	public static String toString(KeggApiGene[] genes, String separator) {
		StringBuilder sb = new StringBuilder(genes[0].toString());
		for (int i=1; i<genes.length; i++) {
			sb.append(separator);
			sb.append(genes[i].toString());
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
