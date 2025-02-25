package edu.utah.kegg.pathway;


public class GeneCount implements Comparable<GeneCount>{
	double fraction = 0;
	String name = null;
	
	public GeneCount(String name, double fraction) {
		this.name = name;
		this.fraction = fraction;
	}
	public int compareTo(GeneCount o) {
		if (this.fraction> o.fraction) return -1;
		if (this.fraction< o.fraction) return 1;
		return this.name.compareTo(o.name);
	}

}
