package edu.utah.kegg.pathway;



public class SelectGene implements Comparable<SelectGene> {

	//fields
	private String geneSymbol;
	private String log2Rto;

	public SelectGene(String[] geneRto){
		geneSymbol = geneRto[0];
		log2Rto = geneRto[1];
	}
	public SelectGene(String geneSymbol, String log2Rto){
		this.geneSymbol = geneSymbol;
		this.log2Rto = log2Rto;
	}
	public int compareTo(SelectGene o) {
		return geneSymbol.compareTo(o.geneSymbol);
	}
	public String getGeneSymbol() {
		return geneSymbol;
	}
	public String getLog2Rto() {
		return log2Rto;
	}

}
