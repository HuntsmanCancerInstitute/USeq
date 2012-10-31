package edu.expr;

public class CorrelationWindow {
	
	//fields
	private float score;
	private float pvalue = 0;	//default is zero
	private String chromosome;
	private int start;
	private int stop;
	private String[] names;
	
	//constructor
	public CorrelationWindow(float score, String chromosome, int start, int stop, String[] names){
		this.score=score;
		this.chromosome = chromosome;
		this.start=start;
		this.stop=stop;
		this.names=names;
	}
	/**Genes NumGenes Chromosome Start Stop CorrelationScore -10Log10(pvalue)*/
	public String toString(){
		StringBuffer sb = new StringBuffer(names[0]);
		for (int i=1; i<names.length; i++){
			sb.append(",");
			sb.append(names[i]);
		}
		
		sb.append("\t");
		sb.append(chromosome);
		sb.append("\t");
		sb.append(start);
		sb.append("\t");
		sb.append(stop);
		sb.append("\t");
		sb.append(names.length);
		sb.append("\t");
		sb.append(score);
		sb.append("\t");
		sb.append(pvalue);
		return sb.toString();
	}
	
	public float getPvalue() {
		return pvalue;
	}

	public void setPvalue(float pvalue) {
		this.pvalue = pvalue;
	}

	public String[] getNames() {
		return names;
	}

	public float getScore() {
		return score;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	public String getChromosome() {
		return chromosome;
	}
}
