package edu.utah.seq.vcf.splice;

public class SpliceJunction {
	
	//fields
	
	/**Three letter code: Gain or Damaged; 5' or 3' relative to gene not genomic; Exonic or Intronic or Splice
	 * e.g. G5E - gain 5' splice in exonic; D3S - damaged 3' splice in splice */
	private String type = null;
	//G or D
	private char gainOrDamage;
	//5 or 3
	private char fiveOrThreePrime;
	//E, I, or S
	private char annotation;
	//for snp, this is the interbase coor for the splice site; for INDEL then it is the variant VCF position
	private int position;
	
	private double referenceScore;
	private String referenceSequence;
	private double alternateScore;
	private String alternateSequence;
	
	/**-10Log10(pvalue)*/
	private double transPValue;
	
	public SpliceJunction(char gainOrDamage, char fiveOrThreePrime, char annotation){
		this.gainOrDamage = gainOrDamage ;
		this.fiveOrThreePrime = fiveOrThreePrime;
		this.annotation = annotation;
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append(getType()); sb.append("\n");
		sb.append(position); sb.append("\tPosition\n");
		sb.append(referenceScore); sb.append("\t"); sb.append(referenceSequence); sb.append("\tRef\n");
		sb.append(alternateScore); sb.append("\t"); sb.append(alternateSequence); sb.append("\tAlt\n");
		sb.append(transPValue); sb.append("\t-10Log10(PValue)\n");
		return sb.toString();
	}
	
	//getters and setters
	public char getFiveOrThreePrime() {
		return fiveOrThreePrime;
	}
	public void setFiveOrThreePrime(char fiveOrThreePrime) {
		this.fiveOrThreePrime = fiveOrThreePrime;
	}
	public int getPosition() {
		return position;
	}
	public void setPosition(int position) {
		this.position = position;
	}
	public String getType() {
		if (type == null) type = new String (gainOrDamage+""+ fiveOrThreePrime +""+ annotation);
		return type;
	}
	public double getReferenceScore() {
		return referenceScore;
	}
	public void setReferenceScore(double referenceScore) {
		this.referenceScore = referenceScore;
	}
	public String getReferenceSequence() {
		return referenceSequence;
	}
	public void setReferenceSequence(String referenceSequence) {
		this.referenceSequence = new String (referenceSequence);
	}
	public double getAlternateScore() {
		return alternateScore;
	}
	public void setAlternateScore(double alternateScore) {
		this.alternateScore = alternateScore;
	}
	public String getAlternateSequence() {
		return alternateSequence;
	}
	public void setAlternateSequence(String alternateSequence) {
		this.alternateSequence = new String(alternateSequence);
	}
	public double getTransPValue() {
		return transPValue;
	}
	public void setTransPValue(double transPValue) {
		this.transPValue = transPValue;
	}
	public char getGainOrDamage() {
		return gainOrDamage;
	}
	public void setGainOrDamage(char gainOrDamage) {
		this.gainOrDamage = gainOrDamage;
	}
	
}
