package edu.utah.seq.vcf;


public class VCFMatch implements Comparable<VCFMatch>{
	private float score;
	private VCFRecord key;
	private VCFRecord test;

	public VCFMatch (VCFRecord key, VCFRecord test){
		this.key = key;
		this.test = test;
		score = test.getScore();
	}
	
	/**Returns VCFRecord[2][#matches] of key and test*/
	public static VCFRecord[][] split(VCFMatch[] matches){
		VCFRecord[] keys = new VCFRecord[matches.length];
		VCFRecord[] tests = new VCFRecord[matches.length];
		for (int i=0; i< matches.length; i++){
			keys[i] = matches[i].key;
			tests[i] = matches[i].test;
		}
		return new VCFRecord[][]{keys, tests};
	}
	
	/**Sorts by score, smallest to largest*/
	public int compareTo(VCFMatch second) {
		if (this.score < second.score) return -1;
		if (this.score > second.score) return 1;
		return 0;
	}
	public float getScore() {
		return score;
	}
	public void setScore(float score) {
		this.score = score;
	}
	public VCFRecord getKey() {
		return key;
	}
	public void setKey(VCFRecord key) {
		this.key = key;
	}
	public VCFRecord getTest() {
		return test;
	}
	public void setTest(VCFRecord test) {
		this.test = test;
	}
}
