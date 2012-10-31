package util.gen;

public class RankedFloatArray {
	
	boolean tiesFound = false;
	float[] ranks;
	
	public RankedFloatArray (boolean tiesFound, float[] ranks){
		this.tiesFound = tiesFound;
		this.ranks = ranks;
	}

	public float[] getRanks() {
		return ranks;
	}

	public boolean isTiesFound() {
		return tiesFound;
	}
}
