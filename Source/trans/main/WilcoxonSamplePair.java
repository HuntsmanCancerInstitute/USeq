package trans.main;

/**
 * Used by the {@link WilcoxonSignedRankTest}.
 */
public class WilcoxonSamplePair implements Comparable{

	//fields
	private float tValue;
	private float cValue;
	private float difference;
	private float absoluteDifference;
	private float rank;

	public WilcoxonSamplePair (float tValue, float cValue){
		this.tValue= tValue;
		this.cValue= cValue;
		difference = tValue-cValue;
		if (difference<0) absoluteDifference = -1*difference;
		else absoluteDifference = difference;
	}
	
	public WilcoxonSamplePair (float absoluteDifference){
		this.absoluteDifference = absoluteDifference;
	}
	
	public int compareTo(Object obj){
		WilcoxonSamplePair other = (WilcoxonSamplePair)obj;
		if (other.absoluteDifference> absoluteDifference) return -1;
		if (other.absoluteDifference<absoluteDifference) return 1;
		return 0;
	}

	public String toString(){
		StringBuffer sb = new StringBuffer("tValue: ");
		sb.append(tValue);
		sb.append(" \tcValue: ");
		sb.append(cValue);
		sb.append(" \tdiff: ");
		sb.append(difference);
		sb.append(" \tabsDiff: ");
		sb.append(absoluteDifference);
		sb.append(" \trank: ");
		sb.append(rank);
		return sb.toString();
	}
	public float getAbsoluteDifference() {
		return absoluteDifference;
	}
	public float getCValue() {
		return cValue;
	}
	public float getDifference() {
		return difference;
	}
	public float getRank() {
		return rank;
	}
	public float getTValue() {
		return tValue;
	}
	public void setAbsoluteDifference(int absoluteDifference) {
		this.absoluteDifference = absoluteDifference;
	}
	public void setRank(float rank) {
		this.rank = rank;
	}
}
