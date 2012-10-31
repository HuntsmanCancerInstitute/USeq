package trans.main;

/**
 * Used by the {@link WilcoxonRankSumTest}.
 */
public class WilcoxonSample implements Comparable{

	//fields
	float value;
	boolean treatment;
	float rank;

	public WilcoxonSample (float value, boolean treatment){
		this.value= value;
		this.treatment = treatment;
	}
	
	public int compareTo(Object obj){
		WilcoxonSample other = (WilcoxonSample)obj;
		if (other.value> value) return -1;
		if (other.value<value) return 1;
		return 0;
	}
}
