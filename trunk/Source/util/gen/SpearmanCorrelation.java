package util.gen;
import java.util.*;

/**Calculates a Spearman correlation coefficient (rho), in the case of ties, averages the ranks 
 * and returns a Pearson correlation coefficient on the ranks.
 * @author Nix
 */
public class SpearmanCorrelation {
	
	//fields
	private RankSampleIndexComparator indexComp;
	private RankSampleValueComparator valueComp;
	private RankedFloatArray aRFA;
	private RankedFloatArray bRFA;
	
	//constructor
	public SpearmanCorrelation(){
		indexComp = new RankSampleIndexComparator();
		valueComp = new RankSampleValueComparator();
	}
	
	/**Returns a Spearman Rank Correlation Coefficient (rho) between the 
	 * two arrays.  If ties are found, their ranks are averaged and a 
	 * Pearson correlation is calculated on the ranks.*/
	public double spearmanCorrelationCoefficient (float[] a, float[] b){
		aRFA = rank (a);
		bRFA = rank (b);
		if (aRFA.tiesFound || bRFA.tiesFound) return PearsonCorrelation.correlationCoefficient(aRFA.ranks, bRFA.ranks);
		else return corrCoeff (aRFA.ranks, bRFA.ranks);
	}
	
	/**Returns a Spearman Rank Correlation Coefficient (rho) between the 
	 * two arrays.  If ties are found, their ranks are averaged and a 
	 * Pearson correlation is calculated on the ranks.*/
	public double spearmanCorrelationCoefficient (RankedFloatArray aRFA, RankedFloatArray bRFA){
		if (aRFA.tiesFound || bRFA.tiesFound) return PearsonCorrelation.correlationCoefficient(aRFA.ranks, bRFA.ranks);
		else return corrCoeff (aRFA.ranks, bRFA.ranks);
	}
	
	/**This is the final step in calculating a Spearman Correlation coeff.
	 * Don't use it unless you know what you are doing
	 * @param a - a values converted to ranked and sorted by original index
	 * @param b - ditto
	 * Call correlate() if you don't have the ranks.*/
	public double corrCoeff (float[] a, float[] b){
		//calculate sum of square differences
		double sumSqrDiffs = 0;
		for (int i=0; i< a.length; i++) {
			double diff = a[i]-b[i];
			sumSqrDiffs += (diff * diff);
		}
		//numerator 6(sumSqrDiffs)
		double numer = 6* sumSqrDiffs;
		//denominator n(n^2-1)
		double n = a.length;
		double sqrN = n*n;
		//denomenator
		double denom = n * (sqrN-1.0);
		//final rho
		return 1-(numer/denom);
	}
	
	public RankedFloatArray rank(float[] f){
		RankSample[] rs = new RankSample[f.length];
		for (int i=0; i< f.length; i++) {
			rs[i] = new RankSample(i, f[i]);
		}
		//sort by value
		Arrays.sort(rs, valueComp);
		//rank
		return rankSamples(rs);
	}
	
	/**Ranks a sorted array of RankSample based on value.
	 * If no ties are found then this is simply their array index number+1.
	 * (ie 1,2,3,4...)
	 * If ties are encountered, ties are assigned the average of their index
	 * positions+1. (ie if index+1's: 2,3,4 have the same absolute difference, all are assigned a
	 * rank of 3).
	 */
	public RankedFloatArray rankSamples(RankSample[] rs){
		int num = rs.length;
		boolean tiesFound = false;
		//assign ranks as index+1
		for (int i=0; i<num; i++) {
			rs[i].rank=i+1;
		}
		//check for ties
		int start;
		int end;
		for (int i=0; i<num; i++){
			start = i;
			end = i;
			//advance stop until the former and latter don't have the same value and the stop
			//	of the array hasn't been reached
			while (++end < num && rs[start].value==rs[end].value){}
			//check if i was advanced
			if (end-start!=1){// ties found
				tiesFound = true;
				//get average of ranks		
				float ave = Num.getAverageInts((int)rs[start].rank, (int)rs[end-1].rank);
				//assign averages
				for (int x=start; x<end; x++) rs[x].rank = ave;
				//reset i
				i=end-1;
			}		
		}
		
		//sort by original position
		Arrays.sort(rs, indexComp);
		//make float[] of ranks
		float[] ranks = new float[rs.length];
		for (int i=0; i< rs.length; i++) {
			ranks[i] = rs[i].rank;
			//System.out.println("Fin "+rs[i].index+" "+rs[i].value+" "+rs[i].rank+" "+tiesFound);			
		}
		rs = null;
		return new RankedFloatArray(tiesFound, ranks);

	}	
	
	//for testing
	public static void main(String[] args) {
		//float[] a = {106 ,86 ,100 ,101 ,99 ,103 ,97 ,113 ,112 ,110};
		//float[] b = {7,0,27,50,28,29,20,12,6,17};
		float[] a = {1,2,9,11};
		float[] b = {5,3,6,2};
		SpearmanCorrelation sp = new SpearmanCorrelation();
		System.out.println(sp.spearmanCorrelationCoefficient(a, b));
	}

	public RankedFloatArray getARFA() {
		return aRFA;
	}

	public void setARFA(RankedFloatArray arfa) {
		aRFA = arfa;
	}

	public RankedFloatArray getBRFA() {
		return bRFA;
	}

	public void setBRFA(RankedFloatArray brfa) {
		bRFA = brfa;
	}

}
