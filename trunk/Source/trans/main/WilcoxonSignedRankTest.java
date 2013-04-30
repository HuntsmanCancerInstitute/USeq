package trans.main;
import java.util.*;

import util.bio.calc.*;
import util.gen.*;

/**
 *Implementation of a Wilcoxon Signed Rank Test.
 *
 *Many issues to consider: what to do about ties, what to do about zero diffs, what to use as the output 
 *(WPlus, WMinus, WCombine, which z-value (Zw or Zs), a p-value, or Affys -10*LogBase10(p-value) transformation, 
 *lastly how to handle situations where N<10. 
 *
 *The minimal constructor here averages ties, throws out zero diffs, uses the Lowry z-value calculation method 
 *and the Affy p-value transformation (transformedPValueFloat) as the final output. 
 * 
 *Note: this transformedPValueFloat may be 0 indicating no significance (or N<5) (p-value > 0.05), 90 "very" significant (p-value < 1E-11)
 *or something measurable in between.  Affys favored 10-5 p-value cutoff corresponds to a LogTrans= 49.4, with bigger numbers
 *"more significant" although in reality, anything bigger than LogTrans 13.1 is, by WSRT, significant, assuming a normal distribution 
 *and >10 non zero difference sample pairs.
 *
 *Uncomment the main() to test and tweak.
 */

public class WilcoxonSignedRankTest {
	
	//fields
	private WilcoxonSamplePair[] wsps;	//array of sample pairs, may or may not contain those with a difference of zero
	private float WPlus; 				//+ number
	private float WMinus;				//- number
	private float W;					//WPlus+WMinus
	private int numberWPlusRanks;		//+ number
	private int numberWMinusRanks;		//+ number
	private int N; 						//+ number
	private float Zw;					//test statistic based on Richard Lowry, Vasser
	private double[] averageIntensityValues;	//difference, ratio
	private double pValue;
	
	/**-10*LogBase10(p-value) transformation, might be "not significant" or ">90"
	 * Some cutoffs (z-score: p-value: logTrans) 1.65: 0.05: 13.1; 4.26: 1.02E-05: 49.4; 6.00: 9.90E-10: 90.0
	 * Will not be calculated if N<10 (N=numberWPlusRanks+numberWMinusRanks)*/
	private String transformedPValueString;
	private double transformedPValue = 0; //0, 90, or the actual transformation, bigger the better, 49.4 is the Affy cutoff
	
	//constructors
	/**Minimal constructor, fastest way to calculate a transformedPValueFloat*/
	public WilcoxonSignedRankTest(float[] treatment, float[] control){
		//create SamplePairs throwing out zeros
		createSamplePairs(treatment,control, true);
 		test();
	}
	
	public double test(){
		//sort based on absolute difference
		Arrays.sort(wsps);	
		
		//assign ranks paying attention to ties
		rankSamplePairs();		
		
		//calculate signed-rank statistics WPlus, WMinus, numberWPlusRanks, numberWMinusRanks, W
		calculateSignedRankStatistics();
		
		//calculate test statistic Zw ala Lowry
		calculateTestStatisticZw();

		//fetch estimated p-value (actually -10LogBase10(p-value) transformation
		transformedPValue = 0;
		pValue = -1;
		if (N>=10 && Zw>0) {
			pValue = Alignment.estimatePValue(Zw);
			if (pValue==0) transformedPValue = 90;
			else transformedPValue= Alignment.transform10Log10(pValue);	
		}
		else if (N< 10 && N>4) transformedPValue = calcSmallSample10Log10PValue(W,N);
		
		return transformedPValue;
	}
	
	public WilcoxonSignedRankTest(){}
	

	//methods
	
	/**Calculates test statistic Z from Lowry .*/ 
	public void calculateTestStatisticZw(){	
		float N = (float)(numberWPlusRanks + numberWMinusRanks);
		float sd = (float)Math.sqrt(((N*(N+1)) *  ((2*N)+1))/6);
		Zw= (W-0.5f)/sd;
	}
	
	/**Calculates signed rank statistic S = sum of the ranks (aka WPlus) for WilcoxonSamplePairs with a 
	 * positive difference and WMinus for those with a negative difference.*/
	public void calculateSignedRankStatistics(){
		numberWPlusRanks = 0;
		numberWMinusRanks =0;
		WPlus =0;
		WMinus =0; 
		for (int i= wsps.length-1; i>=0; i--){
			if (wsps[i].getDifference()>0) {
				WPlus += wsps[i].getRank();
				numberWPlusRanks++;
			} 
			else if (wsps[i].getDifference()<0){
				numberWMinusRanks++;
				WMinus -= wsps[i].getRank();
			}
		}
		W = Math.abs(WPlus+WMinus);
		N = numberWPlusRanks + numberWMinusRanks;
	}
	
	/**Ranks a sorted array of WilcoxonSamplePairs based on absoluteDifference.
	 * If no ties are found then this is simply their array index number+1.
	 * (ie 1,2,3,4...)
	 * If ties are encountered, ties are assigned the average of their index
	 * positions+1. (ie if index+1's: 2,3,4 have the same absolute difference, all are assigned a
	 * rank of 3).*/
	public void rankSamplePairs(){
		int num = wsps.length;
		
		//assign ranks as index+1
		for (int i=0; i<num; i++)wsps[i].setRank(i+1);

		//check for ties
		int start;
		int end;
		
		for (int i=0; i<num; i++){
			start = i;
			end = i;
			//advance stop until the former and latter don't have the same value and the stop
			//	of the array hasn't been reached
			while (++end < num && wsps[start].getAbsoluteDifference() == wsps[end].getAbsoluteDifference()){}
			
			//check if i was advanced
			if (end-start!=1){// ties found
				//get average of ranks		
				float ave = Num.getAverageInts((int)wsps[start].getRank(), (int)wsps[end-1].getRank());
				//assign averages
				for (int x=start; x<end; x++) wsps[x].setRank(ave);
				//reset i
				i=end-1;
			}		
		}
	}	
	
	/**Creates WilcoxonSamplePair objects from two float[]s.
	 * Indicate if you want to throw out pairs with a difference of zero.
	 * (true for no zeros, false for keep zeros*/
	public void createSamplePairs(float[] treatment, float[] control, boolean noZeros){
		int len = treatment.length;
		if (noZeros){
			ArrayList al = new ArrayList(len);
			for (int i=len-1; i>=0; i--) {
				WilcoxonSamplePair sp = new WilcoxonSamplePair(treatment[i], control[i]);
				if (sp.getDifference()!=0) al.add(sp);
			}
			wsps = new WilcoxonSamplePair[al.size()];
			al.toArray(wsps);
		}
		else {
			wsps = new WilcoxonSamplePair[len];
			for (int i=len-1; i>=0; i--) wsps[i]= new WilcoxonSamplePair(treatment[i], control[i]);
		}
	}
	
	/**Returns signed difference between the scores and a median.*/
	public static float[] signedDifference(float[] scores, float median){
		int len = scores.length;
		float[] sd = new float[len];
		for (int i=len-1; i>=0; i--) sd[i]= scores[i]-median;	
		return sd;		
	} 
	
	/**For calculating a -10Log10(pValue) for small sample sizes (N 5-9).*/
	public double calcSmallSample10Log10PValue(float W, int N){
		pValue = -1;
		if (N==5){
			if (W >= 15) pValue = 0.05;
		}
		else if (N==6){
			if (W>=17 && W<21) pValue = 0.05;
			else if (W>21) pValue = 0.025;
		}
		else if (N==7){
			if (W>=22 && W<24) pValue = 0.05;
			else if (W>=24 && W<28) pValue = 0.025;
			else if (W>28) pValue = 0.01;
		}	
		else if (N==8){
			if (W<26) return 0;
			else if (W>=26 && W<30) pValue = 0.05;
			else if (W>=30 && W<34) pValue = 0.025;
			else if (W>=34 && W<36) pValue = 0.01;
			else if (W>36) pValue = 0.005;
		}
		else if (N==9){
			if (W<29) return 0;
			else if (W>=29 && W<35) pValue = 0.05;
			else if (W>=35 && W<39) pValue = 0.025;
			else if (W>=39 && W<43) pValue = 0.01;
			else if (W>=43 && W<50) pValue = 0.005; //b.s.!
			else if (W>50) pValue = 0.0025;//b.s.!
		}
		if (pValue > 0) return Alignment.transform10Log10(pValue);
		return 0;	
	}	

	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		int num = wsps.length;
		for (int i=0; i<num; i++) {
			sb.append(wsps[i].toString());
			sb.append("\n");
		} 
		return sb.toString();
	}
	public WilcoxonSamplePair[] getWsps() {
		return wsps;
	}
	public void setWsps(WilcoxonSamplePair[] pairs) {
		wsps = pairs;
	}
	public float getZw() {
		return Zw;
	}
	public int getNumberWMinusRanks() {
		return numberWMinusRanks;
	}
	public int getNumberWPlusRanks() {
		return numberWPlusRanks;
	}
	public float getW() {
		return W;
	}
	public float getWMinus() {
		return WMinus;
	}
	public float getWPlus() {
		return WPlus;
	}
	public double getTransformedPValue() {
		return transformedPValue;
	}
	public String getTransformedPValueString() {
		return transformedPValueString;
	}
	
	//main for testing  see http://faculty.vassar.edu/lowry/ch12a.html  
	//	and http://www.nist.gov/speech/tests/sigtests/wilcoxon.htmank+test+example&hl=en
	public static void main(String[] args) {
		System.out.println("\n\nUncomment the main method to get this to run.\n");
		
		WilcoxonSignedRankTest test = new WilcoxonSignedRankTest();
		//transformedPvalue 77.3363772668118
		//float[] c = {68,72,51,29,32,28,32,30,48,41,45,38,31,46,43,35,35,83,37,30,32,45,41,39,33,31,30,42,59,60,46,118,77,59,62,53,57,78,90,67,66,99,79,58,67,60,44,74,89};
		//float[] t = {64,61,49,29,26,23,29,27,34,33,42,28,29,37,40,29,27,78,33,24,29,31,38,34,29,31,29,45,77,38,41,79,52,44,43,35,45,67,68,56,47,52,48,43,40,43,42,48,44};
		float[] t = {0.1027545f, 0.09475391f, 0.05565953f, 0.03977493f, 0.551383f, 0.7156163f, 0.1411698f, 0.1365448f, 0.02192605f, 0.7951f, 0.7040997f, 0.3597389f, 0.1528459f, 0.05561343f, 0.0367802f, 0.06106699f, 0.03326163f, 0.03611488f};
		float[] c = {0.1114038f, 0.09000587f, 0.05665263f, 0.04049822f, 0.5642294f, 0.7080179f, 0.1422857f, 0.1016717f, 0.01717411f, 0.7772832f, 0.6866685f, 0.3185638f, 0.1899087f, 0.05370715f, 0.03483688f, 0.0392044f, 0.0435489f, 0.03652769f};
		test.createSamplePairs(t,c, true);
		double transT = test.test();
		
		System.out.println(test.Zw+" "+test.W+ " test transformedPvalue "+transT+" pValue "+test.getPValue());
	}
	public int getN() {
		return N;
	}
	public double[] getAverageIntensityValues() {
		return averageIntensityValues;
	}
	/**A p-value of -1 is the default insignificant number. Check the sample size, must have > 4 pairs to calculate wilcoxon.*/
	public double getPValue() {
		return pValue;
	}


}
