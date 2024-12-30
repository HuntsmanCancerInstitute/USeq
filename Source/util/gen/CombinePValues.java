package util.gen;

import java.util.HashMap;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

/**Combines p-values using Fisher's method. https://en.wikipedia.org/wiki/Fisher%27s_method*/
public class CombinePValues {

	//fields
	HashMap<Integer, ChiSquaredDistribution> dfChiDist = new HashMap<Integer, ChiSquaredDistribution>();

	public double calculateCombinePValues(double[] pvalues) {
		
		double sumNatLgs = 0;
		for (double p: pvalues) sumNatLgs+= Math.log(p);
		double chiSquare = sumNatLgs * -2;
		Integer df = pvalues.length * 2;
		
		ChiSquaredDistribution dist = dfChiDist.get(df);
		if (dist == null) {
			dist = new ChiSquaredDistribution(df);
			dfChiDist.put(df, dist);
		}
		
		return 1.0 - dist.cumulativeProbability(chiSquare);
	}
	
	public static void main(String[] args) {
		double[] pvals = {0.11, 0.12, 0.21, 0.08};
		CombinePValues cp = new CombinePValues();
		IO.pl(cp.calculateCombinePValues(pvals)); //0.03195265009643111
		
		pvals = new double[]{0.06277, 0.020672779};
		IO.pl(cp.calculateCombinePValues(pvals)); //0.009923258826133319
		
		pvals = new double[] {0.891393610977685,0.0234808596059804,4.72855153537908E-5,0.00461861158332152,0.531980303255415,0.00216895843142386,0.0932441525314288,8.24650130277569E-4,6.28071788668598E-4,1.82162997940538E-4,2.22043351387029E-5,2.02021979044664E-4,0.0180229511139275,0.020700987940993,0.0455229853323915,0.00154371013521027,0.947783589224271,0.210981720127937,0.979155031324106,0.00571173100657494,0.734798124036884,9.60790034439448E-4,6.15992769580889E-11};
		IO.pl(cp.calculateCombinePValues(pvals));
	}

}
