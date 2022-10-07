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
	}

}
