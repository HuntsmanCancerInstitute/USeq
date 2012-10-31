package util.bio.calc;
import java.text.*;
import java.util.*;


/**
 * Methods associated with alignments, pvalues, blast.
 *
 */
public class Alignment {
	
	/**Calculates the standard normal cumulative distribution given a z-ratio.
	 * AKA The area under the curve in the tail.
	 * Returns 0 or 1 if z>6 or Z<-6, respectively.
	 * Taken from http://finance.bi.no/~bernt/gcc_prog/recipes/recipes/node23.html */
	public static double estimatePValue(double z){
		if (z>6) return 0;
		if (z<-6) return 1;
		
		double b1 = 0.31938153;
		double b2 = -0.356563782;
		double b3 = 1.781477937;
		double b4 = -1.821255978;
		double b5 = 1.330274429;
		double p = 0.2316419;
		double c2 = 0.3989423;

		if (z<0) z *=-1;
		double t = 1.0/(1.0+z*p);
		double b = c2*Math.exp((-z)*(z/2.0));
		double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
		n= 1.0-b*n;
		return 1.0-n;	
	}
	
	/**Calculates a -10Log10(number) transformation*/
	public static double transform10Log10(double num){
		return -10.0*(Math.log(num)/Math.log(10.0));
	}

	/**Converts a raw S score to a bit score given lambda (in nats) and K*/
	public static double convertRawScoreToBit(
		double lambdaNats,
		double K,
		double rawS) {
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(1);
		double bs = (lambdaNats * rawS - Math.log(K)) / Math.log(2);
		String bsS = f.format(bs);
		bsS = bsS.replaceAll(",", "");
		return Double.parseDouble(bsS);
	}
	
	/**Converts a nat bit score to raw S.*/
	public static double convertBitScoreToRaw(
		double lambdaNats,
		double K,
		double bitS) {
		NumberFormat f = NumberFormat.getNumberInstance();
		f.setMaximumFractionDigits(1);
		double bs = Math.round((bitS * Math.log(2) + Math.log(K)) / lambdaNats);
		String bsS = f.format(bs);
		bsS = bsS.replaceAll(",", "");
		return Double.parseDouble(bsS);
	}
	
	/**Calculates an E value given a variety of alignment parameters*/
	public static String calculateEValue(
		double lambdaNats,
		double K,
		int rawS,
		double effectiveN,
		double effectiveM) {
		double num =
			K * effectiveN
				* effectiveM
				* Math.exp(-1 * lambdaNats * (double) rawS);
		DecimalFormat f = new DecimalFormat("0E0");
		String numS = f.format(num);
		if (numS.equals("0E0"))
			return "<1E-99";
		return numS;

	}
	
	/**Calculates effective m and n using BLAST parameters, returning int[m', n']*/
	public static int[] calculateEffectiveMandN (double m, double n, double K, double H){
		double l= (Math.log(K*m*n))/H;
		return new int[]{(int)Math.round(m-l), (int)Math.round(n-l)};
	}
	
	/**An array of Lambda, K, H, match, misMatch parameters derived from BLASTN.
	 * Use the fetchBLASTParams(int match, int misMatch) method to get at these.*/
	final static String[][] blastParams =
		{
			{
				" 1.10 0.333 0.549 1 -1 ",  //lambda K H for match misMatch
				" 0.264 0.0532 0.0722 2 -1 ",
				" ",
				" ",
				" ",
				" ",
				" ",
				" ",
				" ",
				" " },
			{
			" 1.33 0.621 1.12 1 -2 ",
				" 0.549 0.333 0.549 2 -2 ",
				" 0.271 0.130 0.222 3 -2 ",
				" 0.132 0.0532 0.0722 4 -2 ",
				" 0.0514 0.0250 0.0135 5 -2 ",
				" ",
				" ",
				" ",
				" ",
				" " },
				{
			" 1.37 0.711 1.31 1 -3 ",
				" 0.634 0.408 0.912 2 -3 ",
				" 0.366 0.333 0.549 3 -3 ",
				" 0.227 0.161 0.306 4 -3 ",
				" 0.144 0.0952 0.158 5 -3 ",
				" 0.0882 0.0532 0.0722 6 -3 ",
				" 0.0494 0.0291 0.0264 7 -3 ",
				" 0.0212 0.0268 0.00550 8 -3 ",
				" ",
				" " },
				{
			" 1.38 0.738 1.36 1 -4 ",
				" 0.666 0.621 1.12 2 -4 ",
				" 0.410 0.340 0.811 3 -4 ",
				" 0.275 0.333 0.549 4 -4 ",
				" 0.192 0.176 0.357 5 -4 ",
				" 0.136 1.63 0.222 6 -4 ",
				" 0.0957 0.0809 0.132 7 -4 ",
				" 0.0661 0.0532 0.0722 8 -4 ",
				" 0.0435 0.0327 0.0350 9 -4 ",
				" 0.0257 7.93 0.0135 10 -4 " },
				{
			" 1.39 0.747 1.38 1 -5 ",
				" 0.681 0.515 1.24 2 -5 ",
				" 0.432 0.396 0.997 3 -5 ",
				" 0.301 0.306 0.753 4 -5 ",
				" 0.220 0.333 0.549 5 -5 ",
				" 0.164 0.184 0.390 6 -5 ",
				" 0.124 0.140 0.270 7 -5 ",
				" 0.0945 0.103 0.182 8 -5 ",
				" 0.0713 0.0733 0.118 9 -5 ",
				" 0.0529 0.0532 0.0722 10 -5 " },
				{
			" 1.39 0.749 1.38 1 -6 ",
				" 0.687 0.711 1.31 2 -6 ",
				" 0.444 0.621 1.12 3 -6 ",
				" 0.317 1.17 0.912 4 -6 ",
				" 0.237 0.286 0.716 5 -6 ",
				" 0.183 0.333 0.549 6 -6 ",
				" 0.144 0.189 0.414 7 -6 ",
				" 0.114 1.41 0.306 8 -6 ",
				" 0.0904 1.63 0.222 9 -6 ",
				" 0.0718 1.79 0.158 10 -6 " },
				{
			" 1.39 0.750 1.39 1 -7 ",
				" 0.690 0.548 1.34 2 -7 ",
				" 0.451 0.455 1.21 3 -7 ",
				" 0.327 0.386 1.03 4 -7 ",
				" 0.249 0.327 0.854 5 -7 ",
				" 0.196 0.273 0.690 6 -7 ",
				" 0.157 0.333 0.549 7 -7 ",
				" 0.127 0.192 0.431 8 -7 ",
				" 0.104 0.161 0.334 9 -7 ",
				" 0.0855 0.132 0.256 10 -7 " },
				{
			" 1.39 0.750 1.39 1 -8 ",
				" 0.692 0.738 1.36 2 -8 ",
				" 0.455 0.472 1.27 3 -8 ",
				" 0.333 0.621 1.12 4 -8 ",
				" 0.257 0.357 0.965 5 -8 ",
				" 0.205 1.10 0.811 6 -8 ",
				" 0.167 0.263 0.671 7 -8 ",
				" 0.137 0.333 0.549 8 -8 ",
				" 0.114 0.194 0.445 9 -8 ",
				" 0.0958 1.31 0.357 10 -8 " },
				{
			" 1.39 0.750 1.39 1 -9 ",
				" 0.692 0.558 1.37 2 -9 ",
				" 0.458 0.711 1.31 3 -9 ",
				" 0.337 0.427 1.19 4 -9 ",
				" 0.263 0.378 1.05 5 -9 ",
				" 0.211 1.17 0.912 6 -9 ",
				" 0.174 0.295 0.779 7 -9 ",
				" 0.145 0.255 0.657 8 -9 ",
				" 0.122 0.333 0.549 9 -9 ",
				" 0.104 0.196 0.456 10 -9 " },
				{
			" 1.39 0.750 1.39 1 -10 ",
				" 0.693 0.747 1.38 2 -10 ",
				" 0.460 0.491 1.33 3 -10 ",
				" 0.340 1.06 1.24 4 -10 ",
				" 0.267 0.621 1.12 5 -10 ",
				" 0.216 1.03 0.997 6 -10 ",
				" 0.179 0.321 0.871 7 -10 ",
				" 0.151 1.07 0.753 8 -10 ",
				" 0.128 0.250 0.646 9 -10 ",
				" 0.110 0.333 0.549 10 -10 " }
	};

	/**Returns Lambda, K, H or null for a given match and misMatch from BLASTN.
	Some combinations of match/misMatch are not computed by BLASTN!
	*/
	public static double[] fetchBLASTParams(int match, int misMatch) {
		String x = (blastParams[(misMatch * -1) - 1][match - 1]).trim();
		if (x.equals(""))
			return null;
		String[] nums = x.split("\\s");
		return new double[] {
			Double.parseDouble(nums[0]),
			Double.parseDouble(nums[1]),
			Double.parseDouble(nums[2])};
	}
	
	/**Will score a DNA alignment given the two sequences, and scoring parameters.*/
	public static int ScoreDNAAlignment(
		String seq1,
		String seq2,
		int match,
		int misMatch,
		int gapCreate,
		int gapExtend) {
		/*scores a DNA alignment, returns the score */
		Character dash = new Character('-');
		int score = 0;
		boolean newGap = true;
		int matches = 0;
		int len = seq1.length();
		seq1 = seq1.toUpperCase();
		seq2 = seq2.toUpperCase();
		for (int i = 0; i < len; i++) {
			Character s1 = new Character(seq1.charAt(i));
			Character s2 = new Character(seq2.charAt(i));
			int comp = s1.compareTo(s2);
			//check to see if they are the same
			if (comp == 0) {
				score += match;
				newGap = true;
				matches++;
			}
			//if they are not the same check to see if it is a gap
			else {
				int comp1 = s1.compareTo(dash);
				int comp2 = s2.compareTo(dash);
				if (comp1 == 0 || comp2 == 0) { //there is a gap
					//check to see if this is a new gap
					if (newGap) {
						score += gapCreate;
						newGap = false;
					} else
						score += gapExtend;
				} else { //it's a mismatch
					score += misMatch;
					newGap = true;
				}
			}
		}
		return score;
	}
}
