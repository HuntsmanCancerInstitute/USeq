package util.gen;

//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

/**
 * This does a Fisher Exact test.  The Fisher's Exact test procedure calculates an exact probability value
 * for the relationship between two dichotomous variables, as found in a two by two crosstable. The program
 * calculates the difference between the data observed and the data expected, considering the given marginal
 * and the assumptions of the model of independence. It works in exactly the same way as the Chi-square test
 * for independence; however, the Chi-square gives only an estimate of the true probability value, an estimate
 * which might not be very accurate if the marginal is very uneven or if there is a small value (less than five)
 * in one of the cells.
 *
 * It uses an array of factorials initialized at the beginning to provide speed.
 * There could be better ways to do this.
 *
 * @author Ed Buckler
 * @version $Id: FisherExact.java,v 1
 */

public class FisherExact {
    private double[] f;
    int maxSize;


    /**
     * constructor for FisherExact table
     *
     * @param maxSize is the maximum sum that will be encountered by the table (a+b+c+d)
     */
    public FisherExact(int maxSize) {
        this.maxSize = maxSize;
        f = new double[maxSize + 1];
        f[0] = 0.0;
        for (int i = 1; i <= this.maxSize; i++) {
            f[i] = f[i - 1] + Math.log(i);
        }
    }

    /**
     * calculates the P-value for this specific state
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return the P-value
     */
    public final double getP(int a, int b, int c, int d) {
        int n = a + b + c + d;
        if (n > maxSize) {
            return Double.NaN;
        }
        double p;
        p = (f[a + b] + f[c + d] + f[a + c] + f[b + d]) - (f[a] + f[b] + f[c] + f[d] + f[n]);
        return Math.exp(p);
    }

    /**
     * Calculates the one-tail P-value for the Fisher Exact test.  Determines whether to calculate the right- or left-
     * tail, thereby always returning the smallest p-value.
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return one-tailed P-value (right or left, whichever is smallest)
     */
    public final double getCumlativeP(int a, int b, int c, int d) {
        int min, i;
        int n = a + b + c + d;
        if (n > maxSize) {
            return Double.NaN;
        }
        double p = 0;
        p += getP(a, b, c, d);
        if ((a * d) >= (b * c)) {
            min = (c < b) ? c : b;
            for (i = 0; i < min; i++) {
                p += getP(++a, --b, --c, ++d);
            }
        }
        if ((a * d) < (b * c)) {
            min = (a < d) ? a : d;
            for (i = 0; i < min; i++) {
                double pTemp = getP(--a, ++b, ++c, --d);
                p += pTemp;
            }
        }
        return p;
    }

    /**
     * Calculates the right-tail P-value for the Fisher Exact test.
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return one-tailed P-value (right-tail)
     */
    public final double getRightTailedP(int a, int b, int c, int d) {
        int min, i;
        int n = a + b + c + d;
        if (n > maxSize) {
            return Double.NaN;
        }
        double p = 0;

        p += getP(a, b, c, d);
        min = (c < b) ? c : b;
        for (i = 0; i < min; i++) {
            p += getP(++a, --b, --c, ++d);

        }
        return p;
    }

    /**
     * Calculates the left-tail P-value for the Fisher Exact test.
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return one-tailed P-value (left-tail)
     */
    public final double getLeftTailedP(int a, int b, int c, int d) {
        int min, i;
        int n = a + b + c + d;
        if (n > maxSize) {
            return Double.NaN;
        }
        double p = 0;

        p += getP(a, b, c, d);
        min = (a < d) ? a : d;
        for (i = 0; i < min; i++) {
            double pTemp = getP(--a, ++b, ++c, --d);
            p += pTemp;
        }


        return p;
    }


    /**
     *   Calculates the two-tailed P-value for the Fisher Exact test.
     *
     *   In order for a table under consideration to have its p-value included
     *   in the final result, it must have a p-value less than the original table's P-value, i.e.
     *   Fisher's exact test computes the probability, given the observed marginal
     *   frequencies, of obtaining exactly the frequencies observed and any configuration more extreme.
     *   By "more extreme," we mean any configuration (given observed marginals) with a smaller probability of
     *   occurrence in the same direction (one-tailed) or in both directions (two-tailed).
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return two-tailed P-value or NaN if the table sum exceeds the maxSize
     */
    public final double getTwoTailedP(int a, int b, int c, int d) {
        int min, i;
        int n = a + b + c + d;
        if (n > maxSize) {
            return Double.NaN;
        }
        double p = 0;

        double baseP = getP(a, b, c, d);

        int initialA = a, initialB = b, initialC = c, initialD = d;
        p += baseP;
        min = (c < b) ? c : b;
        for (i = 0; i < min; i++) {
            double tempP = getP(++a, --b, --c, ++d);
            if (tempP <= baseP) {
                p += tempP;
            }
        }

        // reset the values to their original so we can repeat this process for the other side
        a = initialA;
        b = initialB;
        c = initialC;
        d = initialD;

        min = (a < d) ? a : d;
        for (i = 0; i < min; i++) {
            double pTemp = getP(--a, ++b, ++c, --d);
            if (pTemp <= baseP) {
                p += pTemp;
            }
        }
        return p;
    }
    
    /**Returns OR, lower, upper*/
    public static double[] getOddsRatioAnd95thConfidenceInterval(double a, double b, double c, double d) {
    	
    	//Where zeros cause problems with computation of the odds ratio or its standard error, 0.5 is added to all cells (a, b, c, d) (Pagano & Gauvreau, 2000; Deeks & Higgins, 2010).
    	if (a==0.0 || b==0.0 || c==0.0 || d==0.0) {
    		a+=0.5;
    		b+=0.5;
    		c+=0.5;
    		d+=0.5;
    	}

    	double oddsRatio = (a * d) / (b * c);
    	
    	double inner = 1.96 * Math.sqrt(1.0/a + 1.0/b + 1.0/c + 1.0/d);
    	double lnOR = Math.log(oddsRatio);

    	//Upper 95% CI = e ^ [ln(OR) + 1.96 sqrt(1/a + 1/b + 1/c + 1/d)] 
    	double upper = Math.pow(Math.E, (lnOR + inner));
    			
    	//Lower 95% CI = e ^ [ln(OR) - 1.96 sqrt(1/a + 1/b + 1/c + 1/d)] 
    	double lower = Math.pow(Math.E, (lnOR - inner));
    	
    	return new double[] {oddsRatio, lower, upper};
    }
    
    /**Measure of how far from independence the 2x2 table is. 1= independent. Ratio of ratios. Will return Infinity or 0 if cells are zero.*/
    public static double getOddsRatio(double a, double b, double c, double d) {
    	return ((a * d) / (b * c));
    }


    public static void main(String[] args) {

        int[][] argInts = new int[3][4];
        argInts[0] = new int[]{1, 2, 1, 10};
        argInts[1] = new int[]{10, 20, 10, 100};
        argInts[2] = new int[]{8, 2, 4, 6};
        
        //argInts[0] = new int[]{17, 83, 1, 99};
        /*
        argInts[1] = new int[]{3011, 30, 2020, 5};
        argInts[2] = new int[]{3011, 30, 2020, 60};
        argInts[3] = new int[]{1, 2, 0, 3};
        argInts[4] = new int[]{3, 1, 1, 3};
        argInts[5] = new int[]{1, 3, 3, 1};
        argInts[6] = new int[]{0, 1, 1, 0};
        argInts[7] = new int[]{1, 0, 0, 1};
        argInts[8] = new int[]{11, 0, 0, 6};
        argInts[9] = new int[]{10, 1, 1, 5};
        argInts[10] = new int[]{5, 6, 6, 0};
        argInts[11] = new int[]{9, 2, 2, 4};
        argInts[12] = new int[]{6, 5, 5, 1};
        argInts[13] = new int[]{8, 3, 3, 3};
        argInts[14] = new int[]{7, 4, 4, 2};*/

        FisherExact fe = new FisherExact(200);

        for (int i = 0; i < argInts.length; i++) {
            System.out.println("\na=" + argInts[i][0] + " b=" + argInts[i][1] + " c=" + argInts[i][2] + " d=" + argInts[i][3]);
            
           // double cumulativeP = fe.getCumlativeP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
            //System.out.println("\ncumulativeP = " + cumulativeP);
            
            double twoTailedP = fe.getTwoTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
            System.out.println("twoTailedP = " + twoTailedP);

            double leftTailedP = fe.getLeftTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
            System.out.println("leftTailedP = " + leftTailedP);

            double rightTailedP = fe.getRightTailedP(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
            System.out.println("rightTailedP = " + rightTailedP);
       
            
            double[] orLowerUpper = fe.getOddsRatioAnd95thConfidenceInterval(argInts[i][0], argInts[i][1], argInts[i][2], argInts[i][3]);
            System.out.println("OR = " + orLowerUpper[0]);
            System.out.println("Lower = " + orLowerUpper[1]);
            System.out.println("Upper = " + orLowerUpper[2]);
        }
    }
}
