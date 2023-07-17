package util.gen;

import org.apache.commons.math3.distribution.NormalDistribution;

public class ConIntChatGPT {

    public static void main(String[] args) {
        // Example contingency table values
        int a = 10;
        int b = 20;
        int c = 5;
        int d = 15;
        
        //Where zeros cause problems with computation of the odds ratio or its standard error, 0.5 is added to all cells

        // Calculate the odds ratio
        double oddsRatio = calculateOddsRatio(a, b, c, d);

        // Calculate the standard error
        double standardError = calculateStandardError(a, b, c, d);

        // Calculate the z-value for a 95% confidence interval (two-tailed test), 0.975 isn't a mistake
        double zValue = 1.96; new NormalDistribution().inverseCumulativeProbability(0.975);


        // Calculate the lower and upper bounds of the confidence interval
        double lowerBound = Math.exp(Math.log(oddsRatio) - (zValue * standardError));
        double upperBound = Math.exp(Math.log(oddsRatio) + (zValue * standardError));

        // Print the results
        System.out.println("Odds Ratio: " + oddsRatio);
        System.out.println("95% Confidence Interval: [" + lowerBound + ", " + upperBound + "]");
    }

    private static double calculateOddsRatio(double a, double b, double c, double d) {
        return (a * d) / (b * c);
    }

    private static double calculateStandardError(double a, double b, double c, double d) {
        double se = 1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d;
        return Math.sqrt(se);
    }
}