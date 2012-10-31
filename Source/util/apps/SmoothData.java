package util.apps;

import java.io.File;
import util.gen.*;

/**Pulls the first two columns of double from a white space delimited file.
 * Bins the first and tracks the second.  
 * Prints the bin labels and the mean and median of the second scores that fell within the bin.
 * Arg[0] is the matrix file, arg[1] the number of bins.*/
public class SmoothData {
	
	public static void main (String[] args){
		//Any arguments?
		if (args.length != 2) Misc.printExit("\nEnter the file text of a tab delimited text file containing " +
				"two columns of numbers and the number of bins to divide the first column into. A histogram will be" +
				" made using the first column and used to calculate the mean and median, for each bin, of the " +
				"associated values from the second column.\n");
		
		//load matrix
		double[][] matrix = Num.loadDoubleMatrix(new File (args[0]));
		
		//pull first two columns and min and max
		double[] values = new double[matrix.length];
		double[] scores = new double[matrix.length];
		double min = matrix[0][0];
		double max = matrix[0][0];
		for (int i=0; i< matrix.length; i++){
			values[i] = matrix[i][0];
			scores[i] = matrix[i][1];
			//System.out.println(values[i]+"\t"+scores[i]);
			if (values[i]< min) min = values[i];
			else if (values[i]> max) max = values[i];
		}
		
		//how many bins?
		int numBins = Integer.parseInt(args[1]);
		
		Histogram h = new Histogram (min, max+max*0.01, numBins);
		h.countAll(values, scores);
		h.printBins();
	}
}
