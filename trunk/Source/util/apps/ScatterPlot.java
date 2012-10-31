package util.apps;
import javax.swing.*;


/**
 * Draws a simple scatter plot and calculates a Pearson correlation coefficient 
 * for two serialized int[] or float[] arrays.
 * 
 * @author nix
 *
 */
public class ScatterPlot {

	public static void main(String[] args) {
		if (args.length < 2){
			printDocs();
			System.exit(0);
		}
		
		ScatterDrawFrame frame = new ScatterDrawFrame(args);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.show();
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Scatter Plot: Aug 2005                                **\n" +
				"**************************************************************************************\n" +
				"To draw a simple scatter plot and calculates a Pearson correlation coefficient for two\n" +
				"serialized int[] or float[] arrays, enter full path names on the command line. If you\n" +
				"wish to skip zero values in the analysis, type 'skip' after the second file text.\n\n" +

				"Example: java -Xmx750M -jar pathTo/T2/Apps/ScatterPlot /my/int/array1 /my/int/array2 \n\n" +
				
				"**************************************************************************************\n");
	}	
}	
	


