package util.apps;

import java.io.File;

import util.gen.Histogram;
import util.gen.IO;

public class HistogramData {


	public static void main(String[] args) {
		String[] headings = {"GV","MI","MII","PN","CL","MOR","ICM","TROPH"};
		//read in file, double[column][row]
		double[][] colRow = IO.loadTableOfDoubles(new File (args[0]));
		//for each column make a histogram and print
		//for each column
		for (int i=0; i< colRow.length; i++){
			System.out.println("\nCol "+headings[i]);
			//Histogram( double minimum, double maximum, int numberOfBins)
			Histogram h = new Histogram(0,30,50);
			//add row data
			for (int j=0; j<colRow[i].length; j++){
				h.count(colRow[i][j]);
			}
			h.printScaledHistogram();
			
		}
	}

}
