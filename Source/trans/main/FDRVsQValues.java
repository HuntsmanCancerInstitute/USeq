package trans.main;
import util.gen.*;
import java.io.*;
import java.util.*;

public class FDRVsQValues {
	public static void main (String[] args){
		Window[] windows = (Window[])IO.fetchObject(new File(args[0]));
		int numWindows = windows.length;

		//extract log2ratio values
		double[] scores = new double[numWindows];
		for (int i=0; i<numWindows; i++){
			scores[i] = windows[i].getScores()[1];
		}
		
		
		
		//convert to z scores
		double[] zScores = scores;//Num.convertToZScores(scores);
		
		//write to file and histogram
		Histogram h = new Histogram (-2,2,100);
			for (int i=0; i<numWindows; i++){
				
				h.count(zScores[i]);
			}
		
		//print histogram
		System.out.println("\nHistogram of z scores scaled at 0.01\n");
		
		
	}

}
