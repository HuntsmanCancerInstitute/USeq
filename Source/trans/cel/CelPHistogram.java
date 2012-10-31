package trans.cel;
import util.gen.*;
import java.io.*;


public class CelPHistogram {
	
	
	
	public static void main (String[] args){
		File[] celps = IO.extractFiles(new File(args[0]), "celp");
		
		for (int i=0; i< celps.length; i++){
			//load array
			System.out.println("\n"+celps[i]);
			float[] values = (float[]) IO.fetchObject(celps[i]);
			Num.statFloatArray(values, false);
			
			//print histogram (double[] scores, int numberOfBins, double tailTrimFraction)
			//Histogram h = new Histogram (scores, 100, 0.001);
			//public Histogram( double minimum, double maximum, int numberOfBins){
			Histogram h = new Histogram (0, 150, 40);
			h.countAll(values);
			h.printScaledHistogram();
			//h.printCenteredForwardAndReverseCounts();
		}
		
	}
	
}
