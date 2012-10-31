package trans.misc;
import util.gen.*;
import java.io.*;
import trans.main.*;

/**Makes a histogram of scores from a Window[].*/
public class WindowHistogram {
	public static void main (String[] args){
		if ( args.length==0 ){
			Misc.printExit("\nEnter the full path text for the serialized Window[] file and a score index.\n");
		}
		Window[] windows = (Window[])IO.fetchObject(new File(args[0]));
		int scoreIndex = Integer.parseInt(args[1]);
		
		//extract log2ratio values
		double[] scores = new double[windows.length];
		System.out.println("Total number windows "+windows.length);
		for (int i=0; i<windows.length; i++){
			scores[i] = windows[i].getScores()[scoreIndex];
		}

		//write to file and histogram (double[] scores, int numberOfBins, double tailTrimFraction)
		Histogram h = new Histogram (scores, 250, 0);
		//Histogram h = new Histogram (-.5, 13.5, 30);
		//h.countAll(scores);
		h.printScaledHistogram();
		h.printCenteredForwardAndReverseCounts();
	}

}
