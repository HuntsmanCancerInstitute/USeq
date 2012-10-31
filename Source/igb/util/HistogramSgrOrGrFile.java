package igb.util;
import util.gen.*;
import java.io.*;

import trans.misc.Util;
import trans.roc.*;

public class HistogramSgrOrGrFile {
	public static void main (String[] args){
		if ( args.length==0 ){
			Misc.printExit("\nEnter the full path text for an sgr file (tab delimited: chrom, position, score).\n" +
					"Alternatively give full path text and text of chromosome to load a gr file (tab delimited: position, score).\n");
		}
		Sgr[] lines;
		System.out.println("Loading "+args[0]);
		if (args.length==1) lines = trans.misc.Util.loadSgrFile(new File (args[0]));
		else lines = trans.misc.Util.loadGrFile(new File (args[0]), args[1]);
		//extract log2ratio values
		double[] scores = new double[lines.length];
		for (int i=0; i<lines.length; i++){
			scores[i] = lines[i].getScore();
		}
		//System.out.println("Converting to z-scores");
		//scores = Num.convertToZScores(scores);
		//print histogram (double[] scores, int numberOfBins, double tailTrimFraction)
		Histogram h = new Histogram (scores, 500, 0);
		h.printScaledHistogram();
		h.printCenteredForwardAndReverseCounts();
		//write file of z-scores
		File zScoreFile = new File (args[0]+"Zs");
		try {
			PrintWriter out = new PrintWriter( new FileWriter (zScoreFile));
			for (int i=0; i< scores.length; i++){
				out.println(scores[i]);
				out.println("\n"); 
			}
			out.close();
		} catch (IOException e){
			e.printStackTrace();
		}
	}

}
