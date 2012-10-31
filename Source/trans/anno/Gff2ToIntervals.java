package trans.anno;
import java.io.*;

import trans.main.*;
import util.bio.parsers.gff.*;
import util.gen.*;

/**
 * Converts a gff2 file to an interval object file.
 */
public class Gff2ToIntervals {
	static String chromNamePrefix = "chr";
	private static final int sizeOfOligo= 25;
	
	public static void main(String[] args) {
		//load up PCRMs
		File gffFile = new File(args[0]);
		GffParser anno = new GffParser(gffFile, 0, 100000000);
		GffFeature[] features = anno.getGffFeatures();
		
		//make Intervals
		int num = features.length;
		Interval[] intervals = new Interval[num];
		Window window;
		for (int i=0; i<num; i++){
			//String chromosome, int start, int stop, int N, float score, double[] aveIntValues
			window = new Window(chromNamePrefix + features[i].getSeqName(), features[i].getStart(), features[i].getEnd(), 0, new double[]{10,10});
			intervals[i] = new Interval(window, sizeOfOligo);
		}
		
		//save Intervals
		IO.saveObject(new File(IO.getFullPathName(gffFile)+"Int"), intervals);
		
		System.out.println("\nSaved "+num+" Intervals");
		
		
	}
}
