package trans.misc;
import java.io.*;
import trans.graphics.*;
import trans.main.*;
import util.gen.*;

/**
 * Converts binding regions to intervals setting the binding region as the best window and adds 2kb on either side.
 */
public class Regions2Intervals {
	private static int sizeOfOligo = 25;
	
	public static void main(String[] args) {
		if (args.length == 0){
			System.out.println("\nEnter the full path region text file text (tab delimited: " +
					"chrom start stop score notes) and a bp offset to subtract from the start and add to the stop. Set to zero for no offset.\n");
			System.exit(0);
		}
		int offset = Integer.parseInt(args[1]);

		
		//fetch Binding Regions
		File bindingRegionFile = new File(args[0]);
		GenomicRegion[] brs = RankedSetAnalysis.parseGenomicRegionsFile(bindingRegionFile);
		
		//make array of Intervals
		int num = brs.length;
		Interval[] intervals = new Interval[num];
		Window window;
		for (int i=0; i<num; i++){
			//String chromosome, int start, int stop, int N, float score, double[] aveIntValues
			window = new Window(brs[i].getChromosome(), brs[i].getStart()- offset, brs[i].getStop()+offset, 0, new double[]{brs[i].getScore()});
			intervals[i] = new Interval(window, sizeOfOligo);
			intervals[i].setStart1stOligo(intervals[i].getStart1stOligo());
			intervals[i].setStartLastOligo(intervals[i].getStartLastOligo());
		}
		
		//save Intervals
		IO.saveObject(new File(IO.getFullPathName(bindingRegionFile)+"Int"), intervals);
		System.out.println("\nSaved "+num+" Intervals");
		
		
	}
}