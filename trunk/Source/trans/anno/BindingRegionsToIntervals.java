package trans.anno;
import java.io.*;

import trans.main.*;
import util.gen.*;

/**
 * Converts binding regions to intervals setting the binding region as the best window and adds 2kb on either side.
 */
public class BindingRegionsToIntervals {
	private static final int sizeOfOligo= 25;
	
	public static void main(String[] args) {
		if (args.length == 0){
			System.out.println("\nEnter the full path binding region text file text (tab delimited: " +
					"chrom start stop seq).\n");
			System.exit(0);
		}
		
		//fetch Binding Regions
		File bindingRegionFile = new File(args[0]);
		BindingRegion[] brs = IntersectRegions.parseRegionsFile(bindingRegionFile);
		
		//make array of Intervals
		int num = brs.length;
		Interval[] intervals = new Interval[num];
		Window window;
		for (int i=0; i<num; i++){
			//String chromosome, int start, int stop, int N, float score, double[] aveIntValues
			window = new Window(brs[i].getChromosome(), brs[i].getStart(), brs[i].getEnd(), 0, new double[]{0});
			intervals[i] = new Interval(window,sizeOfOligo);
			intervals[i].setStart1stOligo(intervals[i].getStart1stOligo());
			intervals[i].setStartLastOligo(intervals[i].getStartLastOligo()-intervals[i].getSizeOfOligoMinusOne());
		}
		
		//save Intervals
		IO.saveObject(new File(IO.getFullPathName(bindingRegionFile)+"Int"), intervals);
		System.out.println("\nSaved "+num+" Intervals");
		
		
	}
}
