package trans.anno;
import java.io.*;
import java.util.*;

import util.gen.*;
/**
 * Performs an intersection analysis on two peak pick files (chr peakPosition score).
 *
 */
public class IntersectBindingPeaks {
	//fields
	private BindingRegion[] one;
	private BindingRegion[] two;
	private Histogram histogram;

	
	public IntersectBindingPeaks(String[] args){
		System.out.println("\nLaunching...");
		
		//fetch two arrays of BindingRegion to represent the binding peaks
		one = AnnotateRegionsWithGeneList.parsePicksFile(new File(args[0]),0);
		two = AnnotateRegionsWithGeneList.parsePicksFile(new File(args[1]),0);
		
		//find closest peaks, set distance in neighborhood
		findClosestPeak(one, two);
		
		//make histogram object, min value = 0, max value = 10kb, 100 bins
		histogram = new Histogram (0, 10000, 100);
		
		//print
		System.out.println("\n"+one.length+" Num peaks in "+args[0]+"\nChrom\tPeak\tScore\tDist To\tChrom\tPeak\tScore\n");
		for (int i=0; i<one.length; i++){
			System.out.println(one[i].simpleSummaryLine()+"\t"+one[i].getNeighborhood()+"\t"+one[i].getClosestBindingRegion().simpleSummaryLine());
			histogram.count(one[i].getNeighborhood());
		}
		System.out.println();
		histogram.printScaledHistogram();
		
		//Random trials
	     HashMap chromLengths = new HashMap(6);
	     chromLengths.put("2L",new Integer(22407834));
	     chromLengths.put("2R",new Integer(20766785));
	     chromLengths.put("3L",new Integer(23771897));
	     chromLengths.put("3R",new Integer(27905053));
	     chromLengths.put("4",new Integer(1281640));
	     chromLengths.put("X",new Integer(22224390));
	     BindingRegion[] random;
	     histogram = new Histogram (0, 10000, 100);
	     int numRandom = one.length;
	     for (int i=0; i<1000; i++){
	     	//make random set of binding regions
	     	random = AnnotateRegions.makeRandomBindingRegions(one, chromLengths, 1);
	     	//map to two
	     	findClosestPeak(random,two);
	     	//add counts to histogram
			for (int j=0; j<numRandom; j++){
				histogram.count(random[j].getNeighborhood());
			}
	     }
	     //print histogram
	     histogram.printScaledHistogram();
	}
	
	public static void findClosestPeak(BindingRegion[] one, BindingRegion[] two){
		int numOne = one.length;
		int numTwo = two.length;
		for (int i=0; i<numOne; i++){
			int gap = 1000000000;
			BindingRegion closest = null;
			for (int j=0; j<numTwo; j++){
				//check chromosome
				if (one[i].getChromosome().equals(two[j].getChromosome())){
					//calc distance apart
					int dist = Math.abs(one[i].getStart()-two[j].getStart());
					//check gap
					if (dist < gap){
						gap = dist;
						closest = two[j];
						if (gap == 0 ) break;
					}
				}
			}
			//set distance and br
			one[i].setClosestBindingRegion(closest);
			one[i].setNeighborhood(gap);
		}
	}
	
	
	
	public static void main(String[] args){
		if (args.length != 2){
			System.out.println("\n Enter full path file names for two peak picks files (chrom, bp peak, score)\n");
			System.exit(0);
		}
		new IntersectBindingPeaks(args);
	}
	
}
