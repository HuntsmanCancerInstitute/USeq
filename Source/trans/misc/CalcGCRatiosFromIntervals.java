package trans.misc;
import java.util.*;

import java.io.*;

import trans.main.*;
import util.bio.calc.*;
import util.gen.*;

/**
 * Parse out treatment and control oligo info from intervals, calculate GC and ratios to look for bias.
 */
public class CalcGCRatiosFromIntervals {
	
	
	public static void main(String[] args) {
		Interval[] intervals = (Interval[])IO.fetchObject(new File(args[0]));
		
		//for each interval
		for (int i=0; i<intervals.length; i++){
			System.out.println("\n"+ intervals[i].getChromosome()+":" +intervals[i].getStart1stOligo());
			Oligo[] oligos = intervals[i].getOligos();
			int numOligos = oligos.length;
			ArrayList al = new ArrayList();
			
			//load oligo gc and ratio info into GCRatio[]
			for (int x=0; x< numOligos; x++){
				double gc = NucleicAcid.calculateFractionGC(oligos[x].getSequence());
				double aveT = Num.averageFloatArray(oligos[x].getTreatmentIntensities(intervals[i].getNumberTreatmentIntensities()));
				double aveC = Num.averageFloatArray(oligos[x].getControlIntensities(intervals[i].getNumberControlIntensities()));
				double ratio = aveT/aveC;
				al.add(new GCRatio(ratio, (int)Math.round(gc*100)));
			}
			GCRatio[] gcRatios = new GCRatio[al.size()];
			al.toArray(gcRatios);
			
			//for each different GC calculate an average and a standard deviation
			Arrays.sort(gcRatios);
			GCRatio current = gcRatios[0];
			ArrayList ratios = new ArrayList();
			
			for (int j=1; j< gcRatios.length; j++){		
				//last one?
				if (j == (gcRatios.length-1) ){
					if (current.getGc() == gcRatios[j].getGc()){
						ratios.add(new Double(current.getRatio()));
						ratios.add(new Double(gcRatios[j].getRatio()));
						printCalculations(ratios, current.getGc());
					}
					else {
						ratios.add(new Double(current.getRatio()));
						printCalculations(ratios, current.getGc());
						ratios.clear();
						ratios.add(new Double(gcRatios[j].getRatio()));
						printCalculations(ratios, gcRatios[j].getGc());
					}
				}
				//if same gc add
				else if (current.getGc() == gcRatios[j].getGc()){
					ratios.add(new Double(current.getRatio()));
					current = gcRatios[j];
				}
				//print and reset
				else{
					ratios.add(new Double(current.getRatio()));
					printCalculations(ratios, current.getGc());
					//reset
					ratios.clear();
					current = gcRatios[j];
				}
			}
			
		}
	}
	/**Prints %GC, Mean, 3xStdErr*/
	public static void printCalculations(ArrayList doubleAL, int gc){
		double[] d = Num.arrayListOfDoubleToArray(doubleAL);
		double mean = Num.mean(d);
		double stdev = 0;
		if (d.length > 1) stdev = 3* Num.standardError(d, mean);
		System.out.println(gc+"\t"+mean+"\t"+stdev);
	}
	
	
}
