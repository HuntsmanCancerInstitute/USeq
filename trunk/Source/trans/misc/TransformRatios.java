package trans.misc;

import java.io.File;

import trans.main.Interval;
import trans.main.Oligo;
import util.bio.calc.NucleicAcid;
import util.gen.IO;
import util.gen.Num;


/**
 * Transforms ratio scores in Intervals based on the BAC GC Corrector.
 * */
public class TransformRatios {

	public static void main(String[] args) {
		Interval[] intervals = (Interval[])IO.fetchObject(new File(args[0]));
		double[][] mb = makeTransformer();
		
		//for each interval
		for (int i=0; i<intervals.length; i++){
			System.out.println("\n"+ intervals[i].getChromosome()+":" +intervals[i].getStart1stOligo());
			Oligo[] oligos = intervals[i].getOligos();
			int numOligos = oligos.length;
			double[] ratios = new double[numOligos];
			double[] tnsRatios = new double[numOligos];
			for (int x=0; x< numOligos; x++){
				int gc = (int)Math.round(100* NucleicAcid.calculateFractionGC(oligos[x].getSequence()));
				double aveT = Num.averageFloatArray(oligos[x].getTreatmentIntensities(intervals[i].getNumberTreatmentIntensities()));
				double aveC = Num.averageFloatArray(oligos[x].getControlIntensities(intervals[i].getNumberControlIntensities()));
				double ratio = aveT/aveC;
				//y= mX + b
				double trans = mb[gc][0]* ratio + mb[gc][1];
				//System.out.println(intervals[i].getChromosome()+"\t" +oligos[x].getStart()+"\t"+(ratio+trans));
				ratios[x] = ratio;
				tnsRatios[x] = ratio + trans;
			}
			double stnd = Num.standardDeviation(ratios);
			double mean = Num.mean(ratios);
			System.out.println(mean+" "+stnd+" "+(100*stnd/mean));
			stnd = Num.standardDeviation(tnsRatios);
			mean = Num.mean(tnsRatios);
			System.out.println(mean+" "+stnd+" "+(100*stnd/mean));
		}

	}
	
	public static double[][] makeTransformer(){
		double[][] d = new double[100][2];
		//load double array with non transformers m, b
		for (int i=0; i<d.length; i++){
			d[i] = new double[]{0,0};
		}
		//load actual transformers
		d[28]  = new double[] {1.3974, -1.3108};
		d[32]  = new double[] {0.5793, -0.5267};
		d[36]  = new double[] {0.2209, -0.1696};
		d[40]  = new double[] {0.0298, 0.0047};
		d[44]  = new double[] {-0.0068, 0.0082};
		d[48]  = new double[] {-0.0041, -0.0014};
		d[52]  = new double[] {0.0103, -0.0042};
		d[56]  = new double[] {0.1088, -0.0964};
		d[60]  = new double[] {0.2265, -0.1904};
		d[64]  = new double[] {0.4389, -0.3245};
		d[68]  = new double[] {0.8019, -0.7117};
		d[72]  = new double[] {0.9631, -0.4619};
		return d;
	}

}
