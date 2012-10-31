package trans.bkgrnd;
import java.util.*;
import java.io.*;

import util.gen.*;
/**
 * Calculates a correlation coefficient between ratio scores and tm or gc 
 * for a set of CelMapper converted cel/bpmap files.
 * This is a huge memory pig.
 */
public class CalcBkgrndRatios {
	
	//constructor
	public CalcBkgrndRatios(String[] args){
		
			//make array of LoadedOligos
			LoadedOligo[] t1 = fetchLoadedOligos(new File(args[0]));
			int num = t1.length;
			
			//median scale pm values to 50
			medianScale(t1);
			//resort by coordinates
			Arrays.sort(t1);
			
			//make array of LoadedOligos
			LoadedOligo[] t2 = fetchLoadedOligos(new File(args[1]));
			//median scale pm values to 50
			medianScale(t2);
			//resort by coordinates
			Arrays.sort(t2);
			
			//make average
			double[] aveT = new double[num];
			for (int i=0; i<num; i++) aveT[i] = (t1[i].getIntensityPM()+t1[i].getIntensityPM())/2;
			
			t2=null;
			
			//make array of LoadedOligos
			System.out.println("p3...");
			t1 = fetchLoadedOligos(new File(args[2]));
			//median scale pm values to 50
			medianScale(t1);
			//resort by coordinates
			Arrays.sort(t1);
			
			//make array of LoadedOligos
			System.out.println("p4...");
			t2 = fetchLoadedOligos(new File(args[3]));
			//median scale pm values to 50
			medianScale(t2);
			//resort by coordinates
			Arrays.sort(t2);
			
			//make average and get tm and gc
			double[] aveRatios = new double[num];
			double[] tm = new double[num];
			double[] gc = new double[num];
			for (int i=0; i<num; i++) {
				aveRatios[i] = t1[i].getIntensityPM()/t2[i].getIntensityPM();
				tm[i] = t1[i].getTm();
				gc[i] = t1[i].getGc();
			}
			
			t1 =null;
			t2 = null;
	
			//calc corr coefs
			System.out.println(PearsonCorrelation.correlationCoefficient(aveRatios, tm)+" PMRatio:Tm");
			System.out.println(PearsonCorrelation.correlationCoefficient(aveRatios, gc)+" PMRatio:GC");

	}
	
	
	//methods
	public static void medianScale(LoadedOligo[] los){
		int num = los.length;
		double[] med = new double[num];
		for (int i=0; i<num; i++){
			med[i]= los[i].getIntensityPM();
		}
		Arrays.sort(med);
		double scaler = 50.0/Num.median(med);
		for (int i=0; i<num; i++){
			los[i].setIntensityPM(med[i]*scaler);
		}
	}
	
	public static LoadedOligo[] fetchLoadedOligos(File loadedOligoFile){
		ArrayList loAL = new ArrayList();
		try{
			BufferedReader in = new BufferedReader(new FileReader(loadedOligoFile));
			String line;
			LoadedOligo lo;
			while ((line = in.readLine())!= null) {
				lo = new LoadedOligo(line);
				loAL.add(lo);
			}
			in.close();
		} catch (IOException e){
			e.printStackTrace();
		}
		
		LoadedOligo[] los = new LoadedOligo[loAL.size()];
		loAL.toArray(los);
		loAL = null;
		return los;
	}
	
	public static void main(String[] args) {
		if (args.length==0){
			System.out.println("\n Enter a directory containing .celInt files, from CelMapper run using -s flag.\n");
			System.exit(0);
		}
		new CalcBkgrndRatios(args);
	}
}
