package trans.bkgrnd;
import java.util.*;
import java.io.*;

import util.gen.*;
/**
 * For each CelMapper converted cel/bpmap files, this calculates a correlation coeff between pm:mm, pm:tm, pm:gc.
 */
public class CalculateBackgroundStats {
	
	//constructor
	public CalculateBackgroundStats(String[] args){
		//fetch files
		File[] files = IO.extractFiles(new File(args[0]), ".celInt");
		
		for (int j=0; j<files.length; j++){
			System.out.println("\n"+files[j]);
			//make array of LoadedOligos
			LoadedOligo[] loadedOligos = fetchLoadedOligos(files[j]);
			//make double arrays
			int num = loadedOligos.length;
			double[] pm = new double[num];
			double[] mm = new double[num];
			double[] tm = new double[num];
			double[] gc = new double[num];
			for (int i=0; i<num; i++){
				pm[i] = loadedOligos[i].getIntensityPM();
				mm[i] = pm[i] - loadedOligos[i].getIntensityMM();
				tm[i] = loadedOligos[i].getTm();
				gc[i] = loadedOligos[i].getGc();
			}
			//calc corr coefs
			//PM : MM
			System.out.println(PearsonCorrelation.correlationCoefficient(pm, mm)+" PM:MM");
			//PM : Tm
			System.out.println(PearsonCorrelation.correlationCoefficient(pm, tm)+" PM:Tm");
			//PM : GC
			System.out.println(PearsonCorrelation.correlationCoefficient(pm, gc)+" PM:GC");
		}
	}
	
	
	//methods
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
		new CalculateBackgroundStats(args);
	}
}
