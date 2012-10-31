package trans.cel;
import java.io.*;
import java.util.*;
import trans.tpmap.IntensityFeature;
import util.gen.*;

public class StatIntensityFeatureArray {
	
	public static void main (String[] args){
		
		//for each file
		File[] ifFiles = IO.extractFiles(new File(args[0]));
		Arrays.sort(ifFiles);
		float[][] values = new float[ifFiles.length][];
		for (int x=0; x< ifFiles.length; x++){
			
			//get object array
			File ifFile = ifFiles[x];	
			IntensityFeature[] ifs = (IntensityFeature[])IO.fetchObject(ifFile);
			
			//extract intensity values
			float[] intensities = new float[ifs.length];
			for (int i=0; i<ifs.length; i++){
				intensities[i] = ifs[i].intensity;
			}
			values[x] = intensities;
			
			//stat float array
			System.out.println("\n"+ifFiles[x].getName());
			Num.statFloatArray(intensities, false);
			//Arrays.sort(values[x]);
			//System.out.println(ifFiles[x].getName()+"\t"+values[x].length+"\t"+Num.trimmedMean(values[x], 0.1));
		}
		
		//print all pair cc
		System.out.println("\nAll Pair CC:");
		for (int x=0; x< values.length; x++){
			for (int y=x+1; y<values.length; y++){
				System.out.println(ifFiles[x].getName()+" vs "+ifFiles[y].getName()+"\t "+PearsonCorrelation.correlationCoefficient(values[x], values[y]));
			}
		}
		
		//print each text
		if (false){
			System.out.println("\nIntensity dump");
			for (int x=0; x< ifFiles.length; x++){
				System.out.print(ifFiles[x].getName()+"\t");
			}
			System.out.println();
			
			//print values
			for (int x=0; x<values[0].length; x++){
				for (int y=0; y<values.length; y++){
					System.out.print(values[y][x]+"\t");
				}
				System.out.println();
			}
		}
		
	}
	
}
