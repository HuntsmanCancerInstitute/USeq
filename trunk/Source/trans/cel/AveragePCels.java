package trans.cel;
import util.gen.*;
import java.io.*;

/**Averages serialized float[]s, like xxx.celp files, from the CelProcessor app.*/
public class AveragePCels {
	public static void main(String[] args){
		if (args.length ==0) Misc.printExit("\nEnter the directory text containing serialized float[]s to average\n");
		
		System.out.println("\nLaunching...");
		File directory = new File (args[0]);
		
		File average = new File (directory, "ave.celp");
		
		File[] floatFiles = IO.extractFiles(directory);
		
		float[][] floats = new float[floatFiles.length][];
		for (int i=0; i< floats.length; i++){
			floats[i] = (float[]) IO.fetchObject(floatFiles[i]);
		}
		
		float[] ave = Num.averageFloatArraysFlippedToFloat(floats);
		
		IO.saveObject(average, ave);
		System.out.println("Done!");
	}
}
