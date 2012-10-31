package trans.cel;
import util.gen.*;
import java.io.*;

/**Averages two serialized float[][]s, like xxx.cela files, from the CelFileConverter app. Assumes square arrays.*/
public class AverageACels {
	public static void main(String[] args){
		if (args.length!=3) Misc.printExit("\nEnter two full path names for serialized float[]s to average as " +
				"well as the text to give the averaged.\n\tExample: 'java -Xmx512M /my/tech1.cela /my/tech2.cela /my/ave.cela\n");
		
		System.out.println("\nLaunching...");
		//make files
		File f1 = new File(args[0]);
		File f2 = new File(args[1]);
		File fAve = new File(args[2]);
		//fetch floats
		float[][] float1 = (float[][])IO.fetchObject(f1);
		float[][] float2 = (float[][])IO.fetchObject(f2);
		//ave
		float[][] floatAve = Num.mean(float1, float2);
		//save
		IO.saveObject(fAve, floatAve);
		System.out.println("Done!");
	}
}
