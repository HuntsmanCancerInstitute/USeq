package trans.misc;
import util.gen.*;
import java.io.*;
import java.util.*;

public class CorrelateGrFiles {

	public static void main(String[] args) {
		if (args.length == 0){
			Misc.printExit("\nEnter a directory containing xxx.gr files to correlate");
		}
		
		//get files
		File[] grFiles = IO.extractFiles(new File(args[0]), ".gr");
		Arrays.sort(grFiles);
		for (int i=0; i< grFiles.length; i++){
			float[] first = loadGrFile(grFiles[i]);
			for (int j=i+1; j< grFiles.length; j++){
				float[] second = loadGrFile(grFiles[j]);
				double cc = PearsonCorrelation.correlationCoefficient(first, second);
				System.out.println(cc+"\t"+grFiles[i].getName()+"\t"+grFiles[j].getName());
			}
		}
	}

	/**Reads in a gr file returning it's float values.*/
	public static float[] loadGrFile(File grFile){
		ArrayList floats = new ArrayList(100000);
		try{
			String line;
			BufferedReader in = new BufferedReader(new FileReader(grFile));
			while ((line=in.readLine())!=null){
				if (line.trim().length()==0) continue;
				String[] tokens = line.split("\\s+");
				Float x = new Float(tokens[1]);
				if (x.isNaN()) Misc.printExit("\nNaN found! "+grFile);
				floats.add(x);
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
		return Num.arrayListOfFloatToArray(floats);
	}
}
