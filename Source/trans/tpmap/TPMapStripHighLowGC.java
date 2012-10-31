package trans.tpmap;
import java.io.*;
import java.util.*;

import util.bio.calc.*;

/**
 * Strips a tpmap file of lines where the %GC is > 0.6 or < 0.4.
 */
public class TPMapStripHighLowGC {
	//fields
	private static double high = 0.6;
	private static double low = 0.4;
	
	public static void main(String[] args) {
		File bpmap = new File(args[0]);
		File strippedBpmap = new File (args[0]+"strp");
		try {
			BufferedReader in  = new BufferedReader(new FileReader(bpmap));
			PrintWriter out = new PrintWriter(new FileWriter(strippedBpmap));
			String line;
			int pass =0;
			int fail =0;
			while ((line = in.readLine()) !=null){
				String[] tokens = line.split("\\s+");
				if (tokens.length != 8){
					System.out.println("Bad line! "+line);
					System.exit(1);
				}
				double gc = NucleicAcid.calculateFractionGC(tokens[0]);
				if (gc <= high && gc >= low) {
					out.println(line);
					pass++;
				}
				else fail++;
			}
			
			out.close();
			in.close();
			System.out.println("\nPass: "+pass+"\nFail: "+fail);
		} catch (Exception e){
			e.printStackTrace();
		}
	}

}
