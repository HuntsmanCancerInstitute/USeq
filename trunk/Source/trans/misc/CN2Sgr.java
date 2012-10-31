package trans.misc;

import java.io.*;
import java.util.regex.Pattern;

import util.gen.*;

/**Converts xxx.cn files to sgr graphs.  These contain snp probes and multiple graph values.
 * */
public class CN2Sgr {

	public static final Pattern TAB = Pattern.compile("\\t");

	public static void main (String[] args){

		try {
			
			File[] cnFiles = IO.extractFiles(new File (args[0]), ".cn");
			
			for (File f: cnFiles){
				
				BufferedReader in = IO.fetchBufferedReader(f);
				String parentName = f.getName();

				//skip header
				in.readLine();
				
				//read in names line
				String line = in.readLine();
				String[] tokens = TAB.split(line);
				int numFiles = tokens.length -3;
				
				//make print writers
				PrintWriter[] pws = new PrintWriter[numFiles];
				for (int i=0; i< pws.length; i++){
					File r = new File (f.getParentFile(), parentName+"_"+tokens[i+3]+".sgr");
					pws[i] = new PrintWriter( new FileWriter(r));
				}
				
				//print em
				while ((line = in.readLine()) != null){
					tokens = TAB.split(line);
					String chrPos = "chr"+tokens[1]+"\t"+tokens[2]+"\t";
					for (int i=0; i< pws.length; i++){
						if (tokens[i+3].equals("null") == false) pws[i].println(chrPos+tokens[i+3]);
					}
				}

				//close em
				in.close();
				for (PrintWriter pw: pws) pw.close();
			}

		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
