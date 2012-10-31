package trans.misc;

import java.io.*;
import util.gen.*;

public class PrependColumn {
	
	/**
	 * Prepends a column on columns of data.  Useful for converting .gr to .sgr.
	 * Args[0] = full path text to file, [1]= chrom text, [2]= type interval or leave blank
	 */
	public static void main(String[] args) {
		try{
			File file = new File (args[0]);
			BufferedReader in = new BufferedReader(new FileReader(file));
			PrintWriter out = new PrintWriter(new FileWriter(file.getCanonicalPath()+".sgr"));
			String line;
			String[] tokens;
			if (args.length == 2){
				while ((line = in.readLine()) !=null){
					out.println(args[1]+"\t"+line);
				}
			}
			else {
				while ((line = in.readLine()) !=null){
					tokens = line.split("\\s+"); 
					out.println(args[1]+"\t" +tokens[0]+ "\t100");
					out.println(args[1]+"\t" +tokens[1]+ "\t100");
				}
			}
			in.close();
			out.close();
		}catch (Exception e){
			e.printStackTrace();
		}
	}
	
}
