package util.apps;
import java.io.*;

import util.gen.*;

/**
 * For parsing a text file into a tab delimited file, working template.
 */
public class TextParser {
	
	
	public static void main(String[] args) {
		try{
			//read in file
			File file = new File(args[0]);
			String[] lines = IO.loadFileIntoStringArray(file);
			
			//parse and print
			PrintWriter out = new PrintWriter(new FileWriter(file.getCanonicalPath()+".prsd"));
			double[] x = {0,1,2};
			for (int i=0; i< lines.length; i++){
				String[] tokens = lines[i].split("\\s+");
				//print second, first, and last
				if (tokens.length == 3){
					//out.println("ori\t"+tokens[0]+"\t"+tokens[1]+"\t"+tokens[2]);
					double[] y = {
							Double.parseDouble(tokens[0]),
							Double.parseDouble(tokens[1]),
							Double.parseDouble(tokens[2]),
					};
					//make converter
					LinearRegression lr = new LinearRegression (x,y);
					//print converted
					out.println("fix\t"+lr.calculateX(y[0])+"\t"+
							lr.calculateX(y[1])+"\t"+
							lr.calculateX(y[2]));
					
				}
			}
			out.close();
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
}
