package edu.expr;
import java.io.*;
import util.gen.*;


public class StripZeros {


	
	public StripZeros (String[] args){
		File file = new File (args[0]);
		File filtered = new File (file.getParentFile(), Misc.removeExtension(file.getName())+"_Stripped.txt");
		int numStripped = 0;
		try {
			BufferedReader in = new BufferedReader ( new FileReader( file));
			PrintWriter out = new PrintWriter (new FileWriter (filtered));
			String line;
			String[] tokens;
			line = in.readLine();
			out.println(line);
			while ((line=in.readLine()) != null){
				line = line.trim();
				if (line.length() ==0) continue;
				tokens = line.split("\\t");
				boolean noZeros = true;
				for (int i=1; i< tokens.length; i++){
					float value = Float.parseFloat(tokens[i]);
					if (value == 0){
						noZeros = false;
						numStripped ++;
						break;
					}
				}
				if (noZeros) out.println(line);
			}
			System.out.println (numStripped +" Number of zero lines stripped.\n");
			out.close();
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	
	
	public static void main(String[] args) {
		new StripZeros(args);

	}

}
