package util.apps;
import java.io.*;
import util.gen.*;

public class ParserForGff {

	public static void main (String[] args){
		File toParse = new File ("/Users/nix/Desktop/repeatMasker.bed");
		String line;
		String[] tokens;

		try{
			BufferedReader in = new BufferedReader (new FileReader (toParse));
			int counter =0;
			while ((line = in.readLine()) !=null){
				tokens = line.split("\\t+");
				StringBuilder sb = new StringBuilder(tokens[4]);
				
				for (int i=5; i< tokens.length; i++) {
					sb.append("_");
					sb.append(tokens[i]);
				}
				counter++;
				//0	     1          2   3    4
				//chr1	8388322	8388809	+	Dr000331
				System.out.println(tokens[0]+
						"\tRepeatMasker\tRepeatMasker\t"+
						tokens[1]+"\t"+tokens[2]+"\t.\t"+
						tokens[3]+"\t.\t"+sb+"_"+counter);
				
			}
			
			in.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

}
