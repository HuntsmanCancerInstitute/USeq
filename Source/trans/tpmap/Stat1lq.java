package trans.tpmap;
import java.util.*;
import java.io.*;

import util.gen.*;
public class Stat1lq {
	
	public static void main (String[] args){
		if (args.length !=1) Misc.printExit("\nEnter the full path for a 1lq file.\n");
		
		File text1lqFile = new File (args[0]);
		String line;
		String[] tokens;
		HashMap hash =new HashMap();
		try{
			//make files
			BufferedReader in = new BufferedReader(new FileReader(text1lqFile));
			//skip header
			while (true){
				line = in.readLine();
				if (line.startsWith("X")) break;
			}
			
			while ((line = in.readLine()) !=null) {
				//X       Y       Seq     DESTYPE
				//0       1        2       3
				tokens = line.split("\\s+");
				Object obj = hash.get(tokens[3]);
				if (obj != null){
					Integer i = (Integer)obj;
					int num = i.intValue() + 1;
					hash.put(tokens[3],new Integer(num));
				}
				else hash.put(tokens[3], new Integer(1));
			}	
			in.close();
		}catch(Exception e){
			e.printStackTrace();
		}
		System.out.println("DESTYPES with counts: "+hash);
		
		
		
	}
}
