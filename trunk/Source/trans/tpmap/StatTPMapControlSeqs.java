package trans.tpmap;
import java.io.*;
import java.util.*;
import util.gen.*;

public class StatTPMapControlSeqs {
		
		public static void main(String[] args) {
			File bpmap = new File(args[0]);
			try {
				BufferedReader in  = new BufferedReader(new FileReader(bpmap));
				String line;
				HashMap hash = new HashMap();
				HashMap clones = new HashMap();
				while ((line = in.readLine()) !=null){
					//TCTAGGGTTTGGGTGTGTCTGTCTA       t       AT3G12600.1     5       2444    649     2444    650     1
					String[] tokens = line.split("\\s+");
					if (tokens.length != 9){
						System.out.println("Bad line! "+line);
						continue;
					}
					String combo = tokens[0]+"_"+tokens[2];
					if (hash.containsKey(combo)){
						int num = ((Integer)hash.get(combo)).intValue() + 1;
						hash.put(combo, new Integer(num));
					}
					else {
						hash.put(combo, new Integer(1));
					}
					if (clones.containsKey(tokens[2])){
						int num = ((Integer)clones.get(tokens[2])).intValue() + 1;
						clones.put(tokens[2], new Integer(num));
					}
					else {
						clones.put(tokens[2], new Integer(1));
					}
					

				}
				System.out.println(hash);
				System.out.println();
				System.out.println(clones);
				in.close();
			} catch (Exception e){
				e.printStackTrace();
			}
		}

	}
