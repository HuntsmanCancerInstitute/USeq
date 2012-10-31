package trans.tpmap;
import java.io.*;
import java.util.*;
import util.bio.calc.*;
import util.gen.*;

/**
 * Calculates simple stats on a tpmap file, number unique oligos, number of lines, gc, tm...
 */
public class StatTPMap {

	public static void main(String[] args){
		File bpmap = new File(args[0]);
		try {
			String line;
			String[] tokens;

			BufferedReader in = new BufferedReader(new FileReader(bpmap));
			HashSet uniqueOligos = new HashSet(3000000);
			HashSet chromosomes = new HashSet();
			//run thru tpmap
			int counter =0;
			
			while ((line = in.readLine()) !=null) {           
				tokens = line.split("\\s+");
				if (tokens.length < 7) continue;
				chromosomes.add(tokens[2]);
				uniqueOligos.add(tokens[0]);
				counter++;
			}
			
			//calc gc and tm
			float[] gc = new float[uniqueOligos.size()];
			float[] tm = new float[gc.length];
			int index = 0;
			Iterator it = uniqueOligos.iterator();
			while (it.hasNext()){
				String seq = (String)it.next();
				double tmD = NucleicAcid.calcDefaultNearestNeighborTm(seq);
				tm[index] = new Double(tmD).floatValue();
				double gcD = NucleicAcid.calculateFractionGC(seq);
				gc[index] = new Double(gcD).floatValue();
				index++;
			}
			
			/*
			//to print out line
			while ((line = in.readLine()) !=null) { 
				line = line.trim();
				tokens = line.split("\\s+");
				if (tokens.length < 5)continue;
				if(tokens[2].startsWith("chr")) System.out.println(line);
			}
			*/
			in.close();
			
			System.out.println("Number unique oligos "+uniqueOligos.size());
			System.out.println("Number tpmap lines   "+counter);
			System.out.println("Diff Chromosomes: "+chromosomes);
			System.out.println("\nStats on GC:");
			Num.statFloatArray(gc, false);
			Histogram h = new Histogram (0,0.8,20);
			h.countAll(gc);
			h.printScaledHistogram();
			System.out.println("\nStats on Nearest Neighbor TM:");
			Num.statFloatArray(tm, false);
			h = new Histogram (54,86,20);
			h.countAll(tm);
			h.printScaledHistogram();
		} catch (Exception e) {e.printStackTrace();}
		
	}
	
	

	
	
}
