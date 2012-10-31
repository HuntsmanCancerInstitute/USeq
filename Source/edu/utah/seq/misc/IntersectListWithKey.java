  package edu.utah.seq.misc;

import java.io.*;
import java.util.*;
import util.gen.*;

public class IntersectListWithKey {

	//fields
	HashSet key;
	String[] picks;
	ArrayList<String> subPicks;
	double numInKey;
	double priorFDR = -1;
	PrintWriter out;

	public IntersectListWithKey(String[] args) {
		//parse files
		key = IO.loadFileIntoHashSet(new File (args[0]));
		File[] pickFiles = IO.extractFiles(new File(args[1]));

		for (int x=0; x< pickFiles.length; x++){
			priorFDR = -1;
			picks = IO.loadFile(pickFiles[x]);
			numInKey = key.size();
			File results = new File(Misc.removeExtension(pickFiles[x].toString())+"_Int.xls");
			try {
				out = new PrintWriter( new FileWriter(results));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println("# in key "+key.size());
			System.out.println("# in picks "+picks.length);
			out.println("NumSubPicks\tNumHits\tFDR\tRachetFDR\tTPR");

			//for each pick intersect with key and calculate FDRs
			subPicks = new ArrayList<String>();
			for (int i=0; i< picks.length; i++){
				subPicks.add(picks[i]);
				intersect();
			}
			out.close();
			printIntersection();
		}
	}

	public void intersect(){
		double numSubPicks = subPicks.size();
		double numHits = 0;
		for (int i=0; i< numSubPicks; i++){
			if (key.contains(subPicks.get(i))) numHits++;
		}
		double fdr = (numSubPicks - numHits) / numSubPicks;
		if (fdr > priorFDR) priorFDR = fdr;
		double tpr = numHits/numInKey;
		out.println (numSubPicks + "\t" + numHits +"\t"+fdr+"\t"+priorFDR+"\t"+tpr);
	}

	public void printIntersection(){
		System.out.println("\nHits to key");
		for (int i=0; i< picks.length; i++){
			if (key.contains(subPicks.get(i))) System.out.println("Y\t"+picks[i]);
			else System.out.println("N\t"+picks[i]);
		}

	}

	public static void main(String[] args) {
		new IntersectListWithKey(args);
	}

}
