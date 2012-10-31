package trans.misc;
import java.io.*;
import util.bio.annotation.Bed;
import java.util.*;

import edu.utah.seq.useq.data.RegionScoreText;

/**Converts a bed file of say CTCF binding sites into blocks.*/
public class Bed2Blocks {


	public static void main(String[] args) {
		//load bed file
		File bedFile = new File(args[0]);
		HashMap<String,RegionScoreText[]> split = Bed.parseBedFile(bedFile, true);
		
		
		Iterator<String> it = split.keySet().iterator();
		while (it.hasNext()){
			//for each chromosome
			String chrom = it.next();
			RegionScoreText[] bindingSites = split.get(chrom);
			chrom = chrom+"\t";
			
			//convert to center positions
			int[] centerPositions = new int[bindingSites.length];
			for (int i=0; i< bindingSites.length; i++) centerPositions[i] = bindingSites[i].getMiddle();
			
			//write out blocks
			int start = 0;
			for (int i=0; i< bindingSites.length; i++){
				System.out.println(chrom + start + "\t"+ centerPositions[i]);
				start = centerPositions[i];
			}
		}
	}

}
