package trans.anno;
import java.io.*;

import trans.main.*;
import util.bio.parsers.gff.*;
import util.gen.*;

/**
 * Converts a gff2 file to a binding region file (chr start stop mockSequence).
 */
public class Gff2ToBindingRegions {
	static String chromNamePrefix = "chr";

	public static void main(String[] args) {
		//load up PCRMs
		File gffFile = new File(args[0]);
		GffParser anno = new GffParser(gffFile, 0, 100000000);
		GffFeature[] features = anno.getGffFeatures();
		
		//make Intervals
		int num = features.length;
		StringBuffer sb = new StringBuffer();
		for (int i=0; i<num; i++){
			//rank
			sb.append(10); sb.append("\t");
			//chromosome
			sb.append(chromNamePrefix); sb.append(features[i].getSeqName()); sb.append("\t");
			//start
			sb.append(features[i].getStart()); sb.append("\t");
			//stop
			sb.append(features[i].getEnd()); sb.append("\t");
			//mock sequence
			sb.append("gggggggggggggg\n");
			
		}
		
		IO.writeString(sb.toString(), IO.getFullPathName(gffFile)+".br");
		
		System.out.println("\nSaved "+num+" Binding GenomicRegion lines");
		
		
	}
}
