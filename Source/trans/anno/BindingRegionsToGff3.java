package trans.anno;
import java.io.*;

import util.bio.parsers.gff.*;
import util.gen.*;
/**
 * Converts a Binding GenomicRegion file into GFF3 format.
 */
public class BindingRegionsToGff3 {
	static String chromNamePrefix = "";
	
	public static void main(String[] args) {
		if (args.length == 0){
			System.out.println("\nEnter the full path binding region text file text (tab delimited: " +
					"rank subWinMedianScore trimmedMeanPickScore chrom start stop seq).\n");
			System.exit(0);
		}
		//fetch Binding Regions
		File bindingRegionFile = new File(args[0]);
		BindingRegion[] brs = AnnotateRegions.parsePicksFile(bindingRegionFile, 10000);
		
		//make array of GFF3Features
		int num = brs.length;
		Gff3Feature[] gffs = new Gff3Feature[num];
		
		for (int i=0; i<num; i++){
			gffs[i] = new Gff3Feature();
			gffs[i].setSeqId(chromNamePrefix + brs[i].getChromosome());
			gffs[i].setSource("BDTNP");
			gffs[i].setType("CRM");
			gffs[i].setStart(brs[i].getStart());
			gffs[i].setEnd(brs[i].getEnd());
			gffs[i].setScore(brs[i].getRank());
			gffs[i].setId(brs[i].getRank()+"");
		}
		
		//save gff3 file
		Gff3Feature.saveGFF3TextFile(gffs, new File(IO.getFullPathName(bindingRegionFile)+".gff3"));
		
		System.out.println("\nSaved "+num+" GFF3 lines to file");
		
		
	}
}
