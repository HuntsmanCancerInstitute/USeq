package trans.anno;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import trans.main.*;
import util.bio.annotation.*;
import util.bio.parsers.gff.*;
import util.gen.*;


/**
 * Prints genome coordinates for Dmel Rel 4 GFF3 CG names.
 */
public class DmelCGNameCoordinatePrinter {
	
	//fields
	private File gffFile;
	
	public DmelCGNameCoordinatePrinter(String[] args){
		gffFile = new File(args[0]);
		//process gff file generating a CGName: GFF3Feature HashMap
		System.out.println("\tProcessing GFF file...");
		DmelRel4Extractor ex = new DmelRel4Extractor();
		ex.extract(gffFile, true);
		ArrayList geneGrpsAL = ex.getGeneGroupArrayList();
		//make hash
		int num = geneGrpsAL.size();
		HashMap ggs = new HashMap(num);
		for (int i=0; i<num; i++){
			GeneGroup g = (GeneGroup)geneGrpsAL.get(i);
			ggs.put(g.getName(), g);
		}
		//run thru CGNames
		num = args.length;
		ArrayList noMatches = new ArrayList();
		for (int i=1; i<num; i++){
			GeneGroup g = (GeneGroup)ggs.get(args[i]);
			if (g != null) System.out.println(args[i]+"\t"+g.getChromosome()+"\t"+g.getStart()+"\t"+g.getEnd());
			else {
				System.out.println(args[i]+"\tNo Match");
				noMatches.add(args[i]);
			}
		}
		
		if (noMatches.size() !=0) System.out.println("No matches for the following: "+noMatches);
	
		
	}


	

	


	
	public static void main(String[] args) {
		if (args.length ==0){
			System.out.println("\nEnter full path to the DmelRel4.0 gff3 file and then each CGName, space delineated.\n");
			System.exit(0);
		}
		new DmelCGNameCoordinatePrinter (args);
	}
}
