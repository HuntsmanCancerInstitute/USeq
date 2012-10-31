package edu.expr;
import java.io.*;
import java.util.*;
import util.gen.*;

public class ParseWormMart {

	public static void main(String[] args) {
		
		//load worm mart file
		File wormMart = new File (args[1]);
		String[] lines = IO.loadFile(wormMart);
		
		//make HashMap of names from WormMart, skip first line
		HashMap hash = new HashMap();
		String line = null;
		try{
		for (int i=1; i< lines.length; i++){
			line = lines[i];
			//V	-1	9244214	9246162	WBGene00000003 | F07C3.7 | WP:CE32853 | CE32853 | aat-2 | 5J726 | NM_072993
			String[] tokens = lines[i].split("\\t");
			String[] genes = tokens[4].split(" \\| ");
			for (int j=0; j< genes.length; j++){
				hash.put(genes[j], lines[i]);
			}
			
		}
		}catch (Exception e){
			System.err.println(line);
			e.printStackTrace();
		}
		
		//get list of gene names
		File names = new File (args[0]);
		String[] geneNames = IO.loadFileIntoStringArray(names);
		
		//for each gene attempt to find in worm mart hash
		for (int i=0; i< geneNames.length; i++){
			System.out.print(geneNames[i]+"\t");
			if (hash.containsKey(geneNames[i])) System.out.println(hash.get(geneNames[i]));
			else System.out.println();
		}


	}

}
