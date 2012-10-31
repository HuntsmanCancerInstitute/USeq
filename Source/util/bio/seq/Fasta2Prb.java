package util.bio.seq;
import util.bio.parsers.*;
import java.io.*;
import util.gen.*;
import java.util.*;

public class Fasta2Prb {


	public static void main(String[] args) {
		if (args.length == 0) Misc.printExit("\nProvide a multi xxx.fasta file to convert to a Illumina xxx.prb file, base order is ACGT.\n");
		
		//make hash to convert bases
		String baseConversions = "A,40 -40 -40 -40,C,-40 40 -40 -40,G,-40 -40 40 -40,T,-40 -40 -40 40";
		HashMap<String,String> base2Prb = Misc.createHashMap(baseConversions.split(","));

		//fetch seqs to convert
		File fastaFile = new File (args[0]);
		MultiFastaParser mfp = new MultiFastaParser(fastaFile);
		String[] seqs = mfp.getSeqs();
		
		
		StringBuilder sb = new StringBuilder();
		for (int i=0; i< seqs.length; i++){
			char[] bases = seqs[i].toUpperCase().toCharArray();
			//add first
			sb.append(base2Prb.get(new Character(bases[0]).toString()));
			for (int j=1; j< bases.length; j++){
				sb.append("\t");
				sb.append(base2Prb.get(new Character(bases[j]).toString()));
			}
			sb.append("\n");
		}
		
		System.out.println(sb);

	}

}
