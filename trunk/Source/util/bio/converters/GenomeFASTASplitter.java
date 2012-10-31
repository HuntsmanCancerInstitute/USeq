package util.bio.converters;
import java.io.*;

import util.bio.parsers.*;
import util.gen.*;

/**
 * Splits a whole genome fasta into different files, also rewrites the fasta header to just the first word following the >
 */
public class GenomeFASTASplitter {
	
	public static void main(String[] args) {
		if (args.length ==0) {
			System.out.println("Enter the text of a multi FASTA genome file.");
			System.exit(0);
		}
		
		File genomeFile = new File(args[0]);
		
		// read in seqs to memory
		MultiFastaParser p = new MultiFastaParser(genomeFile);
		String[] seqs = p.getSeqs();
		String[] names = p.getNames();
		int num = seqs.length;
		
		//print headers
		System.out.println("FASTA headers:");
		for (int i=0; i< num; i++){
			System.out.println("\t"+names[i]);
		}
		
		//write individual files
		String firstWord;
		for (int i=0; i< num; i++){
			firstWord = names[i].split("\\s+")[0];
			File f = new File(genomeFile.getParentFile(), firstWord+".fasta");
			System.out.println ("Writing "+f);
			IO.writeString(">"+firstWord+"\n"+seqs[i], f);
		}
		
		//write composite
		try{
			PrintWriter out = new PrintWriter(new FileWriter(new File(genomeFile.getParentFile(),"stripped.fasta")));
			System.out.print("Writing stripped.fasta file");
			for (int i=0; i< num; i++){
				System.out.print(".");
				firstWord = names[i].split("\\s+")[0];
				out.println(">"+firstWord+"\n"+seqs[i]);
			}
			out.close();
		} catch (IOException e){
			e.printStackTrace();
		}
		
		System.out.println("\nDone!\n");
		
	}
}
