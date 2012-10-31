package util.bio.seq;

import util.gen.*;
import util.bio.parsers.*;
import java.io.*;
import java.util.*;
import java.util.regex.*;

/**Merges many fasta's into one fasta with a defined number of X spacers.*/
public class ConcatinateFastas {

	//fields
	private File[] fastaFiles;
	private File directory;
	private int numberXs = 1000;
	private String spacer = "";
	private Fasta[] fastas;
	private String chromName = "chrConcat";
	private File merged;
	private File bed;
	private File bounds;
	private File shifter;

	public ConcatinateFastas (String[] args){
		processArgs(args);
		
		//make spacer
		char[] x = new char[numberXs];
		Arrays.fill(x, 'N');
		spacer = new String(x);

		//load fastas and sort by last number in text
		fastas = new MultiFastaParser(fastaFiles).getFastas();
		Arrays.sort(fastas);

		//print merged fasta to file and stats to screen
		//also print a bed file of regions
		try {
			PrintWriter out = new PrintWriter (new FileWriter ( merged));
			PrintWriter outBed = new PrintWriter (new FileWriter ( bed));
			PrintWriter outBounds = new PrintWriter (new FileWriter ( bounds));
			PrintWriter outShifter = new PrintWriter (new FileWriter ( shifter));
			int cumLength = 0;
			//print first seq
			out.println(">"+chromName);
			out.println(fastas[0].getSeq());
			//print to screen
			outShifter.println(fastas[0].getName()+"\t0");
			//print bed
			int newCumLength = fastas[0].getSeq().length();
			outBed.println(chromName+"\t"+cumLength+"\t"+newCumLength+"\t"+fastas[0].getName());
			cumLength = newCumLength+ numberXs;
			//print bounds bed
			outBounds.println(chromName+"\t"+newCumLength+"\t"+cumLength);
			//for each subsequent fasta
			for (int i=1; i< fastas.length; i++){
				//print to file
				out.println(spacer);
				out.println(fastas[i].getSeq());
				//print to screen
				outShifter.println(fastas[i].getName()+"\t"+cumLength);
				//print to bed and set cumLength
				newCumLength = fastas[i].getSeq().length() + cumLength;
				outBed.println(chromName+"\t"+cumLength+"\t"+newCumLength +"\t"+fastas[i].getName());
				cumLength = newCumLength + numberXs;
				//print bounds bed
				outBounds.println(chromName+"\t"+newCumLength+"\t"+cumLength);
			}
			out.close();
			outBed.close();
			outBounds.close();
			outShifter.close();
		} catch (Exception e){
			e.printStackTrace();
		}

	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ConcatinateFastas(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': directory = new File(args[++i]); break;
					case 'f': fastaFiles = IO.extractFiles(new File(args[++i])); break;
					case 'n': numberXs = Integer.parseInt(args[++i]); break;
					case 'c': chromName = args[++i]; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//look for fasta files
		if (fastaFiles == null || fastaFiles.length ==0) Misc.printErrAndExit("\nError: cannot find any fasta sequence files?\n");
		
		//make files
		directory.mkdirs();
		merged = new File(directory,chromName+".fasta");
		bed = new File(directory,chromName+".bed");
		bounds = new File(directory,chromName+"_bounds.bed");
		shifter = new File(directory,chromName+".shifter.txt");
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                              Concatinate Fastas: Oct 2010                        **\n" +
				"**************************************************************************************\n" +
				"Concatinates a directory of fasta files into a single sequence seperated by a defined\n" +
				"number of Ns.  Outputs the merged fasta as well as bed files for the junctions and\n" +
				"spacers as well as a file to be used to shift UCSC gene table annotations. Use this\n" +
				"app to create artificial chromosomes for poorly assembled genomes. \n\n" +

				"Options:\n"+
				"-d Full path directory for saving the results.\n"+
				"-f Full path directory containing fasta files to concatinate.\n"+
				"-n Number of Ns to use as a spacer, defaults to 1000.\n"+
				"-c Name to give the concatinate, defaults to chrConcat .\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/ConcatinateFastas -n 2000 -d\n" +
				"    /zv8/MergedNA_Scaffolds -f /zv8/BadFastas/ -c chrNA_Scaffold \n\n" +

		"**************************************************************************************\n");

	}
}
