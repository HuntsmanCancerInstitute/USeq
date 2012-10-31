package util.apps;
import meme.*;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.parsers.*;

/**
 * Scores sequences for the presence of transcription factor binding sites.
 *
 */
public class ScoreSequences {

	//fields
	private File seqFile;
	private File bindingSiteFile;
	private String[] bindingSites;		//aligned and trimmed
	private MotifScanner motifScanner;
	private double cutOff = 0;

	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new ScoreSequences(args);
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Score Sequences: July 2007                             **\n" +
				"**************************************************************************************\n" +
				"SS scores sequences for the presence of transcription factor binding sites. Use the\n" +
				"following options:\n\n" +
				
				"-g The full path FASTA formatted file text for the sequence(s) to scan.\n" +
				"-t Full path file text for the FASTA file containing aligned trimmed examples of\n"+
				"      transcription factor binding sites.  A log likelihood position specific\n"+
				"      probability matrix will be generated from these sequences and used to scan the\n"+
				"      sequences for hits to the matrix.\n"+
				"-s Score cut off for the matrix. Defaults to zero.\n"+
				"\n" +
				"Example: java -Xmx500M -jar pathTo/T2/Apps/ScoreSequences -g /my/affy/DmelSeqs.fasta\n" +
				"      -t /my/affy/zeste.fasta\n\n" +

		"**************************************************************************************\n");
	}
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'g': seqFile = new File(args[i+1]); i++; break;
					case 't': bindingSiteFile = new File(args[i+1]); i++; break;
					case 's': cutOff = Double.parseDouble(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}

		if (seqFile ==null){
			System.out.println("\nEnter a full path file text for the fasta file containing sequence(s) to scan.\n");
			System.exit(0);
		}

		if (bindingSiteFile==null){
			System.out.println("\nEnter a full path file text for the fasta file containing aligned trimmed examples of transcription factor binding sites.\n");
			System.exit(0);
		}
	}
	

	
	public ScoreSequences(String[] args){
		processArgs(args);
	
		//get binding sites
		MultiFastaParser bsFasta = new MultiFastaParser(bindingSiteFile);
		bindingSites = bsFasta.getSeqs();
		
		//convert binding sites to upper case
		for (int i=0; i< bindingSites.length; i++) bindingSites[i] = bindingSites[i].toUpperCase();

		//make a MotifScanner
		double[][] PSPM = MemeMotif.makePSPM(bindingSites);
		double[][] llPSPM = MemeMotif.makeLLPSPM(PSPM, 0.25, 0.25, 0.25, 0.25);
		motifScanner = new MotifScanner(llPSPM);
		
		//fetch sequences to scan
		MultiFastaParser toScanFasta = new MultiFastaParser(seqFile);
		//motifScanner.scanPrintAllSequences(fastaParser.getSeqs(), fastaParser.getNames());
		System.out.println();
		
		motifScanner.scanPrintNumberHits(cutOff, toScanFasta.getSeqs(), toScanFasta.getNames());
		
		
	}
	
}
