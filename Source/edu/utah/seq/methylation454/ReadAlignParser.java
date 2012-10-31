package edu.utah.seq.methylation454;
import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;

public class ReadAlignParser {
	//fields
	private Reference reference;
	private Read[] reads;
	private HashMap indexCoordinates = null;
	private String sampleId;


	public static void main (String[] args){
		if (args.length == 0) Misc.printExit("\nEnter full path file names for the indexSequenceGenomeCoordinates and readAlign files.\n");
		new ReadAlignParser(new File (args[1]), new File (args[0]));
	}


	public ReadAlignParser (File readAlignFile, File indexFile){
		//make hash map containing the index aligmnemt sequence text (I1, I2, I3...) and the start genomic coordinates
		indexCoordinates = IO.loadFileIntoHashMap(indexFile);
		//for each file
		File[] toParse = IO.extractFiles(readAlignFile);
		for (int i=0; i< toParse.length; i++){
			//parse sampleID
			sampleId = toParse[i].getName().split("_")[0];
			//parse aligned read file
			parseIt(toParse[i]);
			//print lines
			printIt();
		}
	}

	public void printIt(){
		int[] positions = reference.getCpGPositions();
		int[] realPositions = reference.getCpGGenomicPositions();
		for (int x=0; x< reads.length; x++ ){
			//for each position find base in read
			for (int i=0; i< positions.length; i++){
				//is position after start?
				if (positions[i] >= reads[x].getAlignmentStart()){
					//is position before stop?
					if (positions[i] < reads[x].getAlignmentEnd()){
						//get base!
						System.out.println(
								sampleId+"\t"+
								reference.getName()+"\t"+
								reference.getChromosome()+":"+(realPositions[i]+reference.getStartPosition())+"\t"+
								x+"\t"+
								reads[x].getStrand()+"\t"+
								reads[x].getSequence().charAt(positions[i]-reads[x].getAlignmentStart()));
					}
					else break;
				}
			}
		}




	}

	public void parseIt (File file){
		String[] lines = IO.loadFileIntoStringArray(file);
		//make reference
		makeReference(lines);
		//make reads
		int number = (lines.length/2) - 1;
		reads = new Read[number];
		int counter = 0;
		for (int i=2; i<lines.length; i++){
			String header = lines[i++];
			String seq = lines[i];
			reads[counter++] = new Read(header, seq);
		}
	}

	public void makeReference (String[] lines){
		//is it a reference seq?
		if (lines[0].indexOf("Reference") == -1) Misc.printExit("\nError: cannot find reference sequence in read align file "+lines[0]);
		//parse text
		Pattern p = Pattern.compile(".*rNm=\"I(\\d+)\".*");
		Matcher mat = p.matcher(lines[0]);
		String indexSeqName = null;
		if (mat.matches()) indexSeqName = "I"+mat.group(1);
		else Misc.printExit("\nError: cannot extract index sequence text from -> "+lines[0]);
		reference = new Reference (lines[1]);
		//set coordinates
		String chromStart = (String)indexCoordinates.get(indexSeqName);
		String[] tokens = chromStart.split(":");
		reference.setChromosome(tokens[0]);
		reference.setStartPosition(Integer.parseInt(tokens[1]));
		//set text
		reference.setName(indexSeqName+"_"+tokens[2]);

	}

}
