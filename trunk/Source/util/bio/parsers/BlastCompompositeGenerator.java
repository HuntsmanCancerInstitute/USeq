package util.bio.parsers;
import java.io.*;
import java.util.*;
import util.gen.*;
import util.bio.parsers.*;

public class BlastCompompositeGenerator {
	//fields
	private File refSeq;
	private File[] blastResults;
	private BaseCount[] baseCounts;

	//constructor
	public BlastCompompositeGenerator(String[] args){
		System.out.println("\nLaunching ...");
		//assign files
		refSeq = new File (args[0]);
		blastResults = IO.extractFiles(new File(args[1]), ".blast");
		//instantiate BaseCount[]
		makeBaseCounts();
		//count blast results
		scoreBlastResults();
		//print results to screen
		System.out.println("Position\tBase\tG\tA\tT\tC\tN\tDel\tInt");
		Misc.printArray(baseCounts);
	}

	public void scoreGroup(int bcIndex, char[] query, char[] align, char[] sbjct){
		//for each base in query

		for (int i=0; i<query.length; i++){
			//any gaps?
			if (align[i] == '|' || (query[i] != '-' && sbjct[i] != '-') ) baseCounts[bcIndex++].countBase(sbjct[i]);
			//gap in subject, therefore deletion in subject
			else if (sbjct[i] == '-') {
				baseCounts[bcIndex].setDeletions(baseCounts[bcIndex].getDeletions() +1);
				bcIndex++;
			}
			//gap in query, therefore insertion
			else {
				//attribute gap to preceeding base count
				bcIndex--;
				//count number of insertions
				int numGaps = 1;
				while (true){
					//at stop? or not a gap
					int nextIndex = i+1;
					if (nextIndex == query.length || query[nextIndex] != '-') break;
					i++;
					numGaps++;
				}
				baseCounts[bcIndex++].getInsertionLengths().add(new Integer(numGaps));
			}
		}
	}
	
	public char[] align(char[] q, char[] s){
		char[] align = new char[q.length];
		for (int i=0; i< q.length; i++){
			if (q[i] == s[i]) align[i] = '|';
			else align[i] = 'x';
		}
		return align;
	}

	public void scoreLines(String[] lines){
		for (int i=0; i< lines.length; i++){
			//look for alignment
			if (lines[i].startsWith("Query:")){
				//parse start index
				String[] tokens = lines[i].split("\\s+");
				int startIndex = Integer.parseInt(tokens[1]) -1;
				//assign lines
				char[] query = tokens[2].toCharArray();
				i+=2;
				tokens = lines[i].split("\\s+");
				char[] sbjct = tokens[2].toCharArray();
				char[] align = align(query,sbjct);
				//is first query a gap?
				//score group
				scoreGroup(startIndex, query, align, sbjct);

			}
		}
	}

	public void scoreBlastResults() {
		//for each file of blast results 
		for (int i=0; i< blastResults.length; i++){
			//load text
			String[] lines = IO.loadFileIntoStringArray(blastResults[i]);
			//score lines
			scoreLines(lines);
		}
	}

	public void makeBaseCounts(){
		MultiFastaParser mfp = new MultiFastaParser (refSeq);
		String seq = mfp.getSeqs()[0].toLowerCase();
		baseCounts = new BaseCount[seq.length()];
		for (int i=0; i< baseCounts.length; i++) baseCounts[i] = new BaseCount(i+1, seq.charAt(i));
	}
	
	public static void main(String[] args){
		if (args.length ==0) Misc.printExit("\nEnter the full path file names for the reference fasta sequence file and the directory/ file containing xxx.blast file results.\n");
		new BlastCompompositeGenerator (args);
	}
}


