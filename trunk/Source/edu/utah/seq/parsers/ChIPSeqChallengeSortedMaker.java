package edu.utah.seq.parsers;
import java.io.*;
import edu.utah.seq.data.*;
import util.gen.*;

public class ChIPSeqChallengeSortedMaker {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		File pointDataDir = new File(args[0]);
		File[] pointDataFiles = IO.extractFiles(pointDataDir, ".bar.zip");
		File sortedFile = new File(pointDataDir,"fa_sorted.txt");
		//File sortedFile = new File(pointDataDir,"export.txt");
		//File sortedFile = new File(pointDataDir,"bed.txt");
		
		try{
			PrintWriter out = new PrintWriter (new FileWriter(sortedFile));
			for (int i=0; i< pointDataFiles.length; i++){
				System.out.println(pointDataFiles[i]);
				PointData pd = new PointData (pointDataFiles[i], true);
				String chromSpace = pd.getInfo().getChromosome()+" ";
				String spaceStrand = " "+pd.getInfo().getStrand();
				//String strand = "F";
				//if (pd.getInfo().getStrand().equals("-")) strand = "R";
				//String preSortedLine = "HWI\t\t1\t1\t1\t1\t\t\tAAAAAAAAAAAAAAAAAAAAAAAAAA\tZZZZZZZZZZZZZZZZZZZZZZZZZZ\t"+chrom+"\t\t";
				//String postSortedLine = "\t"+strand+"\t26\t74\t\t\t\t\t";
				//String postExportLine = "\t"+strand+"\t26\t74\t\t\t\t\t\tY";
				//String preBed = chrom+"\t";
				//String postBed = "\tp\t74\t"+pd.getInfo().getStrand();
				int[] positions = pd.getPositions();
				float[] scores = pd.getScores();
				for (int j=0; j< positions.length; j++){
					//bed chr, start, stop, text, score, strand
					//out.print(preBed);
					//out.print(positions[j]-13);
					//out.print("\t");
					//out.print(positions[j]+13);
					//out.println(postBed);
					
					//out.print(preSortedLine);
					//out.print(positions[j]-12);
					//out.println(postSortedLine);
					//out.println(postExportLine);
					
					//out.print(chromSpace);
					//out.print(positions[j]-13);
					//out.println(spaceStrand);
					out.println(scores[j]);
				}
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}

	}

}
