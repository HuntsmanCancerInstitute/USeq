package edu.utah.seq.data;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.annotation.ExportIntergenicRegions;
import util.gen.*;

public class Graph2Bed {

	//fields
	private File[] pointDataDirectories;
	private PrintWriter out = null;
	private float threshold = 0;


	//constructors
	public Graph2Bed(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		processArgs(args);

		System.out.println("Parsing:");
		thresholdData();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	//methods
	public void thresholdData(){
		try {
			//for each directory of stair-step point data
			for (File directory: pointDataDirectories){
				System.out.print("\t"+directory.getName());

				//open print writer for results
				File bedFile = new File(Misc.removeExtension(directory.toString())+"_"+threshold+".bed");
				out = new PrintWriter (new FileWriter (bedFile));

				//fetch PointData, using this craziness to remove the reference
				HashMap<String, PointData> pd = PointData.fetchPointData(directory, null, false);
				String[] chromosomes = new String[pd.size()];
				int index = 0;
				for (String c: pd.keySet()) chromosomes[index++] = c;
				for (String chromosome: chromosomes) {
					System.out.println(" "+chromosome);
					parsePointData(pd.remove(chromosome));
				}
				System.out.println();
				
				//cleanup
				out.close();

			}
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}


	public void parsePointData(PointData pd){
		String chromosome = pd.getInfo().getChromosome() +"\t";

		//get positions and scores
		int[] positions = pd.getPositions();
		float[] scores = pd.getScores();

		int lastBase = positions[positions.length-1];
		boolean[] falseMask = new boolean[lastBase+100];
		Arrays.fill(falseMask, true);

		//for each position
		for (int i=0; i< positions.length; i++){
			//good score?
			if (scores[i] > threshold) {
				int startBase = positions[i];
				if (startBase < 0) startBase = 0;
				int stopBase = i+1;
				if (stopBase >= positions.length) {
					stopBase = lastBase;
				}
				else stopBase = positions[stopBase];
				if (stopBase < 0) continue;
				//fill the block
				Arrays.fill(falseMask, startBase, stopBase, false);
			}
		}

		//fetch blocks and print
		int[][] blocks = ExportIntergenicRegions.fetchFalseBlocks(falseMask, 0, 0);
		for (int j=0; j< blocks.length; j++) out.println(chromosome+blocks[j][0]+"\t"+(blocks[j][1]+1));
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new Graph2Bed(args);
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
					case 'p': pointDataDirectories = IO.extractFiles(args[++i]); break;
					case 't': threshold = Float.parseFloat(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for point directories
		if (pointDataDirectories == null || pointDataDirectories[0].isDirectory() == false) Misc.printExit("\nError: cannot find your PointData directories(s)!\n");
		//only one directory look deeper
		if (pointDataDirectories.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(pointDataDirectories[0]);
			if (otherDirs != null && otherDirs.length > 0) pointDataDirectories = otherDirs;
		}


	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                 Graph 2 Bed: Feb 2011                            **\n" +
				"**************************************************************************************\n" +
				"Converts USeq stair step and heat map graphs into region bed files using a threshold.\n" +
				"Do not use this with non USeq generated graphs. Won't work with bar or point graphs.\n\n" +

				"Options:\n"+
				"-p Point Data directories, full path, comma delimited. Should contain chromosome\n" +
				"       specific xxx.bar.zip or xxx_-_.bar files. May point this to a single directory\n" +
				"       of such too.\n"+
				"-t Threshold, regions exceeding it will be saved, defaults to 0.\n"+

				"\nExample: java -Xmx1500M -jar pathTo/USeq/Apps/Graph2Bed -t 9 -p /data/ReadCoverage\n\n"+

		"**************************************************************************************\n");

	}	

}


