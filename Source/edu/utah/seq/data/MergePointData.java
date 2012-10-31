package edu.utah.seq.data;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.annotation.*;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.BarParser;
import trans.tpmap.*;

/**Merges, strips scores, sums identical position, finally merges PointData
 * @author Nix
 * */
public class MergePointData {

	//fields
	private File[] pointDataDirectories;
	private HashMap<String, ArrayList<PointData>>[] splitPointData;
	private File saveDirectory;
	private boolean mergeStrands = false;
	private boolean replaceScoresWithHitCount = true;
	private String chromosome;
	private long totalNumberReadsStarting = 0;
	private long totalNumberReadsFiltered = 0;
	private boolean verbose = true;



	/**For stand alone app.*/
	public MergePointData(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		//process args
		processArgs(args);
		run();
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}


	public void run(){
		//fetch names
		splitPointData = PointData.fetchStrandedPointDataNoMerge (pointDataDirectories);

		//fetch all chromosomes
		HashSet<String> chromosomes = new HashSet<String>();
		chromosomes.addAll(splitPointData[0].keySet());
		chromosomes.addAll(splitPointData[1].keySet());

		//for each chromosome
		Iterator<String> it = chromosomes.iterator();

		while (it.hasNext()){
			chromosome = it.next();
			if (verbose) System.out.println("\t"+chromosome);
			//count starting number of reads
			totalNumberReadsStarting += countNumberReads(splitPointData[0].get(chromosome));
			totalNumberReadsStarting += countNumberReads(splitPointData[1].get(chromosome));
			
			//merge strands?
			if (mergeStrands){
				ArrayList<PointData> all = new ArrayList<PointData>();
				if (splitPointData[0].get(chromosome) !=null) all.addAll(splitPointData[0].get(chromosome));
				if (splitPointData[1].get(chromosome) !=null) all.addAll(splitPointData[1].get(chromosome));
				//PointData p = PointData.efficientMergeSummingScores(all,replaceScoresWithHitCount);
				PointData p = PointData.mergePointData(all,replaceScoresWithHitCount, true);
				p.writePointData(saveDirectory);
				totalNumberReadsFiltered += p.getInfo().getNumberObservations();
			}
			else {
				ArrayList<PointData> pd = splitPointData[0].get(chromosome);
				if (pd != null) {
					//PointData p = PointData.efficientMergeSummingScores(pd,replaceScoresWithHitCount);
					PointData p = PointData.mergePointData(pd,replaceScoresWithHitCount, true);
					p.writePointData(saveDirectory);
					totalNumberReadsFiltered += p.getInfo().getNumberObservations();
				}
				pd = splitPointData[1].get(chromosome);
				if (pd != null) {
					//PointData p = PointData.efficientMergeSummingScores(pd,replaceScoresWithHitCount);
					PointData p = PointData.mergePointData(pd,replaceScoresWithHitCount, true);
					p.writePointData(saveDirectory);
					totalNumberReadsFiltered += p.getInfo().getNumberObservations();
				}
			}
			
		}
			System.out.println("\n"+totalNumberReadsStarting+"\tNumber starting observations");
			System.out.println(totalNumberReadsFiltered+"\tNumber collapsed observations");
			double fraction = (double)totalNumberReadsFiltered/(double)totalNumberReadsStarting;
			System.out.println(Num.formatNumber(fraction, 3)+"\tFraction unique observations");
		
	}

	

	/**Counts number of reads.*/
	public int countNumberReads (ArrayList<PointData> pdAL){
		if (pdAL == null || pdAL.size() == 0) return 0;
		int num = pdAL.size();
		int numberReads = 0;
		for (int i=0; i< num; i++){
			numberReads += pdAL.get(i).getInfo().getNumberObservations();
		}
		return numberReads;
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MergePointData(args);
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
					case 's': saveDirectory = new File(args[++i]); break;
					case 'm': mergeStrands = true; break;
					case 'c': replaceScoresWithHitCount = false; break;
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

		//create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		if (saveDirectory.exists() == false) saveDirectory.mkdir();
		else System.out.println("\nWARNING: the save directory exists, will over write files contained within.\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Merge Point Data: Jan 2011                           **\n" +
				"**************************************************************************************\n" +
				"Efficiently merges PointData, collapsing by position and possibly strand. Identical\n" +
				"position scores are either summed or converted into counts. DO NOT use this app on\n" +
				"PointData that will be part of a primary chIP/RNA-seq analysis.  It is only for\n" +
				"bis-seq and visualization purposes.\n\n" +

				"Options:\n"+
				"-p Point Data directories, full path, comma delimited. Should contain chromosome\n" +
				"       specific xxx.bar.zip or xxx_-_.bar files. Alternatively, provide one directory\n" +
				"       containing multiple PointData directories.\n"+
				"-s Save directory, full path.\n"+
				"-c Don't replace scores with hit count, just sum existing scores.\n"+
				"-m Merge strands\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/MergePointData -p\n" +
				"      /Data/Ets1Rep1/,/Data/Ets1Rep2/ -s /Data/MergedEts1 -m \n\n"+

		"**************************************************************************************\n");

	}		

}
