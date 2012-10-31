package edu.utah.seq.data;

import java.io.*;
import java.text.NumberFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.annotation.*;
import util.gen.*;
import edu.utah.seq.data.*;
import trans.tpmap.*;

/**Merges, strips scores, sums identical position, finally merges PointData
 * @author Nix
 * */
public class PointDataManipulator {

	//fields
	private File[] pointDataDirectories;
	private HashMap<String, ArrayList<PointData>>[] splitPointData;
	private File saveDirectory;
	private boolean stripScores = false;
	private boolean mergeStrands = false;
	private boolean sumIdenticalPositions = false;
	private int shiftPointData = 0;
	private String chromosome;
	private int totalNumberReadsStarting = 0;
	private int totalNumberReadsFiltered = 0;
	private boolean verbose = true;

	/**For integration with the ChIPSeq app.*/
	public PointDataManipulator (File[] pointDataDirectories, File saveDirectory){
		this.pointDataDirectories = pointDataDirectories;
		this.saveDirectory = saveDirectory;
		stripScores = false;
		mergeStrands = false;
		sumIdenticalPositions = true;
		verbose = false;
		run();
	}

	/**For stand alone app.*/
	public PointDataManipulator(String[] args){
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

		//having problems clearing memory so using lots of extra code to kill the objects and any pointers to the objects
		while (it.hasNext()){
			chromosome = it.next();
			it.remove();
			if (verbose) System.out.println("\t"+chromosome);
			ArrayList<PointData> pdPlus = splitPointData[0].get(chromosome);
			ArrayList<PointData> pdMinus = splitPointData[1].get(chromosome);
			splitPointData[0].remove(chromosome);
			splitPointData[1].remove(chromosome);

			//count starting number of reads
			totalNumberReadsStarting += countNumberReads(pdPlus);
			totalNumberReadsStarting += countNumberReads(pdMinus);

			//replace scores with 1? this will load the data
			if (stripScores && sumIdenticalPositions == false) {
				stripScores(pdPlus);
				stripScores(pdMinus);
			}
			//shift position? this will load the data
			if (shiftPointData !=0){
				shiftPointData(pdPlus);
				shiftPointData(pdMinus);
			}
			//merge strands?
			if (mergeStrands){
				ArrayList<PointData> all = new ArrayList<PointData>();
				if (pdPlus !=null) all.addAll(pdPlus);
				if (pdMinus !=null) all.addAll(pdMinus);

				PointData merged;
				if (sumIdenticalPositions) merged = PointData.mergePointData(all, stripScores, true);
				else merged = PointData.combinePointData(all, true);
				merged.writePointData(saveDirectory);
				totalNumberReadsFiltered += merged.getInfo().getNumberObservations();
				merged = null;
			}
			else {
				if (pdPlus != null){
					PointData plus;
					if (sumIdenticalPositions) plus = PointData.mergePointData(pdPlus, stripScores, true);
					else plus = PointData.combinePointData(pdPlus, true);
					plus.writePointData(saveDirectory);
					totalNumberReadsFiltered += plus.getInfo().getNumberObservations();
					plus = null;
				}
				if (pdMinus != null){
					PointData minus;
					if (sumIdenticalPositions) minus = PointData.mergePointData(pdMinus, stripScores, true);
					else minus = PointData.combinePointData(pdMinus, true);
					minus.writePointData(saveDirectory);
					totalNumberReadsFiltered += minus.getInfo().getNumberObservations();
					minus = null;
				}
			}
			//null arrays
			nullArrays(pdPlus);
			nullArrays(pdMinus);
			pdPlus = null;
			pdMinus = null;
			System.gc();
		}

		if (sumIdenticalPositions && verbose){
			System.out.println("\n"+totalNumberReadsStarting+"\tNumber starting reads");
			System.out.println(totalNumberReadsFiltered+"\tNumber collapsed reads");
			double fraction = (double)totalNumberReadsFiltered/(double)totalNumberReadsStarting;
			System.out.println(Num.formatNumber(fraction, 3)+"\tFraction unique reads");
		}

	}

	public void checkMemory(){
		System.out.println("Checking Memory: ");
		for (int i=0; i< 20; i++) System.gc();
		Runtime runtime = Runtime.getRuntime();
		NumberFormat format = NumberFormat.getInstance();
		StringBuilder sb = new StringBuilder();
		long maxMemory = runtime.maxMemory();
		long allocatedMemory = runtime.totalMemory();
		long freeMemory = runtime.freeMemory();
		sb.append("\tfree memory: " + format.format(freeMemory / 1024) + "\n");
		sb.append("\tallocated memory: " + format.format(allocatedMemory / 1024) + "\n");
		sb.append("\tmax memory: " + format.format(maxMemory / 1024) + "\n");
		sb.append("\ttotal free memory: " + format.format((freeMemory + (maxMemory - allocatedMemory)) / 1024) + "\n");
		System.out.println(sb);
	}

	public void nullArrays(ArrayList<PointData> pd){
		if (pd == null) return;
		for (int i=0; i< pd.size(); i++){
			pd.get(i).nullPositionScoreArrays();
		}
	}

	public double fetchFractionUniqueReads(){
		return (double)totalNumberReadsFiltered/(double)totalNumberReadsStarting;
	}

	/**Sets all scores to 1. Typically these are a score associated with the alignment (a probability).*/
	public void stripScores (ArrayList<PointData> pdAL){
		if (pdAL == null || pdAL.size() == 0) return;
		int num = pdAL.size();
		for (int i=0; i< num; i++){
			pdAL.get(i).stripScores();
		}
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

	/**Shifts the base positions of the PointData.*/
	public void shiftPointData (ArrayList<PointData> pdAL){
		if (pdAL == null || pdAL.size() == 0) return;
		int num = pdAL.size();
		for (int i=0; i< num; i++){
			pdAL.get(i).shiftPositions(shiftPointData);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new PointDataManipulator(args);
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
					case 'i': sumIdenticalPositions = true; break;
					case 'd': shiftPointData = Integer.parseInt(args[++i]); break;
					case 'o': stripScores = true; break;
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
				"**                            Point Data Manipulator: Oct 2010                      **\n" +
				"**************************************************************************************\n" +
				"Manipulates PointData to merge strands, shift base positions, replace scores with 1\n" +
				"and sum identical positions. If multiple PointData directories are given, the data is\n" +
				"merged. \n\n" +

				"Options:\n"+
				"-p Point Data directories, full path, comma delimited. Should contain chromosome\n" +
				"       specific xxx.bar.zip or xxx_-_.bar files. Alternatively, provide one directory\n" +
				"       containing multiple PointData directories.\n"+
				"-s Save directory, full path.\n"+
				"-o Replace PointData scores with 1\n"+
				"-d Shift base position XXX bases 3', defaults to 0\n"+
				"-i Sum identical base position scores\n"+
				"-m Merge strands\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/PointDataManipulator -p\n" +
				"      /Data/Ets1Rep1/,/Data/Ets1Rep2/ -s /Data/MergedEts1 -o -i -m \n\n"+

		"**************************************************************************************\n");

	}		

}
