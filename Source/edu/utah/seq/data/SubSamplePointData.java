package edu.utah.seq.data;
import util.gen.*;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.parsers.Tag2Point;

public class SubSamplePointData {

	//fields
	private int desiredNumberObservations;
	private File saveDirectory;
	private File[] pointDataDirectories;
	private int minimumObservationsPerDirectory = 0;

	
	/**Randomly samples out reads to match a desired number.*/
	public SubSamplePointData(String[] args){
		processArgs(args);
		System.out.println("\nStraight/ Even draw from all PointData directories...");
		//how many to fetch per directory?
		int numberPerDirectory = (int)Math.round((double)desiredNumberObservations/(double)pointDataDirectories.length);
		System.out.println("Drawing "+numberPerDirectory+" from each directory");	
		ArrayList<HashMap<String,PointData[]>> al = new ArrayList<HashMap<String,PointData[]>>();
		int numSkipped =0;
		
		//for each point directory
		for (int i=0; i< pointDataDirectories.length; i++){
			//fetch the data
			HashMap<String,PointData[]> strandedData = PointData.fetchStrandedPointData(pointDataDirectories[i]);
			if (strandedData == null || strandedData.size() ==0) Misc.printExit("\nProblem with loading "+pointDataDirectories[i].getName());
			int totalObs = PointData.totalObservationsMultiPointData(strandedData);
			if (totalObs < minimumObservationsPerDirectory) {
				System.out.println("Skipping "+"\t"+pointDataDirectories[i].getName() +"\t"+ totalObs);
				numSkipped++;
			}
			else {
				System.out.println("\t"+pointDataDirectories[i].getName() +"\t"+ totalObs);
				//fetch desired number
				if (totalObs < numberPerDirectory){
					System.out.println("\t\tWarning! Too few reads, adding only "+totalObs);
					al.add(strandedData);
				}
				else {
					if (PointData.subSampleStranded(strandedData, numberPerDirectory) == false) Misc.printErrAndExit("\nProblem subsampling PointData.");
					al.add(strandedData);
				}
			}
		}
		//merge
		//HashMap<String,PointData[]> merged = PointData.combinePointData(al, true);
		HashMap<String,PointData[]> merged = PointData.combinePointDataToHashMap(al, true); 
		HashMap<String,PointData>[] split = PointData.splitStrandedPointData(merged);

		//save
		System.out.println("Saving "+PointData.totalObservationsMultiPointData(merged) + " Points, skipped "+numSkipped);
		PointData.writePointData(split[0], saveDirectory);
		PointData.writePointData(split[1], saveDirectory);

		System.out.println("Done!");
	}
	
	public static void subSamplePointData(int desiredNumberObservations, File pointDataDirectory, File saveDirectory){
		HashMap<String,PointData[]> strandedData = PointData.fetchStrandedPointData(pointDataDirectory);
		PointData.subSampleStranded(strandedData, desiredNumberObservations);
		HashMap<String,PointData>[] split = PointData.splitStrandedPointData(strandedData);
		PointData.writePointData(split[0], saveDirectory);
		PointData.writePointData(split[1], saveDirectory);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SubSamplePointData(args);
	}		

	/**This method will process each argument and assign new varibles*/
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
					case 'f': pointDataDirectories = IO.extractFiles(args[i+1]); i++; break;
					case 'n': desiredNumberObservations = Integer.parseInt(args[++i]); break;
					case 'm': minimumObservationsPerDirectory = Integer.parseInt(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (pointDataDirectories == null || pointDataDirectories[0].canRead() == false) Misc.printExit("\nError: cannot find/ read your PointData directories(s)!\n");
		if (saveDirectory == null) Misc.printExit("\nPlease enter a directory in which to save the merged sub sampled data.\n");
		saveDirectory.mkdir();

		//sub directories?
		if (pointDataDirectories.length == 1){
			File[] subDirs = IO.extractOnlyDirectories(pointDataDirectories[0]);
			if (subDirs != null) pointDataDirectories = subDirs;
		}

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             SubSamplePointData:  Dec 2008                        **\n" +
				"**************************************************************************************\n" +
				"SSPD takes PointData directories and randomly selects points from each directory and\n" +
				"saves the merge.\n\n" +

				"-f Comma delimited full path PointDataDirectories from which to draw or a single \n" +
				"       directory containing multiple PointDataDirectories.\n" +
				"-n Total number of observations desired.\n"+
				"-s Full path file directory in which to save the results.\n"+
				//"-m Minimum number observations in each directory to include in sub sampling.\n"+

				"\nExample: java -Xmx1500M -jar pathTo/USeq/Apps/SubSamplePointData -n 10000000 -f\n" +
				"    /Data/WCE1_Point,/Data/WCE2_Point,/Data/WCE3_Point -s /Data/Sub/ \n\n" +

		"**************************************************************************************\n");

	}


}
