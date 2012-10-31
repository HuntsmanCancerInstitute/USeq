package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.gen.*;
import edu.utah.seq.data.*;



/**
 * @author Nix
 * */
public class BisStatRegionMaker {

	//user defined fields
	private File[] serializedWindowObjects;
	private float minimumFraction = 0;
	private float maximumFraction = 1;
	private int maximumGap = 0;
	private SmoothingWindow[] smoothingWindow;
	private String chromosome;
	private int scoreIndex = 0;
	private long numberERs = 0;

	//constructors
	/**Stand alone.*/
	public BisStatRegionMaker(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		
		EnrichedRegionMaker erm = new EnrichedRegionMaker(maximumGap);
		try{
			String index = "";
			if (scoreIndex !=0) index = "_Indx"+scoreIndex;
			File bed = new File (serializedWindowObjects[0].getParentFile()+index+"Min"+minimumFraction+"Max"+maximumFraction+".bed");
			PrintWriter out = new PrintWriter(new FileWriter( bed));
			//for each window object
			System.out.print("\nScanning chromosomes:");
			for (File f : serializedWindowObjects){
				chromosome = f.getName();
				System.out.print(" "+chromosome);
				smoothingWindow = (SmoothingWindow[])IO.fetchObject(f);
				thresholdWindows();
				if (smoothingWindow != null){
					EnrichedRegion[] ers = erm.fetchEnrichedRegions(smoothingWindow, chromosome);
					numberERs += ers.length;
					for (EnrichedRegion er: ers) out.println(er.getBed5(scoreIndex));
				}
			}
			System.out.println();
			System.out.println("\n"+numberERs+" Regions found.");
			out.close();
			
		} catch (Exception e){
			e.printStackTrace();
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	public void thresholdWindows(){
		ArrayList<SmoothingWindow> al = new ArrayList<SmoothingWindow>();
		for (SmoothingWindow sw : smoothingWindow){
			float fraction = sw.getScores()[scoreIndex];
			if (fraction >= minimumFraction && fraction <= maximumFraction) al.add(sw);
		}
		if (al.size() == 0) smoothingWindow = null;
		else {
			smoothingWindow = new SmoothingWindow[al.size()];
			al.toArray(smoothingWindow);
		}
	}

	//methods
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BisStatRegionMaker(args);
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
					case 's': serializedWindowObjects = IO.extractOnlyFiles(new File(args[++i])); break;
					case 'm': minimumFraction = Float.parseFloat(args[++i]); break;
					case 'x': maximumFraction = Float.parseFloat(args[++i]); break;
					case 'g': maximumGap = Integer.parseInt(args[++i]); break;
					case 'q': scoreIndex = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (serializedWindowObjects == null || serializedWindowObjects.length ==0) Misc.printErrAndExit("\nCannot find any object files in your SerializedWindowObjects directory?! ");

		System.out.println("Params:");
		System.out.println(maximumGap+"\tMaximum Gap");
		System.out.println(minimumFraction+"\tMinimum fraction");
		System.out.println(maximumFraction+"\tMaximum fraction");
		if (scoreIndex !=0) System.out.println(scoreIndex+"\tDensity quartile index");
		
		if (minimumFraction >1 || maximumFraction >1) Misc.printErrAndExit("\nError, min and or max cannot be greater than one!\n");

	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           BisStat Region Maker: March 2012                       **\n" +
				"**************************************************************************************\n" +
				"Takes serialized window objects from BisStat, thresholds based on the min and max\n" +
				"fraction methylation params and prints regions in bed format meeting the criteria.\n" +
				"May also build regions base on the density of a given fraction methylation quartile.\n" +
				"For example, to identify regions where at least 0.8 of the sequenced Cs are low\n" +
				"methylated (<= 0.25 default settings in BisStat) set -q 1 -m 0.8 . To find regions of\n" +
				"with >= 0.9 of the Cs with high methylation (>= 0.75 default BisStat setting), set\n" +
				"-q 3 -m 0.9  . \n\n" +

				"Options:\n"+
				"-s SerializedWindowObject directory from BisStat, full path.\n"+
				"-m Minimum fraction.\n" +
				"-x Maximum fraction.\n" +
				"-g Maximum gap, defaults to 0.\n"+
				"-q Merge windows based on their quartile density score, not fraction methylation, by\n" +
				"      indicating 1,2,or 3 for 1st, 2nd+3rd, or 4th, respectively.\n"+

				"\n"+
				"Example: java -Xmx4G -jar pathTo/USeq/Apps/BisStatRegionMaker -m 0.8 -x 1.0 -g 100\n" +
				"      -s /Data/BisStat/SerializedWindowObjects  \n\n" +

		"**************************************************************************************\n");

	}
}
