package trans.main;
import java.util.*;
import java.util.regex.*;
import util.bio.annotation.*;


import java.io.*;

import util.gen.*;

/**
 *	Collapses an array of Window into and array of Interval based on scores.
 */
public class IntervalMaker {
	//fields
	private double minScore = 0;
	private int minimumNumberOfOligos = 10;
	private int maxGap = 24;
	private String file;
	private Interval[] intervals;
	private ArrayList intervalsAl;
	private int scoreIndex = 1;
	private boolean lastLineReached;
	private int parPassWindows;
	private Interval interval;
	private Window testWindow;
	private Window[] windows;
	private int numberOfWindows;
	private int windowIndex;
	private int numWindowsInInt;	
	private double[] cutoffs = null;
	private boolean flipScore = false;
	boolean verbose = true;
	private int sizeOfOligo= 25;
	private int sizeOfOligoMinusOne = sizeOfOligo -1;
	private File regionsFile;
	private HashMap filterRegions;

	//contructors
	public IntervalMaker(String[] args){
		processArgs(args);
		//get files to process
		File[] files = IO.extractFiles(new File(file));	
		if(flipScore && verbose)System.out.println("\tMultiplied indexed window scores by -1 prior to building intervals.");

		//regions file?
		if (regionsFile !=null){
			if (verbose) System.out.println("\tFiltering windows by regions.");
			Coordinate[] c = Coordinate.parseFile(regionsFile, 0,0);//180, -180);
			Arrays.sort(c);
			filterRegions = Coordinate.splitByChromosome(c);
		}

		//for each Window[] file
		for (int i=0; i<files.length; i++){
			System.out.println("\nProcessing "+files[i]+"...");
			//fetch windows, watch out for non WindowInterval[] files
			windows = (Window[])IO.fetchObject(files[i]);
			//filter?
			if (regionsFile !=null) filterWindowsByRegions();
			//multiply window score by -1?
			if (flipScore) invertScores(windows);	
			makeIntervals(files[i]);
		}
	}


	//for integrating with SetNumberIntervalMaker
	public IntervalMaker(File windowFile, Window[] windows, double[] cutoffs, int sizeOfOligo, int scoreIndex, int maxGap, int minimumNumberOfOligos){
		this.windows = windows;
		numberOfWindows = windows.length;
		this.cutoffs = cutoffs;
		this.sizeOfOligo = sizeOfOligo;
		sizeOfOligoMinusOne = sizeOfOligo -1;
		this.scoreIndex = scoreIndex;
		this.maxGap = maxGap;
		this.minimumNumberOfOligos = minimumNumberOfOligos;
		verbose = false;
		//make em
		makeIntervals(windowFile);
	}

	public void filterWindowsByRegions(){
		String chrom = "";
		Coordinate[] mask = null;
		int index = 0;
		ArrayList filteredWindows = new ArrayList(windows.length/4);
		//for each window
		for (int i=0; i< windows.length; i++){
			//check if chrom has changed
			if (chrom.equals(windows[i].getChromosome()) == false){
				chrom = windows[i].getChromosome();
				mask = (Coordinate[]) filterRegions.get(chrom);
				index = 0;
			}
			//anything to mask?
			if (mask != null) {
				//make start and stop
				int startWin = windows[i].getStart1stOligo();
				int stopWin = windows[i].getStartLastOligo()+ sizeOfOligoMinusOne;
				//look for intersection starting at index, assumes windows and regions are sorted
				boolean intersection = false;
				for (int x=index; x< mask.length; x++){
					if (mask[x].intersects(startWin, stopWin)){
						//System.out.println(windows[i].getChromosome()+ " " +startWin+" "+stopWin+" : "+mask[x]);
						index = x;
						intersection = true;
						break;
					}
				}
				if (intersection == false) filteredWindows.add(windows[i]);
			}
			else filteredWindows.add(windows[i]);
		}
		//convert window array list 
		windows = new Window[filteredWindows.size()];
		filteredWindows.toArray(windows);
	}

	public void makeIntervals(File windowFile){
		numberOfWindows = windows.length;
		//for each score threshold, typically there is only one
		for (int z=0; z< cutoffs.length; z++){
			//set minScore
			minScore = cutoffs[z];

			//initialize params
			intervalsAl = new ArrayList(1000);
			lastLineReached = false;
			parPassWindows = 1;
			interval = null;
			testWindow = null;
			windowIndex = 0;
			numWindowsInInt = 0;

			if (numberOfWindows >0 && windows[0]!=null){

				//advance till good window found, make this the starting interval
				Window firstWin = fetchGoodWindow();

				if (lastLineReached == false && firstWin!=null){
					interval = new Interval (firstWin, sizeOfOligo);
					boolean addIt;			

					for (; windowIndex<numberOfWindows; windowIndex++){
						testWindow = windows[windowIndex];			
						addIt=false;

						if (checkScoreAndNumOligos(testWindow)){
							parPassWindows++;
							//check chromosome
							if (interval.getChromosome().equals(testWindow.getChromosome())){
								addIt=true;
							}
							else{
								//different chromosome, save old interval, start new
								addIt();
								interval = new Interval (testWindow, sizeOfOligo);
								numWindowsInInt = 1;
							}					
						}
						if (addIt){
							//check length
							if ((testWindow.getStart1stOligo()-interval.getStartLastOligo()) <= maxGap){
								//length, score, look good so meld
								interval.setStartLastOligo(testWindow.getStartLastOligo());
								numWindowsInInt++;						
								//Save window
								if (testWindow.getScores()[scoreIndex]>interval.getBestWindow().getScores()[scoreIndex]){
									interval.setBestWindow(testWindow);
								}
							}
							//length bad but good window, thus
							else {
								//write interval to disk
								addIt();	
								numWindowsInInt = 1;			
								//make new interval
								interval = new Interval (testWindow, sizeOfOligo);
							}
						}
					}	

					//close and save final interval, check to see that last line wasn't reached by fetchGoodInterval()
					if (lastLineReached==false) {
						addIt();		
					} 			
				}

				int numInt = intervalsAl.size();
				if (verbose) System.out.println("\t"+numInt+" Intervals created using a Max Gap of "+maxGap+", Score Index: "+scoreIndex+", Min Score: "+minScore+ ", Min # Oligos: "+minimumNumberOfOligos);

				//convert ArrayList to Interval[] and save to disk
				intervals = new Interval[intervalsAl.size()];
				intervalsAl.toArray(intervals);
				File serWinIntArray = new File(windowFile+""+scoreIndex+"Indx"+intervalsAl.size());
				IO.saveObject(serWinIntArray, intervals);

			}
			else {
				System.out.println("\tNo Windows?! Skipping "+windowFile+"...");
			}

		}
	}

	public void addIt(){
		interval.setNumberOfWindows(numWindowsInInt);
		intervalsAl.add(interval);
	}	

	/**Multiplies all of the Window[] scores by -1*/
	public static void invertScores(Window[] windows){
		for (int j=0; j< windows.length; j++){
			double[] scores = windows[j].getScores();
			for (int k=0; k<scores.length; k++){
				scores[k] = -1 * scores[k];
			}
			windows[j].setScores(scores);
		}
	}

	public Window fetchGoodWindow(){
		for (; windowIndex< numberOfWindows; windowIndex++){
			if (checkScoreAndNumOligos(windows[windowIndex])) {
				numWindowsInInt = 1;	
				return windows[windowIndex++]; 
			} 
		}
		//ran thru stop of array add last if good
		if (interval!=null) {
			addIt();
		} 
		lastLineReached = true;
		return null;	
	}

	public boolean checkScoreAndNumOligos(Window win){
		if (win.getScores()[scoreIndex] >= minScore && win.getNumberOligos() >= minimumNumberOfOligos)return true;
		return false;
	}
	public static void main(String[] args) {
		if (args.length<1){
			printDocs();
			System.exit(0);	
		}
		new IntervalMaker(args);
	}

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String cutoffString = null;
		System.out.println("\nParameters: "+Misc.stringArrayToString(args, " "));
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'i': scoreIndex =Integer.parseInt(args[i+1]);i++; break;
					case 'o': minimumNumberOfOligos =Integer.parseInt(args[i+1]);i++; break;
					case 'z': sizeOfOligo =Integer.parseInt(args[i+1]);i++; break;
					case 'g': maxGap =Integer.parseInt(args[i+1]);i++; break;
					case 'm': flipScore = true; break;
					case 'f': file = args[i+1]; i++; break;
					case 'r': regionsFile = new File(args[i+1]); i++; break;
					case 's': cutoffString = args[i+1]; i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		if (Misc.isEmpty(file)){
			System.out.print("\nPlease enter an TiMAT2 results file or directory.\n");
			System.exit(0);
		}
		if (cutoffString !=null) {
			String[] doubles = cutoffString.split(",");
			cutoffs = Num.parseDoubles(doubles);
		}
		else cutoffs = new double[]{minScore};
		//set size of oligo minus one
		sizeOfOligoMinusOne = sizeOfOligo -1;
	}	


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Interval Maker: June 2007                              **\n" +
				"**************************************************************************************\n" +
				"IM combines Windows into larger Intervals given a minimal score, minimum number of\n" +
				"oligos, and a maximum gap.\n\n" +

				"Options:\n\n" +

				"-s Minimal score, defaults to 0. Can also provide a comma delimited list, no spaces,\n" +
				"     to use in generating multiple Interval arrays.\n" +
				"-i Score index, the score you wish to use in merging, see ScanChip or ScanChromosome.\n"+
				"-m Multiply window scores by -1. Useful for looking for reduced regions.\n"+
				"-o Minimum number of oligo positions in each window, defaults to 10.\n"+
				"-g Maximum allowable bp gap between starts of oligos, defaults to 24.\n" +
				"-z Size of oligo, defaults to 25.\n"+
				"-r Enter a full path tab delimited regions file text (chr start stop) to use in\n" +
				"      removing intersecting windows. Coordinates are assumed to be zero based and stop\n" +
				"      inclusive.\n"+
				"-f Full path file text for the Window[] array, if a directory is specified,\n" +
				"      all files will be processed.\n\n" +

				"Example: java -Xmx500M -jar pathTo/T2/Apps/IntervalMaker -f /affy/res/zeste.res -s 50\n" +
				"      -i 0 -g -200 -o 5\n\n" +
		"**************************************************************************************\n");		
	}

}
