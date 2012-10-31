package trans.main;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import java.io.*;

import util.gen.*;

/**
 * Determines the threshold(s) needed to generate a set number of intervals, multiple score indexes accepted.
 *
 */
public class SetNumberIntervalMaker {
	private String file;
	private File[] files;
	private int scoreIndex = 0;
	private boolean scanOnlyOneIndex = false;
	private boolean flipScore = false;
	private int minimumNumberOfOligos = 10;
	private int sizeOfOligo = 25;
	private Window[] windows;
	private int numWindows;
	private int maxGap = 24;
	private int[] numberIntervals;
	private ArrayList intervalsAl;
	private boolean makeIntervals = false;
	
	//internal fields, don't mess with
	private boolean lastLineReached;
	private int parPassWindows;
	private Interval interval;
	private Window testWindow;
	private int windowIndex;
	private int numWindowsInInt;	
	private int numberOfWindows;
	
	public SetNumberIntervalMaker(String[] args){
		processArgs(args);
		if (flipScore) System.out.println("\nMultiplying all window scores by -1.\n");
		//print header
		System.out.println("File\tIndex\tMin # Oligos\tThreshold\t# Intervals");
		//for every window file
		for (int i=0; i<files.length; i++){
			
			//get windows 
			windows = (Window[])IO.fetchObject(files[i]);
			numWindows = windows.length;
			
			//flip scores?
			if (flipScore) IntervalMaker.invertScores(windows);
			
			//for each score index
			int numIndexes;
			int j;
			if (scanOnlyOneIndex){
				j = scoreIndex;
				numIndexes = j + 1;
			}
			else {
				j = 0;
				numIndexes = windows[0].getScores().length;
			}
			for (; j<numIndexes; j++ ){
				scoreIndex = j;
				//scan to find min max scores for a given index
				double[] startingMinMax = minMaxWindowScores();
				double startingThreshold = ((startingMinMax[1] - startingMinMax[0]) / 2) + startingMinMax[0];
				double[] thresholds = new double[numberIntervals.length];
				
				//for each set number of intervals they want to make
				for (int k=0; k< numberIntervals.length; k++){
					double[] minMax = {startingMinMax[0], startingMinMax[1]};
					double threshold = startingThreshold;
					int rounds = 25;
					int num = 0;
					while( rounds-- !=0){
						num = makeMockIntervals(threshold);
						if (num == numberIntervals[k]) {
							break;
						}
						else if (num < numberIntervals[k]){
							minMax[1] = threshold;
						}
						else {
							minMax[0] = threshold;
						}
						threshold = ((minMax[1] - minMax[0]) / 2) + minMax[0];
					}
					thresholds[k] = threshold;
					System.out.println(files[i].getName()+"\t"+scoreIndex+"\t"+minimumNumberOfOligos+"\t"+threshold+"\t"+num);
				}
				
				//make intervals?
				if (makeIntervals){
					//(File windowFile, double[] cutoffs, int sizeOfOligo, int scoreIndex, boolean flipScore, int maxGap, int minimumNumberOfOligos)
					new IntervalMaker(files[i], windows, thresholds, sizeOfOligo, scoreIndex, maxGap,minimumNumberOfOligos);
				}
			}
		}
	}
	
	/**Finds min max*/
	public double[] minMaxWindowScores(){
		double min=0;
		double max=0;
		for (int i=0; i< numWindows; i++){
			double score = windows[i].getScores()[scoreIndex];
			if (score < min) min = score;
			else if (score > max) max = score;
		}
		return new double[] {min, max};
	}
	
	/**Returns the number of intervals that can be made at a given cutoff*/
	public int makeMockIntervals(double minScoreInt){
		//initialize params
		intervalsAl = new ArrayList(1000);
		lastLineReached = false;
		int numIntervals = 0;
		parPassWindows = 1;
		interval = null;
		testWindow = null;
		windowIndex = 0;
		numWindowsInInt = 0;
		try{
			numberOfWindows = windows.length;
			if (numberOfWindows >0 && windows[0]!=null){
				//set sortBy value to the one they want to build intervals with
				for (int i=0; i< numberOfWindows; i++) windows[i].setSortBy((float)windows[i].getScores()[scoreIndex]);
				//advance till good window found, make this the starting interval
				Window firstWin = fetchGoodWindow(minScoreInt);
				
				if (lastLineReached == false && firstWin!=null){
					interval = new Interval (firstWin,sizeOfOligo);
					boolean addIt;			
					
					for (; windowIndex<numberOfWindows; windowIndex++){
						testWindow = windows[windowIndex];			
						addIt=false;
						
						if (testWindow.getSortBy() >= minScoreInt && testWindow.getNumberOligos() >= minimumNumberOfOligos){
							parPassWindows++;
							//check chromosome
							if (interval.getChromosome().equals(testWindow.getChromosome())){
								addIt=true;
							}
							else{
								//different chromosome, save old interval, start new
								addIt();
								interval = new Interval (testWindow,sizeOfOligo);
								numWindowsInInt = 1;
							}					
						}
						
						if (addIt){
							//check length
							if ((testWindow.getStart1stOligo()-interval.getStartLastOligo()) <= maxGap){
								//length, score, N look good so meld
								interval.setStartLastOligo(testWindow.getStartLastOligo());
								numWindowsInInt++;						
								//Save windows 
								//fetch values for the ratio/ fold difference
								double intensityScoreTW = testWindow.getSortBy();
								double intensityScoreI = interval.getBestWindow().getSortBy();
								
								if (intensityScoreTW > intensityScoreI){ 
									interval.setBestWindow(testWindow);
								}
							}
							//length bad but good window, thus
							else {
								//write interval to disk
								addIt();	
								numWindowsInInt = 1;			
								//make new interval
								interval = new Interval (testWindow,sizeOfOligo);
							}
						}
					}	
					//close and save final interval, check to see that last line wasn't reached by fetchGoodInterval()
					if (lastLineReached==false) {
						addIt();		
					} 			
				}
				numIntervals = intervalsAl.size();
			}
		}catch (Exception e){
			System.out.println("\tSkipping "+file+"...");
			e.printStackTrace();
		}
		return numIntervals;
	}	
	public void addIt(){
		intervalsAl.add(interval);
	}		
	
	public Window fetchGoodWindow(double minScoreInt){
		for (; windowIndex< numberOfWindows; windowIndex++){
			if (windows[windowIndex].getSortBy() >= minScoreInt  && windows[windowIndex].getNumberOligos() >= minimumNumberOfOligos) {
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
	
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		String numInts = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': file = args[i+1]; i++; break;
					case 'o': minimumNumberOfOligos =Integer.parseInt(args[i+1]);i++; break;
					case 'z': sizeOfOligo =Integer.parseInt(args[i+1]);i++; break;
					case 'm': maxGap = Integer.parseInt(args[i+1]); i++; break;
					case 'i': flipScore = true; break;
					case 'a': makeIntervals = true; break;
					case 's': scoreIndex = Integer.parseInt(args[i+1]); i++; scanOnlyOneIndex = true; break;
					case 'n': numInts = args[i+1]; i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		if (Misc.isEmpty(file)){
			System.out.print("\nPlease enter a file or directory.\n");
			System.exit(0);
		}
		//get files to process
		File directory = new File (file);
		if (directory.isDirectory()){
			files = directory.listFiles();
			if 	(files == null || files.length==0){
				System.out.println("\nCannot find the directory or files in the directory?!\n\n");
				System.exit(0);	
			}
		}
		else files = new File[]{directory};
		//parse number of intervals
		String[] parsedInts = numInts.split(",");
		numberIntervals = Num.parseInts(parsedInts);
		
	}	
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Set Number Interval Maker: Feb 2007                        **\n" +
				"**************************************************************************************\n" +
				"SNIM determines the threshold needed to create a set number of intervals for each score\n" +
				"index.\n\n" +
				
				"-n Number of intervals, single value or comma delimited list.\n"+
				"-s Particular score index to use, default is to scan all.\n"+
				"-i Multiply indexed window scores by -1. Useful for looking for reduced regions.\n"+
				"-o Minimum number of oligos in each window, defaults to 10.\n"+
				"-z Size of oligo, defaults to 25.\n"+
				"-m Max gap between starts of oligos for collapsing windows when scanning intervals,\n" +
				"      defaults to 24.\n"+
				"-a Make intervals using found thresholds.\n"+
				"-f Full path file text for the Window[] array, if a directory is specified,\n" +
				"      all files will be processed.\n\n" +
				
				"Example: java -Xmx500M -jar pathTo/T2/Apps/SetNumberIntervalMaker -f\n" +
				"                 /affy/res/zeste.res -n 600,1200,24000 -a\n\n" +
				
		"**************************************************************************************\n");		
	}
	
	public static void main(String[] args) {
		if (args.length<2){
			printDocs();
			System.exit(0);	
		}
		new SetNumberIntervalMaker(args);
	}
	
}
