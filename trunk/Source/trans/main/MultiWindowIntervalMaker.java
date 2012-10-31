package trans.main;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.*;
import util.gen.*;


/**
 *	Collapses arrays of Windows into an array of Interval based on scores and minimal number of windows that must pass.
 */
public class MultiWindowIntervalMaker {
	//fields
	private double[] minScores;
	private int maxGap = 24;
	private int minNumWindows = 0;
	private int minNumOfOligos = 10;
	private Interval[] intervals;
	private ArrayList intervalsAl;
	private int scoreIndex = 1;
	private boolean lastLineReached;
	private int parPassWindows;
	private Interval interval;
	private Window testWindow;
	private Window[][] windows;
	private int numberOfWindows;
	private int windowIndex;
	private int numWindowsInInt;	
	private int numberOfWindowArrays;
	private File intervalDir;
	private String compositeName = "compositeInt";
	private int sizeOfOligo= 25;
	

	
	//contructor
	public MultiWindowIntervalMaker(String[] args){
		processArgs(args);
		
		//initialize params
		intervalsAl = new ArrayList(1000);
		lastLineReached = false;
		parPassWindows = 0;
		interval = null;
		testWindow = null;
		windowIndex = 0;
		numWindowsInInt = 0;
		numberOfWindows = windows[0].length;

		//advance till good representative window is found, make this the starting interval
		Window firstWin = fetchGoodWindow();
		
		if (lastLineReached == false && firstWin!=null){
			interval = new Interval (firstWin, sizeOfOligo);
			boolean addIt;		
			parPassWindows++;
			
			//advance one window index at a time
			for (; windowIndex<numberOfWindows; windowIndex++){	
				addIt=false;
				
				//check scores and min num pass
				testWindow = checkWindowScoresAndNumberOfOligos();
				if (testWindow != null){
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
				
				//scores are good check length.
				if (addIt){
					//check length
					if ((testWindow.getStart1stOligo()-interval.getStartLastOligo()) <= maxGap){
						//length, score, look good so meld
						interval.setStartLastOligo(testWindow.getStartLastOligo());
						numWindowsInInt++;						
						//Save window as best representative window in interval?
						if (testWindow.getScores()[scoreIndex]>interval.getBestWindow().getScores()[scoreIndex]){
							interval.setBestWindow(testWindow);
						}
							
					}
					//length bad but good window, thus
					else {
						//save interval, start new
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
		
		System.out.println("\tMax Gap: "+maxGap+" Min Scores: "+Num.doubleArrayToStringOnlyMax(minScores, 4, " "));
		System.out.println("\tMax Gap: "+maxGap+", Score Index: "+scoreIndex+", Min Scores: "+Num.doubleArrayToStringOnlyMax(minScores, 4, " ")+ ", Min # Oligos: "+minNumOfOligos);

		int numInt = intervalsAl.size();
		System.out.println("\t"+numInt+" Intervals were derived from "+parPassWindows+" windows that passed " +
				"the score cutOffs. "+ windows[0].length+ " windows were screened.");
		if (numInt !=0){
			//convert ArrayList to Interval[] and save to disk
			intervals = new Interval[intervalsAl.size()];
			intervalsAl.toArray(intervals);
			//Save file
			File intervalFile = new File(intervalDir, compositeName+numInt);
			System.out.print ("\tSaving ..."+intervalFile+"\n");
			IO.saveObject(intervalFile, intervals);
		}
	}
	
	public void addIt(){
		interval.setNumberOfWindows(numWindowsInInt);
		intervalsAl.add(interval);
	}		
	
	public Window fetchGoodWindow(){
		//advance one at a time
		for (; windowIndex< numberOfWindows; windowIndex++){
			//check all Window[] at a particular index to see if minNumWindows and minScore cutoffs are reached
			Window w = checkWindowScoresAndNumberOfOligos();
			if (w !=null) {
				numWindowsInInt = 1;
				windowIndex++;
				return w; 
			} 
		}
		//ran thru stop of array add last if good
		if (interval!=null) {
			addIt();
		} 
		lastLineReached = true;
		return null;	
	}
	
	/**Checks the scores and number of oligos in each window. 
	 * Returns best one if they pass the minNumWindows and minScore cutoffs.
	 * Otherwise null is returned.*/
	public Window checkWindowScoresAndNumberOfOligos(){
		//check number of oligos
		if (windows[0][windowIndex].getNumberOligos() < minNumOfOligos) return null;
		
		//check scores of each window
		int numPass = 0;
		Window bestWindow = null;
		double bestWindowScore = -10000;
		//for each window array, not each window
		for (int x=0; x<numberOfWindowArrays; x++){
			double score = windows[x][windowIndex].getScores()[scoreIndex];
			if (score >= minScores[x]){
				numPass++;
				//save best window 
				if (bestWindow == null) {
					bestWindow = windows[x][windowIndex]; 
					bestWindowScore = score;
				}
				//best window already set, check to see if score of test is better
				else {
					if (score> bestWindowScore){
						bestWindow = windows[x][windowIndex]; 
						bestWindowScore = score;
					}
				}
			}
		}
		if (numPass >= minNumWindows)return bestWindow;
		return null;
	}
	
	public Window findWindow(Window[] win){
		//get start bp
		int startBp = windows[0][windowIndex].getStart1stOligo();
		//look up
		System.out.println("\nStartBp "+startBp+"\n\tUp: ");
		for (int i=windowIndex; i>=0; i--){
			int startTestWin = win[i].getStart1stOligo();
			System.out.print(startTestWin+ " ");
			if ( startTestWin == startBp) return win[i];
			if ( startTestWin < startBp) break;
		}
		System.out.print("\n\tDwn: ");
		//does windowIndex exceed length of win?
		if (windowIndex >= win.length) {
			System.out.println("\tExceeds windowIndex "+windowIndex+" win.length "+win.length);
			return null;
		}
		//look down
		for (int i=windowIndex; i<win.length; i++){
			int startTestWin = win[i].getStart1stOligo();
			System.out.print(startTestWin+ " ");
			if ( startTestWin == startBp) return win[i];
			if ( startTestWin > startBp) break;
		}
		System.out.println("\n\tWin not found windowIndex "+windowIndex);
		return null;
	}
	
	public static void main(String[] args) {
		if (args.length<1){
			printDocs();
			System.exit(0);	
		}
		new MultiWindowIntervalMaker(args);
	}
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		String minScoreString = null;
		String filesString = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 's': minScoreString =args[i+1];i++; break;
					case 'n': compositeName =args[i+1];i++; break;
					case 'i': scoreIndex =Integer.parseInt(args[i+1]);i++; break;
					case 'm': minNumWindows = Integer.parseInt(args[i+1]);i++; break;
					case 'o': minNumOfOligos =Integer.parseInt(args[i+1]);i++; break;
					case 'z': sizeOfOligo =Integer.parseInt(args[i+1]);i++; break;
					case 'g': maxGap =Integer.parseInt(args[i+1]);i++; break;
					case 'f': filesString = args[i+1]; i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (minScoreString == null || filesString == null){
			Misc.printExit("\nError: you must enter both a list of scores/thresholds and Window[] files\n");
		}
		
		String[] scores = minScoreString.split(",");
		String[] fileSplit = filesString.split(",");
		if (scores.length != fileSplit.length){
			Misc.printExit("\nError: the number of scores/thresholds and number of files do not match!\n");
		}
		numberOfWindowArrays = fileSplit.length;
		
		//fetch window[]s, parse scores
		windows = new Window[numberOfWindowArrays][];
		minScores = new double[numberOfWindowArrays];
		System.out.println();
		for (int i=0; i<numberOfWindowArrays; i++){
			System.out.println("Loading -> "+fileSplit[i]);
			windows[i] = (Window[])IO.fetchObject(new File (fileSplit[i]));
			minScores[i] = Double.parseDouble(scores[i]);
		}
		
		intervalDir = new File(fileSplit[0]).getParentFile();
		
		//set minimum number of windows that must pass
		if (minNumWindows == 0 || minNumWindows > numberOfWindowArrays) minNumWindows = numberOfWindowArrays;		
		
	}
	
	/**Puts longest window array first and it's associated minScore.*/
	public void putLargestWindowArrayFirst (){
		//make arraylists
		ArrayList wins = new ArrayList(windows.length);
		ArrayList scores = new ArrayList(minScores.length);
		//add arrays, find index of longest array
		int longestIndex = 0;
		int number = 0;
		for (int i=0; i< windows.length; i++){
			wins.add(windows[i]);
			scores.add(new Double(minScores[i]));
			if (windows[i].length > number){
				number = windows[i].length;
				longestIndex = i;
			}
		}
		//grab longest
		Window[] longestWin = windows[longestIndex];
		double longestMinScore = minScores[longestIndex];
		//remove longest from arrayLists
		wins.remove(longestIndex);
		scores.remove(longestIndex);
		//add longest on front
		wins.add(0, longestWin);
		scores.add(0, new Double(longestMinScore));
		//convert to arrays
		windows = new Window[windows.length][];
		minScores = new double[windows.length];
		for (int i=0; i< windows.length; i++){
			windows[i] = (Window[])wins.get(i);
			minScores[i] = ((Double)wins.get(i)).doubleValue();
		}
		
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Multi Window Interval Maker: July 2006                     **\n" +
				"**************************************************************************************\n" +
				"MWIM combines Windows into larger Intervals given minimal score(s), a maximum gap, a\n" +
				"minimum number of windows per index that must pass, and a minimum number of oligos in\n" +
				"each window.  For each index, the best window is used to represent the Interval. Must\n" +
				"have all the windows from ScanChip or ScanChromosome, no pruning.\n\n" +
				
				"Options:\n\n" +
				
				"-i Score index, the score you wish to use in merging, see ScanChipNoPerm or \n" +
				"      ScanChromosome.\n"+
				"-o Minimum number of oligo positions in each window, defaults to 10.\n"+
				"-g Maximum allowable bp gap between starts of oligos, defaults to 24.\n" +
				"-z Size of oligo, defaults to 25.\n"+
				"-s Minimal score(s), comma delimited, no spaces, one for each Window[] file.\n" +
				"-f Full path file names for the Window[]s, comma delimited, no spaces.\n"+
				"-n Composite text to use as a base in saving.\n"+
				"-m For a given window index, minimum number of windows from the different arrays that\n" +
				"      must pass to be included in building an interval. Defaults to all.\n\n" +
				
				"Example: java -Xmx1500M -jar pathTo/T2/Apps/MultiWindowIntervalMaker -f \n" +
				"      /affy/win1,/affy/win2,/affy/win3 -i 1 -g -100 -s 0.5,0.5,0.3 -m 2\n\n"+
				
		"**************************************************************************************\n");		
	}
	
}
