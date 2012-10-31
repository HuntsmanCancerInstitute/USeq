package trans.main;
import java.io.*;
import java.util.*;
import java.util.regex.*;
import trans.tpmap.*;
import util.gen.*;

/**
 * Makes all windows given oligo positions.  Useful for calculating coverage of a genome. 
 */
public class MakeAllWindows {
	
	//fields
	private File resultsDirectory;
	private File[] chromosomeOligoPositions = null;
	private Window[] windows;
	private int numWindows;
	private int minOligos = 10;
	private int maxOligoMeasurements = 50;
	private int windowSize = 675;
	private 	String chromosome = null;
	private 	int[] positions;
	private int[][] testWindows;
	private File windowDirectory;
	private String split = null;
	private int splitUnit = 0;
	private int splitBy = 0;
	
	public MakeAllWindows(String[] args){
		//process args
		processArgs(args);
		
		//make directories to hold point gr files and temp ratios when so indicated
		makeDirectories();
		
		//make a WindowMaker to make windows for each chromosome from oligo positions
		WindowMaker windowMaker = new WindowMaker(windowSize, minOligos, maxOligoMeasurements);
		
		//for each chromosome	
		System.out.println();
		for (int i=0; i<chromosomeOligoPositions.length; i++){
			ArrayList windowsAl = new ArrayList(10000);
			chromosome = chromosomeOligoPositions[i].getName();
			System.out.println("Testing chromosome: "+chromosome);
			
			//get positions
			positions = (int[])IO.fetchObject(chromosomeOligoPositions[i]);
			
			//make windows int[windowNumber][start index, stop index]
			testWindows = windowMaker.makeWindows(positions);
			if (testWindows.length == 0) {
				System.out.println("No Windows found, skipping.");
				continue;
			}
			else {
				System.out.println("\t"+testWindows.length+" "+windowSize+"bp windows found with a minimum of "+minOligos+" unique oligo positions and maximum of "+maxOligoMeasurements+" oligo measurements.");
			}
			
			//run the window tests
			numWindows = testWindows.length;
			
			for (int j=0; j<numWindows; j++){
				//get transformed values to test
				int[] startStop = testWindows[j];
				int windowSize = 1+startStop[1]-startStop[0];
				
				//pseudo median ratio?
				double score = 0;
				
				//chromosome, bpStart, bpStop, numberOligos, scores
				Window win = new Window(chromosome, positions[startStop[0]],positions[startStop[1]], windowSize, new double[]{score});
				windowsAl.add(win);
			}
			windows = new Window[windowsAl.size()];
			windowsAl.toArray(windows);
			//save Window[]
			File serWinIntArray = new File(windowDirectory, chromosome);
			System.out.print ("Saving "+numWindows+" windows "+ serWinIntArray+"\n");
			IO.saveObject(serWinIntArray, windows);
		}
	}
	/**Makes directories and files to hold point gr files and tempRatios*/
	public void makeDirectories(){
		windowDirectory = new File(resultsDirectory, "Win");
		if (windowDirectory.exists() == false) windowDirectory.mkdir();
	}
	
	
	
	public static void main(String[] args) {
		if (args.length<4){
			printDocs();
			System.exit(0);
		}	
		new MakeAllWindows(args);
	}		
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		File oligoPositions = null;
		String resultsDirectoryString = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'r': resultsDirectoryString = args[i+1]; i++; break;
					case 'm': minOligos = Integer.parseInt(args[i+1]); i++; break;
					case 'x': maxOligoMeasurements = Integer.parseInt(args[i+1]); i++; break;
					case 'w': windowSize = Integer.parseInt(args[i+1]); i++; break;
					case 'o': oligoPositions = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//look for required parameters
		if (resultsDirectoryString == null || oligoPositions == null){
			Misc.printExit("\nPlease complete one or more of the following required parameters: -t, -o, or -r .\n");
		}
		
		chromosomeOligoPositions = IO.extractFiles(oligoPositions);
		if (chromosomeOligoPositions == null) Misc.printExit("\nProblem parsing oligo positions -> "+oligoPositions);
		
		//make and check results directory
		if (split!=null) {
			String[] x = split.split("-");
			if (x.length !=2) Misc.printExit("\nCannot parse your split chromosome(s) processing request?! Check -> "+split+"\n");
			splitBy = Integer.parseInt(x[0]);
			splitUnit = Integer.parseInt(x[1]);
			resultsDirectory = new File(resultsDirectoryString+"_"+splitUnit+"of"+splitBy);
			System.out.println("Splitting chromosome(s) by "+splitBy+", processing part "+splitUnit+".");
		}
		else resultsDirectory = new File(resultsDirectoryString);
		if (resultsDirectory.exists() == false) resultsDirectory.mkdir();
		
	}	
	
	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Make All Windows: May 2006                             **\n" +
				"**************************************************************************************\n" +
				
				"-o The 'OligoPositions' directory, full path, generated by ChipSetCelProcessor.\n" +
				"-r The full path text to use in creating a directory to save the results.\n" +
				"-m Minimal number of unique oligo positions required in each window. Defaults to\n" +
				"      10. Set to 0 to save all windows.\n" +
				"-x Max number of oligo measurements in each window. Defaults to 50. Guards against\n" +
				"      repetatively tiled oligos stalling the pseudo median calculation.\n"+
				"-w Window size, defaults to 675bp.\n"+

		"**************************************************************************************\n");		
	}
}
