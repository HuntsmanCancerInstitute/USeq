package trans.main;
import java.util.regex.*;
import java.util.*;
import util.gen.*;
import java.io.*;

public class MergeWindowArrays {
	
	private File[] windowFiles;
	private int thresholdScoreIndex = 1;
	private double threshold = 0.2;
	private boolean flipScores = false;
	private File mergedFile = null;
	private boolean restrictToWin = false;
	
	public MergeWindowArrays(String[] args){
		processArgs(args);
		System.out.println("\nLaunched:");
		if (flipScores) System.out.println("\tMultiplying all window scores by -1");
	
		//fetch Window[]s
		Window[][] windows = new Window[windowFiles.length][];
		for (int i=0; i< windowFiles.length; i++){
			System.out.println("\tLoading "+windowFiles[i].getName());
			windows[i] = (Window[])IO.fetchObject(windowFiles[i]);
			if (flipScores) IntervalMaker.invertScores(windows[i]);
			windows[i] = filter( windows[i] );
		}
		
		//collapse
		int num = Num.countObjects(windows);
		Window[] win = new Window[num];
		int k=0;
		//for each file
		for (int i=0; i< windows.length; i++){
			//for each value
			for (int j=0; j< windows[i].length; j++){
				win[k++] = windows[i][j];
			}
		}
		
		//sort!
		System.out.println("\tSorting...");
		Arrays.sort(win, new WindowComparator());
		
		//save
		if (mergedFile == null) mergedFile = new File (windowFiles[0].getParentFile(), "mergedWindows");
		System.out.println("\n\tSaving "+num+ " windows -> "+mergedFile);
		IO.saveObject(mergedFile, win); 
	}
	
	public Window[] filter(Window[] win){
		ArrayList al = new ArrayList(win.length/2);
		for (int i=0; i< win.length; i++){
			double score = win[i].getScores()[thresholdScoreIndex];
			if (score >= threshold) al.add(win[i]);
		}
		win = new Window[al.size()];
		al.toArray(win);
		return win;
	}

	public static void main(String[] args) {
		if (args.length<1){
			printDocs();
			System.exit(0);	
		}
		new MergeWindowArrays(args);
	}
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		File file = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': file = new File(args[i+1]); i++; break;
					case 'i': thresholdScoreIndex = Integer.parseInt(args[i+1]); i++; break;
					case 't': threshold = Double.parseDouble(args[i+1]); i++; break;
					case 'n': mergedFile = new File(args[i+1]); i++; break;
					case 'm': flipScores = true; break;
					case 'r': restrictToWin = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (file ==null || file.isDirectory() == false){
			Misc.printExit("\nPlease enter a directory containing window arrays.\n");
		}
		if (restrictToWin) windowFiles = IO.extractFiles(file, "_Win");
		else windowFiles = IO.extractFiles(file);
	}	
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         Merge Window Arrays: March 2007                          **\n" +
				"**************************************************************************************\n" +
				"Concatinates serialized Window[]s into one.\n" +
				
				"-f Full path directory text for the folder containing eserialized Window[]s\n" +
				"      from the ScanChromosomesCNV or ScanChip app.\n" +
				"-r Restrict files for merging to those ending in '_Win', defaults to all.\n"+
				"-n (Optional) Full path file text to use in saving the merged windows array.\n"+
				"-m Multiply window scores by -1. Useful for looking for reduced regions.\n"+
				"-t Score threshold, tosses any windows with a score < threshold, defaults to 0.2\n"+
				"-i Score index, defaults to 1.  See ScanChromosome for index descriptions.\n\n"+
				
				"Example: java -Xmx5000M -jar /YourPathTo/T2/Apps/MergeWindowArrays -f /affy/win/ \n" +

		"**************************************************************************************\n");		
	}

}
