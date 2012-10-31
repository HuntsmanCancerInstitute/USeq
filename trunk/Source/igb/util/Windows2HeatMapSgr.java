package igb.util;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.*;

import util.gen.IO;


/**Converts a text file (columns: chrom, start, stop, score) into a heat/block map specific 
 * sgr file for import into IGB.*/
public class Windows2HeatMapSgr {
	//fields
	private Window[] windows;
	private int numWindows; 
	private double lowestScore;
	private double score;
	private boolean setCommonScore = false;
	private File[] windowFiles;
	
	public Windows2HeatMapSgr(String[] args){
		processArgs(args);
		
		for (int i=0; i< windowFiles.length; i++){
		//parse txt file building Window[]
		System.out.println("\nCreating Window[]...");
		parseWindowFile(windowFiles[i]);
		
		//sort windows by chrom then position then length, shortest first.
		Arrays.sort(windows);
		
		//make heat map blocks
		System.out.println ("Making blocks...");
		File sgr = new File (windowFiles[i]+".sgr");
		makeBlockMapSgrFile(sgr);
		}
		
		System.out.println("Done!\n");
	}
	
	/**Takes an array of Window */
	public void makeBlockMapSgrFile(File file){
		numWindows = windows.length;	
		
		//find start stop indexes for chromosomes, also lowest score
		ArrayList startStopIndexes = new ArrayList();
		ArrayList chromosomes = new ArrayList();
		int start = 0;
		String chromosome = windows[0].chromosome;
		lowestScore = windows[0].score;
		chromosomes.add(chromosome);
		for (int x=1; x<numWindows; x++){
			if (chromosome.equals(windows[x].chromosome) == false){
				startStopIndexes.add(new int[]{start, x});
				start =x;
				chromosome = windows[x].chromosome;
				chromosomes.add(chromosome);
			}
			if (windows[x].score <lowestScore)lowestScore = windows[x].score; 
		}
		//if setting a common score or all values are positive
		if (setCommonScore || lowestScore > 0) lowestScore = 0;
		
		//add last one
		startStopIndexes.add(new int[]{start, numWindows});
		
		//for each chromosome
		int numChromosomes = startStopIndexes.size();
		try {
			PrintWriter out = new PrintWriter( new FileWriter(file));
			for (int x=0; x<numChromosomes; x++){
				chromosome = (String)chromosomes.get(x);
				int[] startStop = (int[])startStopIndexes.get(x);
				//assign 1st score and 1st position
				int base = windows[startStop[0]].start;
				out.println(chromosome+"\t0\t" +lowestScore);
				out.println(chromosome+"\t"+(base -1)+"\t" + lowestScore);
				double score = windows[startStop[0]].score;
				out.println(chromosome+"\t"+ base+ "\t"+ score);
				
				//for each window in the chromosome 
				for (int y = startStop[0]; y< startStop[1]; y++){
					//for each base in the window
					int endBase = windows[y].stop + 1;
					for (int z=base; z<endBase; z++){
						double highestScore = highestScoringWindow(z, y);
						if (highestScore != score){
							//close old block start new block
							out.println(chromosome+"\t"+ (z-1)+ "\t"+ score);
							out.println(chromosome+"\t"+ z + "\t"+ highestScore);
							score = highestScore;
						}
					}
					//is this the last window in the chromosome?
					if (startStop[1] == (y+1)){
						//close old block
						out.println(chromosome+"\t"+ (endBase-1) + "\t"+ score);
						out.println(chromosome+"\t"+ endBase + "\t"+ lowestScore);
					}
					//Does the next base overlap next window?
					else if (baseOverlapsWindow(endBase, windows[y+1]) == false){
						//close old block
						out.println(chromosome+"\t"+ (endBase-1) + "\t"+ score);
						out.println(chromosome+"\t"+ endBase + "\t" +lowestScore);
						//reset score and base
						score = windows[y+1].score;
						base = windows[y+1].start;
						//start new block
						out.println(chromosome+"\t"+ (base-1) + "\t"+ lowestScore);
						out.println(chromosome+"\t"+ base + "\t"+ score);
					}
					else{
						base = endBase;
					}
				}
			}	
			out.close();
		}catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/**Looks forward at overlapping windows, returns highest score for the given bp position.*/
	public double highestScoringWindow (int base, int currentWindowIndex){
		int num = windows.length;
		double score = windows[currentWindowIndex].score;
		for (int i=currentWindowIndex; i<num; i++){
			if (baseOverlapsWindow(base,windows[i])) {
				//check score
				double testScore = windows[i].score;
				if (testScore> score) score = testScore;				
			}
			else return score;
		}
		return score;
	}
	
	/**Is a given base contained within the window?*/
	public static boolean baseOverlapsWindow(int bp, Window win){
		if (bp >= win.start && bp <= (win.stop)) return true;
		return false;
	}
	public static void main(String[] args){
		//look for file
		if (args.length == 0) {
			printDocs();
			System.exit(0);
		}
		new Windows2HeatMapSgr(args);
	}
	
	public void parseWindowFile(File file){
		try{
			BufferedReader in = new BufferedReader(new FileReader(file));
			String line;
			String[] tokens;
			ArrayList regionsAL = new ArrayList();
			String chromosome= null;
			int start =0;
			int stop =0;
			double parsedScore =0;
			//chrom, start, stop, score
			//  0      1      2     3
			while ((line = in.readLine()) !=null) {
				line = line.trim();
				if (line.equals("")) continue;
				tokens = line.split("\\s+");
				if (tokens.length < 4){
					System.out.println("\nError: the following line does not contain four columns -> "+line+"\n" +
					"Be sure every line contains chrom, start, stop, score.\n");
					System.exit(0);
				}
				chromosome = tokens[0];
				try {
					start = Integer.parseInt(tokens[1]);
					stop = Integer.parseInt(tokens[2]);
					if (setCommonScore) parsedScore = score;
					else parsedScore = Double.parseDouble(tokens[3]);
				} catch (Exception e){
					System.out.println("\nError: the following line contains a malformed number -> "+line+"\n" +
					"Be sure every line contains chrom, start(integer), stop(integer), score(double).\n");
					System.exit(0);
				}
				//check start and stop
				if ((stop-start) <=0){
					System.out.println("\nError: the start must be bigger than the stop position -> "+line+"\n");
					System.exit(0);
				}
				
				regionsAL.add(new Window(chromosome, start, stop, parsedScore));
			}
			windows = new Window[regionsAL.size()];
			regionsAL.toArray(windows);
		}catch (IOException e){
			e.printStackTrace();
		}		
	}
	
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		File dirFile = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': dirFile = new File(args[i+1]); i++; break;
					case 's': score =Double.parseDouble(args[i+1]); i++; setCommonScore = true; break;
					case 'h': printDocs(); System.exit(0);
					default: {
						System.out.println("\nError: unknown option! -> " + mat.group()+"\n");
						System.exit(0);
					}
					}
				}
				catch (Exception e){
					System.out.print("\nError: something doesn't look right with this parameter request: -"+test+"\n");
					System.exit(0);
				}
			}
		}
		//check for tempDirectory
		if (dirFile == null || dirFile.exists()== false){
			System.out.println("\nError: cannot find your window file(s) or directory?! -> "+dirFile+"\n");
			System.exit(0);
		}
		windowFiles = IO.extractFiles(dirFile);
	}
	
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Windows 2 Heat Map: Jan  2006                           **\n" +
				"**************************************************************************************\n" +
				"W2HM converts a list of potentially overlapping windows into an sgr file for\n" +
				"visualization in IGB.\n\n"+
				
				"-f Full path to tab delimited text file or directory containing windows\n" +
				"      (chrom start stop score).\n" +
				"-s Score to assign to all windows, default is to use the given score.\n\n"+
				
				"Example: java -Xmx512M -jar pathTo/T2/Apps/Windows2HeatMap -f /data/winFiles/ -s 100\n\n"+
				
				"Questions, comments, suggestions? Contact Gingeras Group or David_Nix@Affymetrix.com\n" +
		"**************************************************************************************\n");
	}
	
	private class Window implements Comparable{
		//fields
		String chromosome;
		int start;
		int stop;
		double score;
		//constructor
		public Window(String chromosome, int start, int stop, double score){
			this.chromosome = chromosome;
			this.start = start;
			this.stop = stop;
			this.score = score;
		}
		//methods
		/**Compares by chromosome-> start position -> shortest length.*/
		public int compareTo(Object obj){
			Window other = (Window)obj;
			int i = other.chromosome.compareTo(chromosome);
			if (i !=0 ) return -1*i;
			if (other.start< start) return 1;
			if (other.start> start) return -1;
			if ( (other.stop - other.start) < (stop - start) ) return 1;
			if ( (other.stop - other.start) > (stop - start) ) return -1;
			return 0;
		}
	}
}

