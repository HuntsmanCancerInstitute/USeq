package trans.main;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;
import trans.misc.*;
import java.io.*;

/**Converts a Window[] into a heat map specific sgr file for import into IGB. 
 * Assumes windows are sorted by chromosome and start position
 * Score index is the index to use when fetching a score for a window from it's double[] of scores
 * see ScanChip or ScanChromosomesCNV for current score indexes.*/
public class WindowBlockMaker {
	//fields
	private int scoreIndex = 0;	
	private Window[] windows;
	private int numWindows; 
	private float baseScore = 0;
	private boolean flipScore = false;
	private File windowFile;
	private int sizeOfOligo = 25;
	private int sizeOfOligoMinusOne = 24;
	private String genomeVersion = null;
	
	//constructors
	public WindowBlockMaker(int sizeOfOligo){
		sizeOfOligoMinusOne = sizeOfOligo - 1;
	}
	
	public WindowBlockMaker(int sizeOfOligo, String genomeVersion){
		this.genomeVersion = genomeVersion;
		sizeOfOligoMinusOne = sizeOfOligo - 1;
	}
	
	public WindowBlockMaker(String[] args){
		processArgs(args);
		sizeOfOligoMinusOne = sizeOfOligo - 1;
		
		//load Window[]
		System.out.println("\nReading in Window[]...");
		windows = (Window[])IO.fetchObject(windowFile);
		numWindows = windows.length;

		//multiply window score by -1?
		if (flipScore){
			System.out.println("Multiplying indexed scores by -1...");
			IntervalMaker.invertScores(windows);
		}
		System.out.println ("Making blocks...");
		makeHeatMapSgrFile(windows, new File(windowFile+"HM.sgr.zip"), true);
		System.out.println("Done!\n");
		
		/*Testing
		Window[] windows = new Window[5];
		windows[0] = new Window("chr1", 10, 20, new double[]{20}); //chromosome, int start, int stop, double[] scores
		windows[1] = new Window("chr1", 25, 35, new double[]{15});
		windows[2] = new Window("chr1", 30, 40, new double[]{0});
		windows[3] = new Window("chr1", 50, 60, new double[]{5});
		windows[4] = new Window("chr1", 70, 80, new double[]{10});
		bm.makeHeatMapSgrFile(windows, new File("delMeBM.sgr"));*/
	}
	
	/**Be sure to set the scoreIndex, defaults to 0. Automatically zips them*/
	public void makeHeatMapGrFiles(Window[] windows, File directory){
		this.windows = windows;
		numWindows = windows.length;	
		if (numWindows < 1) return;
		//find start stop indexes for chromosomes
		ArrayList startStopIndexes = new ArrayList();
		ArrayList chromosomes = new ArrayList();
		int start = 0;
		String chromosome = windows[0].getChromosome();
		chromosomes.add(chromosome);
		for (int x=1; x<numWindows; x++){
			if (chromosome.equals(windows[x].getChromosome()) == false){
				startStopIndexes.add(new int[]{start, x});
				start =x;
				chromosome = windows[x].getChromosome();
				chromosomes.add(chromosome);
			}
		}
		//add last one
		startStopIndexes.add(new int[]{start, numWindows});
		
		//for each chromosome
		int numChromosomes = startStopIndexes.size();
		try {
			
			for (int x=0; x<numChromosomes; x++){
				chromosome = (String)chromosomes.get(x);
				File gr = new File(directory,chromosome+".gr");
				PrintWriter out = new PrintWriter( new FileWriter(gr));
				int[] startStop = (int[])startStopIndexes.get(x);
				//assign 1st score and 1st position
				int base = windows[startStop[0]].getStart1stOligo();
				out.println("0\t" +baseScore);
				out.println((base -1)+"\t" + baseScore);
				double score = windows[startStop[0]].getScores()[scoreIndex];
				out.println(base+ "\t"+ score);
				
				//for each window in the chromosome 
				for (int y = startStop[0]; y< startStop[1]; y++){
					//for each base in the window
					int endBase = windows[y].getStartLastOligo()+sizeOfOligo;
					for (int z=base; z<endBase; z++){
						double highestScore = highestScoringWindow(z, y);
						if (highestScore != score){
							//close old block start new block
							out.println((z-1)+ "\t"+ score);
							out.println(z + "\t"+ highestScore);
							score = highestScore;
						}
					}
					//is this the last window in the chromosome?
					if (startStop[1] == (y+1)){
						//close old block
						out.println((endBase-1) + "\t"+ score);
						out.println(endBase + "\t"+ baseScore);
					}
					//Does the next base overlap next window?
					else if (baseOverlapsWindow(endBase, windows[y+1]) == false){
						//close old block
						out.println((endBase-1) + "\t"+ score);
						out.println(endBase + "\t" +baseScore);
						//reset score and base
						score = windows[y+1].getScores()[scoreIndex];
						base = windows[y+1].getStart1stOligo();
						//start new block
						out.println((base-1) + "\t"+ baseScore);
						out.println(base + "\t"+ score);
					}
					else{
						base = endBase;
					}
				}
				out.close();
				IO.zipAndDelete(gr);
			}	

		}catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/**Be sure to set the scoreIndex, defaults to 0.*/
	public void makeHeatMapSgrFile(Window[] windows, File file, boolean zipCompressIt){
		this.windows = windows;
		numWindows = windows.length;	
		if (numWindows < 1) return;
		//make tempFile?
		File zipFile = null;
		if (zipCompressIt){
			zipFile = file;
			file = new File(IO.getFullPathName(zipFile)+"tmpFile");
		}
		
		//find start stop indexes for chromosomes, also lowest score
		ArrayList startStopIndexes = new ArrayList();
		ArrayList chromosomes = new ArrayList();
		int start = 0;
		String chromosome = windows[0].getChromosome();
		chromosomes.add(chromosome);
		for (int x=1; x<numWindows; x++){
			if (chromosome.equals(windows[x].getChromosome()) == false){
				startStopIndexes.add(new int[]{start, x});
				start =x;
				chromosome = windows[x].getChromosome();
				chromosomes.add(chromosome);
			}
		}

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
				int base = windows[startStop[0]].getStart1stOligo();
				out.println(chromosome+"\t0\t" +baseScore);
				out.println(chromosome+"\t"+(base -1)+"\t" + baseScore);
				double score = windows[startStop[0]].getScores()[scoreIndex];
				out.println(chromosome+"\t"+ base+ "\t"+ score);
				
				//for each window in the chromosome 
				for (int y = startStop[0]; y< startStop[1]; y++){
					//for each base in the window
					int endBase = windows[y].getStartLastOligo()+sizeOfOligo;
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
						out.println(chromosome+"\t"+ endBase + "\t"+ baseScore);
					}
					//Does the next base overlap next window?
					else if (baseOverlapsWindow(endBase, windows[y+1]) == false){
						//close old block
						out.println(chromosome+"\t"+ (endBase-1) + "\t"+ score);
						out.println(chromosome+"\t"+ endBase + "\t" +baseScore);
						//reset score and base
						score = windows[y+1].getScores()[scoreIndex];
						base = windows[y+1].getStart1stOligo();
						//start new block
						out.println(chromosome+"\t"+ (base-1) + "\t"+ baseScore);
						out.println(chromosome+"\t"+ base + "\t"+ score);
					}
					else{
						base = endBase;
					}
				}
			}	
			out.close();
			if (zipCompressIt){
				IO.zip(file, zipFile);
				file.delete();
			}
		}catch (Exception e){
			e.printStackTrace();
		}
	}
	
	/**Be sure to set the scoreIndex, defaults to 0. Automatically zips them*/
	public void makeHeatMapBarFiles(Window[] windows, File directory){
		if (genomeVersion == null) Misc.printExit("\nError: genomeVersion is not set, cannot launch WindowBlockMaker.makeHeatMapBarFiles()\n");
		this.windows = windows;
		numWindows = windows.length;	
		if (numWindows < 1) return;
		//find start stop indexes for chromosomes
		ArrayList startStopIndexes = new ArrayList();
		ArrayList chromosomes = new ArrayList();
		int start = 0;
		String chromosome = windows[0].getChromosome();
		chromosomes.add(chromosome);
		for (int x=1; x<numWindows; x++){
			if (chromosome.equals(windows[x].getChromosome()) == false){
				startStopIndexes.add(new int[]{start, x});
				start =x;
				chromosome = windows[x].getChromosome();
				chromosomes.add(chromosome);
			}
		}
		//add last one
		startStopIndexes.add(new int[]{start, numWindows});
		
		//for each chromosome
		int numChromosomes = startStopIndexes.size();
		
		for (int x=0; x<numChromosomes; x++){
			chromosome = (String)chromosomes.get(x);
			File barFile = new File(directory,chromosome+".bar");
			int[] startStop = (int[])startStopIndexes.get(x);
			int numberWindows = (startStop[1]-startStop[0])/4;
			ArrayList positions = new ArrayList(numberWindows);
			ArrayList values = new ArrayList(numberWindows);
			//assign 1st score and 1st position
			positions.add(new Integer(0));
			values.add(new Float(baseScore));
			int base = windows[startStop[0]].getStart1stOligo();
			if (base == 0) base = 1;
			positions.add(new Integer((base -1)));
			values.add(new Float(baseScore));
			double score = windows[startStop[0]].getScores()[scoreIndex];
			positions.add(new Integer((base)));
			values.add(new Float(score));
			
			//for each window in the chromosome 
			for (int y = startStop[0]; y< startStop[1]; y++){
				//for each base in the window
				int endBase = windows[y].getStartLastOligo()+sizeOfOligo;
				for (int z=base; z<endBase; z++){
					double highestScore = highestScoringWindow(z, y);
					if (highestScore != score){
						//close old block start new block
						positions.add(new Integer((z-1)));
						values.add(new Float(score));
						positions.add(new Integer((z)));
						values.add(new Float(highestScore));
						score = highestScore;
					}
				}
				//is this the last window in the chromosome?
				if (startStop[1] == (y+1)){
					//close old block
					positions.add(new Integer((endBase-1)));
					values.add(new Float(score));
					positions.add(new Integer((endBase)));
					values.add(new Float(baseScore));
				}
				//Does the next base overlap next window?
				else if (baseOverlapsWindow(endBase, windows[y+1]) == false){
					//close old block
					positions.add(new Integer((endBase-1)));
					values.add(new Float(score));
					positions.add(new Integer((endBase)));
					values.add(new Float(baseScore));
					//reset score and base
					score = windows[y+1].getScores()[scoreIndex];
					base = windows[y+1].getStart1stOligo();
					//start new block
					positions.add(new Integer((base-1)));
					values.add(new Float(baseScore));
					positions.add(new Integer((base)));
					values.add(new Float(score));
				}
				else{
					base = endBase;
				}
			}
			//write bar
			Util.writeSimpleBarFile(chromosome, genomeVersion, ".", Num.arrayListOfIntegerToInts(positions), Num.arrayListOfFloatToArray(values), barFile);
		}	
	}
	
	
	/**Looks forward at overlapping windows, returns highest score for the given bp position.*/
	public double highestScoringWindow (int base, int currentWindowIndex){
		int num = windows.length;
		double score = windows[currentWindowIndex].getScores()[scoreIndex];
		for (int i=currentWindowIndex; i<num; i++){
			if (baseOverlapsWindow(base,windows[i])) {
				//check score
				double testScore = windows[i].getScores()[scoreIndex];
				if (testScore> score) score = testScore;				
			}
			else return score;
		}
		return score;
	}
	
	/**Is a given base contained within the window?*/
	public boolean baseOverlapsWindow(int bp, Window win){
		if (bp>= win.getStart1stOligo() && bp <= (win.getStartLastOligo()+sizeOfOligoMinusOne)) return true;
		return false;
	}
	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);	
		}
		new WindowBlockMaker(args);
	}
	

	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		String file = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': file = args[i+1]; i++; break;
					case 's': scoreIndex = Integer.parseInt(args[i+1]); i++; break;
					case 'i': flipScore=true; break;
					case 'z': sizeOfOligo =Integer.parseInt(args[i+1]);i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		windowFile = new File(file);
		if (windowFile == null || windowFile.canRead() == false){
			Misc.printExit("\nError: cannot read your Window array? "+file);
		}
		
	}	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Window Region Maker: Oct 2006                           **\n" +
				"**************************************************************************************\n" +
				"WBM builds heat maps from overlaping windows for display in IGB.\n\n" +
				
				
				"-f Full path file text for the Window[] array.\n" +
				"-s Score index to use in building heatmap.\n"+
				"-z Size of oligo, defaults to 25.\n"+
				"-i Multiply all scores by -1 prior to processing, good for reduced regions.\n\n"+
				
				"Example: java -Xmx500M trans/main/WindowBlockMaker -f /affy/res/zeste -s 1\n\n" +
				
				"The -Xmx flag increases the allowable heap size for the java virtual machine.  Adjust\n" +
				"accordingly in response to out of memory errors and the available resources.\n\n"+
				
		"**************************************************************************************\n");		
	}

	public int getScoreIndex() {
		return scoreIndex;
	}
	public void setScoreIndex(int scoreIndex) {
		this.scoreIndex = scoreIndex;
	}
	public float getBaseScore() {
		return baseScore;
	}
	public void setBaseScore(float baseScore) {
		this.baseScore = baseScore;
	}
}
