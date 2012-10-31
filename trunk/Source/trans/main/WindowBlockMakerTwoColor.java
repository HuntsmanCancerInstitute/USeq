package trans.main;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;
import trans.misc.*;
import java.io.*;

import edu.utah.seq.parsers.BarParser;

/**Converts a Window[] into a heat map specific bar file for import into IGB. 
 * Assumes windows are sorted by chromosome and start position
 * Score index is the index to use when fetching a score for a window from it's double[] of scores
 * see ScanChip or ScanChromosomesCNV for current score indexes.
 * Note, this will place at the 0 and 1 base two values to maintain symetry in the values.*/
public class WindowBlockMakerTwoColor {
	//fields
	private int scoreIndex = 1;	
	private Window[][] chromosomeWindows; 
	private int chrWinIndex;
	private String chromosome;
	private File[] windowFiles;
	private File barDirectory;
	private float maxValue;
	private float minValue;
	private double minWindowFilter = 0;
	private double maxWindowFilter = 0;
	private int sizeOfOligo = 1;
	private int sizeOfOligoMinusOne = sizeOfOligo-1;
	private String genomeVersion = null;
	private ArrayList bases = new ArrayList(); //Integers
	private ArrayList values = new ArrayList(); //Floats
	private String strand = "."; //+,-,or the default .
	
	//constructors
	public WindowBlockMakerTwoColor(int sizeOfOligo, String genomeVersion, String strand, double minWindowFilter, double maxWindowFilter){
		this.genomeVersion = genomeVersion;
		this.sizeOfOligo = sizeOfOligo;
		sizeOfOligoMinusOne = sizeOfOligo - 1;
		if (strand != null) strand = strand;
		this.minWindowFilter = minWindowFilter;
		this.maxWindowFilter = maxWindowFilter;
	}

	public WindowBlockMakerTwoColor(String[] args){
		processArgs(args);
		sizeOfOligoMinusOne = sizeOfOligo - 1;
		
		if (maxWindowFilter != 0 || minWindowFilter != 0) {
			System.out.println("\nRemoving windows with scores falling within "+minWindowFilter+" - "+maxWindowFilter);
		}
		
		//load Window[]
		for (int i=0; i< windowFiles.length; i++){
			System.out.println("\nReading in Window[] "+windowFiles[i].getName());
			Window[] windows = (Window[])IO.fetchObject(windowFiles[i]);
			barDirectory = new File (windowFiles[i].getParent(), "BarFiles_"+windowFiles[i].getName());
			if (barDirectory.exists() == false ) barDirectory.mkdir();
			if (maxWindowFilter != 0 || minWindowFilter != 0)  windows = filterWindows(windows);
			makeHeatMapBarFiles(windows);
		}
		System.out.println("\nDone!\n");
		

	}
	
	/*public static void main( String[] args){
		//Testing
		Window[] windows = new Window[5];
		windows[0] = new Window("chr1", 10, 20, 10, new double[]{5}); //String chromosome, int start, int stop, int numberOligos, double[] scores
		windows[1] = new Window("chr1", 15, 30, 10, new double[]{10});
		windows[2] = new Window("chr1", 25, 50, 10, new double[]{5});
		windows[3] = new Window("chr1", 45, 75, 10, new double[]{-5});
		windows[4] = new Window("chr1", 70, 80, 10, new double[]{-8});
		WindowBlockMakerTwoColor bm = new WindowBlockMakerTwoColor (1, "hg17", 0);
		bm.makeHeatMapBarFiles(windows);
	}*/
	
	public void makeHeatMapBarFiles(Window[] windows, File barDirectory, int scoreIndex){
		this.barDirectory = barDirectory;
		this.scoreIndex = scoreIndex;
		Window[] filtered = filterWindows(windows);
		if (filtered != null && filtered.length !=0) makeHeatMapBarFiles(filtered);
	}
	
	/**Tosses any windows within the min max filter range.*/
	public Window[] filterWindows(Window[] windows){
		if (maxWindowFilter == 0 && minWindowFilter == 0) return windows;
		ArrayList winAL = new ArrayList (windows.length/2);
		for (int i=0; i< windows.length; i++){
			double score = windows[i].getScores()[scoreIndex];
			if (score <= 0 && score <= minWindowFilter) winAL.add(windows[i]);
			else if (score >= 0 && score >= maxWindowFilter) winAL.add(windows[i]);
		}
		Window[] filtered = new Window[winAL.size()];
		winAL.toArray(filtered);
		return filtered;
	}
		
	/**Makes a stair step heat map from an array of windows in bar format. One per chromosome.
	 * Don't forget to set the barDirectory and score Index!!!!!!!*/
	public void makeHeatMapBarFiles(Window[] windows){
		//make bar parser
		BarParser bp = new BarParser();
		bp.setZipCompress(true);
		HashMap<String,String> tagVals = new HashMap<String,String>();
		tagVals.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		
		//split into different chromosomes
		chromosomeWindows = Util.splitWindowsByChromosome(windows);
		
		//for each chromosome
		System.out.print("Printing... ");
		for (chrWinIndex=0; chrWinIndex<chromosomeWindows.length; chrWinIndex++){
			chromosome = chromosomeWindows[chrWinIndex][0].getChromosome();
			System.out.print(chromosome+" ");
			//add blocks
			assembleHeatMapBlocks();
			//balance by adding max or min at zero base
			balanceValues();
			//write bar file
			File barFile = new File (barDirectory,chromosome+ ".bar");
			bp.writeBarFile(barFile, chromosome, genomeVersion, strand.charAt(0), Num.arrayListOfIntegerToInts(bases), Num.arrayListOfFloatToArray(values), tagVals);
			//clear ArrayLists
			bases.clear();
			values.clear();
		}
		System.out.println();
	}
	
	public void balanceValues(){
		float value2Balance = 0;
		if ((-1*minValue) > maxValue) value2Balance = -1* minValue;
		else value2Balance = -1 * maxValue;
		bases.add(0, new Integer(1));
		values.add(0, new Float(0));
		bases.add(0, new Integer(0));
		values.add(0, new Float(value2Balance));
	}
	
	public void printSgrLines(){
		int num = bases.size();
		for (int i= 0; i< num; i++){
			System.out.println(chromosome +"\t"+bases.get(i)+"\t"+values.get(i));
		}
	}
	
	public void assembleHeatMapBlocks(){
		//reset max and min
		maxValue = 0;
		minValue = 0;
		//find windowed block
		int startIndex = 0;
		int endIndex = 0;
		Window leftWindow = chromosomeWindows[chrWinIndex][startIndex];
		for (int i=1; i< chromosomeWindows[chrWinIndex].length; i++){
			Window rightWindow = chromosomeWindows[chrWinIndex][i];
			//do they overlap?
			//no block found!
			if (overlapOrAbut(leftWindow, rightWindow) == false){			
				endIndex = i;
				//note endIndex is not included when making heatMap blocks
				addHeatMapBlocks(startIndex, endIndex);
				//start new block
				startIndex = i;
				
			}
			leftWindow = rightWindow;
		}
		//find last
		addHeatMapBlocks(startIndex, chromosomeWindows[chrWinIndex].length);
	}
	
	/**Make a heatmap blocks from given window indexes, stop not included.*/
	public void addHeatMapBlocks(int startIndex, int stopIndex){
		
		//create arrays of float, one per base to hold max positive and max negative scores
		int startBase = chromosomeWindows[chrWinIndex][startIndex].getStart1stOligo();
		int stopBase = chromosomeWindows[chrWinIndex][stopIndex-1].getStartLastOligo()+sizeOfOligoMinusOne;
		int numBases = 1+ stopBase- startBase;
		float[] maxPosValues = new float[numBases];
		float[] maxNegValues = new float[numBases];
		
		//load max arrays with max scores
		//for each window
		for (int i=startIndex; i< stopIndex; i++){
			float score = new Double (chromosomeWindows[chrWinIndex][i].getScores()[scoreIndex]).floatValue();
			int baseIndex = chromosomeWindows[chrWinIndex][i].getStart1stOligo() - startBase;
			int length = (chromosomeWindows[chrWinIndex][i].getStartLastOligo() + sizeOfOligoMinusOne) - 
				chromosomeWindows[chrWinIndex][i].getStart1stOligo() + baseIndex+ 1;
			//is the score positive 
			if (score > 0.0f ){
				//It's positive, attempt to increase max
				for (int j=baseIndex; j< length; j++){
					if (score > maxPosValues[j]){ 
						maxPosValues[j] = score;
						if (score > maxValue) maxValue = score;
					}
				}
			}
			//is it negative
			else if (score < 0.0f){
				for (int j=baseIndex; j< length; j++){
					if (score < maxNegValues[j]) {
						maxNegValues[j] = score;
						if (score < minValue) minValue = score;
					}
				}
			}
		}
		
		//average overlapping bases setting in maxPosValues[]
		for (int i=0; i<numBases; i++){
			if (maxPosValues[i] == 0 && maxNegValues[i] != 0) maxPosValues[i] = maxNegValues[i];
			else if (maxPosValues[i] != 0 && maxNegValues[i] != 0 ) maxPosValues[i] = (maxPosValues[i]+maxNegValues[i])/2.0f;
		}		
		
		//build blocks		
		//open first block
		//set zero mark
		int previousBase = startBase -1;
		if (previousBase < 0) previousBase = 0;
		add(previousBase, 0);
		//set block value
		float blockValue = maxPosValues[0];
		add(startBase, blockValue);
		
		//advance each base opening and closing blocks
		for (int i=1; i<numBases; i++){
			float testValue = maxPosValues[i];
			if (testValue != blockValue){
				//close old
				add(i-1+startBase,blockValue);
				//open new
				blockValue = testValue;
				add(i+startBase,blockValue);
			}
		}
		//close last block
		add(numBases-1+startBase,blockValue);
		add(numBases+startBase,0);

		
	}
	
	/**Adds a gr line to the global ArrayLists for later burning to bar.*/
	public void add(int base, float value){
		bases.add(new Integer(base));
		values.add(new Float(value));
	}
	

	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);	
		}
		new WindowBlockMakerTwoColor(args);
	}
	
	/**Assumes right window is to right of left or abuts or doesn't overlap.*/
	public boolean overlapOrAbut (Window left, Window right){
		//overlap or abut, note not using sizeOfOligoMinusOne to include abuts
		if ((left.getStartLastOligo()+sizeOfOligo) >= right.getStart1stOligo()) return true;
		return false;
	}

	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		String file = null;
		String range = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': file = args[i+1]; i++; break;
					case 'r': range = args[i+1]; i++; break;
					case 's': scoreIndex = Integer.parseInt(args[i+1]); i++; break;
					case 'z': sizeOfOligo =Integer.parseInt(args[i+1]);i++; break;
					case 'v': genomeVersion = args[i+1]; i++; break;
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
		windowFiles = IO.extractFiles(new File(file));
		if (windowFiles == null || windowFiles[0].canRead() == false){
			Misc.printExit("\nError: cannot read/ find your Window array file(s)? "+file);
		}
		
		if (range !=null){
			String[] minMax = range.split(",");
			minWindowFilter = Double.parseDouble(minMax[0]);
			maxWindowFilter = Double.parseDouble(minMax[1]);
		}
		
	}	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                     Window Region Maker Two Color: Feb 2007                       **\n" +
				"**************************************************************************************\n" +
				"WBM builds two color heat maps from overlaping windows for display in IGB. A value is\n" +
				"placed at the zero base to keep symmetry between the min and max scores.\n\n" +

				"-f Full path file/ directory text for the Window[] array(s).\n" +
				"-s Score index to use in building heatmap, defaults to 1.\n"+
				"-z Size of oligo, defaults to 25.\n"+
				"-e Exclude windows falling within this range (ie -0.2,0.2).\n"+
				"-v Genome version (ie hg17, dm2, ce2, mm8), get from UCSC Browser,\n" +
				"      http://genome.ucsc.edu/FAQ/FAQreleases, for bar files.\n" +
				
				"Example: java -Xmx500M trans/main/WindowBlockMakerTwoColor -f /affy/res/zeste -s 1\n" +
				"      -v dm2 -e -0.25,0.25\n\n" +
				
		"**************************************************************************************\n");		
	}

	public void setBarDirectory(File barDirectory) {
		this.barDirectory = barDirectory;
	}

	public double getMaxWindowFilter() {
		return maxWindowFilter;
	}

	public void setMaxWindowFilter(double maxWindowFilter) {
		this.maxWindowFilter = maxWindowFilter;
	}

	public double getMinWindowFilter() {
		return minWindowFilter;
	}

	public void setMinWindowFilter(double minWindowFilter) {
		this.minWindowFilter = minWindowFilter;
	}
}
