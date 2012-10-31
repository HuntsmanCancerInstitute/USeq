package trans.misc;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.annotation.*;
import util.gen.*;
import trans.misc.*;
import java.io.*;


import edu.utah.seq.parsers.BarParser;
import edu.utah.seq.useq.data.RegionScoreText;

/**Converts bed files into chromosome specific stair-step bar files for import into IGB. 
 * Assumes interbase coordinates. */
public class Bed2Bar {
	//fields
	private File[] bedFiles;
	private HashMap<String,RegionScoreText[]> bedLinesHash;
	private String chromosome;
	private File barDirectory;
	private File bedFile;
	private RegionScoreText[] windows;
	private float maxValue;
	private float minValue;
	private String genomeVersion = null;
	private ArrayList<Integer> bases = new ArrayList<Integer>(); //Integers
	private ArrayList<Float> values = new ArrayList<Float>(); //Floats
	private boolean sumScores = false;
	private float threshold = 0;
	private float maxGap = 0;
	private PrintWriter bedOut;

	public Bed2Bar(String[] args){
		try {
			processArgs(args);
			//load Window[]
			for (int i=0; i< bedFiles.length; i++){
				bedFile = bedFiles[i];
				System.out.println("Parsing "+ bedFile.getName());
				bedLinesHash = Bed.parseBedFile(bedFile, true);
				if (bedLinesHash ==  null || bedLinesHash.size() ==0) {
					System.out.println("Problem parsing bed file, skipping!");
					continue;
				}
				barDirectory = IO.makeDirectory(bedFile,"");
				File bedOutFile = new File(Misc.removeExtension(bedFile.toString())+"_"+threshold+"_Filt.bed");
				bedOut = new PrintWriter (new FileWriter (bedOutFile));
				makeStairStepBarFiles();
			}
			bedOut.close();
			System.out.println("\nDone!\n");
		} catch (IOException e){
			e.printStackTrace();
		}
	}

	/**Makes a stair step heat map from an array of windows in bar format. One per chromosome.
	 * Don't forget to set the barDirectory and score Index!!!!!!!*/
	public void makeStairStepBarFiles(){
		//make bar parser
		BarParser bp = new BarParser();
		bp.setZipCompress(true);
		HashMap<String,String> tagVals = new HashMap<String,String>();
		tagVals.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		tagVals.put(BarParser.GRAPH_TYPE_COLOR_TAG, "#FF00FF"); //fusha
		tagVals.put(BarParser.SOURCE_TAG, bedFile.toString());

		//for each chromosome
		System.out.print("Printing... ");
		Iterator<String> it = bedLinesHash.keySet().iterator();
		while (it.hasNext()){
			chromosome = it.next();
			System.out.print(chromosome+" ");
			windows = bedLinesHash.get(chromosome);
			//add blocks
			assembleBlocks();
			//balance by adding max or min at zero base
			balanceValues();
			//write bar file
			File barFile = new File (barDirectory,chromosome+ ".bar");
			bp.writeBarFile(barFile, chromosome, genomeVersion, '.', Num.arrayListOfIntegerToInts(bases), Num.arrayListOfFloatToArray(values), tagVals);
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
		//set zero as base value
		bases.add(0, new Integer(1));
		values.add(0, new Float(0));
		//any negative values? Don't set unless so.
		if (minValue < 0){
			bases.add(0, new Integer(0));
			values.add(0, new Float(value2Balance));
		}
	}

	public void assembleBlocks(){
		//reset max and min
		maxValue = 0;
		minValue = 0;
		//find windowed block
		int startIndex = 0;
		int endIndex = 0;
		RegionScoreText leftWindow = windows[startIndex];
		for (int i=1; i< windows.length; i++){
			RegionScoreText rightWindow = windows[i];
			//do they overlap?
			//no block found!
			if (overlapOrAbut(leftWindow, rightWindow) == false){			
				endIndex = i;
				//note endIndex is not included when making heatMap blocks
				if (sumScores) sumHeatMapBlocks(startIndex, endIndex);
				else addHeatMapBlocks(startIndex, endIndex);
				//start new block
				startIndex = i;

			}
			leftWindow = rightWindow;
		}
		//find last
		if (sumScores) sumHeatMapBlocks(startIndex, windows.length);
		else addHeatMapBlocks(startIndex, windows.length);
	}

	public int findMaxStop(int startIndex, int stopIndex){
		int end = stopIndex -1;
		int max = windows[end].getStop();
		for (int i=startIndex; i< end; i++){
			if (windows[i].getStop() > max) max = windows[i].getStop(); 
		}
		return max;
	}

	/**Make a heatmap blocks from given window indexes, stop not included.*/
	public void addHeatMapBlocks(int startIndex, int stopIndex){

		//create arrays of float, one per base to hold max positive and max negative scores
		//windows are sorted by start and length
		int startBase = windows[startIndex].getStart();
		int stopBase = findMaxStop (startIndex, stopIndex);
		int numBases = 1+ stopBase- startBase;
		float[] maxPosValues = new float[numBases];
		float[] maxNegValues = new float[numBases];

		//load max arrays with max scores
		//for each window
		for (int i=startIndex; i< stopIndex; i++){
			float score = windows[i].getScore();
			int baseIndex = windows[i].getStart() - startBase;
			int length = windows[i].getStop() - windows[i].getStart() + baseIndex+ 1;	

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

	/**Make a heatmap blocks from given window indexes, stop not included.*/
	public void sumHeatMapBlocks(int startIndex, int stopIndex){
		//create arrays of float, one per base to hold sums
		//windows are sorted by start and length
		int startBase = windows[startIndex].getStart();
		int stopBase = findMaxStop (startIndex, stopIndex);

		int numBases = 1+ stopBase- startBase;
		float[] sums = new float[numBases];

		//load max arrays with max scores
		//for each window
		for (int i=startIndex; i< stopIndex; i++){
			float score = windows[i].getScore();
			int baseIndex = windows[i].getStart() - startBase;
			int length = windows[i].getStop() - windows[i].getStart() + baseIndex+ 1;	
			for (int j=baseIndex; j< length; j++) sums[j]+= score; 

		}	

		//build blocks		
		//open first block
		//set zero mark
		int previousBase = startBase -1;
		if (previousBase < 0) previousBase = 0;
		add(previousBase, 0);
		//set block value
		float blockValue = sums[0];
		add(startBase, blockValue);

		//advance each base opening and closing blocks
		for (int i=1; i<numBases; i++){
			float testValue = sums[i];
			if (testValue != blockValue){
				//close old
				add(i-1+startBase,blockValue);
				//open new
				blockValue = testValue;
				add(i+startBase,blockValue);
			}
		}
		//close last block
		add(numBases-2+startBase,blockValue);
		add(numBases-1+startBase,0);

		//make merged filtered bed file?
		if (bedOut != null){
			boolean[] falseMask = new boolean[sums.length];
			Arrays.fill(falseMask, true);
			for (int i=0; i< sums.length; i++) if (sums[i]>= threshold) falseMask[i] = false;
			int[][] blocks = ExportIntergenicRegions.fetchFalseBlocks(falseMask, 0, 0);
			//print bed file
			for (int i=0; i< blocks.length; i++){
				bedOut.println(chromosome +"\t"+ (startBase +blocks[i][0])+ "\t" + (startBase+blocks[i][1]));
			}
		}

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
		new Bed2Bar(args);
	}

	/**Assumes right window is to right of left or abuts or doesn't overlap.*/
	public boolean overlapOrAbut (RegionScoreText left, RegionScoreText right){
		//overlap or abut
		//if ((left.getStop()) >= right.getStart()) return true;
		//check distance
		int dist = right.getStart() - left.getStop();
		if (dist <= maxGap) return true;
		return false;
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
					case 'v': genomeVersion = args[i+1]; i++; break;
					case 's': sumScores = true; break;
					case 't': threshold = Float.parseFloat(args[++i]); break;
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
		if (file == null || file.exists()== false) Misc.printErrAndExit("Problem finding your bed files!\n");
		//pull files
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(file,".bed");
		tot[1] = IO.extractFiles(file,".bed.zip");
		tot[2] = IO.extractFiles(file,".bed.gz");
		bedFiles = IO.collapseFileArray(tot);
		if (bedFiles == null || bedFiles.length ==0) Misc.printErrAndExit("Problem finding your xxx.bed(.zip/.gz OK) files!\n");

		//genome version
		if (genomeVersion == null) Misc.printErrAndExit("Please enter a genome version (e.g. H_sapiens_Mar_2006, see http://genome.ucsc.edu/FAQ/FAQreleases\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                  Bed2Bar: June 2010                              **\n" +
				"**************************************************************************************\n" +
				"Bed2Bar builds stair step graphs from bed files for display in IGB. Strands are merged\n" +
				"and text information removed. Will also generate a merged bed file thresholding the \n" +
				"graph at that level. \n\n" +

				"-f Full path file or directory containing xxx.bed(.zip/.gz OK) files\n" +
				"-v Genome version (eg H_sapiens_Mar_2006), get from UCSC Browser,\n" +
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-s Sum bed scores for overlapping regions, defaults to assigning the highest score.\n"+
				"-t Threshold, defaults to 0.\n"+
				"-g Maximum gap, defaults to 0.\n"+

				"\nExample: java -Xmx4G pathTo/Apps/Bed2Bar -f /affy/res/zeste.bed.gz -v \n" +
				"      M_musculus_Jul_2007 -g 1000 -s -t 100 \n\n" +

		"**************************************************************************************\n");		
	}

	public void setBarDirectory(File barDirectory) {
		this.barDirectory = barDirectory;
	}
}
