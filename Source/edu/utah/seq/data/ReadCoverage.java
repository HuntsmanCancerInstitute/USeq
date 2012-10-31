package edu.utah.seq.data;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.annotation.ExportIntergenicRegions;
import util.gen.*;
import edu.utah.seq.parsers.BarParser;
import edu.utah.seq.useq.data.Region;

/**Per base read coverage
 * @author Nix
 * */
public class ReadCoverage {

	//fields
	private File[] pointDataDirectories;
	private HashMap<String, ArrayList<PointData>>[] splitPointData;
	private File saveDirectory;
	private String chromosome;
	private int readLength = 0;
	private int negativeHalfReadLength = 0;
	private HashMap<String,Region[]> regions;
	private Histogram histogram;
	private boolean saveGraphs = true;
	private boolean verbose = true;
	private boolean makeRelativeTracks = true;
	private boolean keepStranded = false;
	private boolean replaceScoresWithHitCount = true;
	private float plusScalar = 0;
	private float minusScalar = 0;
	private float combineScalar = 0;
	private int minimumCounts = 8;
	private PrintWriter out = null;
	private boolean printScalars = true;

	/**For integration with the RNA-Seq app. Makes a read coverage track.*/
	public ReadCoverage (File saveDirectory, File[] pointDataDirectories, boolean keepStranded){
		verbose = false;
		printScalars = false;
		this.keepStranded = keepStranded;
		this.saveDirectory = saveDirectory;
		this.pointDataDirectories = pointDataDirectories;
		splitPointData = PointData.fetchStrandedPointDataNoMerge (pointDataDirectories);		

		//calculate total observations for relative tracks
		if (makeRelativeTracks) setScalars();

		//fetch all chromosomes
		HashSet<String> chromosomes = new HashSet<String>();
		chromosomes.addAll(splitPointData[0].keySet());
		chromosomes.addAll(splitPointData[1].keySet());

		//for each chromosome
		Iterator<String> it = chromosomes.iterator();
		while (it.hasNext()){
			chromosome = it.next();

			//merge strands summing scores of identical positions
			ArrayList<PointData> allPlus = splitPointData[0].get(chromosome);
			ArrayList<PointData> allMinus = splitPointData[1].get(chromosome);
			if (keepStranded == false){
				ArrayList<PointData> all = new ArrayList<PointData>();
				if (allPlus != null && allPlus.size() !=0) all.addAll(allPlus);
				if (allMinus != null&& allMinus.size() !=0) all.addAll(allMinus);
				if (all.size()==0) continue;
				if (readLength == 0) calculateReadLength(all);
				PointData merged = PointData.mergePointData(all, true, true);
				merged.shiftPositionsUnstranded(negativeHalfReadLength);

				//make coverage track
				makeCoverageTrack(merged, combineScalar);
			}
			else{
				if (allPlus != null && allPlus.size() !=0) {
					if (readLength ==0) calculateReadLength(allPlus);
					PointData merged = PointData.mergePointData(allPlus, true, true);
					merged.shiftPositionsUnstranded(negativeHalfReadLength);
					makeCoverageTrack(merged, plusScalar);
				}
				if (allMinus != null && allMinus.size() !=0) {
					if (readLength ==0) calculateReadLength(allMinus);
					PointData merged = PointData.mergePointData(allMinus, true, true);
					merged.shiftPositionsUnstranded(negativeHalfReadLength);
					makeCoverageTrack(merged, minusScalar);
				}
			}
		}
	}

	/**For stand alone app.*/
	public ReadCoverage(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		//process args
		processArgs(args);
		//fetch names
		splitPointData = PointData.fetchStrandedPointDataNoMerge (pointDataDirectories);

		//calculate total observations for relative tracks
		if (makeRelativeTracks) setScalars();

		//fetch all chromosomes
		HashSet<String> chromosomes = new HashSet<String>();
		chromosomes.addAll(splitPointData[0].keySet());
		chromosomes.addAll(splitPointData[1].keySet());

		//for each chromosome
		Iterator<String> it = chromosomes.iterator();
		while (it.hasNext()){
			chromosome = it.next();
			if (chromosome.startsWith("chrAdap")) continue;

			//need to process?
			if (saveGraphs == false && regions != null && regions.containsKey(chromosome) == false) continue;

			System.out.println("\t"+chromosome);

			//merge strands summing scores of identical positions
			ArrayList<PointData> allPlus = splitPointData[0].get(chromosome);
			ArrayList<PointData> allMinus = splitPointData[1].get(chromosome);
			if (keepStranded == false){
				ArrayList<PointData> all = new ArrayList<PointData>();
				if (allPlus != null && allPlus.size() !=0) all.addAll(allPlus);
				if (allMinus != null && allMinus.size() !=0) all.addAll(allMinus);
				if (all.size()==0) continue;
				if (readLength == 0) calculateReadLength(all);
				PointData merged = PointData.mergePointData(all, replaceScoresWithHitCount, true);
				merged.shiftPositionsUnstranded(negativeHalfReadLength);

				//make coverage track
				makeCoverageTrack(merged, combineScalar);
			}
			else{
				if (allPlus != null && allPlus.size() !=0) {
					if (readLength ==0) calculateReadLength(allPlus);
					PointData merged = PointData.mergePointData(allPlus, replaceScoresWithHitCount, true);
					merged.shiftPositionsUnstranded(negativeHalfReadLength);
					makeCoverageTrack(merged, plusScalar);
				}
				if (allMinus != null && allMinus.size() !=0) {
					if (readLength ==0) calculateReadLength(allMinus);
					PointData merged = PointData.mergePointData(allMinus, replaceScoresWithHitCount, true);
					merged.shiftPositionsUnstranded(negativeHalfReadLength);
					makeCoverageTrack(merged, minusScalar);
				}
			}

		}

		//print histogram?
		if (histogram !=null) {
			System.out.println("\nHistogram of the number of interrogated bps with a given read depth:");
			histogram.printScaledHistogram();

			int[] counts = histogram.getBinCounts();
			double total = histogram.getTotalBinCounts();
			System.out.println("\n"+(int)total+"\tTotal interrogated bps\n");

			//calculate % coverage
			System.out.println("Fraction interrogated bps with X or more reads:");
			double numCounts = 0;
			for (int i=0; i< counts.length; i++){
				numCounts += counts[i];
				double fraction = (total-numCounts)/ total;
				if (fraction <= 0) break;
				String formattedFraction = Num.formatNumber(fraction, 5);
				System.out.println(formattedFraction+"\t"+(i+1));
			}

			/*for (int i=1; i< counts.length; i++){
				numCounts += counts[i];
				fraction = (total-numCounts)/ total;
				formattedFraction = Num.formatNumber(fraction, 5);
				System.out.println(i+"\t"+formattedFraction);
			}*/

		}

		if (out != null) out.close();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	public void setScalars(){
		HashMap<String,PointData[]> plusPointData = PointData.convertArrayList2Array(splitPointData[0]);
		HashMap<String,PointData[]> minusPointData = PointData.convertArrayList2Array(splitPointData[1]);
		float totalNumberPlusObservations = PointData.totalObservationsMultiPointData(plusPointData);
		float totalNumberMinusObservations = PointData.totalObservationsMultiPointData(minusPointData);
		if (plusScalar == 0f) plusScalar = totalNumberPlusObservations/ 1000000f;
		if (minusScalar == 0f) minusScalar = totalNumberMinusObservations/ 1000000f;
		if (combineScalar == 0f) combineScalar = (totalNumberPlusObservations + totalNumberMinusObservations)/ 1000000f;
		if (printScalars){
			System.out.println("\nCounts and Scalars (#obs/1,000,000 used to divide the hit count by the scalar):");
			System.out.println("\t"+(int)totalNumberPlusObservations+"\tNumber Plus Observations");
			System.out.println("\t"+(int)totalNumberMinusObservations+"\tNumber Minus Observations");
			System.out.println("\t"+(int)(totalNumberPlusObservations + totalNumberMinusObservations)+"\tNumber Combine Observations");
			if (keepStranded){
				System.out.println("\t"+plusScalar+"\tPlus Scalar");
				System.out.println("\t"+minusScalar+"\tMinus Scalar");
			}
			else System.out.println("\t"+combineScalar+"\tCombine Scalar");
			System.out.println();
		}

	}

	public void makeCoverageTrack(PointData pd, float scalar){		
		//fetch data
		int[] positions = pd.getPositions();
		float[] counts = pd.getScores();

		ArrayList<Point> points = new ArrayList<Point>();
		//for each position
		for (int i=0; i<positions.length; i++){

			int startIndex = i;
			int position = positions[startIndex];
			int stopIndex = i;
			//look ahead to see if any overlapping reads
			i++;
			for (; i< positions.length; i++){
				//check size
				int diff = positions[i] - position;			
				//overlap or adjacent
				if (diff <= readLength) {
					stopIndex = i;
					position = positions[i];
				}
				//no overlap
				else break;
			}
			i--;

			//make base hit count array
			int basesCovered = positions[stopIndex] - positions[startIndex] + readLength;
			int[] baseCount = new int[basesCovered];
			int basePositionStartBlock = positions[startIndex];

			//for each read in block
			for (int j=startIndex; j<=stopIndex; j++){
				//number of overlapping reads
				int numberReads = (int)counts[j];
				//zeroed position
				int pos = positions[j] - basePositionStartBlock;
				//increment counts
				for (int k=0; k< readLength; k++){
					baseCount[pos+k] += numberReads;
				}
			}
			//save block, last base is included!
			//set beginning point
			if (basePositionStartBlock !=0) {
				points.add(new Point(basePositionStartBlock-1, 0));
			}

			//add first point
			int basePosition = basePositionStartBlock;
			int hits = baseCount[0];
			points.add( new Point(basePosition,hits));

			//for each subsequent
			for (int m=1; m<baseCount.length; m++){
				//add Point only if hits differs
				if (baseCount[m] != hits){
					//add old
					points.add( new Point(basePosition,hits));				
					//set new
					basePosition = m + basePositionStartBlock;
					hits = baseCount[m];
					points.add( new Point(basePosition,hits));	
				}
				else basePosition = m + basePositionStartBlock;
			}
			//add last
			points.add( new Point(basePosition,hits));			
			//set ending point
			points.add(new Point(basePosition+1, 0));
		}

		//write into PointData
		Point[] pts = new Point[points.size()];
		points.toArray(pts);

		PointData hitTrack = Point.extractPositionScores(pts);
		Info info = pd.getInfo();
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = new HashMap<String,String>();		
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		//color red
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, "#0066FF"); //light blue
		//what's the source
		String fileNames = Misc.stringArrayToString(IO.fetchFileNames(pointDataDirectories),",");
		map.put(BarParser.SOURCE_TAG, fileNames);
		//description
		String stranded = "";
		if (keepStranded) stranded = " Stranded.";
		if (makeRelativeTracks) map.put(BarParser.DESCRIPTION_TAG, "Relative read coverage track for "+Num.sumArray(counts)+" alignments. (e.g. #Reads/TotalReads/1M)."+stranded);
		else map.put(BarParser.DESCRIPTION_TAG, "Read coverage track for "+Num.sumArray(counts)+" alignments."+stranded);
		//save in info
		info.setNotes(map);
		hitTrack.setInfo(info);
		//calculate interrogated region hits?
		if (regions != null) calculateInterrogatedRegionCoverage(hitTrack);
		//scale values?
		if (scalar !=0){
			for (int z=0; z< pts.length; z++){
				pts[z].setScore(pts[z].getScore()/scalar);
			}
			hitTrack.setScores(Point.extractScores(pts));
		}
		if (saveGraphs) hitTrack.writePointData(saveDirectory);	

	}

	public void calculateInterrogatedRegionCoverage(PointData pd){
		//get start stops
		Region[] ss = regions.get(chromosome);
		if (ss == null) return;
		//get positions and scores
		int[] positions = pd.getPositions();
		float[] scores = pd.getScores();

		//make array of counts
		int lastBase = positions[positions.length-1];
		int[] baseCounts = new int[lastBase+100];		
		for (int i=0; i<positions.length; i++){
			int count = (int)scores[i];
			int startBase = positions[i];
			int stopIndex = i+1;
			int stopBase;
			if (stopIndex >= positions.length) stopBase =lastBase;
			else stopBase = positions[i+1];
			//fill the block
			Arrays.fill(baseCounts, startBase, stopBase, count);
		}

		//for each region fetch and sum scores
		int hits;
		for (int i=0; i< ss.length; i++){
			int[] startStop = ss[i].getStartStop();
			boolean[] falseMask = new boolean[ss[i].getLength()];
			Arrays.fill(falseMask, true);
			int index = 0;
			boolean badBPs = false;
			for (int j=startStop[0]; j < startStop[1]; j++){
				hits = 0;
				if (j < baseCounts.length) hits = baseCounts[j];
				histogram.count(hits);
				if (hits < minimumCounts) {
					falseMask[index] = false;
					badBPs = true;
				}
				index++;
			}
			if (badBPs){
				int[][] blocks = ExportIntergenicRegions.fetchFalseBlocks(falseMask, 0, 0);
				for (int j=0; j< blocks.length; j++){
					out.println(chromosome+"\t"+(blocks[j][0]+startStop[0])+"\t"+(blocks[j][1]+startStop[0]+1));
				}
			}
		}
	}

	/**Shifts the base position from center to the start of the read.*/
	public void calculateReadLength (ArrayList<PointData> pdAL){
		if (pdAL == null || pdAL.size() == 0) return;
		int num = pdAL.size();
		for (int i=0; i< num; i++){
			PointData pd = pdAL.get(i);
			//find size of read
			int readLengthLocal = pd.getInfo().getReadLength();	
			if (readLength !=0){
				if (readLength != readLengthLocal) {
					if(verbose)System.out.println("\t\tWarning: reads of mixed length found! Old length "+readLength+" New length "+readLengthLocal);
					readLength = readLengthLocal;
					negativeHalfReadLength = -1* (int)Math.round(  (((double)readLength) /2));				
				}
			}
			else{
				readLength = readLengthLocal;
				negativeHalfReadLength = -1* (int)Math.round(  (((double)readLength) /2));
			}
		}
		if (readLength ==0) Misc.printErrAndExit("\nError: read length is zero! Aborting.\n");

	}

	/**Sets all scores to 1. Typically these are a score associated with the alignment (a probability).*/
	public void stripScores (ArrayList<PointData> pdAL){
		if (pdAL == null || pdAL.size() == 0) return;
		int num = pdAL.size();
		for (int i=0; i< num; i++){
			pdAL.get(i).stripScores();
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ReadCoverage(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File bedFile = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': pointDataDirectories = IO.extractFiles(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'i': bedFile = new File(args[++i]); break;
					case 'b': saveGraphs = false; break;
					case 'k': keepStranded = true; break;
					case 'm': minimumCounts = Integer.parseInt(args[++i]); break;
					case 'l': plusScalar = Float.parseFloat(args[++i]); break;
					case 'n': minusScalar = Float.parseFloat(args[++i]); break;
					case 'c': combineScalar = Float.parseFloat(args[++i]); break;
					case 'r': makeRelativeTracks = false; break;
					case 'a': replaceScoresWithHitCount = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for point directories
		if (pointDataDirectories == null || pointDataDirectories[0].isDirectory() == false) Misc.printExit("\nError: cannot find your PointData directories(s)!\n");
		//only one directory look deeper
		if (pointDataDirectories.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(pointDataDirectories[0]);
			if (otherDirs != null && otherDirs.length > 0) pointDataDirectories = otherDirs;
		}

		//create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		saveDirectory.mkdirs();

		//looking at interrogated regions?
		if (bedFile != null) {
			if (bedFile.exists() == false) Misc.printExit("\nError: your interrogated region file doesn't exist?! \n");
			regions = Region.parseStartStops(bedFile, 0, 0, 0);
			histogram = new Histogram(0, 151, 151);
			File lowCoverage = new File(saveDirectory, Misc.removeExtension(bedFile.getName())+"_LowCoverage.bed");
			try {
				out = new PrintWriter (new FileWriter (lowCoverage));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Read Coverage: Feb 2012                            **\n" +
				"**************************************************************************************\n" +
				"Generates read coverage stair-step xxx.bar graph files for visualization in IGB. Will\n" +
				"also calculate per base coverage stats for a given file of interrogated regions and\n" +
				"create a bed file of regions with low coverage based on the minimum number of reads.\n" +
				"By default, graph values are scaled per million mapped reads.\n\n" +

				"Options:\n"+
				"-p Point Data directories, full path, comma delimited. Should contain chromosome\n" +
				"       specific xxx.bar.zip or xxx_-_.bar files. Can also provide one dir containing\n" +
				"       PointData dirs.\n"+
				"-s Save directory, full path.\n"+
				"-k Data is stranded, defaults to merging strands while generating graphs.\n"+
				"-a Data contains hit counts due to running it through the MergePointData app.\n"+
				"-r Don't scale graph values. Leave as actual read counts. \n"+
				"-i (Optional) Full path file text for a tab delimited bed file (chr start stop ...)\n" +
				"       containing interrogated regions to use in calculating a per base coverage\n" +
				"       statistics. Interbase coordinates assumed.\n"+
				"-m Minimum number reads for defining good coverage, defaults to 8. Use this in combo\n"+
				"       with the interrogated regions file to identify poor coverage regions.\n"+
				"-b Just calculate stats, skip coverage graph generation.\n"+
				"-l Plus scalar, for stranded RC output, defaults to # plus observations/1000000\n"+
				"-n Minus scalar, for stranded RC output, defaults to # minus observations/1000000\n"+
				"-c Combine scaler, defaults to # observations/1000000\n"+


				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/ReadCoverage -p\n" +
				"      /Data/Ets1Rep1/,/Data/Ets1Rep2/ -s /Data/MergedHitTrckEts1 -i \n" +
				"      /CapSeqDesign/interrogatedExonsChrX.bed\n\n"+

		"**************************************************************************************\n");

	}		

}
