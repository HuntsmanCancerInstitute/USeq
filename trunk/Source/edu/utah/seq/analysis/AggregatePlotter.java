package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.bio.annotation.*;
import util.bio.parsers.*;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import trans.tpmap.*;

/**Generates aggregate plots/ class averages given a list of regions and point data.
 * @author nix
 * */
public class AggregatePlotter {

	//user defined fields
	private File[] pointDataDirs;
	private int peakShift = 0;
	private int halfPeakShift = 0;
	private File bedFile;
	private int scaleSize = 0;
	private boolean scaleScores = false;
	private boolean delogScores = false;
	private boolean replaceScoresWithOne = false;
	private boolean normalize = false;
	private boolean averageRegionScores = false;
	
	//strand usage
	private static final int COMBINE = 0; //draw from both, invert - strand
	private static final int SAME = 1; //draw only from same strand
	private static final int OPPOSITE = 2; //draw only from opp strand
	private int strandUsage = COMBINE; 

	//internal fields
	private String chromosome;
	private String strand;
	private HashMap<String, PointData[]> chromPointData;
	private HashMap<String,ArrayListStartStop[]> chromRegions;
	private int maxSize = 0;
	private float[][] regionScores = null;
	private float[] mergedRegionScores = null;
	float numberOfRegions = 0;


	//constructors
	public AggregatePlotter(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//load gene models from refFlat for refSeq UCSC gene table
		chromRegions = ArrayListStartStop.parseArrayListStartStops(bedFile, 0, 0, true);
		if (chromRegions == null | chromRegions.size() == 0) Misc.printExit("\nProblem loading your bed file.\n");

		//fetch data but don't load it.
		System.out.println("Fetching data...");
		chromPointData = PointData.fetchStrandedCombinePointData(pointDataDirs);
		System.out.println("\t"+PointData.totalObservationsMultiPointData(chromPointData)+" reads");

		//for each region load data
		System.out.print("Loading regions...");
		Iterator<String> it = chromRegions.keySet().iterator();
		while (it.hasNext()){
			//fetch chromosome and strand
			String chromStrand = it.next();
			chromosome = chromStrand.substring(0,chromStrand.length()-1);
			strand = chromStrand.substring(chromStrand.length()-1);
			if (strand.equals(".") && strandUsage != COMBINE) Misc.printErrAndExit("\nStrand usage has been indicated yet regions aren't stranded! Is your bed file properly formated (5 columns)? Aborting.\n"+chromStrand);

			//fetch and shift positions or delog scores
			if (halfPeakShift !=0 || delogScores || replaceScoresWithOne) shiftPositions();

			//load regions
			loadRegions();
			System.out.print(".");
		}
		System.out.println();

		//zero position regions
		resetPositions();
		
		//scale?
		if (scaleSize != 0){
			System.out.println("Scaling regions to "+scaleSize+"bp");
			scaleToSize(scaleSize);
		}
		else if (sameSizes() != true) {
			System.out.println("Regions of different size found. Scaling to max size "+maxSize);
			scaleToSize(maxSize);
			scaleSize = maxSize;
		}
		else scaleSize = maxSize;
		
		//make float[], where each index is a bp
		loadRegionScores();

		//sum columns
		if (averageRegionScores) averageRegionScoreColumns();
		else sumRegionScoreColumns();
		
		//normalize by dividing against the number of regions?
		if(normalize) normalizeSummedRegionScores();

		//print arrays
		printSmoothedArrays();
	
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	/**Prints smoothed arrays of 5, 10, and 20% of region size.*/
	private void printSmoothedArrays(){
		
		int fivePercent = Math.round((float)scaleSize * 0.05f);
		float[] smoothedFive = Num.windowAverageScoresNoShrink(mergedRegionScores, fivePercent);
		int tenPercent = Math.round((float)scaleSize * 0.10f);
		float[] smoothedTen = Num.windowAverageScoresNoShrink(mergedRegionScores, tenPercent);
		int twentyFivePercent = Math.round((float)scaleSize * 0.25f);
		float[] smoothedTwentyFive = Num.windowAverageScoresNoShrink(mergedRegionScores, twentyFivePercent);
		
		System.out.println("\nPrinting region size window smoothed data:\n");
		System.out.println("Index\t0%\t5%("+fivePercent+")\t10%("+tenPercent+")\t25%("+twentyFivePercent+")");
		for (int i=0; i< mergedRegionScores.length; i++){
			System.out.println((i+1)+"\t"+mergedRegionScores[i]+"\t"+smoothedFive[i]+"\t"+smoothedTen[i]+"\t"+smoothedTwentyFive[i]);
		}
		
		
	}
	
	private void normalizeSummedRegionScores(){
		for (int i=0; i< mergedRegionScores.length; i++) mergedRegionScores[i] = mergedRegionScores[i]/numberOfRegions;
	}
	
	/**Sums the columns/ bp for all the regions*/
	private void sumRegionScoreColumns(){
		mergedRegionScores = new float[regionScores[0].length];
		//for each bp
		for (int i=0; i< regionScores[0].length; i++){
			//for each region
			for (int j=0; j< regionScores.length; j++){
				mergedRegionScores[i] += regionScores[j][i];
			}
		}
	}
	
	/**Averages the columns/ bp for all the regions ignoring zeros*/
	private void averageRegionScoreColumns(){
		mergedRegionScores = new float[regionScores[0].length];
		//for each bp
		for (int i=0; i< regionScores[0].length; i++){
			float total = 0;
			float nonZeroObservations = 0;
			//for each region
			for (int j=0; j< regionScores.length; j++){
				if (regionScores[j][i] !=0){
					total += regionScores[j][i];
					nonZeroObservations++;
				}
			}
			if (nonZeroObservations != 0) mergedRegionScores[i] = total/nonZeroObservations;
		}
	}

	/**Loads the scores, inverts the Points for - stranded regions.*/
	private void loadRegionScores(){
		ArrayList<float[]> scoresAL = new ArrayList<float[]>();
		Iterator<String> it = chromRegions.keySet().iterator();	
		//for each chromosome
		while (it.hasNext()){
			//fetch chromosome and strand
			String chromStrand = it.next();
			boolean antiSense = chromStrand.endsWith("-");
			ArrayListStartStop[] alss = chromRegions.get(chromStrand);
			//for each region
			for (int i=0; i< alss.length; i++){
				ArrayList al = alss[i].getArrayList();
				if (al.size()==0) continue;
				Point[] pts = (Point[])al.get(0);
				if (pts!=null){
					float[] s = new float[scaleSize];
					for (int j=0; j<pts.length; j++){
						int pos = pts[j].getPosition();
						if (pos >= scaleSize) pos = scaleSize-1;
						s[pos] = pts[j].getScore();
					}
					if (antiSense) Num.invertArray(s);
					scoresAL.add(s);
				}
			}
		}
		//make array
		int numScores = scoresAL.size();
		if (numScores == 0) Misc.printExit("\nNo regions with scores found?\n");
		regionScores = new float[numScores][scaleSize];
		for (int i=0; i < numScores; i++){
			regionScores[i] = scoresAL.get(i);
		}
	}
	
	/**Subtracts the start position of the region from each associated Point.*/
	private void resetPositions(){		
		Iterator<String> it = chromRegions.keySet().iterator();
		while (it.hasNext()){
			//fetch chromosome and strand
			String chromStrand = it.next();
			ArrayListStartStop[] alss = chromRegions.get(chromStrand);
			for (int i=0; i< alss.length; i++){
				int start = alss[i].getStart();
				ArrayList al = alss[i].getArrayList();
				if (al.size() ==0) continue;
				Point[] pts = (Point[])al.get(0);
				if (pts == null) continue;
				for (int j=0; j<pts.length; j++){
					int pos = pts[j].getPosition();
					pts[j].setPosition(pos-start);					
				}
				
			}
		}		
	}

	/**Scales the Point positions to the defined size.*/
	private void scaleToSize(int size){
		Iterator<String> it = chromRegions.keySet().iterator();
		double sizeMinOne = size -1;
		while (it.hasNext()){
			//fetch chromosome and strand
			String chromStrand = it.next();
			ArrayListStartStop[] alss = chromRegions.get(chromStrand);
			for (int i=0; i< alss.length; i++){
				double length = alss[i].getStop() - alss[i].getStart();
				if (length == sizeMinOne) continue;
				double scalar = sizeMinOne/length;
				ArrayList al = alss[i].getArrayList();
				Point[] pts = (Point[])al.get(0);
				if (pts == null) continue;
				for (int j=0; j<pts.length; j++){
					double pos = pts[j].getPosition();
					int newPos = (int)Math.round(pos*scalar);
					pts[j].setPosition(newPos);					
				}
				//sum identical positions
				al.remove(0);
				al.add(Point.sumIdenticalPositionScores(pts));
			}
		}		
	}

	private boolean sameSizes(){
		Iterator<String> it = chromRegions.keySet().iterator();
		while (it.hasNext()){
			//fetch chromosome and strand
			String chromStrand = it.next();
			ArrayListStartStop[] alss = chromRegions.get(chromStrand);
			for (int i=0; i< alss.length; i++){
				int length = alss[i].getStop() - alss[i].getStart();
				if (length != maxSize) return false;
			}
		}
		return true;
	}


	/**Load Point[] into each region paying attention to strand.*/
	private void loadRegions(){
		//fetch PointData
		PointData[] pd = chromPointData.get(chromosome);
		if (pd == null) return;
		PointData plusPD = pd[0];
		PointData minusPD = pd[1];

		//get regions
		ArrayListStartStop[] regions = chromRegions.get(chromosome+strand);

		numberOfRegions += regions.length;
		
		//for each region
		for (int i=0; i< regions.length; i++){
			//fetch Point[], may be null
			int start = regions[i].getStart();
			int stop = regions[i].getStop();
			int length = stop-start;
			if (length > maxSize) maxSize = length;
			Point[] ptsPlus = null;
			Point[] ptsMinus = null;
			//what to fetch?
			Point[] combo = null;
			//fetch and combine both
			if (strandUsage == COMBINE){
				if (plusPD !=null) ptsPlus = plusPD.fetchPoints(start, stop);
				if (minusPD !=null) ptsMinus = minusPD.fetchPoints(start, stop);
				if (ptsPlus != null){
					if (ptsMinus !=null) combo = Point.mergePoints(ptsPlus, ptsMinus);
					else combo = ptsPlus;
				}
				else if (ptsMinus != null) combo = ptsMinus;
			}

			//fetch same strand as region
			else if (strandUsage == SAME){
				if (strand.equals("+") && plusPD != null){
					combo = plusPD.fetchPoints(start, stop);
				}
				else if (minusPD != null) {				
					combo = minusPD.fetchPoints(start, stop);
				}
			}
			//defaults to opposite
			else {
				if (strand.equals("+") && minusPD != null) combo = minusPD.fetchPoints(start, stop);
				else if (plusPD != null) combo = plusPD.fetchPoints(start, stop);
			}
			//sum scores
			if (combo!= null) {
				combo = Point.sumIdenticalPositionScores(combo);
				//scale scores?
				if (scaleScores) {
					//find total
					float total = 0;
					for (int x=0; x < combo.length; x++) total+= combo[x].getScore();
					float scalar = total/100;
					for (int x=0; x < combo.length; x++) combo[x].setScore(combo[x].getScore()/scalar);
					
				}
			}
			regions[i].getArrayList().add(combo);
		}
	}



	/**Shifts the positions halfPeakShift (+ for sense, - for antisense) sets the positions into the data*/
	private void shiftPositions(){
		//fetch data 
		PointData[] t = chromPointData.get(chromosome);
		if (t== null) return;
		if (t[0]!=null) {
			int[] p = t[0].getPositions();
			if (halfPeakShift !=0) addShift(p,halfPeakShift);
			if (replaceScoresWithOne) t[0].stripScores();
			else if (delogScores) t[0].setScores(Num.antiLog(t[0].getScores(), 2));
		}
		if (t[1]!=null) {
			int[] p = t[1].getPositions();
			if (halfPeakShift !=0) addShift(p, -1*halfPeakShift);
			if (replaceScoresWithOne) t[1].stripScores();
			else if (delogScores) t[1].setScores(Num.antiLog(t[1].getScores(), 2));
		}

	}

	/**Adds the toAdd to each int.*/
	public static void addShift(int[] positions, int toAdd){
		for (int i=0; i< positions.length; i++){
			positions[i] += toAdd;
			if (positions[i]<0) positions[i] = 0;
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AggregatePlotter(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': pointDataDirs = IO.extractFiles(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'p': peakShift = Integer.parseInt(args[++i]); break;
					case 's': scaleSize = Integer.parseInt(args[++i]); break;
					case 'v': scaleScores = true; break;
					case 'd': delogScores = true; break;
					case 'n': normalize = true; break;
					case 'a': averageRegionScores = true; break;
					case 'r': replaceScoresWithOne = true; break;
					case 'u': strandUsage = Integer.parseInt(args[++i]); break;
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
		if (pointDataDirs == null || pointDataDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your PointData directories(s)!\n");
		//only one directory look deeper
		if (pointDataDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(pointDataDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) pointDataDirs = otherDirs;
		}

		//set half peak shift
		halfPeakShift = (int)Math.round( ((double)peakShift)/2 );

		//look for bed file
		if (bedFile == null) Misc.printExit("\nPlease enter a regions file to use in scoring regions.\n");
		
		//print params
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " "));
		System.out.println("\tBed region file\t"+bedFile);
		System.out.println("\tPeak shift\t"+peakShift);
		System.out.println("\tStrand usage\t"+strandUsage);
		if (scaleSize !=0 ) System.out.println("\tScale size\t"+scaleSize);
		System.out.println("\tScale scores\t"+scaleScores);
		System.out.println("\tDivide scores by # regions\t"+normalize);
		System.out.println("\tDelog2 scores\t"+delogScores);
		System.out.println("\tReplace scores with 1\t"+replaceScoresWithOne);
		System.out.println("\tAverage scores instead of summing\t"+averageRegionScores);
		System.out.println();

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Aggregate Plotter:  August 2012                       **\n" +
				"**************************************************************************************\n" +
				"Fetches point data contained within each region, inverts - stranded annotation, zeros\n" +
				"the coordinates, sums, and window averages the values.  Usefull for generating\n" +
				"class averages from a list of annotated regions. Use a spreadsheet app to graph the\n" +
				"results.\n\n" +

				"Options:\n"+
				"-t PointData directories, full path, comma delimited. These should contain chromosome\n" +
				"       specific xxx.bar.zip files.\n" +
				"-b Bed file (chr, start, stop, text, score, strand(+/-/.), full path, containing\n" +
				"       regions to stack. Must be all the same size.\n"+
				"-p Peak shift, average distance between + and - strand peaks. Will be used to shift\n" +
				"       the PointData by 1/2 the peak shift, defaults to 0. \n"+
				"-u Strand usage, defaults to 0 (combine), 1 (use only same strand), 2 (opposite\n" +
				"       strand), or 3 (ignore).\n"+
				"       this option to select particular stranded data to aggregate.\n"+
				"-r Replace scores with 1.\n"+
				"-d Delog2 scores. Do it if your data is in log2 space.\n"+
				"-v Convert each region scores to % of total.\n"+
				"-n Divide scores by the number of regions.\n"+
				"-s Scale all regions to a particular size. Defaults to max region size.\n"+
				"-a Average region scores instead of summing.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/AgregatePlotter -t\n" +
				"      /Data/PolIIRep1/,/Data/PolIIRep2/ -b /Anno/tssSites.bed -p 73 -u 1\n\n" +

		"**************************************************************************************\n");

	}
}
