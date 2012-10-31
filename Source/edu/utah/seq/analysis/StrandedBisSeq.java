package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import util.bio.cluster.HierarchicalClustering;
import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import trans.tpmap.*;

/**For processing bisulfite sequencing data
 * @author Nix
 * */
public class StrandedBisSeq {

	//user defined fields
	private File[] convertedPointDirs;
	private File[] nonConvertedPointDirs;
	private File saveDirectory;
	private int windowSize = 500;
	private int minimumNumberCsInWindow = 4;
	private HashMap<String, File> chromosomeFastaFiles = new HashMap<String, File>();
	private boolean printGraphs = false;
	private File fullPathToR = new File ("/usr/bin/R");
	private int minimumCoverage = 2;

	//internal fields
	private String genomeVersion;
	private String[] chromosomes;
	private String genomicSequence = null;
	private HashMap<String,PointData[]> convertedPlusPointData;
	private HashMap<String,PointData[]> convertedMinusPointData;
	private HashMap<String,PointData[]> nonConvertedPlusPointData;
	private HashMap<String,PointData[]> nonConvertedMinusPointData;
	private File pvalueDirectory;
	private File log2RatioDirectory;
	private File windowsDirectory;
	private int numberSkippedWindows = 0;
	private int numberPassingWindows = 0;
	private float fdrThreshold = 30;
	private float log2RatioThreshold = 1.585f;
	private int maxGap = 500;
	private EnrichedRegion[] enrichedRegions = null;
	private int log2RatioIndex = 0;
	private int pValIndex = 1;

	//window scanning
	private WindowMaker windowMaker; 
	private int[][] windows;
	private SmoothingWindow[] smoothingWindow;
	private SmoothingWindowIndex[] smoothingWindowIndex;
	private float windowMergeMultiplier = 0.80f;
	private ArrayList<File> swiFile = new ArrayList<File>();
	private double cacheNumber = 1000;
	private FisherExact fischerExact = new FisherExact((int)cacheNumber);
	private ChiSquareTest chiSquare = new ChiSquareTest();
	private boolean nonCGFound = false;

	//by chromosome
	private String chromosome;
	private PointData convertedMergedChromPlus = null;
	private PointData convertedMergedChromMinus = null;	
	private PointData nonConvertedMergedChromPlus = null;
	private PointData nonConvertedMergedChromMinus = null;	


	//constructors
	/**Stand alone.*/
	public StrandedBisSeq(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//scan each chromosome for particular mCG
		bpScan();

		//window scan chromosomes for larger regions with biased mCG
		if (nonCGFound){
			System.out.println("\nSkipping window scanning since nonCG contexts found.  Remove these using the ParsePointDataContexts and restart.");
		}
		else if (windowSize !=0) windowScan();

		IO.deleteDirectory(windowsDirectory);
		windowsDirectory.deleteOnExit();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	//methods

	public void windowScan(){
		//reset counters
		numberSkippedWindows = 0;
		numberPassingWindows = 0;
		swiFile.clear();

		//make dirs to hold graph tracks?
		if (printGraphs){
			pvalueDirectory = new File (saveDirectory, "WinPValues");
			pvalueDirectory.mkdir();
			log2RatioDirectory = new File (saveDirectory, "WinLog2Ratio");
			log2RatioDirectory.mkdir();
		}

		//load hashes
		loadDataHashes();

		//scan for pvalues
		System.out.print("\nScanning windows of mCG for stranded methylation... ");

		for (int i=0; i< chromosomes.length; i++){		
			chromosome = chromosomes[i];	

			//fetch data
			fetchDataAndRemove();

			//check
			if (convertedMergedChromPlus == null || convertedMergedChromMinus == null || nonConvertedMergedChromPlus == null || nonConvertedMergedChromMinus == null) {
				System.out.print(" - Couldn't find all 4 datasets, skipping!");
				continue;
			}

			//generate windows
			windowScanChromosome();

			//print graphs?
			if (printGraphs && smoothingWindow != null && smoothingWindow.length !=0) writeBarFileGraphsOnCurrentSM();

			//merge windows		
			mergeWindows();

			//filter windows for log2ratio and pval thresholds
			filterWindows();

			//save windows
			File obFile = new File(windowsDirectory,chromosome);
			IO.saveObject(obFile, smoothingWindow);
			swiFile.add(obFile);
			obFile.deleteOnExit();
		}
		System.out.println();

		//B & H correct pvalues and make EnrichedRegions
		System.out.println("\nApplying a Benjamini and Hochberg FDR correction to the p-values...\n");
		correctPValues();

		//any enriched regions?
		if (enrichedRegions == null) System.out.println("No Enriched/Reduced Regions found that pass the log2Ratio ("+ log2RatioThreshold +") and FDR ("+ fdrThreshold +") thresholds.");
		else {
			System.out.print("Loading "+enrichedRegions.length+" Enriched/Reduced Regions with count data... ");
			loadEnrichedRegions();
			//sort by abs log2ratio
			Arrays.sort(enrichedRegions, new ComparatorEnrichedRegionAbsoluteScore(log2RatioIndex));
			printEnrichedRegions("Win");
		}
	}

	private void mergeWindows(){
		//make array of SmoothingWindowIndex and sort by score largest to smallest
		smoothingWindowIndex = SmoothingWindowIndex.makeSmoothingWindowIndex(smoothingWindow);
		Arrays.sort(smoothingWindowIndex, new ComparatorSmoothingWindowIndexScore(pValIndex));
		//find max window and expand
		for (int i=0; i< smoothingWindowIndex.length; i++){
			//is SmoothingWindow object null?
			if (smoothingWindowIndex[i].getSmoothingWindow().getScores() != null){
				expand(smoothingWindowIndex[i].getIndex());
			}
		}
		//reset smoothingWindow[]
		ArrayList<SmoothingWindow> alSW = new ArrayList<SmoothingWindow>();
		for (int i=0; i< smoothingWindow.length; i++){
			float[] scores = smoothingWindow[i].getScores();
			if (scores != null) alSW.add(smoothingWindow[i]);
		}
		smoothingWindow = new SmoothingWindow[alSW.size()];
		alSW.toArray(smoothingWindow);
	}


	private void expand(int indexWindow){
		SmoothingWindow win = smoothingWindow[indexWindow];
		float[] scores = win.getScores();
		float minPVal = scores[pValIndex] * windowMergeMultiplier;
		float numWindows = 1;
		//expand left
		int index = indexWindow;
		while (true){
			//fetch left window
			index--;
			if (index < 0) break;
			SmoothingWindow leftWin = smoothingWindow[index];
			//does it overlap
			if (overlap(win, leftWin) == false) break;
			//expand?
			if (checkScore(minPVal, leftWin)) {
				//reset left position
				win.setStart(leftWin.getStart());
				//kill the left wins score so that it isn't considered later
				leftWin.setScores(null);
				numWindows++;
			}
			else break;
		}

		//expand right
		index = indexWindow;
		while (true){
			//fetch left window
			index++;
			if (index == smoothingWindow.length) break;
			SmoothingWindow rightWin = smoothingWindow[index];
			//does it overlap
			if (overlap(win, rightWin) == false) break;
			//expand?
			if (checkScore(minPVal, rightWin)) {
				//reset right position
				win.setStop(rightWin.getStop());
				//kill the right wins score so that it isn't considered later
				rightWin.setScores(null);
				numWindows++;
			}
			else break;
		}
		//add on end scores the number merged
		scores = Num.appendFloat(scores, numWindows);
		win.setScores(scores);
	}

	private boolean overlap(SmoothingWindow a, SmoothingWindow b){
		if (a.getStop() < b.getStart() || b.getStop() < a.getStart()) return false;
		return true;
	}

	private boolean checkScore(float minScore, SmoothingWindow win){
		float[] scores = win.getScores();
		if (scores == null || scores.length > 2 || scores[pValIndex] < minScore) return false;
		return true;
	}

	public void bpScan(){
		loadDataHashes();

		//this also removes the phiX and lambda data from that to be scanned if resent
		fetchAllChromosomes();

		//scan for pvalues
		System.out.print("Scanning bp chromosomes... ");
		String oldChrom = "";
		for (int i=0; i< chromosomes.length; i++){		
			chromosome = chromosomes[i];	
			//check if present
			//fetch the chrom sequence?
			if (oldChrom.equals(chromosome) ==false){
				File seqFile = chromosomeFastaFiles.get(chromosome);
				if (seqFile == null) continue;
				else System.out.print(chromosome+" ");
				MultiFastaParser mfp = new MultiFastaParser(seqFile);				
				genomicSequence = mfp.getSeqs()[0].toUpperCase();
				oldChrom = chromosome;
			}
			//fetch data
			fetchDataAndRemove();

			//check
			if (convertedMergedChromPlus == null || convertedMergedChromMinus == null || nonConvertedMergedChromPlus == null || nonConvertedMergedChromMinus == null) {
				System.out.print(" - Couldn't find all 4 datasets, skipping!");
				continue;
			}
			//calculate base level stats
			baseScanCGsForPValues();

			//print graphs?
			if (printGraphs && smoothingWindow != null && smoothingWindow.length !=0) writeBarFileGraphsOnCurrentSM();

			//filter windows for log2ratio and pval thresholds
			filterWindows();

			//save windows
			File obFile = new File(windowsDirectory,chromosome);
			IO.saveObject(obFile, smoothingWindow);
			swiFile.add(obFile);
			obFile.deleteOnExit();
		}
		System.out.println();

		//B & H correct pvalues and make EnrichedRegions
		System.out.println("\nApplying a Benjamini and Hochberg FDR correction to the p-values...\n");
		correctPValues();

		//any enriched regions?
		if (enrichedRegions == null) System.out.println("No Enriched/Reduced Regions found that pass the log2Ratio ("+ log2RatioThreshold +") and FDR ("+ fdrThreshold +") thresholds.");
		else {
			System.out.print("Loading "+enrichedRegions.length+" Enriched/Reduced Regions with count data... ");
			loadEnrichedRegions();
			//sort by abs log2ratio
			Arrays.sort(enrichedRegions, new ComparatorEnrichedRegionAbsoluteScore(log2RatioIndex));
			printEnrichedRegions("Bp");
		}
	}

	/**Returns nonConPlus, conPlus, nonConMinus, conMinus*/
	private double[] scoreRegion (int start, int stop){
		double conPlus = convertedMergedChromPlus.sumScoreBP(start, stop);
		double nonConPlus = nonConvertedMergedChromPlus.sumScoreBP(start, stop);
		double conMinus = convertedMergedChromMinus.sumScoreBP(start, stop); 
		double nonConMinus = nonConvertedMergedChromMinus.sumScoreBP(start, stop);

		return new double[]{nonConPlus, conPlus, nonConMinus, conMinus};
	}

	/**Returns nonConPlus, conPlus, nonConMinus, conMinus for those with min 2 obs*/
	private double[] scoreRegionWithMinObs (int start, int stop){
		//fetch Point[]
		Point[] conPlus = convertedMergedChromPlus.fetchPoints(start, stop);
		Point[] nonConPlus = nonConvertedMergedChromPlus.fetchPoints(start, stop);
		Point[] conMinus = convertedMergedChromMinus.fetchPoints(start, stop+1); 
		Point[] nonConMinus = nonConvertedMergedChromMinus.fetchPoints(start, stop+1);
		//obs
		double numConPlus = 0;
		double numNonConPlus = 0;
		double numConMinus = 0;
		double numNonConMinus = 0;
		//indexes
		int indexConPlus =0;
		int indexNonConPlus=0;
		int indexConMinus =0;
		int indexNonConMinus =0;

		//for each position
		for (int i=start; i< stop; i++){
			int position = i;
			int positionPlusOne = i +1;
			//ConPlus
			float countsConPlus = 0;
			if (conPlus != null){
				for (int x=indexConPlus; x< conPlus.length; x++){
					int testPos = conPlus[x].getPosition();
					//same position
					if (testPos == position) {
						countsConPlus = conPlus[x].getScore();
						indexConPlus = x +1;
						break;
					}
					//greater than
					else if (testPos > position) {
						indexConPlus = x;
						break;
					}
					//less than so do nothing
				}
			}
			//NonConPlus
			float countsNonConPlus = 0;
			if (nonConPlus != null){
				for (int x=indexNonConPlus; x< nonConPlus.length; x++){
					int testPos = nonConPlus[x].getPosition();
					//same position
					if (testPos == position) {
						countsNonConPlus = nonConPlus[x].getScore();
						indexNonConPlus = x +1;
						break;
					}
					//greater than
					else if (testPos > position) {
						indexNonConPlus = x;
						break;
					}
					//less than so do nothing
				}
			}

			//enough obs on plus strand
			if ((countsConPlus+ countsNonConPlus) < minimumCoverage) continue;
			
			//ConMinus
			float countsConMinus = 0;
			if (conMinus != null){
				for (int x=indexConMinus; x< conMinus.length; x++){
					int testPos = conMinus[x].getPosition();
					//same position
					if (testPos == positionPlusOne) {
						countsConMinus = conMinus[x].getScore();
						indexConMinus = x +1;
						break;
					}
					//greater than
					else if (testPos > positionPlusOne) {
						indexConMinus = x;
						break;
					}
					//less than so do nothing
				}
			}
			//NonConMinus
			float countsNonConMinus = 0;
			if (nonConMinus != null){
				for (int x=indexNonConMinus; x< nonConMinus.length; x++){
					int testPos = nonConMinus[x].getPosition();
					//same position
					if (testPos == positionPlusOne) {
						countsNonConMinus = nonConMinus[x].getScore();
						indexNonConMinus = x +1;
						break;
					}
					//greater than
					else if (testPos > positionPlusOne) {
						indexNonConMinus = x;
						break;
					}
					//less than so do nothing
				}
			}
			
			//enough obs on minus strand
			if ((countsConMinus+ countsNonConMinus) < minimumCoverage) continue;

			numConPlus += countsConPlus;
			numNonConPlus += countsNonConPlus;
			numConMinus += countsConMinus;
			numNonConMinus += countsNonConMinus;
			
		}
		return new double[]{numNonConPlus, numConPlus, numNonConMinus, numConMinus};
	}



	/**Assumes nonConPlus, conPlus, nonConMinus, conMinus*/
	private double calculateLog2Ratio(double[] counts){
		double tNumCon = counts[1]+1;
		double tNumNonCon = counts[0]+1;
		double cNumCon = counts[3]+1;
		double cNumNonCon = counts[2]+1;

		//calculate log2Ratio
		double tMethylated = tNumNonCon/(tNumNonCon+tNumCon);
		double cMethylated = cNumNonCon/(cNumNonCon+cNumCon);
		return Num.log2(tMethylated/cMethylated);
	}


	public void loadEnrichedRegions(){
		//load data hashes
		loadDataHashes();
		chromosome = "";
		for (int i=0; i< enrichedRegions.length; i++){
			String erChromosome = enrichedRegions[i].getChromosome();
			//load PointData?
			if (erChromosome.equals(chromosome) == false){
				chromosome = erChromosome;
				System.out.print(chromosome +" ");
				fetchDataAndRemove();
			}
			//load ER
			loadEnrichedRegion(enrichedRegions[i]);
		}
		System.out.println();
	}

	public void loadEnrichedRegion(EnrichedRegion er){
		//fetch er stats
		double[] erCounts = scoreRegionWithMinObs(er.getStart(), er.getStop());
		double erLog2Ratio = calculateLog2Ratio (erCounts);
		er.setLog2Ratio((float)erLog2Ratio);
		er.setCounts(Num.doubleArrayToIntArray(erCounts));
		//fetch best window stats
		SmoothingWindow bestWindow = er.getBestWindow();
		double[] bwCounts = scoreRegionWithMinObs(bestWindow.getStart(), bestWindow.getStop());
		//there are two, log2Ratio and FDR, append on counts numNonConPlus, numConPlus, numNonConMinus, numConMinus
		float[] scores = bestWindow.getScores();
		float[] allScores = new float[]{scores[log2RatioIndex], scores[pValIndex], (float)bwCounts[0], (float)bwCounts[1], (float)bwCounts[2], (float)bwCounts[3]};
		bestWindow.setScores(allScores);
	}


	public void printEnrichedRegions(String winBp){
		printEnrichedRegionSpreadSheet(winBp);
		//split by log2Ratio into enriched and reduced
		EnrichedRegion[] ers = null;
		ArrayList<EnrichedRegion> ersAL = new ArrayList<EnrichedRegion>();
		//enriched
		for (int i=0; i< enrichedRegions.length; i++){
			float[] scores = enrichedRegions[i].getBestWindow().getScores();
			if (scores[log2RatioIndex] > 0) ersAL.add(enrichedRegions[i]);
		}
		if (ersAL.size() != 0){
			ers = new EnrichedRegion[ersAL.size()];
			ersAL.toArray(ers);
			printEGRFile(ers, new File (saveDirectory, "plusStrMethylRegions"+winBp+".egr"));
			printBedFile(ers, new File (saveDirectory, "plusStrMethylRegions"+winBp+".bed"), "+");
			ersAL.clear();
		}
		//reduced
		for (int i=0; i< enrichedRegions.length; i++){
			float[] scores = enrichedRegions[i].getBestWindow().getScores();
			if (scores[log2RatioIndex] < 0) ersAL.add(enrichedRegions[i]);
		}
		if (ersAL.size() != 0){
			ers = new EnrichedRegion[ersAL.size()];
			ersAL.toArray(ers);
			printEGRFile(ers, new File (saveDirectory, "minusStrMethylRegions"+winBp+".egr"));
			printBedFile(ers, new File (saveDirectory, "minusStrMethylRegions"+winBp+".bed"), "-");
			ersAL.clear();
		}

	}

	/**Writes out an egr file, used by IGB for graph display.*/
	public  void printEGRFile (EnrichedRegion[] egrs, File file){
		try{
			PrintWriter out = new PrintWriter (new FileWriter (file));
			//genome version
			out.println("# genome_version = "+genomeVersion);

			//score names
			out.println("# score0 = log2(plusMethyl/minusMethyl)");
			out.println("# score1 = FDR");

			//transform scores to 1:1000
			float[][] tScores = new float[2][egrs.length];
			//for each ER
			for (int i=0; i< egrs.length; i++){
				float[] scores = egrs[i].getBestWindow().getScores();
				for (int j=0; j< 2; j++){
					tScores[j][i] = scores[j];
				}
			}
			//transform
			for (int i=0; i< tScores.length; i++) Num.scale1To1000(tScores[i]);

			//for each ER
			for (int i=0; i< egrs.length; i++){
				out.print(egrs[i].getChromosome());
				out.print("\t");
				out.print(egrs[i].getStart());
				out.print("\t");
				out.print(egrs[i].getStop());
				out.print("\t.\t");
				out.print(tScores[0][i]);
				for (int j=1; j< tScores.length; j++){
					out.print("\t");
					out.print(tScores[j][i]);
				}
				out.println();
			}
			out.close();
		}catch (Exception e){
			e.printStackTrace();
		}
	}

	/**Writes out an egr file, used by IGB for graph display.*/
	public  void printBedFile (EnrichedRegion[] egrs, File file, String strand){
		try{
			PrintWriter out = new PrintWriter (new FileWriter (file));
			//genome version
			out.println("# genome_version = "+genomeVersion);

			//score names
			out.println("# name = log2(tMethyl/cMethyl) : -10Log10(FDR)");
			out.println("# score = scaled 1:1000 log2(plusMethyl/minusMethyl)");

			//transform log2Ratio scores to 1:1000
			float[] tScores = new float[egrs.length];
			//for each ER
			for (int i=0; i< egrs.length; i++){
				float[] scores = egrs[i].getBestWindow().getScores();
				tScores[i] = scores[log2RatioIndex];
			}
			//transform
			Num.scale1To1000(tScores);

			//for each ER
			for (int i=0; i< egrs.length; i++){
				float[] scores = egrs[i].getBestWindow().getScores();
				out.print(egrs[i].getChromosome());
				out.print("\t");
				out.print(egrs[i].getStart());
				out.print("\t");
				out.print(egrs[i].getStop());
				out.print("\t");
				out.print(Num.formatNumber(scores[log2RatioIndex], 2)+":"+(int)scores[pValIndex]);
				out.print("\t");
				out.print(tScores[i]);
				out.println("\t"+strand);
			}
			out.close();
		}catch (Exception e){
			e.printStackTrace();
		}
	}


	/**Writes out an excel compatible tab delimited spreadsheet with hyperlinks for IGB.*/
	public void printEnrichedRegionSpreadSheet(String winBp){
		try{
			File file = new File(saveDirectory, "diffMethylRegions"+winBp+".xls");
			PrintWriter out = new PrintWriter (new FileWriter (file));
			//print header line
			out.println("#"+genomeVersion+"_IGBHyperLinks\tChr\tStart\tStop\t#Windows\tLog2(plusMethyl/minusMethyl)\t#PlusNonCon\t#PlusCon\t#MinusNonCon\t#MinusCon\t"+
			"BW_Start\tBW_Stop\tBW_Log2(plusMethyl/minusMethyl)\tBW_FDR\t#BW_PlusNonCon\t#BW_PlusCon\t#BW_MinusNonCon\t#BW_MinusCon");

			String url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";
			String tab = "\t";
			//for each ER
			for (int i=0; i< enrichedRegions.length; i++){
				//url
				int winStart = enrichedRegions[i].getStart() - 10000;
				if (winStart < 0) winStart = 0;
				int winEnd = enrichedRegions[i].getStop() + 10000;
				out.print(url+enrichedRegions[i].getChromosome()+"&start="+winStart+"&end="+winEnd+"\",\""+(i+1)+"\")\t");
				//chrom
				out.print(enrichedRegions[i].getChromosome()); out.print(tab);
				//start
				out.print(enrichedRegions[i].getStart()); out.print(tab);
				//stop
				out.print(enrichedRegions[i].getStop()); out.print(tab);
				//#Windows
				out.print(enrichedRegions[i].getNumberOfWindows()); out.print(tab);
				//Log2Ratio
				out.print(enrichedRegions[i].getLog2Ratio()); out.print(tab);
				//Counts
				int[] erCounts = enrichedRegions[i].getCounts();
				for (int x=0; x< erCounts.length; x++){
					out.print(erCounts[x]);
					out.print(tab);
				}
				//best window data
				SmoothingWindow bw = enrichedRegions[i].getBestWindow();
				//coordinates
				out.print(bw.getStart()); out.print(tab);
				out.print(bw.getStop()); out.print(tab);
				//scores
				float[] scores = bw.getScores();
				//	log2Ratio
				out.print(scores[0]); out.print(tab);
				//	FDR
				out.print(scores[1]);
				//	counts
				for (int j=2; j< scores.length; j++){
					out.print(tab);
					out.print((int)scores[j]);
				}
				out.println();
			}
			out.close();
		}catch (Exception e){
			System.out.println("\nError: problem printing spreadsheet report");
			e.printStackTrace();
		}
	}

	private void filterWindows(){
		ArrayList<SmoothingWindow> goodWin = new ArrayList<SmoothingWindow>();
		for (int i=0; i< smoothingWindow.length; i++){
			float[] scores = smoothingWindow[i].getScores();
			if (Math.abs(scores[log2RatioIndex]) < log2RatioThreshold || scores[pValIndex]< fdrThreshold) numberSkippedWindows++;
			else goodWin.add(smoothingWindow[i]);
		}
		smoothingWindow = new SmoothingWindow[goodWin.size()];
		goodWin.toArray(smoothingWindow);
		numberPassingWindows+= smoothingWindow.length;
	}

	public void correctPValues(){
		//load Point[] with all of the pvalues
		Point[] point = new Point[numberPassingWindows];
		//for each chromosome of data
		int index = 0;
		for (int i=0; i< swiFile.size(); i++){
			SmoothingWindow[] sw = (SmoothingWindow[])IO.fetchObject(swiFile.get(i));
			for (int j=0; j<sw.length; j++){
				point[index] = new Point(index, sw[j].getScores()[pValIndex]);
				index++;
			}
		}
		//sort by score smallest -10Log10(pval) to largest
		Arrays.sort(point, new ComparatorPointAscendingScore());
		//correct
		Point.benjaminiHochbergCorrect(point, numberSkippedWindows);
		//sort back to original position
		Arrays.sort(point, new ComparatorPointPosition());
		//make ER generator threshold on fdr (already thresholded on log2ratio)
		EnrichedRegionMaker erm = new EnrichedRegionMaker(maxGap, new int[]{pValIndex}, new float[]{fdrThreshold}, pValIndex);
		ComparatorSmoothingWindowPosition comparator = new ComparatorSmoothingWindowPosition();
		boolean noSWs = true;

		//for each chromosome of data
		index = 0;
		for (int i=0; i< swiFile.size(); i++){
			ArrayList<SmoothingWindow> good = new ArrayList<SmoothingWindow>();
			SmoothingWindow[] sw = (SmoothingWindow[])IO.fetchObject(swiFile.get(i));
			chromosome = swiFile.get(i).getName();

			//swap pvals for fdrs
			for (int j=0; j<sw.length; j++){				
				float[] scores = sw[j].getScores();
				scores[pValIndex] = point[index++].getScore();
				//save it? threshold on fdr, already thresholded on log2ratio
				if (scores[pValIndex] >= fdrThreshold) good.add(sw[j]);
				else sw[j] = null;
			}
			//convert als to [] and make ERS
			//	enriched regions
			if (good.size() !=0){
				noSWs = false;
				sw = new SmoothingWindow[good.size()];
				good.toArray(sw);
				Arrays.sort(sw, comparator);
				erm.addEnrichedRegions(sw, chromosome, log2RatioIndex);
			}
		}
		//make enriched regions? otherwise leave null
		if (noSWs == false) {
			enrichedRegions = new EnrichedRegion[erm.getEnrichedRegionsAL().size()];
			erm.getEnrichedRegionsAL().toArray(enrichedRegions);
		}
	}


	/**Fetches the names of all the chromosomes in the data excluding lambda and phiX if present.*/
	public void fetchAllChromosomes(){
		HashSet<String> c = new HashSet<String>();
		Iterator<String> it = convertedPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = convertedMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = nonConvertedPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = nonConvertedMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());

		//any lambda or phiX found?
		it = c.iterator();
		Pattern lambda = Pattern.compile(".*lambda.*", Pattern.CASE_INSENSITIVE);
		Pattern phiX = Pattern.compile(".*phix.*", Pattern.CASE_INSENSITIVE);
		String lambdaChromosome = null;
		String phiXChromosome = null;
		while (it.hasNext()){
			String chr = it.next();
			if (lambda.matcher(chr).matches()) lambdaChromosome = chr;
			else if (phiX.matcher(chr).matches()) phiXChromosome = chr;
		}
		if (lambdaChromosome != null) c.remove(lambdaChromosome);
		if (phiXChromosome != null) c.remove(phiXChromosome);
		chromosomes=  Misc.hashSetToStringArray(c);
	}

	/**Fetchs the data for a particular chromosome.*/
	public void fetchDataAndRemove(){
		ArrayList<PointData> al = null;
		PointData[] pd = null;
		//merge converted
		convertedMergedChromPlus = null;
		if (convertedPlusPointData.containsKey(chromosome)) {
			pd = convertedPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			convertedMergedChromPlus = PointData.mergePointData(al, false, true);
		}
		//set version
		if (genomeVersion == null) genomeVersion = pd[0].getInfo().getVersionedGenome();

		convertedMergedChromMinus = null;
		if (convertedMinusPointData.containsKey(chromosome)) {
			pd = convertedMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			convertedMergedChromMinus = PointData.mergePointData(al, false, true);
		}
		//merge nonConverted
		nonConvertedMergedChromPlus = null;
		if (nonConvertedPlusPointData.containsKey(chromosome)) {
			pd = nonConvertedPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			nonConvertedMergedChromPlus = PointData.mergePointData(al, false, true);
		}
		nonConvertedMergedChromMinus = null;
		if (nonConvertedMinusPointData.containsKey(chromosome)) {
			pd = nonConvertedMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			nonConvertedMergedChromMinus = PointData.mergePointData(al, false, true);
		}
		pd = null;
		al = null;
	}

	/**Window scans a chromosome.*/
	public void windowScanChromosome(){
		//fetch and merge positions 
		int[] positions = fetchMergePositions();
		//make windows using all of the reads
		makeWindows(positions);
		//any windows?
		if (windows.length == 0){
			System.out.println("\n\tSkipping "+chromosome+". No windows found with minimum reads of "+minimumNumberCsInWindow+" within a window size of "+windowSize);
			return;
		}
		//scan
		scoreWindows();
	}

	/**Scores a chromosome for window based standed methylation.*/
	private void scoreWindows(){
		ArrayList<SmoothingWindow> pValueWindows = new ArrayList<SmoothingWindow>();
		ArrayList<SmoothingWindow> forRWindows = new ArrayList<SmoothingWindow>();

		//for each window
		for (int i=0; i< windows.length; i++){
			//score region : nonConPlus, conPlus, nonConMinus, conMinus
			double[] counts = scoreRegionWithMinObs(windows[i][0], windows[i][1]);

			//check for zero obs in either dataset
			if (counts[0] == 0 && counts[1] ==0) continue;
			if (counts[2] == 0 && counts[3] ==0) continue;

			//calculate log2Ratio
			double log2Ratio = calculateLog2Ratio(counts);

			//calculate p-value for differences in methylation
			float[] scores = null;
			double totalObservations = Num.sumArray(counts);

			//use fischer's?
			if (totalObservations < cacheNumber) {
				double pLog = Num.minus10log10(fischerExact.getTwoTailedP((int)counts[1], (int)counts[0], (int)counts[3], (int)counts[2]));
				if (BisSeq.pValCheck(pLog) == false) continue;
				scores = new float[]{(float)log2Ratio, (float)pLog};
			}

			//use chi-square?
			else {
				//any low count cells
				if (Num.findMinMaxDoubleValues(counts)[0] < 5) continue;
				//too many for java chi-square?
				long[][] c = new long[][]{{(long)counts[1], (long)counts[0]},{(long)counts[3], (long)counts[2]}};
				double chiSquareStat = Num.chiSquareTestStatistic(c);
				if (chiSquareStat > 1400) scores = new float[]{(float)log2Ratio, (float)counts[1], (float)counts[0], (float)counts[3], (float)counts[2]};
				//use apache chi-square
				else{
					double pNoLog = -1;
					try {
						pNoLog = chiSquare.chiSquareTest(c);
					} catch (Exception e) {}
					if (BisSeq.pValCheck(pNoLog) == false) scores = new float[]{(float)log2Ratio, (float)counts[1], (float)counts[0], (float)counts[3], (float)counts[2]};
					else scores = new float[]{(float)log2Ratio, Num.minus10log10Float(pNoLog)};
				}
			}

			//make smoothing window
			SmoothingWindow win = new SmoothingWindow (windows[i][0], windows[i][1], scores);
			if (scores.length == 2) {
				if (scores[pValIndex] !=0) pValueWindows.add(win);
			}
			else forRWindows.add(win);			
		}


		//make []
		int numPValWindows = pValueWindows.size();
		int numRWindows = forRWindows.size();
		smoothingWindow = new SmoothingWindow[numPValWindows + numRWindows];
		for (int i=0; i< numPValWindows; i++) smoothingWindow[i] = pValueWindows.get(i);

		//any difficult counts? Returns the window with log2 and pvalue
		if (numRWindows !=0) {				
			//use R to calculate
			//BisSeq.calculateChiSquareInRWithLookup(forRWindows, saveDirectory, fullPathToR, log2RatioIndex);
			//add to []
			int index = 0;
			for (int i=numPValWindows; i< smoothingWindow.length; i++) {
				smoothingWindow[i] = forRWindows.get(index++);
			}				
			//sort			
			Arrays.sort(smoothingWindow, new ComparatorSmoothingWindowPosition());
		}	
	}


	/**Scores a chromosome for non-converted to total at base level.*/
	private void baseScanCGsForPValues(){
		ArrayList<SmoothingWindow> pValueWindows = new ArrayList<SmoothingWindow>();
		ArrayList<SmoothingWindow> forRWindows = new ArrayList<SmoothingWindow>();

		//fetch arrays
		int[] positionsNonConPlus = nonConvertedMergedChromPlus.getPositions();
		float[] readsNonConPlus = nonConvertedMergedChromPlus.getScores();
		int[] positionsConPlus = convertedMergedChromPlus.getPositions();
		float[] readsConPlus = convertedMergedChromPlus.getScores();

		int[] positionsNonConMinus = nonConvertedMergedChromMinus.getPositions();
		float[] readsNonConMinus = nonConvertedMergedChromMinus.getScores();
		int[] positionsConMinus = convertedMergedChromMinus.getPositions();
		float[] readsConMinus = convertedMergedChromMinus.getScores();

		//for each position in the methylated plus strand data, look on minus strand CpG for observations
		int indexConPlus =0;
		int indexNonConMinus =0;
		int indexConMinus =0;
		for (int i=0; i< positionsNonConPlus.length; i++){
			//bp positions
			int testPos = positionsNonConPlus[i];
			int testPosPlusOne = testPos +1;

			//is it a CpG?
			String genSeq = genomicSequence.substring(testPos,testPos+2);
			if (genSeq.equals("CG") == false) {
				nonCGFound = true;
				continue;
			}

			//counts
			float numNonConPlus = readsNonConPlus[i];
			float numConPlus =0;
			float numNonConMinus =0;
			float numConMinus =0;

			//plus strand
			//present in conPlus?
			for (int j=indexConPlus; j < positionsConPlus.length; j++){
				//match!
				if (testPos == positionsConPlus[j]){
					numConPlus = readsConPlus[j];
					indexConPlus++;
					break;
				}
				//less than
				if (testPos < positionsConPlus[j]) break;
				//greater than so keep advancing
				indexConPlus = j;
			}

			//minus strand
			//present in nonConMinus?
			for (int j=indexNonConMinus; j < positionsNonConMinus.length; j++){
				//match!
				if (testPosPlusOne == positionsNonConMinus[j]){
					numNonConMinus = readsNonConMinus[j];
					indexNonConMinus++;
					break;
				}
				//less than
				if (testPosPlusOne < positionsNonConMinus[j]) break;
				//greater than so keep advancing
				indexNonConMinus = j;
			}
			//present in conMinus?
			for (int j=indexConMinus; j < positionsConMinus.length; j++){
				//match!
				if (testPosPlusOne == positionsConMinus[j]){
					numConMinus = readsConMinus[j];
					indexConMinus++;
					break;
				}
				//less than
				if (testPosPlusOne < positionsConMinus[j]) break;
				//greater than so keep advancing
				indexConMinus = j;
			}

			//enough coverage?
			float numPlus = numNonConPlus + numConPlus;
			float numMinus = numNonConMinus + numConMinus;

			if (numPlus < minimumCoverage || numPlus < minimumCoverage) continue;

			//calculate log2Ratio
			double plusMethylated = (numNonConPlus+1)/(numNonConPlus+numConPlus+2);
			double minusMethylated = (numNonConMinus+1)/(numNonConMinus+numConMinus+2);
			double log2Ratio = Num.log2(plusMethylated/minusMethylated);

			//calculate p-value for differences in methylation
			float[] scores = null;
			double totalObservations = numPlus + numMinus;

			//use fischer's?
			if (totalObservations < cacheNumber) {
				double pLog = Num.minus10log10(fischerExact.getTwoTailedP((int)numConPlus, (int)numNonConPlus, (int)numConMinus, (int)numNonConMinus));
				if (BisSeq.pValCheck(pLog) == false) continue;
				scores = new float[]{(float)log2Ratio, (float)pLog};
			}

			//use chi-square?
			else {
				//any low count cells
				if (numConPlus < 5 || numNonConPlus < 5 || numConMinus < 5 || numNonConMinus < 5) continue;
				//too many for java chi-square?
				long[][] counts = new long[][]{{(long)numConPlus, (long)numNonConPlus},{(long)numConMinus, (long)numNonConMinus}};
				double chiSquareStat = Num.chiSquareTestStatistic(counts);
				if (chiSquareStat > 1400) scores = new float[]{(float)log2Ratio, (float)numConPlus, (float)numNonConPlus, (float)numConMinus, (float)numNonConMinus};
				//use apache chi-square
				else{
					double pNoLog = -1;
					try {
						pNoLog = chiSquare.chiSquareTest(counts);
					} catch (Exception e) {}
					if (BisSeq.pValCheck(pNoLog) == false) scores = new float[]{(float)log2Ratio, (float)numConPlus, (float)numNonConPlus, (float)numConMinus, (float)numNonConMinus};
					else scores = new float[]{(float)log2Ratio, Num.minus10log10Float(pNoLog)};
				}
			}

			//make smoothing window
			SmoothingWindow win = new SmoothingWindow (testPos, testPos+2, scores);
			if (scores.length == 2) {
				if (scores[pValIndex] !=0) pValueWindows.add(win);
			}
			else forRWindows.add(win);			
		}

		//make []
		int numPValWindows = pValueWindows.size();
		int numRWindows = forRWindows.size();
		smoothingWindow = new SmoothingWindow[numPValWindows + numRWindows];
		for (int i=0; i< numPValWindows; i++) smoothingWindow[i] = pValueWindows.get(i);

		//any difficult counts? Returns the window with log2 and pvalue
		if (numRWindows !=0) {				
			//use R to calculate
			//BisSeq.calculateChiSquareInRWithLookup(forRWindows, saveDirectory, fullPathToR, log2RatioIndex);
			//add to []
			int index = 0;
			for (int i=numPValWindows; i< smoothingWindow.length; i++) {
				smoothingWindow[i] = forRWindows.get(index++);
			}				
			//sort			
			Arrays.sort(smoothingWindow, new ComparatorSmoothingWindowPosition());
		}
	}

	/**Collects and calculates a bunch of stats re the PointData.*/
	private void loadDataHashes(){
		//fetch converted PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (convertedPointDirs);
		convertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		convertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
		combo = PointData.fetchStrandedPointDataNoMerge (nonConvertedPointDirs);
		nonConvertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		nonConvertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
	}

	/**Fetches all of the positions and merges them.*/
	private int[] fetchMergePositions(){
		ArrayList<int[]> posAL = new ArrayList<int[]>();
		posAL.add(nonConvertedMergedChromPlus.getPositions());
		posAL.add(nonConvertedMergedChromMinus.getPositions());
		posAL.add(convertedMergedChromMinus.getPositions());
		posAL.add(convertedMergedChromMinus.getPositions());

		//merge
		int[][] toMerge = new int[posAL.size()][];
		for (int i=0; i< posAL.size(); i++) toMerge[i] = posAL.get(i);

		return Num.returnUniques(toMerge);
	}

	/**Makes a common set of windows using merged positions.*/
	public void makeWindows(int[] positions){
		windows = windowMaker.makeWindows(positions);
		//assign bp positions
		for (int i=0; i< windows.length; i++){
			windows[i][0] = positions[windows[i][0]];
			windows[i][1] = positions[windows[i][1]]+1;	//last base isn't included
		}
	}

	/**Writes stair step window bar graph files*/
	public void writeBarFileGraphsOnCurrentSM(){

		//build info object
		Info info;
		if (convertedMergedChromPlus != null) info= convertedMergedChromPlus.getInfo();
		else info = convertedMergedChromMinus.getInfo();
		HashMap<String,String> notes = new HashMap<String,String>();
		notes.put(BarParser.DESCRIPTION_TAG, "Log2Ratio (fraction methylated in T/ fraction methylated in C)");
		notes.put(BarParser.UNIT_TAG, "Ratio");
		info.setNotes(notes);
		info.setStrand(".");	

		//save log2 ratio data
		saveSmoothedHeatMapData (0, smoothingWindow, info, log2RatioDirectory, "#FF0000", true); //red

		//save pvalues
		notes.put(BarParser.DESCRIPTION_TAG, "Uncorrected p-values of unstranded differential methylation");
		notes.put(BarParser.UNIT_TAG, "-10Log10(pval)");
		saveSmoothedHeatMapData (1, smoothingWindow, info, pvalueDirectory, "#00FF00", false); //green

	}

	/**Saves bar heatmap/ stairstep graph files*/
	public void saveSmoothedHeatMapData (int scoreIndex, SmoothingWindow[] sm, Info info, File dir, String color, boolean posNeg){
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = info.getNotes();		
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		//color red
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, color);
		//what's window size
		map.put(BarParser.WINDOW_SIZE, windowSize+"");
		//save in info
		info.setNotes(map);
		//get heatmap positions and values
		PointData pd;		
		if (posNeg){
			HeatMapMakerPosNeg hm = new HeatMapMakerPosNeg(scoreIndex, 0, 0);
			pd = hm.makeHeatMapPositionValues(sm);
		}
		else {
			HeatMapMaker hm = new HeatMapMaker(scoreIndex,0);
			pd = hm.makeHeatMapPositionValues(sm, false);
		}

		pd.setInfo(info);
		pd.writePointData(dir);
		//clean up
		pd.nullPositionScoreArrays();
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new StrandedBisSeq(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		File[] fastas = null;
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': fastas = IO.extractFiles(new File(args[++i])); break;
					case 'c': convertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'n': nonConvertedPointDirs = IO.extractFiles(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'w': windowSize = Integer.parseInt(args[++i]); break;
					case 'm': minimumNumberCsInWindow = Integer.parseInt(args[++i]); break;
					case 'l': log2RatioThreshold = Float.parseFloat(args[++i]); break;
					case 'p': fdrThreshold = Float.parseFloat(args[++i]); break;	
					case 'g': printGraphs = true; break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'x': maxGap = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//look for fasta files
		if (fastas == null || fastas.length ==0) Misc.printErrAndExit("\nError: cannot find any fasta sequence files?\n");

		//look for point directories
		if (convertedPointDirs == null || convertedPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your converted PointData directories(s)!\n");
		//only one directory look deeper
		if (convertedPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(convertedPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) convertedPointDirs = otherDirs;
		}

		//nonConverted data
		if (nonConvertedPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your nonConverted PointData directories(s)!\n");
		//only one directory look deeper
		if (nonConvertedPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(nonConvertedPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) nonConvertedPointDirs = otherDirs;
		}

		//make window maker 
		windowMaker = new WindowMaker(windowSize,minimumNumberCsInWindow);

		//look for and or create the save directory
		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		else if (saveDirectory.exists() == false) saveDirectory.mkdirs();
		if (printGraphs) {
			pvalueDirectory = new File (saveDirectory, "BpPValues");
			pvalueDirectory.mkdir();
			log2RatioDirectory = new File (saveDirectory, "BpLog2Ratio");
			log2RatioDirectory.mkdir();
		}
		windowsDirectory = new File (saveDirectory, "BinaryWindowData");
		windowsDirectory.mkdir();


		//load fasta files into hash
		chromosomeFastaFiles = new HashMap<String,File>();
		Pattern chrom = Pattern.compile("(.+)\\.fa.*");
		for (int i=0; i< fastas.length; i++){
			Matcher mat = chrom.matcher(fastas[i].getName());
			if (mat.matches()) chromosomeFastaFiles.put(mat.group(1), fastas[i]);
		}


		//check for R and required libraries
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}

		System.out.println("Params:");
		System.out.println(fdrThreshold+ "\tFDR threshold");
		System.out.println(log2RatioThreshold+ "\tLog2Ratio threshold");
		System.out.println(windowSize+ "\tWindow size");
		System.out.println(minimumNumberCsInWindow+"\tMinimum #C obs in window");
		System.out.println(minimumCoverage + "\tMinimum read coverage over + and - strands to score each CG");
		System.out.println(maxGap+ "\tMax gap between significant CGs to merge");
		System.out.println(windowMergeMultiplier+"\tWindow merge multiplier");
		System.out.println(printGraphs+ "\tPrint graphs");
		System.out.println();

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                 StandedBisSeq: Feb 2011                          **\n" +
				"**************************************************************************************\n" +
				"Looks for strand bias in CG methylation from one dataset using fischer or chi-square\n" +
				"tests followed by a Benjamini and Hochberg FDR correction. Merges significant CGs\n" +
				"within max gap into larger regions. WARNING: many bisulfite datasets display strand\n" +
				"bias due to preferential breakage of C rich regions.  Use this app with caution.\n\n" +

				"Options:\n"+
				"-s Save directory, full path.\n"+
				"-c Converted PointData directories, full path, comma delimited. These should\n" +
				"       contain stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories. These will be merged. Use the ParsePointDataContexts to filter\n" +
				"       for just CG contexts.\n" +
				"-n Non-converted PointData directories, ditto. \n" +
				"-f Fasta files for each chromosome.\n"+
				"\n"+
				"Default Options:\n"+
				"-p Minimimal FDR for stranded methylation, defaults to 30, a -10Log10(FDR = 0.001)\n" +
				"       conversion.\n"+
				"-l Log2Ratio threshold for stranded methylation, defaults to 1.585 (3x).\n"+
				"-w Window size, defaults to 500.\n"+
				"-m Minimum #C obs in window, defaults to 4. \n" +
				"-o Minimum coverage for CG bp methylation scanning, defaults to 2.\n"+
				"-x Max gap between significant CGs to merge, defaults to 500bp.\n"+
				"-g Generate graph files for IGB, defaults to just identifying biased regions.\n"+
				"-r Full path to R, defaults to '/usr/bin/R'\n" +
				"\n"+

				"Example: java -Xmx12G -jar pathTo/USeq/Apps/StandedBisSeq -c /Data/Sperm/Converted -n \n" +
				"      /Data/Sperm/NonConverted -s /Data/Sperm/StrandedBisSeqRes -g -p 20 - l 1 -f\n" +
				"      /Genomes/Hg18/Fastas/ \n\n" +

		"**************************************************************************************\n");

	}
}
