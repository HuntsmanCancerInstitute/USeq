package edu.utah.seq.analysis;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.annotation.*;
import util.gen.*;
import edu.utah.seq.data.*;
import trans.tpmap.*;

/**Scans stranded point data for a skew in the distribution of peaks, a signature of a real chIP-seq peak.  Calculates the bp difference.
 *
 * @author Nix
 * */
public class PeakShiftFinder {

	//user defined fields
	private File[] treatmentPointDirs;
	private File[] controlPointDirs;
	private File resultsDirectory;
	//windowing parameters
	private int windowSize =50;
	private int minimumNumberWindowReads = 10;
	//peak parameters
	private int numberPeaksToMerge = 100;
	private int fivePrimePeakSpan = 500;
	private int threePrimePeakSpan = 1000;
	//thresholds for calling a peak
	private float minDiffScore = 2.5f;
	private double ratioT2CReads = 5;
	private boolean scanReduced = false;

	//internal fields
	private HashMap<String, PointData[]> treatments;
	private HashMap<String, PointData[]> controls;
	private WindowMaker windowMaker;
	private int[][] windows;
	private ArrayList<Point> pointsALPlus;
	private ArrayList<Point> pointsALMinus;
	private PointData[] plusMinusT;
	private PointData[] plusMinusC;
	private ArrayList<StrandedPeak> strandedPeaksAL = new ArrayList<StrandedPeak>();
	private boolean verbose = true;

	//product
	private StrandedPeak[] strandedPeaks;
	private Point[] compositeTPlus;
	private Point[] compositeTMinus;
	private Point[] compositeCPlus;
	private Point[] compositeCMinus;
	private int bpPlusPeak;
	private int bpMinusPeak;
	private double compositePeakShift =0;
	private double medianIndividualPeakShift =0;
	private double stndDevIndividualPeakShift =0;
	private int[] peakDistances;

	//constructors
	
	/**For stand alone app.*/
	public PeakShiftFinder(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		findPeakShift();
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");

	}
	/**For integration with the ChIPSeq app.*/
	public PeakShiftFinder(File[] treatmentPointDirs, File[] controlPointDirs, File resultsDirectory, boolean scanReduced){
		this.treatmentPointDirs = treatmentPointDirs;
		this.controlPointDirs = controlPointDirs;
		this.resultsDirectory = resultsDirectory;
		this.scanReduced = scanReduced;
		verbose = false;
		findPeakShift();
	}
	
	public void findPeakShift(){
			//fetch stranded treatment
			if (verbose) System.out.println("Fetching stranded PointData...");
			treatments = PointData.fetchStrandedCombinePointData(treatmentPointDirs);
			
			//fetch stranded control
			controls = PointData.fetchStrandedCombinePointData(controlPointDirs);

			//instantiate window maker
			windowMaker = new WindowMaker(windowSize,minimumNumberWindowReads);

			//pick peaks 
			if (verbose) System.out.println("\nPicking peaks...");
			pickPeaks();

			//merge peaks into a composite
			int numStrandedPeaks = strandedPeaks.length;
			if (numStrandedPeaks < 10) {
				System.out.println("\nToo few peaks, cannot reliably estimate peak shift.\n");
				return;
			}
			if (numStrandedPeaks < numberPeaksToMerge) {
				if (verbose) System.out.println("Only found "+numStrandedPeaks+" Stranded Peaks, merging...");
				numberPeaksToMerge = numStrandedPeaks;
			}
			else {
				if (verbose) System.out.println("\nMerging top "+numberPeaksToMerge+" Stranded Peaks...");
			}
			mergePeaks();

			//window scan composite
			System.out.println("\nCalculating composite shift...");
			findShift();
			System.out.println("\t"+bpPlusPeak+"\t5' Peak\n\t"+bpMinusPeak+"\t3' Peak\n\t"+compositePeakShift+"\tComposite peak shift");
			
			//calculate median and mean distances
			System.out.println("\nCalculating individual peak shifts for top "+numberPeaksToMerge+"...");
			collectPeakDistances();
			medianIndividualPeakShift = Num.median(peakDistances);
			System.out.println("\t"+Num.formatNumber(medianIndividualPeakShift, 0)+"\tMedian peak shift");
			System.out.println("\t"+Num.formatNumber(Num.mean(peakDistances), 0)+"\tMean peak shift");
			stndDevIndividualPeakShift = Num.standardDeviation(peakDistances);
			System.out.println("\t"+Num.formatNumber(stndDevIndividualPeakShift, 1)+"\tStnd dev peak shifts");
			
			//move dir
			File resultsDirectoryMoved = new File(resultsDirectory+"_"+(int)compositePeakShift+"bp_"+(int)medianIndividualPeakShift+"_bp");
			resultsDirectory.renameTo(resultsDirectoryMoved);
	}
	
	/**Makes a best guess at the peak shift and window size. Returns null if the numbers look unreliable.*/
	public int[] fetchPeakShiftAndWindowSize(){
		if (compositePeakShift < 50 || medianIndividualPeakShift < 50) return null;
		//are the two estimates within 80% of each other?
		double fraction;
		int smaller;
		if (compositePeakShift > medianIndividualPeakShift){
			fraction = medianIndividualPeakShift/compositePeakShift; 
			smaller = (int)medianIndividualPeakShift;
		}
		else {
			fraction = compositePeakShift / medianIndividualPeakShift;
			smaller = (int) compositePeakShift;
		}
		if (fraction < 0.8) return null;
		int windowSize = smaller+ (int)stndDevIndividualPeakShift;
		if (windowSize > (smaller *2)) windowSize = smaller * 2;
		return new int[]{smaller, windowSize};
	}
	
	public void collectPeakDistances(){
		peakDistances = new int[numberPeaksToMerge];
		for (int i=0; i< numberPeaksToMerge; i++){
			peakDistances[i] = strandedPeaks[i].distance;
		}
		Arrays.sort(peakDistances);
	}

	public void findShift(){
		File plusBar = new File (resultsDirectory, "sMergedPeak"+windowSize+"BpWin.bar");
		File minusBar = new File (resultsDirectory, "asMergedPeak"+windowSize+"BpWin.bar");
		//reset window size
		windowMaker = new WindowMaker(windowSize,minimumNumberWindowReads);
		bpPlusPeak = findHighestBPPosition(compositeTPlus, compositeCPlus, plusBar,"+");
		bpMinusPeak = findHighestBPPosition(compositeTMinus, compositeCMinus, minusBar,"-");
		compositePeakShift = bpMinusPeak - bpPlusPeak;
	}

	private int findHighestBPPosition(Point[] t, Point[] c, File barFile, String strand){
		//make windows based on treatment
		PointData tPlus = Point.extractPositionScores(t);
		makeWindows(tPlus.getPositions());

		//score windows returning a list of Point sorted by score with position the index of a window
		PointData cPlus = Point.extractPositionScores(c);
		pointsALPlus = scoreWindows(tPlus, cPlus, false);
		
		//convert point positions to actual bp instead of window indexes
		Point[] points = new Point[pointsALPlus.size()];
		pointsALPlus.toArray(points);	
		 
		for (int i=0; i< points.length; i++){
			int[] win = windows[points[i].getPosition()];
			int center = Num.calculateMiddleIntergenicCoordinates(win[0], win[1]);
			points[i].setPosition(center);
		}
		
		//save to file
		PointData pd = Point.extractPositionScores(points);
		//public Info (String text, String versionedGenome, String chromosome, String strand, int readLength, HashMap<String,String> notes){
		Info info = new Info("USeqPeakComposite", "USeqPeakComposite", "chr1", strand, windowSize, null);
		pd.setInfo(info);
		pd.writePointDataToFile(barFile);		
		

		//find best point and return center position
		Arrays.sort(points, new ComparatorPointDecendingScore());	
		return points[0].getPosition();

	}

	/**Merges top peaks into a single composite.*/
	public void mergePeaks(){
		//composites assign first
		StrandedPeak p = strandedPeaks[0];
		compositeTPlus = p.tPlus;
		compositeTMinus = p.tMinus;
		compositeCPlus = p.cPlus;
		compositeCMinus = p.cMinus;

		for (int i=1; i< numberPeaksToMerge; i++){
			//merge, scores contain the sum of identical positions!
			compositeTPlus = Point.mergePoints(compositeTPlus, strandedPeaks[i].tPlus);
			compositeTMinus = Point.mergePoints(compositeTMinus, strandedPeaks[i].tMinus);
			compositeCPlus = Point.mergePoints(compositeCPlus, strandedPeaks[i].cPlus);
			compositeCMinus = Point.mergePoints(compositeCMinus, strandedPeaks[i].cMinus);
		}
		
		//write composite
		//public Info (String text, String versionedGenome, String chromosome, String strand, int readLength, HashMap<String,String> notes){
		Info info = new Info("USeqPeakComposite", "USeqPeakComposite", "chr1", "+", windowSize, null);
		
		PointData pd = Point.extractPositionScores(compositeTPlus);
		pd.setInfo(info);
		pd.writePointDataToFile(new File (resultsDirectory,"sCompT"+windowSize+"Bp.bar"));
		
		pd = Point.extractPositionScores(compositeCPlus);
		pd.setInfo(info);
		pd.writePointDataToFile(new File (resultsDirectory,"sCompC"+windowSize+"Bp.bar"));
		
		info = new Info("USeqPeakComposite", "USeqPeakComposite", "chr1", "-", windowSize, null);
		pd = Point.extractPositionScores(compositeTMinus);
		pd.setInfo(info);
		pd.writePointDataToFile(new File (resultsDirectory,"asCompT"+windowSize+"Bp.bar"));
		
		pd = Point.extractPositionScores(compositeCMinus);
		pd.setInfo(info);
		pd.writePointDataToFile(new File (resultsDirectory,"asCompC"+windowSize+"Bp.bar"));
		
		
	}

	/**Attempts to pick peaks from each chromosome of stranded point data producing an array of StrandedPeaks[]*/
	public void pickPeaks(){
		//for each chromosome in treatment dataset
		Iterator<String> chroms = treatments.keySet().iterator();
		while (chroms.hasNext()){
			String chrom = chroms.next();
			if (verbose) System.out.println("\t"+chrom);
			//fetch point data
			plusMinusT = treatments.get(chrom);
			plusMinusC = controls.get(chrom);
			//check for plus and minus
			if (plusMinusT== null || plusMinusC == null || plusMinusT[0] == null || plusMinusT[1] == null || plusMinusC[0] == null || plusMinusC[1] == null) {
				System.err.println("Warning: one of your PointData is missing or is not paired with + and - stranded data. Skipping "+chrom);
				continue;
			}

			//subsample so all strands match
			PointData.subSamplePointData(plusMinusT, plusMinusC);
			//make windows on merged treatment strands, convert to base positions
			ArrayList<PointData> mergedAL = new ArrayList<PointData>();
			mergedAL.add(plusMinusT[0]);
			mergedAL.add(plusMinusT[1]);
			PointData merged = PointData.combinePointData(mergedAL, false);
			if (makeWindows(merged.getPositions()) == false){
				if (verbose) System.out.println("\t\tNo windows, skipping!");
				continue;
			}
			
			//assign 1 to all scores
			if (plusMinusT[0] != null) Arrays.fill(plusMinusT[0].getScores(), 1f);
			if (plusMinusT[1] != null) Arrays.fill(plusMinusT[1].getScores(), 1f);
			if (plusMinusC[0] != null) Arrays.fill(plusMinusC[0].getScores(), 1f);
			if (plusMinusC[1] != null) Arrays.fill(plusMinusC[1].getScores(), 1f);

			//score windows from plus strand
			pointsALPlus = scoreWindows(plusMinusT[0],plusMinusC[0], true);
			if (pointsALPlus == null ){
				if (verbose) System.out.println("\t\tNo peaks found, skipping!");
				continue;
			}
			//score windows from minus strand
			pointsALMinus = scoreWindows(plusMinusT[1],plusMinusC[1],true);
			if (pointsALPlus == null){
				if (verbose) System.out.println("\t\tNo peaks found, skipping!");
				continue;
			}
			
			//extract peaks from windows
			int numPeaks = extractPeaks(chrom, false);
			if (verbose) System.out.println("\t\t"+numPeaks+"\t t vs c peaks passing thresholds ");
			
			////// flip t and c and look for reduced peaks?
			if (scanReduced){
			//make windows on merged treatment strands, convert to base positions
			mergedAL = new ArrayList<PointData>();
			mergedAL.add(plusMinusC[0]);
			mergedAL.add(plusMinusC[1]);
			merged = PointData.combinePointData(mergedAL, false);
			if (makeWindows(merged.getPositions()) == false){
				if (verbose) System.out.println("\t\tNo windows, skipping!");
				continue;
			}
			
			//assign 1 to all scores
			if (plusMinusT[0] != null) Arrays.fill(plusMinusT[0].getScores(), 1f);
			if (plusMinusT[1] != null) Arrays.fill(plusMinusT[1].getScores(), 1f);
			if (plusMinusC[0] != null) Arrays.fill(plusMinusC[0].getScores(), 1f);
			if (plusMinusC[1] != null) Arrays.fill(plusMinusC[1].getScores(), 1f);

			//score windows from plus strand
			pointsALPlus = scoreWindows(plusMinusC[0],plusMinusT[0], true);
			if (pointsALPlus == null ){
				if (verbose) System.out.println("\t\tNo peaks found, skipping!");
				continue;
			}
			//score windows from minus strand
			pointsALMinus = scoreWindows(plusMinusC[1],plusMinusT[1],true);
			if (pointsALPlus == null){
				if (verbose) System.out.println("\t\tNo peaks found, skipping!");
				continue;
			}
			
			//extract peaks from windows
			numPeaks = extractPeaks(chrom, true);
			if (verbose) System.out.println("\t\t"+numPeaks+"\t c vs t peaks passing thresholds ");
			}
		}
		
		
		//convert and sort strandedPeaksAL by score, highest to lowest
		strandedPeaks = new StrandedPeak[strandedPeaksAL.size()];
		strandedPeaksAL.toArray(strandedPeaks);
		Arrays.sort(strandedPeaks);
	}

	private boolean makeWindows(int[] bpPositions){
		windows = windowMaker.makeWindows(bpPositions);
		if (windows == null || windows.length == 0){
			return false;
		}
		//assign bp positions
		for (int i=0; i< windows.length; i++){
			windows[i][0] = bpPositions[windows[i][0]];
			windows[i][1] = bpPositions[windows[i][1]]+1;	
		}
		return true;
	}

	private ArrayList<Point> scoreWindows(PointData t, PointData c, boolean applyThresholds){
		//make array list to hold strandedPeaksAL
		ArrayList<Point> pointsAL = new ArrayList<Point>(windows.length/2);

		//for each window, sum the scores, not the positions!
		for (int i=0; i< windows.length; i++){
			//fetch scores
			float tSum = t.sumScoreBP(windows[i][0], windows[i][1]);
			float cSum = c.sumScoreBP(windows[i][0], windows[i][1]);
			//calculate norm diff score
			float score = (tSum- cSum) / new Double(Math.pow(tSum+ cSum, 0.5)).floatValue();
			//save score ?
			if (applyThresholds){
				if (score > minDiffScore && ((cSum * ratioT2CReads) < tSum)) pointsAL.add(new Point(i,score));
			}
			else pointsAL.add(new Point(i,score));
		}
		return pointsAL;
	}
	
	/**Finds 1st highest point. Returns index of point.*/
	private int findHighestPlusStrandPoint(){
		int num = pointsALPlus.size();
		float highestScore = pointsALPlus.get(0).getScore();
		int highestIndex = 0;
		for (int i=1; i< num; i++){
			float testScore = pointsALPlus.get(i).getScore();
			if (testScore > highestScore){
				highestScore = testScore;
				highestIndex = i;
			}
		}
		return highestIndex;
	}
	
	private double[] findHighestNegativeStrandPoint(int centerBP){
		//find 1st highest scoring Point/window just downstream of centerBP on negative strand
		int numPoints = pointsALMinus.size();
		float highestScore = 0;
		int firstBestPosition = 0;
		int firstBestIndex = 0;
		boolean inRange = false;
		for (int i=0; i< numPoints; i++){
			Point testPoint = pointsALMinus.get(i);
			//get associated window (start stop in bps)
			int[] ss = windows[testPoint.getPosition()];
			int testCenter = Num.calculateMiddleIntergenicCoordinates(ss[0], ss[1]);
			int distance = testCenter - centerBP;
			//within range?
			if (distance >=0 && distance <= threePrimePeakSpan) {
				//check score, note <= to maximize distance
				if (highestScore < testPoint.getScore()) {
					highestScore = testPoint.getScore();
					firstBestPosition = testCenter;
					firstBestIndex = i;
				}
				inRange = true;
			}
			else if (inRange) break;
		}	
		
		//now look for additional peaks
		int secondBestIndex = -1;
		int secondBestPosition = 0;
		for (int i=firstBestIndex+1; i< numPoints; i++){
			Point testPoint = pointsALMinus.get(i);
			//get associated window (start stop in bps)
			int[] ss = windows[testPoint.getPosition()];
			int testCenter = Num.calculateMiddleIntergenicCoordinates(ss[0], ss[1]);
			int distance = testCenter - centerBP;
			//within range?
			if (distance <= threePrimePeakSpan) {
				//check score, note <= to maximize distance
				if (highestScore == testPoint.getScore()) {
					secondBestIndex = i;
					secondBestPosition = testCenter;
				}
			}
			else break;
		}
		//calc position;
		if (secondBestIndex != -1){
			firstBestPosition = ((int)Math.round(((double)secondBestPosition-(double)firstBestPosition)/2)) + firstBestPosition;
		}
		return new double[]{highestScore, firstBestPosition};
	}
	
	private int extractPeaks(String chrom, boolean inverted){
		//take top window, make StrandedPeak, mask intersecting windows, repeat to maxNumPeaksPerChrom
		int numPeaks = 0;		
		while (numberPeaksToMerge != numPeaks && pointsALPlus.size() > 0){
			//pull best window from plus strand
			int bestIndex = findHighestPlusStrandPoint();
			Point startPt = pointsALPlus.remove(bestIndex);
			float score = startPt.getScore();

			//look for additional with same score?
			Point stopPt = null;
			while (true){
				//at stop?
				if (bestIndex >= pointsALPlus.size()) break;
				Point pt = pointsALPlus.get(bestIndex);
				float testScore = pt.getScore();
				if (score == testScore) {					
					stopPt = pt;
					pointsALPlus.remove(bestIndex);
				}
				else break;
			}
			//find center position
			int[] ss;
			if (stopPt == null) ss = windows[startPt.getPosition()];
			else {
				int[] start = windows[startPt.getPosition()];
				int[] stop = windows[stopPt.getPosition()];
				ss = new int[]{start[0], stop[1]};
			}
			int centerBP = Num.calculateMiddleIntergenicCoordinates(ss[0], ss[1]); 


			//find highest scoring peak on negative strand within threePrimePeakSpan down stream
			double[] negScorePosition = findHighestNegativeStrandPoint(centerBP);
			
			//is neg strand score within .4/.6 of pos strand?
			double fraction = negScorePosition[0]/(negScorePosition[0]+score);
			if (fraction < 0.4f || fraction > 0.6f)  continue;
			
			int distance = (int)negScorePosition[1]- centerBP;
			
			//make slices
			int start = centerBP - fivePrimePeakSpan;
			if (start < 0) start = 0;
			int negStart = -1 * start;
			int stop = centerBP + threePrimePeakSpan;

			//slice and shift	
			Point[] plusTPoint = fetchSumAndShift(plusMinusT[0], start, stop, negStart);
			Point[] minusTPoint = fetchSumAndShift(plusMinusT[1], start, stop, negStart);
			Point[] plusCPoint = fetchSumAndShift(plusMinusC[0], start, stop, negStart);
			Point[] minusCPoint = fetchSumAndShift(plusMinusC[1], start, stop, negStart);
			
			//save StrandedPeak
			ScoredCoordinate sc = new ScoredCoordinate(chrom, ss[0], ss[1], startPt.getScore());
			if (inverted) strandedPeaksAL.add(new StrandedPeak(sc,plusCPoint, minusCPoint, plusTPoint, minusTPoint, distance));
			else strandedPeaksAL.add(new StrandedPeak(sc,plusTPoint, minusTPoint, plusCPoint, minusCPoint, distance));
			numPeaks++;
			//mask
			start = centerBP - threePrimePeakSpan;
			if (start < 0) start =0;
			for (int i=0; i< pointsALPlus.size(); i++){
				Point test = pointsALPlus.get(i);
				int[] win = windows[test.getPosition()];
				//does the win intersect the peak?			
				if (win[0]>= stop || win[1]<= start) {
					continue;
				}
				//intersects so remove
				pointsALPlus.remove(i);
				i--;
			}

			//stop?
			if (pointsALPlus.size()==0) break;
		}
		return numPeaks;
	}

	public void printPeaks(){
		for (int i=0; i< strandedPeaks.length; i++){
			System.out.println(strandedPeaks[i]);
		}

	}

	public static Point[] fetchSumAndShift(PointData p, int startBP, int stopBP, int shift){
		Point[] pts = p.fetchPoints(startBP, stopBP);
		if (pts == null) return null;
		//strip scores
		for (int i=0; i< pts.length; i++) pts[i].setScore(1f);
		pts = Point.sumIdenticalPositionScores (pts);
		for (int i=0; i< pts.length; i++) pts[i].setPosition(pts[i].getPosition()+shift);
		return pts;
	}

	private class StrandedPeak implements Comparable{

		ScoredCoordinate scoredCoordinate; //of best window
		Point[] tPlus;
		Point[] tMinus;
		Point[] cPlus;
		Point[] cMinus;
		int distance;

		public StrandedPeak(ScoredCoordinate scoredCoordinate, Point[] tPlus, Point[] tMinus, Point[] cPlus, Point[] cMinus, int distance){
			this.scoredCoordinate = scoredCoordinate;
			this.tPlus = tPlus;
			this.tMinus = tMinus;
			this.cPlus = cPlus;
			this.cMinus = cMinus;
			this.distance = distance;
		}

		public int compareTo(Object arg0) {
			double otherScore = ((StrandedPeak)arg0).scoredCoordinate.getScore();
			double thisScore = scoredCoordinate.getScore();
			if (otherScore > thisScore) return 1;
			if (otherScore < thisScore) return -1;
			return 0;
		}

		public String toString(){
			StringBuilder sb = new StringBuilder();
			sb.append("\nScore\t");
			sb.append(scoredCoordinate.toString());
			sb.append("\ntPlus\n");
			sb.append(joinPoints(tPlus));

			sb.append("\ntMinus\n");
			sb.append(joinPoints(tMinus));

			sb.append("\ncPlus\n");
			sb.append(joinPoints(cPlus));

			sb.append("\ncMinus\n");
			sb.append(joinPoints(cMinus));
			
			sb.append("\nDist\t"+distance+"\n");

			return sb.toString();
		}

		public String joinPoints(Point[] p){
			StringBuilder sb = new StringBuilder();
			for (int i=0; i< p.length; i++) {
				sb.append(p[i].toString());
				sb.append("\n");
			}
			return sb.toString();
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new PeakShiftFinder(args);
	}		



	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': treatmentPointDirs = IO.extractFiles(args[++i]); break;
					case 'c': controlPointDirs = IO.extractFiles(args[++i]); break;
					case 's': resultsDirectory = new File(args[++i]); break;
					case 'w': windowSize = Integer.parseInt(args[++i]); break;
					case 'a': minimumNumberWindowReads = Integer.parseInt(args[++i]); break;
					case 'n': numberPeaksToMerge = Integer.parseInt(args[++i]); break;
					case 'p': fivePrimePeakSpan = Integer.parseInt(args[++i]); break;
					case 'm': threePrimePeakSpan = Integer.parseInt(args[++i]); break;
					case 'd': minDiffScore = Float.parseFloat(args[++i]); break;
					case 'r': ratioT2CReads = Double.parseDouble(args[++i]); break;
					case 'e': scanReduced = true; break;
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
		if (treatmentPointDirs == null || treatmentPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your treatment Point Data directories(s)!\n");
		//only one directory look deeper
		if (treatmentPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(treatmentPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) treatmentPointDirs = otherDirs;
		}

		if (controlPointDirs != null){
			if (controlPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your control Point Data directories(s)!\n");
			//only one directory look deeper
			if (controlPointDirs.length == 1){
				File[] otherDirs = IO.extractOnlyDirectories(controlPointDirs[0]);
				if (otherDirs != null && otherDirs.length > 0) controlPointDirs = otherDirs;
			}
		}
		//look for and or create the save directory
		if (resultsDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		if (resultsDirectory.exists() == false) resultsDirectory.mkdir();
		else System.out.println("\nWARNING: the save directory exists, will over write files contained within.\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               PeakShiftFinder: May 2010                          **\n" +
				"**************************************************************************************\n" +
				"PeakShiftFinder estimates the bp difference between sense and antisense proximal chIP-\n" +
				"seq peaks. It calculates the shift int two ways: by generating a composite peak from a\n" +
				"set of the top peaks in a dataset and by taking the median shift for the top peaks.\n" +
				"The latter appears more reliable for some datasets. Inspect the results in IGB by\n" +
				"loading the xxx.bar graphs. When in doubt, run ScanSeqs with just your\n" +
				"treatment data setting the peak shift to 0 and window size to 50 and manually inspect\n" +
				"the shift in IGB.\n\n" +

				"Options:\n"+
				"-t Treatment Point Data directories, full path, comma delimited. These should\n" +
				"       contain stranded chromosome specific xxx_+_.bar.zip and xxx_-_.bar.zip  files.\n"+
				"-c Control Point Data directories, ditto. \n" +
				"-s Save directory, full path.\n\n"+
				
				"Advanced Options:\n"+
				"-e Two chIP samples are provided, no input, scan for reduced peaks too.\n"+
				"-w Window size in bps, defaults to 50.\n"+
				"-a Minimum number window reads, defaults to 10\n" +
				"-d Minimum normalized window score, defaults to 2.5\n"+
				"-r Minimum fold of treatment to control window reads, defaults to 5\n"+
				"-n Number of peaks to merge for composite, defaults to 100\n"+
				"-p Distance off peak center to collect from 5' stop, defaults to 500\n"+
				"-m Distance off peak center to collect from 3' stop, defaults to 1000\n"+
				
				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/PeakShiftFinder -t\n" +
				"      /Data/Ets1Rep1/,/Data/Ets1Rep2/ -c /Data/Input1/,Data/Input2/ -s\n" +
				"      /Results/Ets1PeakShiftResults -w 25 -d 5\n\n" +

		"**************************************************************************************\n");

	}

	public double getCompositePeakShift() {
		return compositePeakShift;
	}

	public double getMedianIndividualPeakShift() {
		return medianIndividualPeakShift;
	}

	public double getStndDevIndividualPeakShift() {
		return stndDevIndividualPeakShift;
	}		

}
