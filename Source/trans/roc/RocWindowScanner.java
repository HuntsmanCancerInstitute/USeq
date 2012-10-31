package trans.roc;
import java.util.regex.*;
import trans.anno.*;
import trans.main.Interval;
import trans.main.Window;
import trans.misc.*;
import java.util.*;

import java.io.*;
import util.gen.*;


/**
 *Scans a window array over a range of thresholds and intersects with a list of positives to 
 *generate the false positive rate, true positive rate, and false discovery rates.
 */
public class RocWindowScanner {
	//fields
	private String file;
	private int scoreIndex = 0;
	private Window[] windows;
	private double[] maxMinWindowScores;
	private double numberOfWindows;
	private double totalIntersectingWindows = 0;
	private double totalNonIntersectingWindows = 0;
	private int minNumberOligos = 10;
	private int sizeOfOligoMinusOne= 24;

	/**When a window partially overlaps a region, the fraction of the region that must be covered
	 * by a window to return true when intersecting the two.  Returns true if contained within or by.
	 * 0.5 = 50% converage of the region by the window.*/
	private double fractionAcceptibleCoverage = 0.75;
	private File positiveRegionsFile;
	private GenomicRegion[] posRegions;
	private File maskedRegionsFile;
	private GenomicRegion[] maskedRegions;
	private double targetMaxFDR = 0.25;
	private double targetFDR = 0.05;
	private double targetMinFDR = 0.008;
	private StringBuffer targetLines = new StringBuffer(); 
	
	//contructor
	public RocWindowScanner(String[] args){
		//process args
		processArgs(args);
		
		System.out.println("At a given threshold:\n\tFDR = #FP/(#FP+#TP)\n\tPVal = #FP/Total # Windows\n");
		
		//print params
		System.out.println("\nFraction Acceptible Coverage: "+fractionAcceptibleCoverage);
		System.out.println("\nMin Number Probes: "+minNumberOligos);
		
		//get files to process
		File[] files = IO.extractFiles(new File(file));
		
		//load positive regions
		posRegions = GenomicRegion.parseRegions(positiveRegionsFile);
		System.out.println("Positive regions parsed: "+posRegions.length);
		
		//load regions to skip (ie those repeat masked)
		if (maskedRegionsFile != null) {
			maskedRegions = GenomicRegion.parseRegions(maskedRegionsFile);
			System.out.println("Masked regions parsed: "+maskedRegions.length);
		}
		
		
		//for each Window[] file
		for (int i=0; i<files.length; i++){
			//fetch windows and strip of windows that don't contain the min number of oligos.
			windows = (Window[])IO.fetchObject(files[i]);
			System.out.println("\nProcessing "+files[i]+"...");
			windows = Util.removeLowOligoWindows(windows, minNumberOligos);
			//strip masked windows?
			if (maskedRegionsFile != null) removeMaskedRegions();
			//count window intersection stats
			calculateTotalWindowIntersectionStats();

			//for each score index
			for (int z=0; z<windows[0].getScores().length; z++){	
			//for (int z=15; z<17; z++){
				scoreIndex = z;
				System.out.println("Threshold\t#Windows Passing Threshold\t#Pos\t#Neg\tFPR\tTPR\tFDR\tPVal\tScore index="+scoreIndex);
				
				//set initial start stop thresholds to min max
				maxMinWindowScores = Util.minMaxWindowArray(windows, scoreIndex);
				double startThreshold = maxMinWindowScores[0];
				double stopThreshold = maxMinWindowScores[1];
				
				//find start and stop thresholds for ROC curve that bracket the maxFDR and minFDR
				double[] startStop = findStartStopThresholds(startThreshold, stopThreshold);
				startThreshold = startStop[0];
				stopThreshold = startStop[1];
				
				//check for problems
				if (startThreshold == stopThreshold){
					targetLines.append(files[i].getName());
					targetLines.append("\t");
					targetLines.append(scoreIndex);
					targetLines.append("\t");
					targetLines.append("Could not set initial thresholds?!  Score bound Min: "+maxMinWindowScores[0]+"  Max: "+maxMinWindowScores[1]+"\n");
					continue;
				}

				//reset thresholds for ROC curve
				double increment = (stopThreshold- startThreshold) / 500;
				double threshold = startThreshold;	

				//print ROC lines, for each threshold scan windows and calculate fpr fdr ..., also find thresholds that bracket FDR
				boolean fDRStopFound = false;
				boolean fDRStartFound = false;
				double startFDR = 0;
				double stopFDR = 0;
				while (threshold <= stopThreshold){
					//intersect
					double numPositivesHit = 0;
					double numWindowsPassingThreshold = 0;
					//for each window
					for (int k=0; k<numberOfWindows; k++){
						boolean positiveHitCounted = false;
						//check threshold
						if (windows[k].getScores()[scoreIndex] >= threshold){ 
							numWindowsPassingThreshold++;
							//for each region
							for (int j=0; j<posRegions.length; j++){
								if (positiveHitCounted == false && posRegions[j].intersect(windows[k], sizeOfOligoMinusOne, fractionAcceptibleCoverage)) {
									numPositivesHit++;
									positiveHitCounted = true;
								}
							}
						}
					}
					double numNegWindows = numWindowsPassingThreshold - numPositivesHit;
					double tpr = numPositivesHit/totalIntersectingWindows;
					double fpr = numNegWindows/totalNonIntersectingWindows;
					double fdr = numNegWindows/(numNegWindows+numPositivesHit);
					double pVal = numNegWindows/numberOfWindows;
					//print threshold, #pass, #pos, #neg, fpr, tpr
					StringBuffer sb = new StringBuffer();
					sb.append(threshold);
					sb.append("\t");
					sb.append(numWindowsPassingThreshold);
					sb.append("\t");
					sb.append(numPositivesHit);
					sb.append("\t");
					sb.append(numNegWindows);
					sb.append("\t");
					sb.append(fpr);
					sb.append("\t");
					sb.append(tpr);
					sb.append("\t");
					sb.append(fdr);
					sb.append("\t");
					sb.append(pVal);
					System.out.println(sb);
					
					if (fDRStopFound == false && fdr > targetFDR){
						startFDR = threshold;					
					}
					if (fDRStartFound == false && fdr < targetFDR){
						fDRStopFound = true;
						fDRStartFound = true;
						stopFDR = threshold;
					}
					if (fdr <= 0.005) break;
					//increment threshold
					threshold += increment;
				}
				
				//check if fdr found
				if (fDRStartFound == false || fDRStopFound == false){
					targetLines.append(files[i].getName());
					targetLines.append("\t");
					targetLines.append(scoreIndex);
					targetLines.append("\t");
					targetLines.append("Could not set FDR start and stop.\n");
					//Misc.printExit("Could not set FDR start and stop.\n");	
				}
				
				
				//reset thresholds for FDR determination
				stopThreshold = stopFDR;
				startThreshold = startFDR;
				increment = (stopThreshold- startThreshold) / 100;
				threshold = startThreshold;	
				
				boolean targetFDRFound = false;
				StringBuffer sb = null;
				fDRStopFound = false;
				fDRStartFound = false;
				double beginningFDR = 0;
				double endingFDR = 0;
				double beginningTPR = 0;
				double endingTPR = 0;

				//for each threshold scan windows and count whether they intersect
				while (threshold <= stopThreshold){
					//intersect
					double numPositivesHit = 0;
					double numWindowsPassingThreshold = 0;
					//for each window
					for (int k=0; k<numberOfWindows; k++){
						boolean positiveHitCounted = false;
						//check threshold
						if (windows[k].getScores()[scoreIndex] >= threshold){ 
							numWindowsPassingThreshold++;
							//for each region
							for (int j=0; j<posRegions.length; j++){
								if (positiveHitCounted == false && posRegions[j].intersect(windows[k], sizeOfOligoMinusOne, fractionAcceptibleCoverage)) {
									numPositivesHit++;
									positiveHitCounted = true;
								}
							}
						}
					}
					double numNegWindows = numWindowsPassingThreshold - numPositivesHit;
					double fdr = numNegWindows/(numNegWindows+numPositivesHit);
					
					if (fDRStopFound == false && fdr > targetFDR){
						endingFDR = fdr;
						endingTPR = numPositivesHit/totalIntersectingWindows;
					}
					if (fDRStartFound == false && fdr <= targetFDR){
						fDRStopFound = true;
						fDRStartFound = true;
						beginningFDR = fdr;
						beginningTPR = numPositivesHit/totalIntersectingWindows;
					}
					//look for target and save if found
					if (targetFDRFound == false && fdr <= targetFDR){
						targetFDRFound = true;
						double tpr = numPositivesHit/totalIntersectingWindows;
						double fpr = numNegWindows/totalNonIntersectingWindows;
						double pVal = numNegWindows/numberOfWindows;
						//print threshold, #pass, #pos, #neg, fpr, tpr
						sb = new StringBuffer();
						sb.append(threshold);
						sb.append("\t");
						sb.append(numWindowsPassingThreshold);
						sb.append("\t");
						sb.append(numPositivesHit);
						sb.append("\t");
						sb.append(numNegWindows);
						sb.append("\t");
						sb.append(fpr);
						sb.append("\t");
						sb.append(tpr);
						sb.append("\t");
						sb.append(fdr);	
						sb.append("\t");
						sb.append(pVal);	
						sb.append("\t");
						//calc regressed TPR
						double rTPR = LinearRegression.interpolateY(beginningFDR, beginningTPR, endingFDR, endingTPR, targetFDR);
						sb.append(rTPR);
						//prepend file and index info to summary targetLines
						targetLines.append(files[i].getName());
						targetLines.append("\t");
						targetLines.append(scoreIndex);
						targetLines.append("\t");
						targetLines.append(sb.toString());
						targetLines.append("\n");
						break;
					}
					//increment threshold
					threshold += increment;
				}
				//check to see if target FDR was reached
				if (targetFDRFound == false){
					targetFDRFound = true;
					targetLines.append(files[i].getName());
					targetLines.append("\t");
					targetLines.append(scoreIndex);
					targetLines.append("\t");
					targetLines.append("Failed to reach FDR target!!!\n");
				}
			}
		}
		//print target lines
		System.out.println("\nName\tScoreIndex\tThreshold\t#Windows Pass Threshold\t#Pos\t#Neg\tFPR\tTPR\tFDR\tPVal\tInterpolated TPR at "+targetFDR+" FDR");
		System.out.println(targetLines);
	}
	
	public static void main(String[] args) {
		if (args.length<1){
			printDocs();
			System.exit(0);	
		}
		new RocWindowScanner(args);
	}
	
	/**Calculate total # windows that intersect pos regions, no counting, ie one window overlaps two pos regions is counted once.
	 * Also strip out windows that partially overlap a positive region but don't meet the fractionAcceptibleCoverage.*/
	public void calculateTotalWindowIntersectionStats(){
		numberOfWindows = windows.length;
		totalIntersectingWindows = 0;
		totalNonIntersectingWindows = 0;
		ArrayList al = new ArrayList((int)numberOfWindows);
		//for each window
		for (int k=0; k<numberOfWindows; k++){
			//for each region
			boolean found = false;
			boolean partialIntersection = false;
			for (int j=0; j<posRegions.length; j++){
				//does it partially intersect
				if (posRegions[j].intersect(windows[k], sizeOfOligoMinusOne, 0.001)){
					//does it intersect a positive region?
					if (posRegions[j].intersect(windows[k], sizeOfOligoMinusOne, fractionAcceptibleCoverage)) {
						if (found == false) totalIntersectingWindows++;
						found = true;
					}
					else partialIntersection = true;
				}
			}
			if (partialIntersection == false){
				found = true;
				totalNonIntersectingWindows++;
			}
			if (found) al.add(windows[k]);
		}
		numberOfWindows = al.size();
		windows = new Window[(int)numberOfWindows]; 
		al.toArray(windows);			
		System.out.println("Total intesecting windows: "+(int)totalIntersectingWindows+" Total non-intersecting windows: "+(int)totalNonIntersectingWindows);
	}
	
	/**Remove windows that intersect the masked regions.*/
	public void removeMaskedRegions(){
		ArrayList al = new ArrayList(windows.length);
		//for each window
		int masked = 0;
		for (int k=0; k<windows.length; k++){
			//for each masked region
			boolean noIntersection = true;
			for (int j=0; j<maskedRegions.length; j++){
				//does it partially intersect
				if (maskedRegions[j].intersect(windows[k], sizeOfOligoMinusOne, 0.00001)){
					masked++;
					noIntersection = false;
					break;
				}
			}
			if (noIntersection) al.add(windows[k]);
		}
		numberOfWindows = al.size();
		windows = new Window[al.size()]; 
		al.toArray(windows);			
		System.out.println("Total masked windows: "+masked);
	}
	
	/**Scans an array of Window intersecting with positives to calculate an FDR at a given threshold.*/
	public double calcFDR(double threshold){
		//intersect
		double numPositivesHit = 0;
		double numWindowsPassingThreshold = 0;
		//for each window
		for (int k=0; k<numberOfWindows; k++){
			boolean positiveHitCounted = false;
			//check threshold
			if (windows[k].getScores()[scoreIndex] >= threshold){ 
				numWindowsPassingThreshold++;
				//for each region
				for (int j=0; j<posRegions.length; j++){
					if (positiveHitCounted == false && posRegions[j].intersect(windows[k], sizeOfOligoMinusOne, fractionAcceptibleCoverage)) {
						numPositivesHit++;
						positiveHitCounted = true;
					}
				}
			}
		}
		double numNegWindows = numWindowsPassingThreshold - numPositivesHit;
		double fdr = numNegWindows/(numNegWindows+numPositivesHit);
		return fdr;
	}
	
	/**Attempts to find a threshold that gives the target FDR, will run through maxCycles and return closest threshold if targetFDR is not reached.*/
	public double findTargetFDRThreshold(double fdrTarget, int maxCycles){
		//use float to give some wiggle room for ==
		float targetFDR = new Double(fdrTarget).floatValue();
		//set starting max min
		double min = maxMinWindowScores[0];
		double max = maxMinWindowScores[1];
		double threshold = 0;
		for (int x=0; x< maxCycles; x++){
			//set starting threshold
			threshold = ((max-min)/2.0) + min;
			//calc FDR
			float testFdr = new Double (calcFDR(threshold)).floatValue();		
			//is it equal to the target fdr
			if (testFdr == targetFDR) {
				return threshold;
			}
			//is it higher 
			if (testFdr > targetFDR) {
				min = threshold;
			}
			//must be lower
			else {
				max = threshold;
			}
		}
		return threshold;
	}
	
	public double[] findStartStopThresholds(double startThreshold, double stopThreshold){
		return new double[] {			
				findTargetFDRThreshold(targetMaxFDR, 50),
				findTargetFDRThreshold(targetMinFDR, 50)
		};
	}
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'r': minNumberOligos =Integer.parseInt(args[i+1]);i++; break;
					case 'z': sizeOfOligoMinusOne =Integer.parseInt(args[i+1])-1;i++; break;
					case 'c': fractionAcceptibleCoverage =Double.parseDouble(args[i+1]);i++; break;
					case 'f': file = args[i+1]; i++; break;
					case 'p': positiveRegionsFile = new File(args[i+1]); i++; break;
					case 'm': maskedRegionsFile = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (Misc.isEmpty(file)){
			Misc.printExit("\nPlease enter an TiMAT results file or directory.\n");
		}
	}	
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Roc Window Scanner: June 2006                           **\n" +
				"**************************************************************************************\n" +
				"RWS scans windows over a range of cutoff scores and intersects\n" +
				"each with a positive regions file generating TPR and FPR for ROC curves.\n\n" +
				
				"Options:\n\n" +
				
				"-c  Minimal fraction of coverage for a region when intersecting with a window, defaults\n" +
				"       to 0.75. Regions containing or contained by a window always return true.\n" +
				"-r  Required number of oligos in each window, defaults to 10.\n"+
				"-z Size of oligo, defaults to 25.\n"+
				"-f  Full path file text for the Window[] array, if a directory is specified,\n" +
				"       all files will be processed.\n\n" +
				"-p  Full path file text for the positive regions text file, tab delimited (chrom,\n"+
				"       start, stop, etc.)\n"+
				"-m  Full path file text for a list of regions to exclude, formatted as above.\n\n"+
				
				"Example: java -Xmx500M trans/main/RocWindowScanner -f /affy/res/zeste -p /affy/pos.txt\n\n"+
				
		"**************************************************************************************\n");		
	}
	
}
