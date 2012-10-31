package trans.main;

import java.io.*;
import java.util.regex.*;


import java.util.*;

import trans.tpmap.WindowMaker;
import util.gen.*;

/**
 * Finds the best sub window and binding peaks (using the {@link PeakPicker}) within each interval.
 */
public class FindSubBindingRegions {
	//fields
	private File[] intervalFiles;
	//sub window
	private int subWindowSize = 350; 
	private int subWindowMinNumOligos = 4;
	//peak picker, peaks
	private int peakPickerSmoothinWindowSize = 230;
	private int maxNumBindingPeaks = 4;
	private double bindingPeakScoreCutOff = 1.3;
	//flat top
	private int maxBpGap = 75;		//max bp gap between oligo starts
	private double maxScoreFractionDrop = 0.80;	//expanded score cannot be less than this * original sub window score
	private boolean makeFlattops = false;

	
	public FindSubBindingRegions(String[] args){
		//process user params
		processArgs(args);
		
		//for each interval file
		int numIntFiles = intervalFiles.length;
		Interval[] intervals;
		int numIntervals;
		Oligo[] oligos;
		int[] starts = null;
		float[][] treatmentIntensities = null;
		float[][] controlIntensities = null;
		int numTreatmentValues;
		int numControlValues;
		int numOligos;
		WindowMaker windowMaker = new WindowMaker(subWindowSize, subWindowMinNumOligos);
		int[][] subWindows = null;
		
		try {
			for (int i=0; i<numIntFiles; i++){
				//attempt to fetch Interval[]
				System.out.println("\tProcessing "+intervalFiles[i]);
				intervals = (Interval[])IO.fetchObject(intervalFiles[i]);
				
				//for each interval get oligo[]
				numIntervals = intervals.length;
				for (int j=0; j<numIntervals; j++){
					oligos = intervals[j].getOligos();
					if (oligos == null){
						System.out.println("\nThis interval file has not been processed by the LoadIntervalOligoInfo program.  Do so and " +
						"rerun this program. \n");
						System.exit(0);
					}
					//create int[] of bp starts and int[][] of treatment intensities and int[][] of control intensities
					numOligos = oligos.length;
					numTreatmentValues = intervals[j].getNumberTreatmentIntensities();
					numControlValues = intervals[j].getNumberControlIntensities();
					starts = new int[numOligos];
					treatmentIntensities = new float[numOligos][];
					controlIntensities = new float[numOligos][];
					for (int k=0; k<numOligos; k++){
						starts[k] = oligos[k].getStart();
						treatmentIntensities[k] = oligos[k].getTreatmentIntensities(numTreatmentValues);
						controlIntensities[k] = oligos[k].getControlIntensities(numControlValues);
					}
					double[] aveT = Num.averageFloatArrays(treatmentIntensities);
					double[] aveC = Num.averageFloatArrays(controlIntensities);
					
					//create smoothed oligo values using a trimmed mean and set in each Oligo
					double[] ratios = Num.ratio(aveT, aveC);
					double[] smoothedScores = Num.smoothScores(ratios, starts, peakPickerSmoothinWindowSize);
					for (int k=0; k<numOligos; k++){
						oligos[k].setSmoothedScore(smoothedScores[k]);
					}
					//create subwindows int[][]
					subWindows = null;
					if (starts.length >1) subWindows = windowMaker.makeWindows(starts);
					//check to see if any windows were found
					if (subWindows != null && subWindows.length!=0){
						//find highes scoring sub window
						SubWindow sw = findHighestScoringWindow(subWindows, aveT, aveC);
						int[] startStopOligoArrayIndex = subWindows[sw.getIndex()];
						//make sub oligos and set in Interval
						int numSubOligos = startStopOligoArrayIndex[1] - startStopOligoArrayIndex[0] + 1;
						Oligo[] subOligos = new Oligo[numSubOligos];
						int counter = 0;
						int numRun = startStopOligoArrayIndex[1]+1;
						for (int k= startStopOligoArrayIndex[0]; k<numRun; k++){
							subOligos[counter++] = oligos[k];
						}
						sw.setOligos(subOligos);
						//save SubWindow in Interval
						intervals[j].setBestSubWindow(sw);
					}
				}
				
				//find and set BindingPeaks
				PeakPicker p = new PeakPicker();
				p.setMaxNumPeaks(maxNumBindingPeaks);
				p.setScoreCutOff(bindingPeakScoreCutOff);
				if (makeFlattops){	
					p.setMakeFlattops(makeFlattops);
					p.setMaxBpGap(maxBpGap);
					p.setMaxScoreFractionDrop(maxScoreFractionDrop);
				}
				p.pickPeaks(intervals);
				
				//save intervals 
				File save = new File (intervalFiles[i].getCanonicalPath()+"Sub");
				System.out.println("\tSaving sub windowed Interval file "+save.getCanonicalPath());
				IO.saveObject(save, intervals);
			}
			
		} catch (Exception e){
			System.out.println("\nSomething is wrong with one of your apparent interval files, skipping");
			e.printStackTrace();
		}
	}	
		
	/**Returns a partially complete SubWindow.*/
	public static SubWindow findHighestScoringWindow(int[][] win, double[] t, double[] c){
		int maxIndex = 0;
		double maxDiff = 0;
		double[] ratios;
		int numWindows = win.length;
		int numRun;
		double diff;
		int counter;
		//for each window, sum t and sum c, calc ratio and compare to max
		for (int i=0; i<numWindows; i++){
			numRun = win[i][1]+1;	
			ratios = new double[numRun-win[i][0]];
			counter =0;
			//calculate ratios
			for (int j=win[i][0]; j<numRun; j++){
				ratios[counter++] = t[j]/c[j];
			}

			//toss any Nan, Infinity
			ratios = Num.removeNaNInfinite(ratios);

			if (ratios.length>1){
				//median
				Arrays.sort(ratios);
				diff = Num.median(ratios);
			}
			else diff = Num.mean(ratios);
			//bigger than max?
			if (diff> maxDiff){
				maxDiff = diff;						
				maxIndex = i;
			}
		}
		//int index, double medianRatio
		return new SubWindow(maxIndex, maxDiff);
	}
	
	
	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}
		new FindSubBindingRegions(args);
	}
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		File directory = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'i': directory = new File (args[i+1]); i++; break;
					case 'w': subWindowSize = Integer.parseInt(args[i+1]); i++; break;
					case 's': peakPickerSmoothinWindowSize = Integer.parseInt(args[i+1]); i++; break;
					case 'n': subWindowMinNumOligos = Integer.parseInt(args[i+1]); i++; break;
					case 'm': maxNumBindingPeaks = Integer.parseInt(args[i+1]); i++; break;
					case 'c': bindingPeakScoreCutOff = Double.parseDouble(args[i+1]); i++; break;
					case 'f': makeFlattops = true; break;
					case 'd': maxScoreFractionDrop = Double.parseDouble(args[i+1]); i++; break;
					case 'g': maxBpGap = Integer.parseInt(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter: -"+test);
					System.out.println();
					System.exit(0);
				}
			}		
		}
		//check directories and files
		if (directory == null){
			System.out.println("\nPlease enter a directory containing serialized Interval[] files or the text of one file.\n");
			System.exit(0);
		}
		intervalFiles = IO.extractFiles(directory);
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Find Sub Binding Regions:    Oct 2006                      **\n" +
				"**************************************************************************************\n" +
				"FSBR takes Interval[] that have been loaded with oligo information and scans each for\n" +
				"the highest median ratio scoring sub window. FSBR also picks binding peaks using\n" +
				"smoothed intensity values (trimmed mean ratio sliding window) .\n\n" +
				
				"Parameters:\n" +
				"-i The full path file text for a serialize Interval[] file, alternatively, give a\n" +
				"       directory, and all the Interval files within will be processed.\n" +
				"-w Sub Window size, default is 350bp.\n" +
				"-n Minimum number of oligos required in sub window, defaults to 4.\n" +
				"-s Peak Picker smoothing window size, default is 230bp.\n" +
				"-m Max number of peaks to find, defaults to 4.\n"+
				"-c Minimum binding peak ratio score, defaults to 1.3, set internally for flattops.\n"+
				"-f Pick flattop peaks, defaults to sloped peaks.\n"+
				"-g Max bp gap for expansion of flattop peaks, defaults to 75bp.\n"+
				"-d Maximum fraction score drop for expanding flattop peaks, defaults to 0.8.\n\n"+
				
				"Example: java -Xmx256M -jar pathTo/T2/Apps/FindSubBindingRegions -i /affy/Ints/ -w\n" +
				"       300 -n 4\n\n" +
				
		"**************************************************************************************\n");		
	}	
}
