package trans.main;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import trans.tpmap.MapSplitter;
import util.gen.*;

/**
 * Scores windows of intensity values using several tests.
 * If the user processed their cel intensities using the MM data then remember the following transform was 
 * performed: intensity = max((pm-mm),1)
 * Rather inefficient.
 */
public class TestScanChip {
	
	//fields
	private String tpmapFile; 
	private String infoFile;
	private String resultsFile; 
	private String treatmentDir; 
	private String controlDir; 
	private boolean writeSgrs = false;
	private boolean multiplyByNumberValues = false;  //multiplies by square root of the number of t and c's
	private static final double log2 = Math.log(2);
	private boolean pairs = false;
	private int sizeOfOligoMinusOne = 24;
	
	
	public TestScanChip(String[] args){
		processArgs(args);
		
		if (multiplyByNumberValues) System.out.println("\nMultiplying scores by the square root of the window size!\n");
		if (pairs) System.out.println("\nMaking paired scores, not layered.\n");
		
		System.out.println ("No filtering, processing all windows...");
		
		//get .tpmap file info
		ArrayList info = (ArrayList)IO.fetchObject(new File(infoFile));
		//fetch float[][] of processed intensity values for treatment and controls
		//break both apart by chromosome using the bpmapInfo ArrayList, save to disk
		float[][] treatment = ScanChip.fetchFloatArrays(treatmentDir,false);
		System.out.println("\nProcessing "+treatment.length+" treatment arrays...");
		String uniqueIdT = MapSplitter.breakSaveIntensityValues(info, treatment,treatmentDir);
		treatment = null;
		Misc.printArray(IO.extractFiles(treatmentDir,".celp"));
		
		float[][] control = ScanChip.fetchFloatArrays(controlDir, false);
		System.out.println("Processing "+control.length+" control arrays...");
		String uniqueIdC = MapSplitter.breakSaveIntensityValues(info, control,controlDir);
		System.out.println();
		control = null;
		Misc.printArray(IO.extractFiles(controlDir,".celp"));
		
		//params
		int num = info.size();
		String chromosome;
		int[][] testWindows;
		int numWindows;
		int[] startStop;
		int[] positions;
		int windowSize;
		ArrayList windowsAl = new ArrayList(100000);
		PrintWriter wilcoxonRankSumLogPValueFile  = null;
		PrintWriter psmMeanRatioFile  = null;
		PrintWriter psmMeanLogRatioFile  = null;
		PrintWriter psmMeanRelDiffFile  = null;
		PrintWriter psmLayeredRatioFile  = null;
		PrintWriter psmLayeredLogRatioFile  = null;
		PrintWriter psmLayeredRelDiffFile  = null;
		PrintWriter trmMeanRatioFile  = null;
		PrintWriter trmMeanLogRatioFile  = null;
		PrintWriter trmMeanRelDiffFile  = null;
		PrintWriter trmLayeredRatioFile  = null;
		PrintWriter trmLayeredLogRatioFile  = null;
		PrintWriter trmLayeredRelDiffFile  = null;
		PrintWriter tkbMeanRatioFile  = null;
		PrintWriter tkbMeanLogRatioFile  = null;
		PrintWriter tkbMeanRelDiffFile  = null;
		PrintWriter tkbLayeredRatioFile  = null;
		PrintWriter tkbLayeredLogRatioFile  = null;
		PrintWriter tkbLayeredRelDiffFile  = null;
		PrintWriter gmMeanRatioFile  = null;
		PrintWriter gmLayeredRatioFile  = null;
		PrintWriter mdMeanRatioFile  = null;
		PrintWriter mdMeanLogRatioFile  = null;
		PrintWriter mdMeanRelDiffFile  = null;
		PrintWriter mdLayeredRatioFile  = null;
		PrintWriter mdLayeredLogRatioFile  = null;
		PrintWriter mdLayeredRelDiffFile  = null;
		PrintWriter mMeanRatioFile  = null;
		PrintWriter mMeanLogRatioFile  = null;
		PrintWriter mMeanRelDiffFile  = null;
		PrintWriter mLayeredRatioFile  = null;
		PrintWriter mLayeredLogRatioFile  = null;
		PrintWriter mLayeredRelDiffFile  = null;
		PrintWriter matFile  = null;
		
		
		try{
			//print writers
			if (writeSgrs){
				wilcoxonRankSumLogPValueFile  = new PrintWriter(new FileWriter(resultsFile+"wcxnRkSumPValue.sgr"));
				psmMeanRatioFile  = new PrintWriter(new FileWriter(resultsFile+"psmMnRto.sgr"));
				psmMeanLogRatioFile  = new PrintWriter(new FileWriter(resultsFile+"psmMnLgRto.sgr"));
				psmMeanRelDiffFile  = new PrintWriter(new FileWriter(resultsFile+"psmMnRD.sgr"));
				psmLayeredRatioFile  = new PrintWriter(new FileWriter(resultsFile+"psmLydRto.sgr"));
				psmLayeredLogRatioFile  = new PrintWriter(new FileWriter(resultsFile+"psmLydLgRto.sgr"));
				psmLayeredRelDiffFile  = new PrintWriter(new FileWriter(resultsFile+"psmLydRD.sgr"));
				trmMeanRatioFile  = new PrintWriter(new FileWriter(resultsFile+"trmMnRto.sgr"));
				trmMeanLogRatioFile  = new PrintWriter(new FileWriter(resultsFile+"trmMnLgRto.sgr"));
				trmMeanRelDiffFile  = new PrintWriter(new FileWriter(resultsFile+"trmMnRD.sgr"));
				trmLayeredRatioFile  = new PrintWriter(new FileWriter(resultsFile+"trmLydRto.sgr"));
				trmLayeredLogRatioFile  = new PrintWriter(new FileWriter(resultsFile+"trmLydLgRto.sgr"));
				trmLayeredRelDiffFile  = new PrintWriter(new FileWriter(resultsFile+"trmLydRD.sgr"));
				tkbMeanRatioFile  = new PrintWriter(new FileWriter(resultsFile+"tkbMnRto.sgr"));
				tkbMeanLogRatioFile  = new PrintWriter(new FileWriter(resultsFile+"tkbMnLgRto.sgr"));
				tkbMeanRelDiffFile  = new PrintWriter(new FileWriter(resultsFile+"tkbMnRD.sgr"));
				tkbLayeredRatioFile  = new PrintWriter(new FileWriter(resultsFile+"tkbLydRto.sgr"));
				tkbLayeredLogRatioFile  = new PrintWriter(new FileWriter(resultsFile+"tkbLydLgRto.sgr"));
				tkbLayeredRelDiffFile  = new PrintWriter(new FileWriter(resultsFile+"tkbLydRD.sgr"));
				gmMeanRatioFile  = new PrintWriter(new FileWriter(resultsFile+"gmMnRto.sgr"));
				gmLayeredRatioFile  = new PrintWriter(new FileWriter(resultsFile+"gmLydRto.sgr"));
				mdMeanRatioFile  = new PrintWriter(new FileWriter(resultsFile+"mdMnRto.sgr"));
				mdMeanLogRatioFile  = new PrintWriter(new FileWriter(resultsFile+"mdMnLgRto.sgr"));
				mdMeanRelDiffFile  = new PrintWriter(new FileWriter(resultsFile+"mdMnRD.sgr"));
				mdLayeredRatioFile  = new PrintWriter(new FileWriter(resultsFile+"mdLydRto.sgr"));
				mdLayeredLogRatioFile  = new PrintWriter(new FileWriter(resultsFile+"mdLydLgRto.sgr"));
				mdLayeredRelDiffFile  = new PrintWriter(new FileWriter(resultsFile+"mdLydRD.sgr"));
				mMeanRatioFile  = new PrintWriter(new FileWriter(resultsFile+"mMnRto.sgr"));
				mMeanLogRatioFile  = new PrintWriter(new FileWriter(resultsFile+"mMnLgRto.sgr"));
				mMeanRelDiffFile  = new PrintWriter(new FileWriter(resultsFile+"mMnRD.sgr"));
				mLayeredRatioFile  = new PrintWriter(new FileWriter(resultsFile+"mLydRto.sgr"));
				mLayeredLogRatioFile  = new PrintWriter(new FileWriter(resultsFile+"mLydLgRto.sgr"));
				mLayeredRelDiffFile  = new PrintWriter(new FileWriter(resultsFile+"mLydRD.sgr"));
				matFile  = new PrintWriter(new FileWriter(resultsFile+"mat.sgr"));
			}
			
			//for each chromosome	
			for (int i=1; i<num; i+=4){
				//initialize tester
				WilcoxonRankSumTest wt = new WilcoxonRankSumTest();
				//WilcoxonSignedRankTest wsrt = new WilcoxonSignedRankTest();
				
				//fetch testWindows to scan, the int[][] made in the WindowMaker
				chromosome = (String)info.get(i);
				System.out.println("Testing chromosome: "+chromosome);
				testWindows = (int[][])IO.fetchObject(new File(tpmapFile+chromosome+"Win"));
				
				//fetch treatment and control normalized, transformed, and chromosome divided int[]'s
				treatment = (float[][])IO.fetchObject(new File(treatmentDir+File.separator+chromosome+uniqueIdT));
				control = (float[][])IO.fetchObject(new File(controlDir+File.separator+chromosome+uniqueIdC));
				double numTreatments = treatment.length;
				double numControls = control.length;
				
				//fetch bp positions to use in converting window indexes to real base pairs
				positions = (int[])IO.fetchObject(new File(tpmapFile+chromosome));
				
				//run the window tests
				numWindows = testWindows.length;
				for (int j=0; j<numWindows; j++){
					
					//get transformed values to test
					startStop = testWindows[j];
					windowSize = 1+startStop[1]-startStop[0];
					
					//build array of sub treatment and sub control windows int[windowNumber][oligo values]
					float[][] subT = new float[(int)numTreatments][windowSize];
					float[][] subC = new float[(int)numControls][windowSize];
					for (int k=0; k<numTreatments; k++){
						System.arraycopy(treatment[k],startStop[0], subT[k], 0, windowSize);
					}
					for (int k=0; k<numControls; k++){
						System.arraycopy(control[k],startStop[0], subC[k], 0, windowSize);
					}	
					//average arrays, if mm data used then max((pm-mm),1) was performed!
					double[] aveT = Num.averageFloatArraysFlipped(subT);
					double[] aveC = Num.averageFloatArraysFlipped(subC);
					int numAves = aveT.length;
					
					//calc square root of the total number of intensities
					double varStbl = 1;
					if (multiplyByNumberValues) {
						varStbl = Math.sqrt(windowSize);
					}
					
					//first check if trimmed mean log ratio > 0 then do other computations if so
					double[] logRatios = new double[numAves];
					for (int x=0; x<numAves; x++){
						logRatios[x] = Math.log(aveT[x]/aveC[x]) /log2;
					}	
					double trmMeanLogRatio = Num.trimmedMean(logRatios,0.1);
					
					
					//make layers or pairs
					double[] layeredRatios;
					double[] layeredLog2Ratios;
					double[] layeredRelDiffs;
					if (pairs){
						layeredRatios = Num.pairedRatios(subT, subC);
						layeredLog2Ratios = Num.pairedLogRatios(subT, subC);
						layeredRelDiffs = Num.pairedRelativeDifferences(subT, subC);
					}
					else {
						layeredRatios = Num.layeredRatios(subT, subC);
						layeredLog2Ratios = Num.layeredLogRatios(subT, subC);
						layeredRelDiffs = Num.layeredRelativeDifferences(subT, subC);
					}
					
					//take geometric mean instead of mean, assumes no zeros or neg numbers.
					double[] geoT = Num.geometricMean(subT);
					double[] geoC = Num.geometricMean(subC);
					
					//for pooled T and C for statistical tests
					float[] allT = Num.collapseFloatArray(subT);
					float[] allC = Num.collapseFloatArray(subC);
					
					//calculate stats
					double[] ratios = new double[numAves];
					double[] relDiffs = new double[numAves];
					double[] geoRatios = new double[numAves];
					double[] geoRelDiffs = new double[numAves];
					double[] geoLogRatios = new double[numAves];
					
					for (int x=0; x<numAves; x++){
						//for means
						ratios[x] = aveT[x]/aveC[x];
						relDiffs[x] = Num.relativeDifference(aveT[x], aveC[x]);
						//for geometric means
						geoRatios[x] = geoT[x]/geoC[x];
						geoLogRatios[x] = Math.log(geoRatios[x]) /log2;
						geoRelDiffs[x] = Num.relativeDifference(geoT[x], geoC[x]);
						
					}	
					
					//sort scores
					Arrays.sort(ratios);
					Arrays.sort(logRatios);
					Arrays.sort(relDiffs);
					Arrays.sort(layeredRatios);
					Arrays.sort(layeredLog2Ratios);
					Arrays.sort(layeredRelDiffs);
					Arrays.sort(geoRatios);
					Arrays.sort(geoLogRatios);
					Arrays.sort(geoRelDiffs);
					
					//stat tests
					double wilcoxonRankSumLogPValue = wt.test(allT, allC);
					//double wilcoxonSignedRankLogPValue = wsrt.test(Num.convertToInt(aveT), Num.convertToInt(aveC));
					
					double matScore = Num.matScore(allT, allC);
					
					//pseudoMedian
					double psmMeanRatio ;
					double psmMeanLogRatio ;
					double psmMeanRelDiff ;
					double psmLayeredRatio ;
					double psmLayeredLogRatio ;
					double psmLayeredRelDiff ;
					double psmGeoMeanRatio ;
					double psmGeoMeanLogRatio ;
					double psmGeoMeanRelDiff ;
					
					if (ratios.length == 1){
						psmMeanRatio = ratios[0] * varStbl;
						psmMeanLogRatio = logRatios[0] * varStbl;
						psmMeanRelDiff = relDiffs[0] * varStbl;
						psmGeoMeanRatio = geoRatios[0] * varStbl;
						psmGeoMeanLogRatio = geoLogRatios[0] * varStbl;
						psmGeoMeanRelDiff = geoRelDiffs[0] * varStbl;
					}
					else {
						psmMeanRatio = Num.pseudoMedian(ratios) * varStbl;
						psmMeanLogRatio = Num.pseudoMedian(logRatios) * varStbl;
						psmMeanRelDiff = Num.pseudoMedian(relDiffs) * varStbl;
						psmGeoMeanRatio = Num.pseudoMedian(geoRatios) * varStbl;
						psmGeoMeanLogRatio = Num.pseudoMedian(geoLogRatios) * varStbl;
						psmGeoMeanRelDiff = Num.pseudoMedian(geoRelDiffs) * varStbl;
					}
					if (layeredRatios.length == 1){
						psmLayeredRatio = layeredRatios[0] * varStbl;
						psmLayeredLogRatio = layeredLog2Ratios[0] * varStbl;
						psmLayeredRelDiff = layeredRelDiffs[0] * varStbl;
					}
					else{
						psmLayeredRatio = Num.pseudoMedian(layeredRatios) * varStbl;
						psmLayeredLogRatio = Num.pseudoMedian(layeredLog2Ratios) * varStbl;
						psmLayeredRelDiff = Num.pseudoMedian(layeredRelDiffs) * varStbl;
					}
					
					
					
					//trimmed mean, if < 3 probes, will just take the mean!
					trmMeanLogRatio = trmMeanLogRatio * varStbl;
					double trmMeanRatio = Num.trimmedMean(ratios,0.1) * varStbl;
					double trmMeanRelDiff = Num.trimmedMean(relDiffs,0.1) * varStbl;
					double trmLayeredRatio = Num.trimmedMean(layeredRatios,0.1) * varStbl;
					double trmLayeredLogRatio = Num.trimmedMean(layeredLog2Ratios,0.1) * varStbl;
					double trmLayeredRelDiff = Num.trimmedMean(layeredRelDiffs,0.1) * varStbl;
					
					double trmGeoMeanRatio = Num.trimmedMean(geoRatios,0.1) * varStbl;
					double trmGeoMeanRelDiff = Num.trimmedMean(geoRelDiffs,0.1) * varStbl;
					double trmGeoMeanLogRatio = Num.trimmedMean(geoLogRatios,0.1) * varStbl;
					
					//Tukeybiweight
					double tkbMeanRatio = Num.tukeyBiWeight(ratios) * varStbl;
					double tkbMeanLogRatio = Num.tukeyBiWeight(logRatios) * varStbl;
					double tkbMeanRelDiff = Num.tukeyBiWeight(relDiffs) * varStbl;
					double tkbLayeredRatio = Num.tukeyBiWeight(layeredRatios) * varStbl;
					double tkbLayeredLogRatio = Num.tukeyBiWeight(layeredLog2Ratios) * varStbl;
					double tkbLayeredRelDiff = Num.tukeyBiWeight(layeredRelDiffs) * varStbl;
					
					double tkbGeoMeanRatio = Num.tukeyBiWeight(geoRatios) * varStbl;
					double tkbGeoMeanLogRatio = Num.tukeyBiWeight(geoLogRatios) * varStbl;
					double tkbGeoMeanRelDiff = Num.tukeyBiWeight(geoRelDiffs) * varStbl;
					
					//geometric mean, the only one that will work is the first!
					double gmMeanRatio = Num.geometricMean(ratios) * varStbl;
					double gmLayeredRatio = Num.geometricMean(layeredRatios) * varStbl;
					double gmGeoMeanRatio = Num.geometricMean(geoRatios) * varStbl;
					
					//median
					double mdMeanRatio = Num.median(ratios) * varStbl;
					double mdMeanLogRatio = Num.median(logRatios) * varStbl;
					double mdMeanRelDiff = Num.median(relDiffs) * varStbl;
					double mdLayeredRatio = Num.median(layeredRatios) * varStbl;
					double mdLayeredLogRatio = Num.median(layeredLog2Ratios) * varStbl;
					double mdLayeredRelDiff = Num.median(layeredRelDiffs) * varStbl;
					
					double mdGeoMeanRatio = Num.median(geoRatios) * varStbl;
					double mdGeoMeanLogRatio = Num.median(geoLogRatios) * varStbl;
					double mdGeoMeanRelDiff = Num.median(geoRelDiffs) * varStbl;
					
					//mean
					double mMeanRatio = Num.mean(ratios) * varStbl;
					double mMeanLogRatio = Num.mean(logRatios) * varStbl;
					double mMeanRelDiff = Num.mean(relDiffs) * varStbl;
					double mLayeredRatio = Num.mean(layeredRatios) * varStbl;
					double mLayeredLogRatio = Num.mean(layeredLog2Ratios) * varStbl;
					double mLayeredRelDiff = Num.mean(layeredRelDiffs) * varStbl;
					
					double mGeoMeanRatio = Num.mean(geoRatios) * varStbl;
					double mGeoMeanLogRatio = Num.mean(geoLogRatios) * varStbl;
					double mGeoMeanRelDiff = Num.mean(geoRelDiffs) * varStbl;
					
					
					//print .sgr files for sum test and ave intensity value? centered in window
					double diff = sizeOfOligoMinusOne + positions[startStop[1]] - positions[startStop[0]];
					String pos = chromosome +"\t"+ ( (int)Math.round(diff/2.0) + positions[startStop[0]] )+"\t";
					
					//print sgr lines	
					if (writeSgrs){
						wilcoxonRankSumLogPValueFile.println(pos+	wilcoxonRankSumLogPValue);
						psmMeanRatioFile.println(pos+	psmMeanRatio );
						psmMeanLogRatioFile.println(pos+ psmMeanLogRatio  );
						psmMeanRelDiffFile.println(pos+ psmMeanRelDiff  );
						psmLayeredRatioFile.println(pos+ psmLayeredRatio  );
						psmLayeredLogRatioFile.println(pos+ psmLayeredLogRatio  );
						psmLayeredRelDiffFile.println(pos+ psmLayeredRelDiff  );
						trmMeanRatioFile.println(pos+ trmMeanRatio  );
						trmMeanLogRatioFile.println(pos+ trmMeanLogRatio  );
						trmMeanRelDiffFile.println(pos+ trmMeanRelDiff  );
						trmLayeredRatioFile.println(pos+ trmLayeredRatio  );
						trmLayeredLogRatioFile.println(pos+ trmLayeredLogRatio  );
						trmLayeredRelDiffFile.println(pos+ trmLayeredRelDiff  );
						tkbMeanRatioFile.println(pos+ tkbMeanRatio  );
						tkbMeanLogRatioFile.println(pos+ tkbMeanLogRatio  );
						tkbMeanRelDiffFile.println(pos+ tkbMeanRelDiff  );
						tkbLayeredRatioFile.println(pos+ tkbLayeredRatio  );
						tkbLayeredLogRatioFile.println(pos+ tkbLayeredLogRatio  );
						tkbLayeredRelDiffFile.println(pos+ tkbLayeredRelDiff  );
						gmMeanRatioFile.println(pos+ gmMeanRatio  );
						gmLayeredRatioFile.println(pos+ gmLayeredRatio  );
						mdMeanRatioFile.println(pos+ mdMeanRatio  );
						mdMeanLogRatioFile.println(pos+ mdMeanLogRatio  );
						mdMeanRelDiffFile.println(pos+ mdMeanRelDiff  );
						mdLayeredRatioFile.println(pos+ mdLayeredRatio  );
						mdLayeredLogRatioFile.println(pos+ mdLayeredLogRatio  );
						mdLayeredRelDiffFile.println(pos+ mdLayeredRelDiff  );
						mMeanRatioFile.println(pos+ mMeanRatio  );
						mMeanLogRatioFile.println(pos+ mMeanLogRatio  );
						mMeanRelDiffFile.println(pos+ mMeanRelDiff  );
						mLayeredRatioFile.println(pos+ mLayeredRatio  );
						mLayeredLogRatioFile.println(pos+ mLayeredLogRatio  );
						mLayeredRelDiffFile.println(pos+ mLayeredRelDiff  );
						matFile.println(pos+ matScore  );
						
					}
					
					//save window
					double[] scores= new double[] {
							wilcoxonRankSumLogPValue,
							psmMeanRatio,
							psmMeanLogRatio,
							psmMeanRelDiff,
							psmGeoMeanRatio,
							psmGeoMeanLogRatio,
							psmGeoMeanRelDiff,
							psmLayeredRatio,
							psmLayeredLogRatio,
							psmLayeredRelDiff,
							trmMeanRatio,
							trmMeanLogRatio,
							trmMeanRelDiff,
							trmGeoMeanRatio,
							trmGeoMeanRelDiff,
							trmGeoMeanLogRatio,
							trmLayeredRatio,
							trmLayeredLogRatio,
							trmLayeredRelDiff,
							tkbMeanRatio,
							tkbMeanLogRatio,
							tkbMeanRelDiff,
							tkbLayeredRatio,
							tkbLayeredLogRatio,
							tkbLayeredRelDiff,
							tkbGeoMeanRatio,
							tkbGeoMeanLogRatio,
							tkbGeoMeanRelDiff,
							gmMeanRatio,
							gmLayeredRatio,
							gmGeoMeanRatio,
							mdMeanRatio,
							mdMeanLogRatio,
							mdMeanRelDiff,
							mdLayeredRatio,
							mdLayeredLogRatio,
							mdLayeredRelDiff,
							mdGeoMeanRatio,
							mdGeoMeanLogRatio,
							mdGeoMeanRelDiff,
							mMeanRatio,
							mMeanLogRatio,
							mMeanRelDiff,
							mGeoMeanRatio,
							mGeoMeanLogRatio,
							mGeoMeanRelDiff,
							mLayeredRatio,
							mLayeredLogRatio,
							mLayeredRelDiff,
							matScore};
					Window win = new Window(chromosome, positions[startStop[0]],positions[startStop[1]], windowSize, scores);
					windowsAl.add(win);
				}
				
				
			}
			
			//close files 
			if (writeSgrs){
				wilcoxonRankSumLogPValueFile.close(); 
				psmMeanRatioFile.close();
				psmMeanLogRatioFile.close();
				psmMeanRelDiffFile.close();
				psmLayeredRatioFile.close();
				psmLayeredLogRatioFile.close();
				psmLayeredRelDiffFile.close();
				trmMeanRatioFile.close();
				trmMeanLogRatioFile.close();
				trmMeanRelDiffFile.close();
				trmLayeredRatioFile.close();
				trmLayeredLogRatioFile.close();
				trmLayeredRelDiffFile.close();
				tkbMeanRatioFile.close();
				tkbMeanLogRatioFile.close();
				tkbMeanRelDiffFile.close();
				tkbLayeredRatioFile.close();
				tkbLayeredLogRatioFile.close();
				tkbLayeredRelDiffFile.close();
				gmMeanRatioFile.close();
				gmLayeredRatioFile.close();
				mdMeanRatioFile.close();
				mdMeanLogRatioFile.close();
				mdMeanRelDiffFile.close();
				mdLayeredRatioFile.close();
				mdLayeredLogRatioFile.close();
				mdLayeredRelDiffFile.close();
				mMeanRatioFile.close();
				mMeanLogRatioFile.close();
				mMeanRelDiffFile.close();
				mLayeredRatioFile.close();
				mLayeredLogRatioFile.close();
				mLayeredRelDiffFile.close();
				matFile.close();
			}
			
			
			//delete Tmp files
			IO.deleteFiles(treatmentDir, uniqueIdT);
			IO.deleteFiles(controlDir, uniqueIdC);
			
			//convert ArrayList to Window[] 
			numWindows = windowsAl.size();
			Window[] windows = new Window[numWindows];
			windowsAl.toArray(windows);
			windowsAl.clear();
			System.out.println("Number of windows screened: "+numWindows+"\n");
			File serWinIntArray = new File(resultsFile);
			System.out.print ("Saving windows "+serWinIntArray+"\n");
			IO.saveObject(serWinIntArray, windows);
			
		}catch (IOException e){
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		if (args.length<4){
			printDocs();
			System.exit(0);
		}	
		new TestScanChip(args);
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
					case 't': treatmentDir = args[i+1]; i++; break;
					case 'c': controlDir = args[i+1]; i++; break;
					case 'r': resultsFile = args[i+1]; i++; break;
					case 'b': tpmapFile = new File(args[i+1],"tpmap.fa").toString(); i++; break;
					case 'p': pairs = true; break;
					case 'n': multiplyByNumberValues = true; break;
					case 'z': sizeOfOligoMinusOne =Integer.parseInt(args[i+1])-1;i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					
					//don't forget to set square root boolean and print sgr boolean
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		//check if directories are directories
		File treat = new File(treatmentDir);
		File control = new File(controlDir);
		if (treat==null || control==null || treat.isDirectory()==false || control.isDirectory()== false){
			System.out.println("\nCheck your directories, one is incorrect!\n");
			System.exit(1);
		}	
		infoFile = tpmapFile+"Info";
		if (new File(infoFile).exists()==false){
			System.out.println("\nCheck your TPMapFiles directory, it appears to be incomplete!\n");
			System.exit(1);
		}
	}	
	
	public static void printDocs(){ 
		System.out.println("\n" 
		);		
	}
	
}
