package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import net.sf.samtools.*;
import util.bio.annotation.Bed;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.gen.*;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.useq.apps.Text2USeq;


/**
 * @author Nix
 * */
public class OverdispersedRegionScanSeqs {

	//user defined fields
	private File[] treatmentBamFiles;
	private File[] controlBamFiles;
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private File bedFile;
	private File refSeqFile;
	private boolean dataIsStranded = false;
	private float minimumSpliceLog2Ratio = 1;
	private float minFDR = 10;
	private float minLog2Ratio = 1f;
	private int minReads = 20;
	private boolean calculateChiSquare = true;
	private boolean scoreIntrons = false;
	private boolean removeOverlappingRegions = true;
	private boolean printExonReplicaMatrix = false;
	private boolean deleteTempFiles = true;
	private boolean usePseudoMedianLog2Ratio = false;

	//internal fields
	private HashMap<String,UCSCGeneLine[]> geneModels;
	private UCSCGeneLine[] allGeneLines;
	private ArrayList<UCSCGeneLine> geneLinesWithReadsAL = new ArrayList<UCSCGeneLine>();
	private ArrayList<UCSCGeneLine> geneLinesWithOutReadsAL = new ArrayList<UCSCGeneLine>();
	private UCSCGeneLine[] geneLinesWithReads;
	private String chromosome;
	private SAMFileReader[] treatmentReaders;
	private SAMFileReader[] controlReaders;
	private boolean reverseStrand;
	private boolean verbose = true;
	private boolean pseudoMedianCalculated = false;
	private boolean varianceOutlierStatsPresent = false;
	private int[][] tObs;
	private int[][] cObs;
	private static Pattern CIGAR_SUB = Pattern.compile("(\\d+)([MSDHN])");

	private double totalInterrogatedTreatmentObservations = 0;
	private double totalInterrogatedControlObservations = 0;
	private int numberTreatmentReplicas;
	private int numberControlReplicas;
	private int totalNumberReplicas;
	private double millionMappedTreatmentReads;
	private double millionMappedControlReads;
	private String genomeVersion;
	private double scalarTC;
	private double scalarCT;
	private int numberGenesPassingThresholds;
	private int log2ScoreIndex = 3;
	private int chiSquarePValIndex = 4;
	private int chiSquareLog2RatioIndex = 5;
	private int numberOfScores = -1;

	//constructors

	/**For integration with RNASeq app. Note processedRefSeqFile should contain non overlapping exons! */
	public OverdispersedRegionScanSeqs(File[] treatmentBamFiles, File[] controlBamFiles, String genomeVersion, File saveDirectory, File fullPathToR, File processedRefSeqFile, boolean dataIsStranded, boolean scoreIntrons, boolean verbose){
		this.treatmentBamFiles = treatmentBamFiles;
		this.controlBamFiles = controlBamFiles;
		this.saveDirectory = saveDirectory;
		this.genomeVersion = genomeVersion;
		this.fullPathToR = fullPathToR;
		this.refSeqFile = processedRefSeqFile;
		removeOverlappingRegions = false;
		this.dataIsStranded = dataIsStranded;
		this.scoreIntrons = scoreIntrons;
		this.verbose = verbose;

		numberTreatmentReplicas = treatmentBamFiles.length;
		numberControlReplicas = controlBamFiles.length;
		totalNumberReplicas = numberTreatmentReplicas + numberControlReplicas;
		run();
	}

	/**Stand alone.*/
	public OverdispersedRegionScanSeqs(String[] args){
		long startTime = System.currentTimeMillis();
		//set fields
		processArgs(args);
		//launch
		run();
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		if (verbose) System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}

	public void makeSamReaders(){
		treatmentReaders = new SAMFileReader[treatmentBamFiles.length];
		for (int i=0; i< treatmentReaders.length; i++) {
			treatmentReaders[i] = new SAMFileReader(treatmentBamFiles[i]);
			treatmentReaders[i].enableIndexMemoryMapping(false);
		}
		controlReaders = new SAMFileReader[controlBamFiles.length];
		for (int i=0; i< controlReaders.length; i++) {
			controlReaders[i] = new SAMFileReader(controlBamFiles[i]);
			controlReaders[i].enableIndexMemoryMapping(false);
		}
	}

	public void closeSamReaders(){
		for (int i=0; i< treatmentReaders.length; i++) treatmentReaders[i].close();
		for (int i=0; i< controlReaders.length; i++) controlReaders[i].close();
	}


	public void run(){
		//load gene models
		if (verbose) System.out.println("Loading regions/ gene models...");
		try {
			loadGeneModels();
		} catch (IOException e1) {
			System.err.println(e1.getMessage());
			//e1.printStackTrace();
			return;
		}

		//make readers on each bam file
		makeSamReaders();

		//for each chromosome of gene models
		if (verbose) System.out.print("\nScanning regions by chromosome (very slow using Picard query lookup)... ");
		Iterator<String> it = geneModels.keySet().iterator();
		while (it.hasNext()){
			chromosome = it.next();
			if (verbose) {
				System.out.print(" ");
				System.out.print(chromosome);
			}
			//if (chromosome.equals("chrM") == false) continue;
			scanGenes();
		}
		if (verbose) System.out.println("\n");

		//close readers
		closeSamReaders();

		//call after running the ScanGenes so that the interrogated counts are accurate
		if (verbose) System.out.println("Calculating read count stats...");
		calculateReadCountStatistics();

		//load geneLines with those that have observations
		geneLinesWithReads = new UCSCGeneLine[geneLinesWithReadsAL.size()];
		geneLinesWithReadsAL.toArray(geneLinesWithReads);

		//convert binomial pvalues to FDRs
		if (verbose) System.out.println("\nCalculating negative binomial p-values and FDRs in R using DESeq (http://www-huber.embl.de/users/anders/DESeq/)...");
		calculateDESeqPValues();

		//append on pseudo median log2 ratios?
		if (numberTreatmentReplicas > 2 && numberControlReplicas > 2){
			if (verbose) {
				if (usePseudoMedianLog2Ratio) System.out.println("\nCalculating pseudo median log2 ratios...");
				else {
					System.out.println("\nCalculating smallest all pair log2 ratio per gene...");
					if (numberTreatmentReplicas >3 || numberControlReplicas >3) System.err.println("\nWARNING! Use the pseudo median log2 ratio calculation for datasets with 4 or more replicas!");
				}
			}
			appendLog2Ratios();
			pseudoMedianCalculated = true;
		}

		//calculate chiSquare test of independence
		if (calculateChiSquare){ 
			if (verbose) System.out.println("\nEstimating alternative splicing read distribution, bonferroni corrected, chi-square p-values in R (very slow)...");
			if (printExonReplicaMatrix) estimateDifferencesInReadDistributionsWithReplicas();
			else estimateDifferencesInReadDistributions();
		}

		//sort by pValue
		if (verbose) System.out.println("\nThresholding and printing results...");
		Arrays.sort(geneLinesWithReads, new UCSCGeneLineComparatorScoreBigToSmall(1));

		//print all gene report
		if (verbose) System.out.println("\tSave dir\t"+saveDirectory);
		File allGenesFile = new File (saveDirectory, "all.xls");
		if (verbose) System.out.println("\t"+allGenesFile.getName()+"\t: all gene/region spreadsheet");
		printGeneModels(geneLinesWithReads, allGenesFile, true);

		//threshold genes
		UCSCGeneLine[] good = thresholdGenes();
		numberGenesPassingThresholds = good.length;

		if (good.length != 0){
			Arrays.sort(good, new UCSCGeneLineComparatorScoreBigToSmall(1));
			String name = "diffExprFDR"+minFDR+"LgRto"+Num.formatNumber(minLog2Ratio,3);
			File partialGenesFile = new File (saveDirectory, name+".xls.gz");
			printGeneModels(good, partialGenesFile, false);

			//print egr file 
			File egrFile = new File (saveDirectory, name+".egr.gz");
			printGeneModelsEgr(good, egrFile);

			//print bed file 
			File bedFileUp = new File (saveDirectory, name+"_Up.bed.gz");
			File bedFileDown = new File (saveDirectory, name+"_Down.bed.gz");
			printGeneModelsBed(good, bedFileUp, bedFileDown);

			//convert to useq
			try {
				if (bedFileUp != null && bedFileUp.exists()) new Text2USeq(bedFileUp, genomeVersion, "#FF0000");
				if (bedFileDown != null && bedFileDown.exists()) new Text2USeq(bedFileDown, genomeVersion, "#00FF00");
			} catch (IOException e) {
				e.printStackTrace();
			}


		}
		if (verbose) System.out.println("\n"+numberGenesPassingThresholds+" genes/regions passed thresholds.");
	}

	/**Replaces the log2 ratio with the pseudo median*/
	public void appendLog2Ratios(){
		//fetch em
		double[][] pseAndSmallest = Num.calculatePseudoMedianLog2Ratios(tObs, cObs);
		double[] log2;
		if (usePseudoMedianLog2Ratio) log2 = pseAndSmallest[0];
		else log2 = pseAndSmallest[1];
		//set em
		for (int x=0; x< geneLinesWithReads.length; x++){
			//get scores
			float[] scores = geneLinesWithReads[x].getScores();
			//replace log2 ratio 
			scores[log2ScoreIndex] = (float)log2[x];
		}
	}

	/**Calculates a log2( (tSum+1)/(cSum+1) ) on linearly scaled tSum and cSum based on the total observations.*/
	public float calculateLog2Ratio( double tSum, double cSum){
		double t;
		double c;
		if (tSum !=0 ) {
			t = tSum * scalarCT;
			c = cSum;
		}
		else {
			c = cSum * scalarTC;
			t = tSum;
		}
		double ratio = (t+1)/(c+1);
		return (float)Num.log2(ratio);
	}


	public UCSCGeneLine[] thresholdGenes(){
		ArrayList<UCSCGeneLine> passing = new ArrayList<UCSCGeneLine>();
		for (int i=0; i< geneLinesWithReads.length; i++){
			//scores = totalExonicBPs, pval, padj, (might have pval2, padj2) log2ratio, chiSquarePVal, chiSquareLog2Ratio, tFPKM, cFPKM, countsT1,2,3, countsC1,2,3, varCorrT1,2,3,  varCorrC1,2,3
			float[] scores = geneLinesWithReads[i].getScores();
			if (scores[2]>= minFDR && Math.abs(scores[log2ScoreIndex])>=minLog2Ratio) passing.add(geneLinesWithReads[i]);
		}
		UCSCGeneLine[] good = new UCSCGeneLine[passing.size()];
		passing.toArray(good);
		return good;
	}

	public File writeOutObservations(String randomWord){
		File matrixFile = null;
		try {
			//load matrix for possible pseudo median calculation
			tObs = new int[geneLinesWithReads.length][numberTreatmentReplicas];
			cObs = new int[geneLinesWithReads.length][numberControlReplicas];

			//write matrix of name, t1,t2,t3...c1,c2,c3 to file for genes with observations
			matrixFile = new File(saveDirectory, randomWord+"_CountMatrix.txt");
			PrintWriter out = new PrintWriter( new FileWriter(matrixFile));
			for (int i=0; i< geneLinesWithReads.length; i++){
				float[] scores = geneLinesWithReads[i].getScores();
				out.print(geneLinesWithReads[i].getDisplayNameThenName());
				//out.print("G"+i);
				for (int j=0; j< scores.length; j++){
					out.print("\t"+(int)scores[j]);
				}
				out.println();
				//load scores for pse
				for (int c=0; c<numberTreatmentReplicas; c++) tObs[i][c]= (int)scores[c]; 
				int index = 0;
				for (int c=numberTreatmentReplicas; c<totalNumberReplicas; c++) cObs[i][index++] = (int)scores[c];
			}
			out.close();
		}
		catch (Exception e){
			System.err.println("Problem writing out observations for R.");
			e.printStackTrace();
		}
		return matrixFile;
	}


	/**Returns the deseq stats file and the variance stabilized data.
	private File[] executeDESeqBlind(File matrixFile, String randomWord){
		File rResultsStats = new File (saveDirectory, randomWord+"_DESeqBlindResultsStats.txt");
		File rResultsVarCorr = new File (saveDirectory, randomWord+"_DESeqBlindResultsVarCorr.txt");
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("library(DESeq)\n");
			sb.append("countsTable = read.delim('"+matrixFile.getCanonicalPath()+"', header=FALSE)\n");
			sb.append("rownames(countsTable) = countsTable$V1\n");
			sb.append("countsTable = countsTable[,-1]\n");
			sb.append("conds = c(");
			for (int i=0; i< numberTreatmentReplicas; i++) sb.append("'T',");
			sb.append("'N'");
			for (int i=1; i< numberControlReplicas; i++) sb.append(",'N'");
			sb.append(")\n");
			sb.append("cds = newCountDataSet( countsTable, conds)\n");
			sb.append("cds = estimateSizeFactors( cds )\n");
			sb.append("cds = estimateDispersions( cds, method='blind', sharingMode='fit-only' )\n");
			sb.append("res = nbinomTest( cds, 'N', 'T')\n");
			//Recalculate log2 ratio with add one, note flip
			sb.append("res[,6] = log2((1+res[,4])/(1+res[,3]))\n");
			//Fred  pvalues
			sb.append("res[,7] = -10 * log10(res[,7])\n");
			sb.append("res[,8] = -10 * log10(res[,8])\n");
			//Parse padj, log2ratio; note flip of A and B back to T and C
			sb.append("res = res[,c(7,8,6)]\n");
			//note, the order of the rows is the same as the input
			sb.append("write.table(res, file = '"+rResultsStats.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");
			//pull variance stabilized data for each replica, note these are not flipped
			sb.append("res = getVarianceStabilizedData( cds )\n");
			sb.append("write.table(res, file = '"+rResultsVarCorr.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");

			//write script to file
			File scriptFile = new File (saveDirectory, randomWord+"_RScript.txt");
			File rOut = new File(saveDirectory, randomWord+"_RScript.txt.Rout");
			IO.writeString(sb.toString(), scriptFile);

			//make command
			String[] command = new String[] {
					fullPathToR.getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					scriptFile.getCanonicalPath(),
					rOut.getCanonicalPath()};

			//execute command
			IO.executeCommandLine(command);

			//read in results watching for Inf values from R
			if (rResultsStats.exists() == false || rResultsVarCorr.exists() == false) throw new IOException("\nR results file doesn't exist. Check temp files in save directory for error.\n");

			//check for warnings
			String[] res = IO.loadFile(rOut);
			for (String s: res) if (s.contains("Warning message")) throw new IOException("\nR outputFile contains a warning message. Check, fix, and restart. See -> "+rOut+"\n");

			//cleanup
			if (deleteTempFiles){
				rOut.deleteOnExit();
				scriptFile.deleteOnExit();
				rResultsStats.deleteOnExit();
				rResultsVarCorr.deleteOnExit();
			}

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
			return null;
		}
		return new File[]{rResultsStats,rResultsVarCorr};
	}*/


	/**Returns the deseq stats file and the variance stabilized data.*/
	private File[] executeDESeq(File matrixFile, String randomWord, boolean filterOutliers, boolean useLocalFitType, boolean useMethodBlind){
		File rResultsStats = new File (saveDirectory, randomWord+"_DESeqResultsStats.txt");
		File rResultsVarCorr = new File (saveDirectory, randomWord+"_DESeqResultsVarCorr.txt");
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("library(DESeq)\n");
			sb.append("countsTable = read.delim('"+matrixFile.getCanonicalPath()+"', header=FALSE)\n");
			sb.append("rownames(countsTable) = countsTable$V1\n");
			sb.append("countsTable = countsTable[,-1]\n");
			sb.append("conds = c(");
			for (int i=0; i< numberTreatmentReplicas; i++) sb.append("'T',");
			sb.append("'N'");
			for (int i=1; i< numberControlReplicas; i++) sb.append(",'N'");
			sb.append(")\n");
			sb.append("cds = newCountDataSet( countsTable, conds)\n");
			sb.append("cds = estimateSizeFactors( cds )\n");

			String outlier = "";
			String local = "";
			String method = "";
			if (filterOutliers == false) outlier= ", sharingMode='fit-only'";
			if (useLocalFitType) local = ", fitType='local'";
			if (useMethodBlind) method = ", method='blind'";

			sb.append("cds = estimateDispersions( cds" +outlier + local + method+")\n");
			sb.append("res = nbinomTest( cds, 'N', 'T')\n");
			//Recalculate log2 ratio with add one, note flip
			sb.append("res[,6] = log2((1+res[,4])/(1+res[,3]))\n");
			//Fred  pvalues
			sb.append("res[,7] = -10 * log10(res[,7])\n");
			sb.append("res[,8] = -10 * log10(res[,8])\n");
			//Parse padj, log2ratio; note flip of A and B back to T and C
			sb.append("res = res[,c(7,8,6)]\n");
			//note, the order of the rows is the same as the input
			sb.append("write.table(res, file = '"+rResultsStats.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");
			//pull variance stabilized data for each replica, note these are not flipped
			sb.append("res = getVarianceStabilizedData( cds )\n");
			sb.append("write.table(res, file = '"+rResultsVarCorr.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");

			//write script to file
			File scriptFile = new File (saveDirectory, randomWord+"_RScript.txt");
			File rOut = new File(saveDirectory, randomWord+"_RScript.txt.Rout");
			IO.writeString(sb.toString(), scriptFile);

			//make command
			String[] command = new String[] {
					fullPathToR.getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					scriptFile.getCanonicalPath(),
					rOut.getCanonicalPath()};

			//execute command
			IO.executeCommandLine(command);

			//check for warnings
			if (rOut.exists()){
				String[] res = IO.loadFile(rOut);
				for (String s: res) {
					if (s.contains("Dispersion fit did not converge") || s.contains("Parametric dispersion fit failed")) {
						System.err.println("\n\tWarning, DESeq's GLM dispersion fit failed. Relaunching using fitType='local'");
						if (deleteTempFiles == false){
							rOut.delete();
							scriptFile.delete();
						}
						return executeDESeq(matrixFile, randomWord+"_Local", filterOutliers, true, useMethodBlind);
					}
				}
			}

			//look for results files
			if (rResultsStats.exists() == false || rResultsVarCorr.exists() == false) throw new IOException("\nDESeq results files don't exist? Check temp files (xxx.outR) in save directory for error.\n");


			//cleanup
			if (deleteTempFiles){
				rOut.deleteOnExit();
				scriptFile.deleteOnExit();
				rResultsStats.deleteOnExit();
				rResultsVarCorr.deleteOnExit();
			}

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
			return null;
		}
		return new File[]{rResultsStats,rResultsVarCorr};
	}



	private void parseDESeqStatResults(File results){
		try {
			BufferedReader in = new BufferedReader (new FileReader (results));
			String line;
			String[] tokens;
			Pattern tab = Pattern.compile("\t");
			float maxPVal = 0;
			float maxAdjPVal = 0;
			numberOfScores = 8+ (totalNumberReplicas*2);
			for (int x=0; x< geneLinesWithReads.length; x++){
				line= in.readLine();
				//parse line: pval, padj, log2ratio
				tokens = tab.split(line);
				if (tokens.length!=3) Misc.printErrAndExit("One of the DESeq stats R results rows is malformed -> "+line);
				//get existing scores: read counts in T's and C's
				float[] scores = geneLinesWithReads[x].getScores();
				//build: totalExonicBPs, pval, padj, log2ratio, chiSquarePval, chiSquareLog2, tFPKM, cFPKM, countsT1,2,3, countsC1,2,3, varCorrT1,2,3,  varCorrC1,2,3
				float[] allScores = new float[numberOfScores];
				//totalExonicBPs
				allScores[0] = geneLinesWithReads[x].getTotalExonicBasePairs();
				//pval
				allScores[1] = parseFloatWatchForBad(tokens[0]);
				if (allScores[1]> maxPVal) maxPVal = allScores[1];
				//adjPval
				allScores[2] = parseFloatWatchForBad(tokens[1]);
				if (allScores[2]> maxAdjPVal) maxAdjPVal = allScores[2];
				//log2ratio
				allScores[log2ScoreIndex] = Float.parseFloat(tokens[2]);
				//chisquare, wait for block calculation
				//allScores[4], allScores[5]
				//tFPKM
				int totalT = 0;
				for (int i=0; i<numberTreatmentReplicas; i++) totalT+= scores[i]; 
				allScores[6] = calculateFPKM(millionMappedTreatmentReads, allScores[0], totalT);
				//cFPKM
				int totalC = 0;
				for (int i=numberTreatmentReplicas; i<totalNumberReplicas; i++) totalC+= scores[i]; 
				allScores[7] = calculateFPKM(millionMappedControlReads, allScores[0], totalC);
				//add counts
				System.arraycopy(scores, 0, allScores, 8, totalNumberReplicas);
				//add
				geneLinesWithReads[x].setScores(allScores);
			}
			in.close();
			//convert Inf to max values * 1%
			maxPVal = maxPVal *1.01f;
			maxAdjPVal = maxAdjPVal * 1.01f;
			for (int i=0; i< geneLinesWithReads.length; i++){
				float[] scores = geneLinesWithReads[i].getScores();
				//check pval
				if (scores[1] == Float.MIN_VALUE) scores[1] = maxPVal;
				//check adjPVal
				if (scores[2] == Float.MIN_VALUE) scores[2] = maxAdjPVal;
			}
		} catch (Exception e){
			System.err.println("Problem parsing DESeq stats results from R.\n");
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void parseDESeqVarianceOutlierStats(File results){
		try {
			BufferedReader in = new BufferedReader (new FileReader (results));
			String line;
			String[] tokens;
			Pattern tab = Pattern.compile("\t");
			float maxPVal = 0;
			float maxAdjPVal = 0;
			for (int x=0; x< geneLinesWithReads.length; x++){
				line= in.readLine();
				//parse line: pval, padj, log2ratio
				tokens = tab.split(line);
				if (tokens.length!=3) Misc.printErrAndExit("One of the DESeq stats R results rows is malformed -> "+line);
				float pval= parseFloatWatchForBad(tokens[0]);
				if (pval > maxPVal) maxPVal = pval;
				float padj = parseFloatWatchForBad(tokens[1]);
				if (padj> maxAdjPVal) maxAdjPVal = padj;

				//get the scores : totalExonicBPs, pval, padj, log2ratio, chiSquarePval, chiSquareLog2, tFPKM, cFPKM, countsT1,2,3, countsC1,2,3, varCorrT1,2,3,  varCorrC1,2,3
				float[] scores = geneLinesWithReads[x].getScores();

				//insert additional pval and padj after first ones: totalExonicBPs, pval, padj, pval2, padj2, log2ratio, chiSquarePval, chiSquareLog2, tFPKM, cFPKM, countsT1,2,3, countsC1,2,3, varCorrT1,2,3,  varCorrC1,2,3
				float[] allScores = new float[scores.length+2];
				allScores[0] = scores[0];
				allScores[1] = scores[1];
				allScores[2] = scores[2];
				allScores[3] = pval;
				allScores[4] = padj;
				System.arraycopy(scores, 3, allScores, 5, scores.length-3);

				//add
				geneLinesWithReads[x].setScores(allScores);
			}
			in.close();
			//convert Inf to max values * 1%
			maxPVal = maxPVal *1.01f;
			maxAdjPVal = maxAdjPVal * 1.01f;
			for (int i=0; i< geneLinesWithReads.length; i++){
				float[] scores = geneLinesWithReads[i].getScores();
				//check pval
				if (scores[3] == Float.MIN_VALUE) scores[3] = maxPVal;
				//check adjPVal
				if (scores[4] == Float.MIN_VALUE) scores[4] = maxAdjPVal;
			}
			//set flags for other methods
			varianceOutlierStatsPresent = true;
			log2ScoreIndex+=2;
			chiSquareLog2RatioIndex+=2;
			chiSquarePValIndex+=2;
			numberOfScores+=2;

		} catch (Exception e){
			System.err.println("Problem parsing DESeq filtered stats results from R.\n");
			e.printStackTrace();
			System.exit(1);
		}
	}



	private void parseDESeqDataResults(File results){
		try {
			BufferedReader in = new BufferedReader (new FileReader (results));
			String line;
			String[] tokens;
			Pattern tab = Pattern.compile("\t");
			for (int x=0; x< geneLinesWithReads.length; x++){
				line= in.readLine();
				//parse line: pval, padj, log2ratio
				tokens = tab.split(line);
				if (tokens.length!=totalNumberReplicas) Misc.printErrAndExit("One of the DESeq data R results rows is malformed -> "+line);
				//parse the scores
				float[] vars = new float[totalNumberReplicas];
				for (int i=0; i< totalNumberReplicas; i++) vars[i] = Float.parseFloat(tokens[i]);
				//get scores
				float[] scores = geneLinesWithReads[x].getScores();
				//copy in
				System.arraycopy(vars, 0, scores, scores.length-totalNumberReplicas, totalNumberReplicas);
				//add
				geneLinesWithReads[x].setScores(scores);
			}
			in.close();
		} catch (Exception e){
			System.err.println("Problem parsing DESeq data results from R.\n");
			e.printStackTrace();
			System.exit(1);
		}
	}

	public float parseFloatWatchForBad(String f){
		if (f.equals("Inf") || f.equals("NA")) return Float.MIN_VALUE;
		else return Float.parseFloat(f);
	}

	public void calculateDESeqPValues(){
		//write matrix of name, t1,t2,t3...c1,c2,c3 to file for genes with observations
		String random = Passwords.createRandowWord(6);
		File matrixFile = writeOutObservations(random);
		//call DESeq
		File[] deseqResults;
		File[] deseqResultsFiltered = null;
		if (totalNumberReplicas ==2){
			//no variance filtering and blind
			deseqResults = executeDESeq(matrixFile, random+"NoVarOutFilt", false, false, true);
		}
		else {
			//without variance outlier filtering
			deseqResults = executeDESeq(matrixFile, random+"NoVarOutFilt", false, false, false);
			//with
			deseqResultsFiltered = executeDESeq(matrixFile, random+"VarOutFilt", true, false, false);
			//note the VarCorr values are the same regardless of the filtering.

		}
		//parse results and populate gene scores
		parseDESeqStatResults (deseqResults[0]);
		parseDESeqDataResults (deseqResults[1]);
		if (deseqResultsFiltered != null) parseDESeqVarianceOutlierStats(deseqResultsFiltered[0]);

		//clean up
		if (deleteTempFiles) matrixFile.delete();
	}

	/**Calculates the reads per kb per million mapped reads 
	 * # Observed reads in the region/ bp size of the region / 1000/ total number reads/ 1000000 */
	public float calculateFPKM(double millionTotalMappedReads, double interrogatedRegionBPSize, double numberObservedReadsInRegion){
		double exonicBasesPerKB = interrogatedRegionBPSize/1000;
		double rpkm = numberObservedReadsInRegion/exonicBasesPerKB/millionTotalMappedReads;
		float val = new Double(rpkm).floatValue();
		return val;
	}

	public void estimateDifferencesInReadDistributions(){
		//find genes with multiple exons and max number exons
		int maxNumberExons = -1;
		ArrayList<UCSCGeneLine> al = new ArrayList<UCSCGeneLine>();
		for (int i=0; i< geneLinesWithReads.length; i++){
			//check exon conts to see if they need to be modified for the chiSquare
			if (checkForMinimums(geneLinesWithReads[i])){
				al.add(geneLinesWithReads[i]);
				int numEx = geneLinesWithReads[i].getExonCounts()[0].length;
				if (numEx > maxNumberExons) maxNumberExons = numEx;
			}
		}
		//any genes?
		if (al.size() == 0) return;

		UCSCGeneLine[] genesWithExonsAndReads = new UCSCGeneLine[al.size()];
		al.toArray(genesWithExonsAndReads);

		//collect counts
		int[][] treatment = new int[genesWithExonsAndReads.length][maxNumberExons];
		int[][] control = new int[genesWithExonsAndReads.length][maxNumberExons];

		for (int i=0; i< genesWithExonsAndReads.length; i++){
			float[][] tc = genesWithExonsAndReads[i].getExonCounts();
			Arrays.fill(treatment[i], -1);
			Arrays.fill(control[i], -1);
			int[] t = Num.convertToInt(tc[0]);
			int[] c = Num.convertToInt(tc[1]);
			System.arraycopy(t, 0, treatment[i], 0, t.length);
			System.arraycopy(c, 0, control[i], 0, c.length);
		}

		//estimate chi-square pvalues using R for resolution of extreemly small p-values, radiculously slow
		double[] pVals = Num.chiSquareIndependenceTest(treatment, control, saveDirectory, fullPathToR, true);

		//bonferroni correction
		float bc = (float)Num.minus10log10(genesWithExonsAndReads.length);

		//add back
		for (int i=0; i< genesWithExonsAndReads.length; i++){
			//set corrected p-value 
			//scores = totalExonicBPs, pval, padj, (might have pval2, padj2), log2ratio, chiSquarePVal, chiSquareLog2Ratio, tFPKM, cFPKM, countsT1,2,3, countsC1,2,3 varT1,2,3  varC1,2,3
			float[] scores = genesWithExonsAndReads[i].getScores();
			scores[chiSquarePValIndex] = (float) pVals[i] + bc;
			if (scores[chiSquarePValIndex] < 0) scores[chiSquarePValIndex] = 0;
			scores[chiSquareLog2RatioIndex] = genesWithExonsAndReads[i].getLog2Ratio();
		}
	}

	public void estimateDifferencesInReadDistributionsWithReplicas(){
		//find genes with multiple exons and max number exons
		int maxNumberExons = -1;
		ArrayList<UCSCGeneLine> al = new ArrayList<UCSCGeneLine>();
		for (int i=0; i< geneLinesWithReads.length; i++){
			int numEx = geneLinesWithReads[i].getExons().length;
			if (numEx > 1 && geneLinesWithReads[i].getTreatmentExonCounts() != null) {
				al.add(geneLinesWithReads[i]);
				if (numEx > maxNumberExons) maxNumberExons = numEx;
			}
		}
		UCSCGeneLine[] genesWithExonsAndReads = new UCSCGeneLine[al.size()];
		al.toArray(genesWithExonsAndReads);

		//just print them to file
		File matrix = new File (this.saveDirectory, "exonMatrixFile.txt");
		PrintWriter out;
		try {
			out = new PrintWriter( new FileWriter (matrix));
			out.println("# Number genes with >1 exons and counts\t"+genesWithExonsAndReads.length);
			out.println("# Number treatment replicas\t"+ numberTreatmentReplicas);
			out.println("# Number control replicas\t"+ numberControlReplicas);
			out.println("## GeneName\t Number exons");
			out.println("## Per exon treatment counts");
			out.println("## Per exon control counts");
			for (UCSCGeneLine gene: genesWithExonsAndReads){
				out.println(gene.getDisplayName()+"\t"+gene.getExons().length);
				printExonCounts(out, gene.getTreatmentExonCounts());
				printExonCounts(out, gene.getControlExonCounts());

			}
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		//collect counts
		int[][] treatment = new int[genesWithExonsAndReads.length][maxNumberExons];
		int[][] control = new int[genesWithExonsAndReads.length][maxNumberExons];

		for (int i=0; i< genesWithExonsAndReads.length; i++){
			float[][] tc = genesWithExonsAndReads[i].getExonCounts();
			Arrays.fill(treatment[i], -1);
			Arrays.fill(control[i], -1);
			int[] t = Num.convertToInt(tc[0]);
			int[] c = Num.convertToInt(tc[1]);
			System.arraycopy(t, 0, treatment[i], 0, t.length);
			System.arraycopy(c, 0, control[i], 0, c.length);
		}

		//estimate chi-square pvalues using R for resolution of extreemly small p-values, radiculously slow
		double[] pVals = Num.chiSquareIndependenceTest(treatment, control, saveDirectory, fullPathToR, true);

		//bonferroni correction
		float bc = (float)Num.minus10log10(genesWithExonsAndReads.length);

		//add back
		for (int i=0; i< genesWithExonsAndReads.length; i++){
			//set corrected p-value 
			//scores = totalExonicBPs, pval, padj, log2ratio, chiSquareLog2Ratio, chiSqrLog2Ratio, tFPKM, cFPKM, countsT1,2,3, countsC1,2,3 varT1,2,3  varC1,2,3
			float[] scores = genesWithExonsAndReads[i].getScores();
			scores[chiSquarePValIndex] = (float) pVals[i] + bc;
			if (scores[chiSquarePValIndex] < 0) scores[chiSquarePValIndex] = 0;
			scores[chiSquareLog2RatioIndex] = genesWithExonsAndReads[i].getLog2Ratio();
		}

	}

	public static void printExonCounts (PrintWriter out, float[][] repsCounts){
		int numReps = repsCounts.length;
		int numExons = repsCounts[0].length;
		for (int i=0; i< numReps; i++){
			float[] counts = repsCounts[i];
			out.print((int)counts[0]);
			for (int j =1; j< numExons; j++){
				out.print("\t");
				out.print((int)counts[j]);
			}
			out.println();
		}
	}

	public void printGeneModelsEgr(UCSCGeneLine[] genes, File results){
		try {
			Gzipper out = new Gzipper( results);
			//genome version
			out.println("# genome_version = "+genomeVersion);
			//score names
			if (usePseudoMedianLog2Ratio && pseudoMedianCalculated) out.println("# score0 = PseudoMedianLog2Ratio");
			else out.println("# score0 = Log2Ratio");
			out.println("# score1 = -10Log10(FDR_NoVarOutFilt)");
			//for each gene
			for (int i=0; i< genes.length; i++){
				out.print(genes[i].getDisplayNameThenName());
				out.print("\t");
				out.print(genes[i].getChrom());
				out.print("\t");
				out.print(genes[i].getTxStart());
				out.print("\t");
				out.print(genes[i].getTxEnd());
				out.print("\t");
				out.print(genes[i].getStrand());
				out.print("\t");
				//scores = totalExonicBPs, pval, padj, (might have pval2, padj2) log2ratio, chiSquarePVal, chiSquareLog2Ratio, tFPKM, cFPKM, countsT1,2,3, countsC1,2,3 varT1,2,3  varC1,2,3 
				float[] s = genes[i].getScores();
				out.print(s[log2ScoreIndex]);
				out.print("\t");
				out.print(s[2]);
				out.println();
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void printGeneModelsBed(UCSCGeneLine[] genes, File resultsUp, File resultsDown){
		try {
			Gzipper up = new Gzipper(resultsUp);
			Gzipper down = new Gzipper(resultsDown);
			boolean upFound = false;
			boolean downFound = false;

			//genome version
			up.println("# genome_version = "+genomeVersion);
			down.println("# genome_version = "+genomeVersion);
			//score names
			up.println("# score = abs(Log2Ratio) scaled 1-1000");
			up.println("# score = abs(Log2Ratio) scaled 1-1000");


			//for each gene fetch log2Ratio
			float[] ratios = new float[genes.length];
			for (int i=0; i< genes.length; i++){
				ratios[i] =genes[i].getScores()[log2ScoreIndex];
			}
			//scale
			Num.scale1To1000(ratios);
			//print
			for (int i=0; i< genes.length; i++){
				if (genes[i].getStrand().equals("-")) {
					down.println(genes[i].toStringBed12((int)ratios[i]));
					downFound = true;
				}
				else {
					up.println(genes[i].toStringBed12((int)ratios[i]));
					upFound = true;
				}
			}

			up.close();
			down.close();

			if (upFound == false) {
				resultsUp.delete();
				resultsUp = null;
			}
			if (downFound == false) {
				resultsDown.delete();
				resultsDown = null;
			}
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void printGeneModels(UCSCGeneLine[] genes, File results, boolean appendFailedGenes){
		try {
			Gzipper out = new Gzipper (results);

			//build header line
			String log2Ratio = "Log2Ratio";
			if (pseudoMedianCalculated && usePseudoMedianLog2Ratio) log2Ratio = "Log2((pseT+1)/(pseC+1))";
			String varOutPVals = "";
			if (varianceOutlierStatsPresent) varOutPVals = "NegBinomPVal_VarOutFilt\tBH_FDR_VarOutFilt\t";

			if (genes[0].getDisplayName() != null) out.print("#DisplayName\tName\t");
			else out.print("#Name\t");
			out.print("Chr\tStrand\tStart\tStop\tTotalRegionBPs\tNegBinomPVal_NoVarOutFilt\tBH_FDR_NoVarOutFilt\t"+ varOutPVals + log2Ratio +"\tAltSpliceChiSqrPVal\tAltSpliceMaxLog2Ratio\ttFPKM\tcFPKM");

			//for each replica
			for (int i=0; i< numberTreatmentReplicas; i++) out.print("\t#T"+(i+1));
			for (int i=0; i< numberControlReplicas; i++) out.print("\t#C"+(i+1));
			for (int i=0; i< numberTreatmentReplicas; i++) out.print("\tvarCorrT"+(i+1));
			for (int i=0; i< numberControlReplicas; i++) out.print("\tvarCorrC"+(i+1));
			out.print("\tGenomeVersion="+genomeVersion+ ", QueriedTreatObs="+(int)totalInterrogatedTreatmentObservations+ ", QueriedCtrlObs="+(int)totalInterrogatedControlObservations);

			out.println();

			//make hot link
			String gv = "";
			if (genomeVersion != null) gv = "version="+genomeVersion+"&";
			String url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?"+ gv +"seqid=";


			//for each gene that had counts
			for (int i=0; i< genes.length; i++){
				String name = genes[i].getDisplayNameThenName();
				//url
				int start = genes[i].getTxStart() - 10000;
				if (start < 0) start = 0;
				int end = genes[i].getTxEnd() + 10000;
				out.print(url+genes[i].getChrom()+"&start="+start+"&end="+end+"\",\""+name+"\")\t");
				//print second text?
				if (genes[i].getDisplayName() !=null) out.print(genes[i].getName()+"\t");
				//scores 
				float[] s = genes[i].getScores();
				out.println(genes[i].coordinates()+"\t"+Misc.floatArrayToString(s, "\t"));
			}

			//append on genes that missed minimum counts?
			if (appendFailedGenes){
				StringBuilder tabs = new StringBuilder();
				for (int x=0; x< this.totalNumberReplicas; x++) tabs.append("\t0");

				for (UCSCGeneLine g: geneLinesWithOutReadsAL){
					String name = g.getDisplayNameThenName();
					//url
					int start = g.getTxStart() - 10000;
					if (start < 0) start = 0;
					int end = g.getTxEnd() + 10000;
					out.print(url+g.getChrom()+"&start="+start+"&end="+end+"\",\""+name+"\")\t");
					//print second text?
					if (g.getDisplayName() !=null) out.print(g.getName()+"\t");
					//scores, just the counts for each of the reps
					out.println(g.coordinates() + "\t"+g.getTotalExonicBasePairs()+"\t0\t0\t0\t0\t0\t0\t0\t0\t0\t"+Misc.floatArrayToString(g.getScores(), "\t")+tabs);

				}
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void loadGeneModels() throws IOException{
		//load gene models from refFlat for refSeq UCSC gene table
		UCSCGeneModelTableReader reader = null;
		if (refSeqFile != null){
			reader = new UCSCGeneModelTableReader(refSeqFile, 0);
			if (scoreIntrons){
				if (verbose) System.out.println("\tParsing/scoring introns instead of exons from gene models.");
				reader.swapIntronsForExons();
			}
			if (removeOverlappingRegions) {
				if (verbose) System.out.print("\tRemoving overlapping regions from gene models");
				String deletedGenes = reader.removeOverlappingExons();
				if (deletedGenes.length() !=0) if (verbose) System.out.println("\t\tWARNING: the following genes had more than 1/2 of their exonic bps removed -> "+deletedGenes);
				//write out gene table sans overlaps
				File f = new File (saveDirectory, Misc.removeExtension(refSeqFile.getName())+"_NoOverlappingExons.txt") ;
				if (verbose) System.out.println("\t\tWrote modified gene table to the save directory.");
				reader.writeGeneTableToFile(f);
			}
			//any genes?
			reader.splitByChromosome();
			geneModels = reader.getChromSpecificGeneLines();
			allGeneLines = reader.getGeneLines();
		}
		//or from bed file
		else if (bedFile != null) {
			calculateChiSquare = false;
			Bed[] bed = Bed.parseFile(bedFile, 0, 0);
			allGeneLines = new UCSCGeneLine[bed.length];
			boolean addName = bed[0].getName().trim().equals("");
			for (int i=0; i< bed.length; i++){
				if (addName) {
					String name = (i+1)+"_"+bed[i].getChromosome()+":"+bed[i].getStart()+"-"+bed[i].getStop();
					bed[i].setName(name);
				}
				allGeneLines[i] = new UCSCGeneLine(bed[i]);
			}
			reader = new UCSCGeneModelTableReader();
			reader.setGeneLines(allGeneLines);
			reader.splitByChromosome();
			geneModels = reader.getChromSpecificGeneLines();
		}
		if (geneModels == null || allGeneLines == null || allGeneLines.length == 0) throw new IOException ("\nProblem loading your USCS gene model table or bed file? No genes/ regions found?\n");
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your regions's coordinates are reversed. Check that each start is less than the stop.\n");
		//check number of regions
		if (reader.getGeneLines().length < 100) {
			System.err.println("\nWarning! Too few regions to scan! DRSS needs at least 100 regions to properly estimate FDRs.\n");
		}
		//zero scores
		for (int i=0; i< allGeneLines.length; i++) allGeneLines[i].setScores(new float[]{0,0});
	}



	public ArrayList<String> fetchAlignmentNames (ExonIntron ei, SAMFileReader reader){
		ArrayList<String> al = new ArrayList<String>();
		SAMRecordIterator i = reader.queryOverlapping(chromosome, (ei.getStart()+1), ei.getEnd());
		if (dataIsStranded){
			while (i.hasNext()){
				SAMRecord sam = i.next();
				if (reverseStrand == sam.getReadNegativeStrandFlag()){
					//fetch blocks of actual alignment
					ArrayList<int[]> blocks = SamAlignment.fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);
					//check to see if any intersect the exon
					for (int[] b : blocks){
						if (ei.intersects(b[0], b[1])){
							al.add(sam.getReadName());
							break;
						}
					}
				}
			}
		}
		else {
			while (i.hasNext()) {
				SAMRecord sam = i.next();
				//fetch blocks of actual alignment
				ArrayList<int[]> blocks = SamAlignment.fetchAlignmentBlocks(sam.getCigarString(), sam.getUnclippedStart()-1);
				//check to see if any intersect the exon
				for (int[] b : blocks){
					//System.out.println(b[0]+" "+b[1]);
					if (ei.intersects(b[0], b[1])){
						al.add(sam.getReadName());
						break;
					}
				}
			}
		}
		i.close();
		i = null;
		return al;
	}

	public ArrayList<String>[] fetchOverlappingAlignmentNames (ExonIntron[] ei, SAMFileReader reader){
		ArrayList<String>[] als = new ArrayList[ei.length];
		for (int i=0; i< ei.length; i++) { 
			als[i] = fetchAlignmentNames( ei[i], reader);
		}
		return als;
	}

	public ArrayList<String>[][] fetchOverlappingAlignmentNames (ExonIntron[] ei, SAMFileReader[] readers){
		ArrayList<String>[][] ohMy = new ArrayList[readers.length][];
		for (int i=0; i< readers.length; i++) {
			ohMy[i] = fetchOverlappingAlignmentNames(ei, readers[i]);
		}
		return ohMy;
	}


	public float countUniqueAlignments (ArrayList<SAMRecord>[] alignments){
		HashSet<String> uniqueReadNames = new HashSet<String>();
		for (ArrayList<SAMRecord> al : alignments){
			for (SAMRecord sam: al) uniqueReadNames.add(new String (sam.getReadName()));
		}
		return uniqueReadNames.size();
	}

	public float countUniqueNames (ArrayList<String>[] alignments){
		HashSet<String> uniqueReadNames = new HashSet<String>();
		for (ArrayList<String> al : alignments){
			for (String sam: al) uniqueReadNames.add(sam);
		}
		return uniqueReadNames.size();
	}

	/**Collects number of observations under each gene's exons. Stranded if so designated.*/
	private void scanGenes(){
		//get UCSCGeneLine[]
		UCSCGeneLine[] genes = geneModels.get(chromosome);

		//for each gene 
		for (int i=0; i< genes.length; i++){
			//set strand
			reverseStrand = genes[i].getStrand().equals("-");

			//get exons
			ExonIntron[] exons = genes[i].getExons();

			//fetch the overlapping alignments [replicas][exons], switched to using names in an attempt to cut memory usage
			ArrayList<String>[][] tReps = fetchOverlappingAlignmentNames (exons, treatmentReaders);


			//counters
			int totalTreatment = 0;
			int totalControl = 0;

			//calculate number of unique alignments for each replica, put into the scores
			float[] scores = new float[totalNumberReplicas];
			for (int j=0; j< numberTreatmentReplicas; j++) {
				scores[j] = countUniqueNames(tReps[j]);
				totalTreatment+= scores[j];
			}

			ArrayList<String>[][] cReps = fetchOverlappingAlignmentNames (exons, controlReaders);
			int index = 0;
			for (int j=numberTreatmentReplicas; j< totalNumberReplicas; j++) {
				scores[j] = countUniqueNames(cReps[index++]);
				totalControl += scores[j];
			}

			genes[i].setScores(scores);

			//enough observations?
			if ((totalTreatment + totalControl) < minReads) {
				geneLinesWithOutReadsAL.add(genes[i]);
				continue;
			}

			//fetch counts for each exon, note some of these might be duplicates where one read spans several exons
			float[] tExonCountsSummed = new float[exons.length];
			float[][] tExonCounts = new float[numberTreatmentReplicas][exons.length];
			//for each replica
			for (int j=0; j< numberTreatmentReplicas; j++){
				//for each exon
				for (int k=0; k< exons.length; k++) {
					tExonCounts[j][k] = tReps[j][k].size();
					tExonCountsSummed[k] += tExonCounts[j][k];
				}
			}

			float[][] cExonCounts = new float[numberControlReplicas][exons.length];
			float[] cExonCountsSummed = new float[exons.length];
			//for each replica
			for (int j=0; j< numberControlReplicas; j++){
				//for each exon
				for (int k=0; k< exons.length; k++) {
					cExonCounts[j][k] = cReps[j][k].size();
					cExonCountsSummed[k] += cExonCounts[j][k];
				}
			}


			geneLinesWithReadsAL.add(genes[i]);

			//these will be for Ken's OD test of independence
			genes[i].setTreatmentExonCounts(tExonCounts);
			genes[i].setControlExonCounts(cExonCounts);

			//save exon counts for possible chi square test
			genes[i].setExonCounts(new float[][]{tExonCountsSummed, cExonCountsSummed});

			totalInterrogatedTreatmentObservations+= totalTreatment;
			totalInterrogatedControlObservations+= totalControl;

			tReps = null;
			cReps = null;
		}

		System.gc();

	}



	/**Looks for minimum number of reads and minimum log2Ratio difference between exon counts.*/
	private boolean checkForMinimums(UCSCGeneLine gene){
		float[] tExonCounts = gene.getExonCounts()[0];
		float[] cExonCounts = gene.getExonCounts()[1];
		//for each exon
		ArrayList<Integer> goodIndexes = new ArrayList<Integer>();
		float maxLogRatio = 0;
		int maxLogRatioIndex = -1;
		float[] log2Ratios = new float[tExonCounts.length];
		for (int i=0; i< tExonCounts.length; i++){
			//enough reads?
			if (tExonCounts[i] < 4 || cExonCounts[i] < 4) continue;
			//check ratio
			log2Ratios[i] = calculateLog2Ratio(tExonCounts[i], cExonCounts[i]);
			float logRatio = Math.abs(log2Ratios[i]);
			if (logRatio > maxLogRatio) {
				maxLogRatio = logRatio;
				maxLogRatioIndex = i;
			}
			goodIndexes.add(new Integer(i));
		}
		//check ratio and good exons, need at least two for alt splicing
		int numGoodExons = goodIndexes.size();
		if (maxLogRatio < minimumSpliceLog2Ratio || numGoodExons < 2) {
			gene.setExonCounts(null);
			return false;
		}
		//set log2Ratio 
		gene.setLog2Ratio(log2Ratios[maxLogRatioIndex]);
		//set exon counts
		if (numGoodExons == tExonCounts.length) return false;
		float[] tExonCountsSub = new float[numGoodExons];
		float[] cExonCountsSub = new float[numGoodExons];

		for (int i=0; i< numGoodExons; i++){
			int goodIndex = goodIndexes.get(i);
			tExonCountsSub[i] = tExonCounts[goodIndex];
			cExonCountsSub[i] = cExonCounts[goodIndex];
		}
		gene.setExonCounts(new float[][]{tExonCountsSub, cExonCountsSub});
		return true;
	}

	/**Collects and calculates a bunch of stats re the PointData.*/
	private void calculateReadCountStatistics(){

		millionMappedTreatmentReads = totalInterrogatedTreatmentObservations/1000000.0;
		if (verbose) System.out.println("\t"+(int)totalInterrogatedTreatmentObservations+" Queried Treatment Observations");
		for (File f : treatmentBamFiles) if (verbose) System.out.println("\t\t"+f.getName());

		if (verbose) System.out.println("\t"+(int)totalInterrogatedControlObservations+" Queried Control Observations");	
		for (File f : controlBamFiles) if (verbose) System.out.println("\t\t"+f.getName());
		millionMappedControlReads = totalInterrogatedControlObservations/1000000.0;

		if (totalInterrogatedTreatmentObservations < 1000 || totalInterrogatedControlObservations < 1000) Misc.printErrAndExit("\nError: too few observations?! Aborting.\n");

		//scalars for log2ratio
		scalarTC = totalInterrogatedTreatmentObservations/ totalInterrogatedControlObservations;
		scalarCT = totalInterrogatedControlObservations/totalInterrogatedTreatmentObservations;
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new OverdispersedRegionScanSeqs(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': treatmentBamFiles = IO.extractFiles(args[++i], ".bam"); break;
					case 'c': controlBamFiles = IO.extractFiles(args[++i], ".bam"); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'u': refSeqFile = new File(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'i': scoreIntrons = true; break;
					case 'p': usePseudoMedianLog2Ratio = true; break;
					case 'a': dataIsStranded = true; break;
					case 'v': genomeVersion = args[++i]; break;
					case 'f': minFDR = Float.parseFloat(args[++i]); break;
					case 'l': minLog2Ratio = Float.parseFloat(args[++i]); break;
					case 'e': minReads = Integer.parseInt(args[++i]); break;
					case 'x': calculateChiSquare = false; break;
					case 'o': removeOverlappingRegions = false; break;
					case 'd': deleteTempFiles = false; break;
					case 'z': printExonReplicaMatrix = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//genome version?
		if (genomeVersion == null) Misc.printErrAndExit("\nPlease provide a versioned genome (e.g. H_sapiens_Mar_2006).\n");

		//look for bam files
		if (treatmentBamFiles == null || treatmentBamFiles.length == 0) Misc.printErrAndExit("\nError: cannot find any treatment xxx.bam files?\n");
		//control data
		if (controlBamFiles == null || controlBamFiles.length == 0) Misc.printErrAndExit("\nError: cannot find any control xxx.bam files?\n");

		//look for bai index files
		lookForBaiIndexes(treatmentBamFiles, false);
		lookForBaiIndexes(controlBamFiles, false);

		//set number of replicas
		numberTreatmentReplicas = treatmentBamFiles.length;
		numberControlReplicas = controlBamFiles.length;
		totalNumberReplicas = numberTreatmentReplicas + numberControlReplicas;
		if (totalNumberReplicas < 3 && verbose) System.out.println("WARNING: DESeq tends to be very conservative when analyzing data with < 3 replicas.  Ideally you should have 3-4 treatment and 3-4 control biological replicas.\n");


		//look for and or create the save directory
		if (saveDirectory == null) Misc.printErrAndExit("\nError: enter a directory text to save results.\n");
		saveDirectory.mkdir();

		//check for R and required libraries
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printErrAndExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}
		else {
			String errors = IO.runRCommandLookForError("library(DESeq)", fullPathToR, saveDirectory);
			if (errors == null || errors.length() !=0){
				Misc.printErrAndExit("\nError: Cannot find the required R library.  Did you install DESeq " +
						"(http://www-huber.embl.de/users/anders/DESeq/)?  See the author's websites for installation instructions. Once installed, " +
						"launch an R terminal and type 'library(DESeq)' to see if it is present. R error message:\n\t\t"+errors+"\n\n");
			}
			//look for estimateDispersions() function
			boolean useEstimateDispersions = DefinedRegionDifferentialSeq.estimateDispersions(fullPathToR, saveDirectory);
			if (useEstimateDispersions == false){
				Misc.printErrAndExit("\nError: Please upgrade DESeq to the latest version, see http://www-huber.embl.de/users/anders/DESeq/ \n");
			}
		}



		//look for bed file
		if (refSeqFile == null && bedFile == null){
			Misc.printErrAndExit("\nPlease enter a regions file to use in scoring regions.\n");
		}
	}	

	/**Looks for xxx.bam.bai and xxx.bai for each bamFile, prints error and exits if missing.*/
	public static void lookForBaiIndexes (File[] bamFiles, boolean onlyResetLastModifiedDate){
		for (File f: bamFiles){
			File index = new File (f+".bai");
			if (index.exists() == false){
				int len = f.toString().length() - 3;
				index = new File(f.toString().substring(0, len) + "bai");
				if (onlyResetLastModifiedDate == false && index.exists() == false) Misc.printErrAndExit("\nError: failed to find a xxx.bai index file for -> "+f);
			}
			//reset date?
			if (index.exists() && index.lastModified() < f.lastModified()) index.setLastModified(f.lastModified()+1);
		}
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Overdispersed Region Scan Seqs: May 2012                  **\n" +
				"**************************************************************************************\n" +
				
				"\nWARNING: this application is depreciated and no longer maintained, use the\n"+
				"DefinedRegionDifferentialSeq app instead!\n\n"+
				
				"ORSS takes bam alignment files and extracts reads under each region or gene's exons to\n" +
				"calculate several statistics. Makes use of Simon Anders' DESeq R package to with its\n" +
				"negative binomial p-value test to control for overdispersion. A Benjamini-Hochberg FDR\n" +
				"correction is used to control for multiple testing. DESeq is run with and without\n" +
				"variance outlier filtering. A chi-square test of independence between the exon read\n" +
				"count distributions is used to score alternative splicing. Several read count measures\n" +
				"are provided including counts for each replica, FPKMs (# frags per kb of int region\n" +
				"per total mill mapped reads) as well as DESeq's variance adjusted counts(use these for\n" +
				"clustering, correlation, and other distance type analysis). If replicas are provided\n" +
				"either the smallest all pair log2Ratio is reported (default) or the pseudomedian.\n" +
				"Several results files are written: two spread sheets containing all of the genes,\n" +
				"those that pass the thresholds, as well as egr, bed12, and useq region files for\n" +
				"visualization in genome browsers.\n\n"+

				"Required Options:\n"+
				"-s Save directory.\n"+
				"-t Treatment directory containing one xxx.bam file with xxx.bai index per biological\n" +
				"       replica. The BAM files should be sorted by coordinate and have passed Picard\n" +
				"       validation. Use the SamTranscriptomeParser to convert your aligned transcriptome\n" +
				"       data to genomic coordinates.\n" +
				"-c Control directory, ditto. \n" +
				"-u UCSC RefFlat or RefSeq Gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (name1 name2(optional) chrom strand\n" +
				"       txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds). WARNING!!!!!!\n"+
				"       This table should contain only one composite transcript per gene. Use the\n" +
				"       MergeUCSCGeneTable app to collapse Ensembl transcripts downloaded from UCSC in\n" +
				"       RefFlat format.\n"+
				"-b (Or) a bed file (chr, start, stop,...), full path, See,\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1\n"+
				"-v Versioned Genome (ie H_sapiens_Mar_2006, D_rerio_Jul_2010), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases. \n" +


				"\nAdvanced/ Default Options:\n"+
				"-o Don't remove overlapping exons, defaults to filtering gene annotation for overlaps.\n"+
				"-i Score introns instead of exons.\n"+
				"-a Data is stranded. Only collect reads from the same strand as the annotation.\n"+
				"-f Minimum FDR threshold, defaults to 10 (-10Log10(FDR=0.1))\n"+
				"-l Minimum absolute log2 ratio threshold, defaults to 1 (2x)\n"+
				"-e Minimum number mapping reads per region, defaults to 20\n"+
				"-d Don't delete temp files used by DESeq\n"+
				"-p Use a pseudo median log2 ratio in place of the smallest all pair log2 ratios for\n"+
				"      scoring the degree of differential expression when replicas are present.\n" +
				"      Recommended for experiments with 4 or more replicas.\n"+
				"-r Full path to R loaded with DESeq library, defaults to '/usr/bin/R' file, see\n" +
				"       http://www-huber.embl.de/users/anders/DESeq/ . Type 'library(DESeq)' in\n" +
				"       an R terminal to see if it is installed.\n"+


				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/OverdispersedRegionScanSeqs -t\n" +
				"      /Data/PolIIRep1/,/Data/PolIIRep2/ -c /Data/Input1/,Data/Input2/ -s\n" +
				"      /Data/PolIIResults/ -f 30 -e 30 -u /Anno/mergedZv9EnsemblGenes.ucsc.gz\n\n" +

		"**************************************************************************************\n");

	}

	public double getTotalInterrogatedTreatmentObservations() {
		return totalInterrogatedTreatmentObservations;
	}

	public double getTotalInterrogatedControlObservations() {
		return totalInterrogatedControlObservations;
	}

	public int getNumberGenesPassingThresholds() {
		return numberGenesPassingThresholds;
	}
}
