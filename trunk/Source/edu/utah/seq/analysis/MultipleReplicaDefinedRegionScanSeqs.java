package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.bio.annotation.Bed;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.gen.*;
import edu.utah.seq.data.*;
import util.bio.cluster.*;


/**
 * @author Nix
 * */
public class MultipleReplicaDefinedRegionScanSeqs {

	//user defined fields
	private File[] treatmentPointDirs;
	private File[] controlPointDirs;
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private File bedFile;
	private File refSeqFile;
	private int peakShift = 0;
	private boolean dataIsStranded = false;
	private int minimumExonicReadsForChiSquare = 10;
	private float minimumSpliceLog2Ratio = 1;
	private float minFDR = 10;
	private float minLog2Ratio = 0.5849625f;
	private boolean calculateChiSquare = true;
	private boolean cluster = true;
	private boolean scoreIntrons = false;
	private boolean removeOverlappingRegions = true;
	private boolean printExonReplicaMatrix = false;

	//internal fields
	private HashMap<String,UCSCGeneLine[]> geneModels;
	private UCSCGeneLine[] allGeneLines;
	private ArrayList<UCSCGeneLine> geneLinesWithReadsAL = new ArrayList<UCSCGeneLine>();
	private UCSCGeneLine[] geneLinesWithReads;
	private int halfPeakShift = 0;
	private int totalInterrogatedTreatmentObservations = 0;
	private int totalInterrogatedControlObservations = 0;
	private double numberTreatmentObservations;
	private double numberControlObservations;
	private double millionMappedTreatmentReads;
	private double millionMappedControlReads;
	private int numberTreatmentReplicas;
	private int numberControlReplicas;
	private int totalNumberReplicas;
	private int[] numberTreatmentReads;
	private int[] numberControlReads;
	private HashMap<String,PointData[]> treatmentPlusPointData;
	private HashMap<String,PointData[]> treatmentMinusPointData;
	private HashMap<String,PointData[]> controlPlusPointData;
	private HashMap<String,PointData[]> controlMinusPointData;
	private String genomeVersion;
	private boolean deleteMatrixFile = false;
	private double scalarTC;
	private double scalarCT;

	
	//by chromosome
	private String chromosome;
	private PointData[] treatmentChromPlus = null;
	private PointData[] treatmentChromMinus = null;	
	private PointData[] controlChromPlus = null;
	private PointData[] controlChromMinus = null;

	//constructor
	/**For integration with RNASeq app.*/
	public MultipleReplicaDefinedRegionScanSeqs(File[] treatmentPointDirs, File[] controlPointDirs, File saveDirectory, File fullPathToR, File refSeqFile, boolean dataIsStranded, boolean scoreIntrons){
		this.treatmentPointDirs = treatmentPointDirs;
		this.controlPointDirs = controlPointDirs;
		this.saveDirectory = saveDirectory;
		this.fullPathToR = fullPathToR;
		this.refSeqFile = refSeqFile;
		this.dataIsStranded = dataIsStranded;
		this.scoreIntrons = scoreIntrons;
		removeOverlappingRegions = false;
		run();
	}
	
	/**Stand alone.*/
	public MultipleReplicaDefinedRegionScanSeqs(String[] args){
		long startTime = System.currentTimeMillis();
		//set fields
		processArgs(args);
		
		System.err.println("\n\nWARNING: This application is depreciated and no longer maintained.  Use the OverdispersedRegionScanSeqs.\n\n");
		
		//launch
		run();
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	public void run(){
		//load gene models
		System.out.println("Loading regions/ gene models...");
		loadGeneModels();

		//fetch counts
		System.out.println("Calculating read count stats...");
		calculateReadCountStatistics();

		//for each chromosome of gene models
		System.out.print("Scanning regions by chromosome");
		Iterator<String> it = geneModels.keySet().iterator();
		while (it.hasNext()){
			chromosome = it.next();
			scanChromosome();
		}
		System.out.println();

		//load geneLines with those that have observations
		geneLinesWithReads = new UCSCGeneLine[geneLinesWithReadsAL.size()];
		geneLinesWithReadsAL.toArray(geneLinesWithReads);

		//convert binomial pvalues to FDRs
		System.out.println("Calculating negative binomial p-values and FDRs in R using DESeq (http://www-huber.embl.de/users/anders/DESeq/)...");
		calculateDESeqPValues();

		//calculate chiSquare test of independence
		if (calculateChiSquare){
			System.out.println("Estimating alternative splicing read distribution, bonferroni corrected, chi-square p-values in R (very slow)...");
			if (printExonReplicaMatrix) estimateDifferencesInReadDistributionsWithReplicas();
			else estimateDifferencesInReadDistributions();
		}
		else System.out.println("Skipping the alternative splicing read distribution calculation ...");
		
		File clusterFile =null;
		if (cluster){
			System.out.println("Hierarchical clustering replicas based on DESeq's variance corrected replica read counts...");
			clusterFile = new File (saveDirectory,"clusterPlot.png");
			cluster(clusterFile);
		}

		//sort by score[0] pValue
		System.out.println("Thresholding and printing results...");
		Arrays.sort(allGeneLines, new UCSCGeneLineComparatorScoreBigToSmall(1));

		//print all gene report
		System.out.println("\tSave dir\t"+saveDirectory);
		if (clusterFile!= null) System.out.println("\t"+clusterFile.getName()+"\t: replica cluster plot");
		File allGenesFile = new File (saveDirectory, "all.xls");
		System.out.println("\t"+allGenesFile.getName()+"\t: all gene/region spreadsheet");
		printGeneModels(allGeneLines, allGenesFile);
		
		//threshold genes
		UCSCGeneLine[] good = thresholdGenes();
		
		if (good.length != 0){
			Arrays.sort(good, new UCSCGeneLineComparatorScoreBigToSmall(1));
			File partialGenesFile = new File (saveDirectory, "diffExprFDR"+minFDR+"LgRto"+Num.formatNumber(minLog2Ratio,3)+".xls");
			System.out.println("\t"+partialGenesFile.getName()+"\t: thresholded gene/region spreadsheet");
			printGeneModels(good, partialGenesFile);
			
			//print egr file 
			File egrFile = new File (saveDirectory, "diffExprFDR"+minFDR+"LgRto"+Num.formatNumber(minLog2Ratio,3)+".egr");
			System.out.println("\t"+egrFile.getName()+"\t: thresholded gene/region graph file");
			printGeneModelsEgr(good, egrFile);
			
		}
		System.out.println("\n"+good.length+" genes/regions passed thresholds.");
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
	
	public void cluster(File clusterFile){
		try {
			//write single float[] files to the save directory
			File[] files = new File[totalNumberReplicas];
			Pattern slash = Pattern.compile("/");
			for (int i=0; i< numberTreatmentReplicas; i++){
				String fullPathName = treatmentPointDirs[i].getCanonicalPath();
				fullPathName = slash.matcher(fullPathName).replaceAll("_");
				files[i] = new File(saveDirectory, fullPathName);
			}
			int index = 0;
			for (int i=numberTreatmentReplicas; i< totalNumberReplicas ; i++){
				String fullPathName = controlPointDirs[index++].getCanonicalPath();
				fullPathName = slash.matcher(fullPathName).replaceAll("_");
				files[i] = new File(saveDirectory, fullPathName);
			}
			//for each gene with reads
			float[][] varCounts = new float[totalNumberReplicas][geneLinesWithReads.length];
			int adder = 7+totalNumberReplicas;
			for (int i=0; i< geneLinesWithReads.length; i++){
				//scores = totalExonicBPs, pval, padj, log2ratio, chiSquare, tRPKM, cRPKM, countsT1,2,3, countsC1,2,3, varCorrT1,2,3,  varCorrC1,2,3
				float[] scores = geneLinesWithReads[i].getScores();
				for (int x=0; x< totalNumberReplicas; x++) varCounts[x][i]= scores[x+adder];
			}
			//save float[] objects
			for (int x=0; x< totalNumberReplicas; x++) IO.saveObject(files[x], varCounts[x]);
			
			//cluster
			HierarchicalClustering hc = new HierarchicalClustering(files, clusterFile);
			
			//delete
			for (int x=0; x< totalNumberReplicas; x++) files[x].delete();
			
		} catch (IOException e) {
			Misc.printErrAndExit("\nProblem clustering data.");
			e.printStackTrace();
		}
		
	}
	
	public UCSCGeneLine[] thresholdGenes(){
		ArrayList<UCSCGeneLine> passing = new ArrayList<UCSCGeneLine>();
		for (int i=0; i< geneLinesWithReads.length; i++){
			//scores = totalExonicBPs, pval, padj, log2ratio, chiSquare, tRPKM, cRPKM, countsT1,2,3, countsC1,2,3, varCorrT1,2,3,  varCorrC1,2,3
			float[] scores = geneLinesWithReads[i].getScores();
			if (scores[2]>= minFDR && Math.abs(scores[3])>=minLog2Ratio) passing.add(geneLinesWithReads[i]);
		}
		UCSCGeneLine[] good = new UCSCGeneLine[passing.size()];
		passing.toArray(good);
		return good;
	}
	
	public File writeOutObservations(String randomWord){
		File matrixFile = null;
		try {
			//write matrix of name, t1,t2,t3...c1,c2,c3 to file for genes with observations
			matrixFile = new File(saveDirectory, randomWord+"_CountMatrix.txt");
			PrintWriter out = new PrintWriter( new FileWriter(matrixFile));
			for (int i=0; i< geneLinesWithReads.length; i++){
				float[] scores = geneLinesWithReads[i].getScores();
				//out.print(geneLinesWithReads[i].getName());
				out.print("G"+i);
				for (int j=0; j< scores.length; j++){
					out.print("\t"+(int)scores[j]);
				}
				out.println();
			}
			out.close();
		}
		catch (Exception e){
			System.err.println("Problem writing out observations for R.");
			e.printStackTrace();
		}
		return matrixFile;
	}

	/**Returns the deseq stats file and the variance stabilized data.*/
	private File[] executeDESeq(File matrixFile, String randomWord){
		File rResultsStats = new File (saveDirectory, randomWord+"_DESeqResults.txt");
		File rResultsData = new File (saveDirectory, randomWord+"_DESeqResultsData.txt");
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("library(DESeq)\n");
			sb.append("countsTable = read.delim('"+matrixFile.getCanonicalPath()+"', header=FALSE)\n");
			sb.append("rownames(countsTable) = countsTable$V1\n");
			sb.append("countsTable = countsTable[,-1]\n");
			sb.append("conds = c(");
			for (int i=0; i< numberTreatmentReplicas; i++){
				sb.append("'T',");
			}
			sb.append("'N'");
			for (int i=1; i< numberControlReplicas; i++){
				sb.append(",'N'");
			}
			sb.append(")\n");
			sb.append("cds = newCountDataSet( countsTable, conds)\n");
			sb.append("cds = estimateSizeFactors( cds )\n");
			if (totalNumberReplicas == 2) sb.append("cds = estimateVarianceFunctions( cds, pool=TRUE )\n");
			else sb.append("cds = estimateVarianceFunctions( cds )\n");
			sb.append("res = nbinomTest( cds, 'N', 'T', pvals_only = FALSE)\n");
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
			sb.append("write.table(res, file = '"+rResultsData.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");

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
			if (rResultsStats.exists() == false || rResultsData.exists() == false) throw new IOException("\nR results file doesn't exist. Check temp files in save directory for error.\n");

			//cleanup
			rOut.deleteOnExit();
			scriptFile.deleteOnExit();

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
			return null;
		}
		return new File[]{rResultsStats,rResultsData};
	}

	private void parseDESeqStatResults(File results){
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
				//get existing scores: read counts in T's and C's
				float[] scores = geneLinesWithReads[x].getScores();
				//build: totalExonicBPs, pval, padj, log2ratio, chiSquare, tRPKM, cRPKM, countsT1,2,3, countsC1,2,3, varCorrT1,2,3,  varCorrC1,2,3
				float[] allScores = new float[7+ (totalNumberReplicas*2)];
				//totalExonicBPs
				allScores[0] = geneLinesWithReads[x].getTotalExonicBasePairs();
				//pval
				allScores[1] = parseFloat(tokens[0]);
				if (allScores[1]> maxPVal) maxPVal = allScores[1];
				//adjPval
				allScores[2] = parseFloat(tokens[1]);
				if (allScores[2]> maxAdjPVal) maxAdjPVal = allScores[2];
				//log2ratio
				allScores[3] = Float.parseFloat(tokens[2]);
				//chisquare, wait for block calculation
				//allScores[4]
				//tRPKM
				int totalT = 0;
				for (int i=0; i<numberTreatmentReplicas; i++) totalT+= scores[i]; 
				allScores[5] = this.calculateRPKM(millionMappedTreatmentReads, allScores[0], totalT);
				//cRPKM
				int totalC = 0;
				for (int i=numberTreatmentReplicas; i<totalNumberReplicas; i++) totalC+= scores[i]; 
				allScores[6] = this.calculateRPKM(millionMappedControlReads, allScores[0], totalC);
				//add counts
				System.arraycopy(scores, 0, allScores, 7, totalNumberReplicas);
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
			//clean up
			results.delete();
		} catch (Exception e){
			System.err.println("Problem parsing DESeq stats results from R.\n");
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
			results.delete();
		} catch (Exception e){
			System.err.println("Problem parsing DESeq data results from R.\n");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public float parseFloat(String f){
		if (f.equals("Inf") || f.equals("NA")) return Float.MIN_VALUE;
		else return Float.parseFloat(f);
	}

	public void calculateDESeqPValues(){
		//write matrix of name, t1,t2,t3...c1,c2,c3 to file for genes with observations
		String random = Passwords.createRandowWord(6);
		File matrixFile = writeOutObservations(random);
		//call DESeq
		File[] deseqResults = executeDESeq(matrixFile, random);
		//parse results and populate gene scores
		parseDESeqStatResults (deseqResults[0]);
		parseDESeqDataResults (deseqResults[1]);
		//clean up
		if (deleteMatrixFile) matrixFile.delete();
	}

	/**Calculates the reads per kb per million mapped reads 
	 * # Observed reads in the region/ bp size of the region / 1000/ total number reads/ 1000000 */
	public float calculateRPKM(double millionTotalMappedReads, double interrogatedRegionBPSize, double numberObservedReadsInRegion){
		double exonicBasesPerKB = interrogatedRegionBPSize/1000;
		double rpkm = numberObservedReadsInRegion/exonicBasesPerKB/millionTotalMappedReads;
		return new Double(rpkm).floatValue();
	}

	public void estimateDifferencesInReadDistributions(){
		//find genes with multiple exons and max number exons
		int maxNumberExons = -1;
		ArrayList<UCSCGeneLine> al = new ArrayList<UCSCGeneLine>();
		for (int i=0; i< geneLinesWithReads.length; i++){
			int numEx = geneLinesWithReads[i].getExons().length;
			if (numEx > 1 && geneLinesWithReads[i].getExonCounts() != null) {
				al.add(geneLinesWithReads[i]);
				if (numEx > maxNumberExons) maxNumberExons = numEx;
			}
		}
		//anything pass?
		if (al.size() ==0){
			System.out.println("\nWARNING: no genes with introns or minimal number reads found, skipping alternative splice detection.\n");
			return;
		}
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
			float[][] tc = genesWithExonsAndReads[i].getExonCounts();
			//set corrected p-value 
			//scores = totalExonicBPs, pval, padj, log2ratio, chiSquare, tRPKM, cRPKM, countsT1,2,3, countsC1,2,3 varT1,2,3  varC1,2,3
			float[] scores = genesWithExonsAndReads[i].getScores();
			scores[4] = (float) pVals[i] + bc;
			if (scores[4] < 0) scores[4] = 0;
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
		
		//anything pass?
		if (al.size() ==0){
			System.out.println("\nWARNING: no genes with introns or minimal number reads found, skipping alternative splice detection.\n");
			return;
		}
		
		UCSCGeneLine[] genesWithExonsAndReads = new UCSCGeneLine[al.size()];
		al.toArray(genesWithExonsAndReads);
		
		//just print them to file
		File matrix = new File (this.saveDirectory, "exonMatrixFile.txt");
		PrintWriter out;
		try {
			out = new PrintWriter( new FileWriter (matrix));
			out.println("# Number genes with >1 exons and counts\t"+genesWithExonsAndReads.length);
			out.println("# Number treatment replicas\t"+this.numberTreatmentReplicas);
			out.println("# Number control replicas\t"+this.numberControlReplicas);
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
		
		System.out.println("Matrix printed exiting...... ");
		System.exit(0);
		
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
			float[][] tc = genesWithExonsAndReads[i].getExonCounts();
			//set corrected p-value 
			//scores = totalExonicBPs, pval, padj, log2ratio, chiSquare, tRPKM, cRPKM, countsT1,2,3, countsC1,2,3 varT1,2,3  varC1,2,3
			float[] scores = genesWithExonsAndReads[i].getScores();
			scores[4] = (float) pVals[i] + bc;
			if (scores[4] < 0) scores[4] = 0;
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
			PrintWriter out = new PrintWriter (new FileWriter( results));
			//genome version
			out.println("# genome_version = "+genomeVersion);
			//score names
			out.println("# score0 = Log2Ratio");
			out.println("# score1 = -10Log10(FDR)");
			//for each gene
			for (int i=0; i< genes.length; i++){
				out.print(genes[i].getName());
				out.print("\t");
				out.print(genes[i].getChrom());
				out.print("\t");
				out.print(genes[i].getTxStart());
				out.print("\t");
				out.print(genes[i].getTxEnd());
				out.print("\t");
				out.print(genes[i].getStrand());
				out.print("\t");
				//scores = totalExonicBPs, pval, padj, log2ratio, chiSquare, tRPKM, cRPKM, countsT1,2,3, countsC1,2,3 varT1,2,3  varC1,2,3 
				float[] s = genes[i].getScores();
				out.print(s[3]);
				out.print("\t");
				out.print(s[2]);
				out.println();
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void printGeneModels(UCSCGeneLine[] genes, File results){
		try {
			PrintWriter out = new PrintWriter (new FileWriter( results));
			
			//build header line
			if (genes[0].getDisplayName() != null) out.print("#DisplayName\tName\t");
			else out.print("#Name\t");
			out.print("Chr\tStrand\tStart\tStop\tTotalRegionBPs\tNegBinomPVal\tBH_FDR\tLog2((sumT+1)/(sumC+1))\tChiSqrPValDiffDist\ttRPKM\tcRPKM");
			
			//for each replica
			for (int i=0; i< numberTreatmentReplicas; i++) out.print("\t#T"+(i+1));
			for (int i=0; i< numberControlReplicas; i++) out.print("\t#C"+(i+1));
			for (int i=0; i< numberTreatmentReplicas; i++) out.print("\tvarCorrT"+(i+1));
			for (int i=0; i< numberControlReplicas; i++) out.print("\tvarCorrC"+(i+1));
			out.println("\tGenomeVersion="+genomeVersion+ ", TotalTreatObs="+(int)numberTreatmentObservations+ ", TotalCtrlObs="+(int)numberControlObservations);
			
			//make hot link
			String url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";
			
			//for each gene
			for (int i=0; i< genes.length; i++){
				String name;
				if (genes[i].getDisplayName() !=null) name = genes[i].getDisplayName();
				else name = genes[i].getName();
				//url
				int start = genes[i].getTxStart() - 10000;
				if (start < 0) start = 0;
				int end = genes[i].getTxEnd() + 10000;
				out.print(url+genes[i].getChrom()+"&start="+start+"&end="+end+"\",\""+name+"\")\t");
				//print second text?
				if (genes[i].getDisplayName() !=null) out.print(genes[i].getName()+"\t");
				//scores 
				float[] s = genes[i].getScores();
				if (s.length == 2) out.print(genes[i].coordinates());
				else out.print(genes[i].coordinates()+"\t"+Misc.floatArrayToString(s, "\t"));
				out.println();
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void loadGeneModels(){
		//load gene models from refFlat for refSeq UCSC gene table
		UCSCGeneModelTableReader reader = null;
		if (refSeqFile != null){
			reader = new UCSCGeneModelTableReader(refSeqFile, 0);
			if (scoreIntrons){
				System.out.println("\tParsing/scoring introns instead of exons from gene models.");
				reader.swapIntronsForExons();
			}
			if (removeOverlappingRegions) {
				System.out.println("\tRemoving overlapping regions from gene models.\n");
				String deletedGenes = reader.removeOverlappingExons();
				if (deletedGenes.length() !=0) System.out.println("\t\tWARNING: the following genes had more than 1/2 of their exonic bps removed -> "+deletedGenes);
			}
			System.out.println();
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
				if (addName) bed[i].setName((i+1)+"");
				allGeneLines[i] = new UCSCGeneLine(bed[i]);
			}
			reader = new UCSCGeneModelTableReader();
			reader.setGeneLines(allGeneLines);
			reader.splitByChromosome();
			geneModels = reader.getChromSpecificGeneLines();
		}
		if (geneModels == null || allGeneLines == null || allGeneLines.length == 0) Misc.printExit("\nProblem loading your USCS gene model table or bed file? No genes/ regions?\n");
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your regions's coordinates are reversed. Check that each start is less than the stop.\n");
		//check number of regions
		if (reader.getGeneLines().length < 100) {
			System.err.println("\nWarning! Too few regions to scan! MRDRSS needs at least 100 regions to properly estimate FDRs.\n");
		}
		//zero scores
		for (int i=0; i< allGeneLines.length; i++) allGeneLines[i].setScores(new float[]{0,0});
	}

	/**Fetchs the data for a particular chromosome.*/
	public boolean fetchData(){
		//fetch treatments
		treatmentChromPlus = treatmentPlusPointData.get(chromosome);
		treatmentChromMinus = treatmentMinusPointData.get(chromosome);
		//fetch controls
		controlChromPlus = controlPlusPointData.get(chromosome);
		controlChromMinus = controlMinusPointData.get(chromosome);
		//check for nulls
		if (treatmentChromPlus == null || treatmentChromMinus == null || controlChromPlus == null || controlChromMinus == null){
			System.err.println("\n\tWarning: failed to find stranded data from all datasets for "+chromosome+". Skipping!");
			return false;
		}
		
		//check number
		if (treatmentChromPlus.length != numberTreatmentReplicas || treatmentChromMinus.length != numberTreatmentReplicas || controlChromPlus.length != numberControlReplicas || controlChromMinus.length != numberControlReplicas){
			System.err.println("\n\tWarning: The number of treatment or control replicas differ for chromosome -> "+chromosome+". Skipping!");
			return false;
		}
		
		
		//increment number of obs
		for (int i=0; i< numberTreatmentReplicas; i++){
			numberTreatmentReads[i] += treatmentChromPlus[i].getInfo().getNumberObservations();
			numberTreatmentReads[i] += treatmentChromMinus[i].getInfo().getNumberObservations();
		}
		for (int i=0; i< numberControlReplicas; i++){
			numberControlReads[i] += controlChromPlus[i].getInfo().getNumberObservations();
			numberControlReads[i] += controlChromMinus[i].getInfo().getNumberObservations();
		}		//assign genome
		if (genomeVersion == null){
			if (treatmentChromPlus != null) genomeVersion = treatmentChromPlus[0].getInfo().getVersionedGenome();
			else if (treatmentChromMinus != null) genomeVersion = treatmentChromMinus[0].getInfo().getVersionedGenome();
		}
		return true;
	}

	/**Scans a chromosome collecting read count data and then collects the associated scores.*/
	public void scanChromosome(){
		System.out.print(".");
		//fetch data
		if (fetchData() == false) return;
		//fetch, shift positions from the treatment and control
		fetchShiftStripPointData();
		//scan
		if (printExonReplicaMatrix) scanGenesWithExonReplicas();
		else scanGenes();
	}

	/**Collects number of observations under each gene's exons. Stranded if so designated.*/
	private void scanGenes(){
		//get UCSCGeneLine[]
		UCSCGeneLine[] genes = geneModels.get(chromosome);

		//for each gene 
		for (int i=0; i< genes.length; i++){

			//get exons
			ExonIntron[] exons = genes[i].getExons();

			//fetch scores for each exon and some summary stats
			float[] scores = new float[totalNumberReplicas];
			float[] tExonCounts = new float[exons.length];
			float[] cExonCounts = new float[exons.length];
			int totalCounts = 0;

			//is data stranded?
			if (dataIsStranded){
				String strand = genes[i].getStrand();
				//plus stranded annotation
				if (strand.equals("+")){
					for (int k=0; k< exons.length; k++){
						int start = exons[k].getStart();
						int stop = exons[k].getEnd();
						for (int j=0; j< numberTreatmentReplicas; j++){
							float plus = treatmentChromPlus[j].sumScoreBP(start, stop);
							scores[j] += plus;
							totalCounts += plus;
							tExonCounts[k] += plus;
						}
						for (int j=0; j< numberControlReplicas; j++){
							float plus = controlChromPlus[j].sumScoreBP(start, stop);
							scores[j+numberTreatmentReplicas] += plus;
							totalCounts += plus;
							cExonCounts[k] += plus;
						}
					}
				}
				//minus stranded
				else {
					for (int k=0; k< exons.length; k++){
						int start = exons[k].getStart();
						int stop = exons[k].getEnd();
						for (int j=0; j< numberTreatmentReplicas; j++){
							float minus = treatmentChromMinus[j].sumScoreBP(start, stop);
							scores[j] += minus;
							totalCounts += minus;
							tExonCounts[k] += minus;
						}
						for (int j=0; j< numberControlReplicas; j++){
							float minus = controlChromMinus[j].sumScoreBP(start, stop);
							scores[j+numberTreatmentReplicas] += minus;
							totalCounts += minus;
							cExonCounts[k] += minus;
						}
					}
				}
			}
			//nope combine strand counts
			else {
				for (int k=0; k< exons.length; k++){
					int start = exons[k].getStart();
					int stop = exons[k].getEnd();
					for (int j=0; j< numberTreatmentReplicas; j++){
						float plus = treatmentChromPlus[j].sumScoreBP(start, stop);
						float minus = treatmentChromMinus[j].sumScoreBP(start, stop);
						float total = plus+minus;
						scores[j] += total;
						totalCounts += total;
						tExonCounts[k] += total;
					}
					for (int j=0; j< numberControlReplicas; j++){
						float plus = controlChromPlus[j].sumScoreBP(start, stop);
						float minus = controlChromMinus[j].sumScoreBP(start, stop);
						float total = plus+minus;
						scores[j+numberTreatmentReplicas] += total;
						totalCounts+= total;
						cExonCounts[k] += total;
					}
				}
			}
			//save scores ?
			if (totalCounts==0) continue;
			genes[i].setScores(scores);
			geneLinesWithReadsAL.add(genes[i]);
			totalInterrogatedTreatmentObservations += Num.sumArray(tExonCounts);
			totalInterrogatedControlObservations += Num.sumArray(cExonCounts);
			//save exon counts? limiter needed since chi square test so slow
			if (checkForMinimums(tExonCounts, cExonCounts)) genes[i].setExonCounts(new float[][]{tExonCounts, cExonCounts});
		}

	}
	
	
	/**Collects number of observations under each gene's exons. Stranded if so designated.*/
	private void scanGenesWithExonReplicas(){
		//get UCSCGeneLine[]
		UCSCGeneLine[] genes = geneModels.get(chromosome);

		//for each gene 
		for (int i=0; i< genes.length; i++){

			//get exons
			ExonIntron[] exons = genes[i].getExons();

			//fetch scores for each exon and some summary stats
			float[] scores = new float[totalNumberReplicas];
			float[][] tExonCounts = new float[numberTreatmentReplicas][exons.length];
			float[][] cExonCounts = new float[numberControlReplicas][exons.length];
			int totalCounts = 0;

			//is data stranded?
			if (dataIsStranded){
				String strand = genes[i].getStrand();
				//plus stranded annotation
				if (strand.equals("+")){
					for (int k=0; k< exons.length; k++){
						int start = exons[k].getStart();
						int stop = exons[k].getEnd();
						for (int j=0; j< numberTreatmentReplicas; j++){
							float plus = treatmentChromPlus[j].sumScoreBP(start, stop);
							scores[j] += plus;
							totalCounts += plus;
							tExonCounts[j][k] = plus;
						}
						for (int j=0; j< numberControlReplicas; j++){
							float plus = controlChromPlus[j].sumScoreBP(start, stop);
							scores[j+numberTreatmentReplicas] += plus;
							totalCounts += plus;
							cExonCounts[j][k] = plus;
						}
					}
				}
				//minus stranded
				else {
					for (int k=0; k< exons.length; k++){
						int start = exons[k].getStart();
						int stop = exons[k].getEnd();
						for (int j=0; j< numberTreatmentReplicas; j++){
							float minus = treatmentChromMinus[j].sumScoreBP(start, stop);
							scores[j] += minus;
							totalCounts += minus;
							tExonCounts[j][k] = minus;
						}
						for (int j=0; j< numberControlReplicas; j++){
							float minus = controlChromMinus[j].sumScoreBP(start, stop);
							scores[j+numberTreatmentReplicas] += minus;
							totalCounts += minus;
							cExonCounts[j][k] = minus;
						}
					}
				}
			}
			//nope combine strand counts
			else {
				for (int k=0; k< exons.length; k++){
					int start = exons[k].getStart();
					int stop = exons[k].getEnd();
					for (int j=0; j< numberTreatmentReplicas; j++){
						float plus = treatmentChromPlus[j].sumScoreBP(start, stop);
						float minus = treatmentChromMinus[j].sumScoreBP(start, stop);
						float total = plus+minus;
						scores[j] += total;
						totalCounts += total;
						tExonCounts[j][k] = total;
					}
					for (int j=0; j< numberControlReplicas; j++){
						float plus = controlChromPlus[j].sumScoreBP(start, stop);
						float minus = controlChromMinus[j].sumScoreBP(start, stop);
						float total = plus+minus;
						scores[j+numberTreatmentReplicas] += total;
						totalCounts+= total;
						cExonCounts[j][k] = total;
					}
				}
			}
			//save scores ?
			if (totalCounts==0) continue;
			
			genes[i].setScores(scores);
			geneLinesWithReadsAL.add(genes[i]);
			totalInterrogatedTreatmentObservations += Num.sumArray(tExonCounts);
			totalInterrogatedControlObservations += Num.sumArray(cExonCounts);
			genes[i].setTreatmentExonCounts(tExonCounts);
			genes[i].setControlExonCounts(cExonCounts);
		}

	}
	
	
	
	/**Looks for minimum number of reads and minimum log2Ratio difference between exon counts.*/
	private boolean checkForMinimums(float[] tExonCounts, float[] cExonCounts){
		for (int i=0; i< tExonCounts.length; i++){
			//check number of reads
			float num = tExonCounts[i] + cExonCounts[i];
			if (num < minimumExonicReadsForChiSquare) continue;
			//check ratio
			float logRatio = Math.abs(calculateLog2Ratio(tExonCounts[i], cExonCounts[i]));
			if (logRatio >= minimumSpliceLog2Ratio) return true;
		}
		return false;
	}

	/**Collects and calculates a bunch of stats re the PointData.*/
	private void calculateReadCountStatistics(){
		//fetch treatment PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (treatmentPointDirs);
		treatmentPlusPointData = PointData.convertArrayList2Array(combo[0]);
		treatmentMinusPointData = PointData.convertArrayList2Array(combo[1]);
		numberTreatmentObservations = PointData.totalObservationsMultiPointData(treatmentPlusPointData);
		numberTreatmentObservations += PointData.totalObservationsMultiPointData(treatmentMinusPointData);
		millionMappedTreatmentReads = numberTreatmentObservations/1000000;
		System.out.println("\t"+(int)numberTreatmentObservations+" Treatment Observations");

		//likewise control data
		combo = PointData.fetchStrandedPointDataNoMerge (controlPointDirs);
		controlPlusPointData = PointData.convertArrayList2Array(combo[0]);
		controlMinusPointData = PointData.convertArrayList2Array(combo[1]);
		numberControlObservations = PointData.totalObservationsMultiPointData(controlPlusPointData);
		numberControlObservations += PointData.totalObservationsMultiPointData(controlMinusPointData);
		System.out.println("\t"+(int)numberControlObservations+" Control Observations");		
		millionMappedControlReads = numberControlObservations/1000000;
		
		//scalars for log2ratio
		scalarTC = numberTreatmentObservations/ numberControlObservations;
		scalarCT = numberControlObservations/numberTreatmentObservations;

		//set number of replicas
		numberTreatmentReplicas = treatmentPointDirs.length;
		numberControlReplicas = controlPointDirs.length;
		totalNumberReplicas = numberTreatmentReplicas + numberControlReplicas;

		numberTreatmentReads = new int[numberTreatmentReplicas];
		numberControlReads = new int[numberControlReplicas];
	}

	/**Shifts the positions halfPeakShift (+ for sense, - for antisense) sets the positions into the data
	 * returns all of the positions after sorting. May replace all scores with 1 if stripScores == true.*/
	private void fetchShiftStripPointData(){
		for (int i=0; i< treatmentChromPlus.length; i++){
			int[] p = treatmentChromPlus[i].getPositions();
			if (halfPeakShift !=0) addShift(p,halfPeakShift);
			treatmentChromPlus[i].stripScores();
		}
		for (int i=0; i< treatmentChromMinus.length; i++){
			int[] p = treatmentChromMinus[i].getPositions();
			if (halfPeakShift !=0) addShift(p, -1*halfPeakShift);
			treatmentChromMinus[i].stripScores();
		}
		for (int i=0; i< controlChromPlus.length; i++){
			int[] p = controlChromPlus[i].getPositions();
			if (halfPeakShift !=0) addShift(p,halfPeakShift);
			controlChromPlus[i].stripScores();
		}
		for (int i=0; i< controlChromMinus.length; i++){
			int[] p = controlChromMinus[i].getPositions();
			if (halfPeakShift !=0) addShift(p, -1*halfPeakShift);
			controlChromMinus[i].stripScores();
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
		new MultipleReplicaDefinedRegionScanSeqs(args);
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
					case 's': saveDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'p': peakShift = Integer.parseInt(args[++i]); break;
					case 'u': refSeqFile = new File(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'i': scoreIntrons = true; break;
					case 'f': minFDR = Float.parseFloat(args[++i]); break;
					case 'l': minLog2Ratio = Float.parseFloat(args[++i]); break;
					case 'x': calculateChiSquare = false; break;
					case 'o': removeOverlappingRegions = false; break;
					case 'z': printExonReplicaMatrix = true; break;
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
		if (treatmentPointDirs == null || treatmentPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your treatment PointData directories(s)!\n");
		//only one directory look deeper
		if (treatmentPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(treatmentPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) treatmentPointDirs = otherDirs;
		}
		//control data
		if (controlPointDirs != null){
			if (controlPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your control PointData directories(s)!\n");
			//only one directory look deeper
			if (controlPointDirs.length == 1){
				File[] otherDirs = IO.extractOnlyDirectories(controlPointDirs[0]);
				if (otherDirs != null && otherDirs.length > 0) controlPointDirs = otherDirs;
			}
		}
		//set half peak shift and windowSize
		halfPeakShift = (int)Math.round( ((double)peakShift)/2 );

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		saveDirectory.mkdir();
		
		//check for R and required libraries
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}
		else {
			String errors = IO.runRCommandLookForError("library(DESeq)", fullPathToR, saveDirectory);
			if (errors == null || errors.length() !=0){
				Misc.printExit("\nError: Cannot find the required R library.  Did you install DESeq " +
						"(http://www-huber.embl.de/users/anders/DESeq/)?  See the author's websites for installation instructions. Once installed, " +
						"launch an R terminal and type 'library(DESeq)' to see if it is present. R error message:\n\t\t"+errors+"\n\n");
			}
		}

		//look for bed file
		if (refSeqFile == null && bedFile == null){
			Misc.printExit("\nPlease enter a regions file to use in scoring regions.\n");
		}
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                   Multiple Replica Defined Region Scan Seqs: April 2010          **\n" +
				"**************************************************************************************\n" +
				
				"\nWARNING: This application is depreciated and no longer maintained. Use the\n"+
				"                       OverdispersedRegionScanSeqs.\n\n"+
				
				"MRDRSS takes chromosome specific PointData xxx.bar.zip files and extracts reads under\n" +
				"each region or gene's exons to calculate several statistics using S. Anders' DESeq\n" +
				"package including a negative binomial p-value and Benjamini-Hochberg FDR for \n" +
				"differential expression. A chi-square test of independence between the\n" +
				"exon read count distributions is used to detect possible alternative splicing. Several\n" +
				"measures of read counts are provided including counts for each replica, a normalized\n" +
				"log2 ratio, and RPKMs (# reads per kb of interrogated region per total million\n" +
				"mapped reads) as well as DESeq's variance adjusted count data (use these values for\n" +
				"clustering, correlation, and other microarray type analysis). Three results files are\n" +
				"written: two spread sheets containing all of the regions/genes and those that pass\n" +
				"the thresholds as well as an egr region graph file for visualization in IGB.\n\n"+

				"Options:\n"+
				"-s Save directory, full path.\n"+
				"-t Treatment replica PointData directories, full path, comma delimited, no spaces,\n" +
				"       one per biological replica. Use the PointDataManipulator app to merge same\n" +
				"       replica and technical replica datasets. Each directory should contain stranded\n" +
				"       chromosome specific xxx_-/+_.bar.zip files. Alternatively, provide one\n" +
				"       directory that contains multiple biological replical PointData directories.\n" +
				"-c Control PointData directories, ditto. \n" +
				"-r Full path to R loaded with DESeq library, defaults to '/usr/bin/R' file, see\n" +
				"       http://www-huber.embl.de/users/anders/DESeq/ . Type 'library(DESeq)' in\n" +
				"       an R terminal to see if it is installed.\n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (name1 name2(optional) chrom strand\n" +
				"       txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds). NOTE:\n"+
				"       this table should contain only one composite transcript per gene (e.g. use\n" +
				"       Ensembl genes NOT transcripts). Otherwise set the -o option.\n"+
				"-b (Or) a bed file (chr, start, stop,...), full path, See,\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1\n"+
				"-d Data is stranded. Only collect reads from the same strand as the annotation.\n"+

				"\nAdvanced Options:\n"+
				"-o Don't remove overlapping exons, defaults to filtering gene annotation for overlaps.\n"+
				"-i Score introns instead of exons.\n"+
				"-f Minimum FDR threshold, defaults to 10 (-10Log10(FDR=0.1))\n"+
				"-l Minimum absolute log2 ratio threshold, defaults to 0.585 (1.5x)\n"+
				"-p Peak shift, average distance between + and - strand peaks for chIP-Seq data, see\n" +
				"       PeakShiftFinder. For RNA-Seq consider setting it to the smallest expected\n" +
				"       fragment size. Will be used to shift the PointData 3' by 1/2 the peak shift.\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/MultipleReplicaDefinedRegionScanSeqs -t\n" +
				"      /Data/PolIIRep1/,/Data/PolIIRep2/ -c /Data/Input1/,Data/Input2/ -s\n" +
				"      /Data/PolIIResults/ -f 5 -l 0.322 \n\n" +

		"**************************************************************************************\n");

	}

	public int getTotalInterrogatedTreatmentObservations() {
		return totalInterrogatedTreatmentObservations;
	}

	public int getTotalInterrogatedControlObservations() {
		return totalInterrogatedControlObservations;
	}

	public double getNumberTreatmentObservations() {
		return numberTreatmentObservations;
	}

	public double getNumberControlObservations() {
		return numberControlObservations;
	}
}
