package edu.utah.kohli;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.CombinePValues;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class SmallPanelCopyAlterationAnalysis {
	
	//user input
	private File patientAnnotationFile = null;
	private File regionReadCountData = null;
	private File sampleNormalizerGenes = null;
	private File resultsDirectory = null;
	private File fullPathToR = new File("/usr/local/bin/R");
	private boolean useBestScalar = true;
	
	private TreeMap<String, KohliPatient> patients = new TreeMap<String,KohliPatient>();
	private TreeMap<String, KohliSample> samples = new TreeMap<String,KohliSample>();
	private ArrayList<CaptureRegion> captureRegions = new ArrayList<CaptureRegion>();
	private LinkedHashMap<String, KohliGene> genes = null;
	private ArrayList<String> germlineSamplesToRemoveFromPoN = null;
	private String testedGenes = null;
	private KohliSample testSample = null;
	private File workingSaveDir = null;
	private TreeSet<String> genesToNormWith = new TreeSet<String>();
	private double bestScalar = -1;
	private CombinePValues combinePValues = new CombinePValues();
	
	//params to pick genes to use in normalizing a test sample
	private double maximumZScore = 3;
	private double minimumFractionPassingGeneCaptureRegions = 0.75;
	private double minPosZScores = 0.4;
	private double maxPosZScores = 0.6;
	private double minimumNumberPassingGenes = 3;
	
	//ttest pval thresholds
	private double minimumAdjustedCRTTestPval = 0.01;
	
	
	public SmallPanelCopyAlterationAnalysis(String[] args) throws Exception {
		processArgs(args);
		
		//load patient sample annotations
		loadPatients();
		
		//load read coverage for each sample, region, and gene
		loadReadCoverageData();
		
		//might be null or just for some samples
		loadSampleNormalizerGenes();

		//mean scale the germline samples to 1000, needed for resetting
		IO.pl("Initial all germline all region normalization to 1000...");
		if (germlineSamplesToRemoveFromPoN.size()!=0) IO.pl("\tExcluding from PoN: "+germlineSamplesToRemoveFromPoN);
		meanScaleGermlineReadCoverageDatasetsAllGenes(1000.0, captureRegions, true);
		
		runCfDNASamples();
		
		//runAllGermlinesAsTest();

		IO.pl("\nCOMPLETE!");

	}
	
	private void runAllGermlinesAsTest() throws IOException {
		IO.pl("Running germline samples...");
		int numNotProc = 0;
		int numProc = 0;
		for (KohliSample ks: samples.values()) {
			if (ks.isGermlineSample()==false) continue;
			testSample = ks;

			workingSaveDir = new File(resultsDirectory, testSample.getPatient().getHciPatientId()+"_"+testSample.getSampleId());
			workingSaveDir.mkdirs();
			IO.pl("Testing Sample "+workingSaveDir.getName()+"...");
			
			//exclude the test from the pon
			germlineSamplesToRemoveFromPoN.add(testSample.getSampleId());
			testSample.setExcludeFromPon(true);
			
			//mean scale the pon minus the "test" sample
			meanScaleGermlineReadCoverageDatasetsAllGenes(1000.0, captureRegions, true);
			
			//any user defined normalizer genes?  From the AccuGenomicGeneSelector app?
			if (testSample.getGenesToUseInNormalization()!=null) {

				IO.pl("\tUser defined genes for normalization: "+testSample.getGenesToUseInNormalization());
				// Collect the Capture Regions to use in normalizing
				ArrayList<CaptureRegion> crsToNorm = new ArrayList<CaptureRegion>();
				genesToNormWith.clear();
				for (String gn: testSample.getGenesToUseInNormalization()) {
					KohliGene kg = genes.get(gn);
					crsToNorm.addAll(kg.getCaptureRegions());
					genesToNormWith.add(gn);
				}

				//Normalize the PoN
				IO.pl("\tNormalizing the PoN with the select genes...");
				meanScaleGermlineReadCoverageDatasetsAllGenes(1000.0, crsToNorm, false);

				//Normalize the test
				IO.pl("\tNormalizing the test with the select genes...");
				meanScaleTestReadCoverageDatasetsWithDefinedCaptureRegions(1000.0, crsToNorm);
				numNotProc++;
			}

			// nope want this app to scan for genes to normalize with
			else {
				IO.pl("\tScanning scalars to identify genes in the test that resemble the PoN for normalization...");
				
				//ZScore gene selection method
				genesToNormWith = maximizeGenesWithInPoNZScores();

				//Too few genes?
				if (genesToNormWith.size()==0) {
					IO.pl("\tERROR: too few genes found to normalize! SKIPPING SAMPLE");
					numNotProc++;
					continue;
				}
				//All on one chrom
				else if (chromosomeCheck(genesToNormWith) == false) {
					IO.pl("\tERROR: all genes to normalize were found on one chromosome! SKIPPING SAMPLE");
					numNotProc++;
					continue;
				}
				numProc++;

				//Do they want to use the best scalar method?
				if (useBestScalar) {
					//normalize the test with the best scalar
					IO.pl("\tNormalizing the test sample with the best scalar "+bestScalar);
					meanScaleTestReadCoverageDatasetsWithDefinedCaptureRegions(bestScalar , captureRegions);
				}
				//Nope, want to normalize using just the select genes
				else {
					IO.pl("\tGenes for normalization: "+genesToNormWith);
					// Collect the Capture Regions to use in normalizing
					ArrayList<CaptureRegion> crsToNorm = new ArrayList<CaptureRegion>();
					for (String gn: genesToNormWith) {
						KohliGene kg = genes.get(gn);
						crsToNorm.addAll(kg.getCaptureRegions());
					}

					//Normalize the PoN
					IO.pl("\tNormalizing the PoN with the select genes...");
					meanScaleGermlineReadCoverageDatasetsAllGenes(1000.0, crsToNorm, false);

					//Normalize the test
					IO.pl("\tNormalizing the test with the select genes...");
					meanScaleTestReadCoverageDatasetsWithDefinedCaptureRegions(1000.0, crsToNorm);
				}
			}

			//Generate some plots for each gene
			generateBoxPlots("("+Misc.treeSetToString(genesToNormWith, ",")+")");
			generateZScoreHistograms("("+Misc.treeSetToString(genesToNormWith, ",")+")");

			//Calculate ttest pvalues
			calculateCaptureRegionTTestPValues();

			//Stat copy changes
			generateGeneStatistics();			
			//remove the test from the PoN blocker
			germlineSamplesToRemoveFromPoN.remove(germlineSamplesToRemoveFromPoN.size()-1);
			testSample.setExcludeFromPon(false);
		}
		IO.pl("\tNum Samples Normalized:\t"+numProc+"\tSkipped:\t"+numNotProc+"\n");
	}

	private void runCfDNASamples() throws IOException {
		IO.pl("\nRunning cfDNA samples...");
		int numNotProc = 0;
		int numProc = 0;
		for (KohliPatient kp: patients.values()) {
			IO.pl("\nTesting Patient "+kp.getHciPatientId()+"...");

			for (KohliSample ks: kp.getNonGermlineSamples()) {
				testSample = ks;
				
				//"22712X28") || id.equals("22712X9") ||				
//String id = testSample.getSampleId();

//if (id.equals("22712X16") || id.equals("22712X28") || id.equals("22712X9")) {  
				
				workingSaveDir = new File(resultsDirectory, kp.getHciPatientId()+"_"+testSample.getSampleId());
				workingSaveDir.mkdirs();
				IO.pl("Testing Sample "+workingSaveDir.getName()+"...");

				//any user defined normalizer genes?  From the AccuGenomicGeneSelector app?
				if (testSample.getGenesToUseInNormalization()!=null) {

					IO.pl("\tUser defined genes for normalization: "+testSample.getGenesToUseInNormalization());
					// Collect the Capture Regions to use in normalizing
					ArrayList<CaptureRegion> crsToNorm = new ArrayList<CaptureRegion>();
					genesToNormWith.clear();
					for (String gn: testSample.getGenesToUseInNormalization()) {
						KohliGene kg = genes.get(gn);
						crsToNorm.addAll(kg.getCaptureRegions());
						genesToNormWith.add(gn);
					}

					//Normalize the PoN
					IO.pl("\tNormalizing the PoN with the select genes...");
					meanScaleGermlineReadCoverageDatasetsAllGenes(1000.0, crsToNorm, false);

					//Normalize the test
					IO.pl("\tNormalizing the test with the select genes...");
					meanScaleTestReadCoverageDatasetsWithDefinedCaptureRegions(1000.0, crsToNorm);
					numProc++;
				}

				// nope want this app to scan for genes to normalize with
				else {
					IO.pl("\tScanning scalars to identify genes in the test that resemble the PoN for normalization...");
					
					//ZScore gene selection method
					genesToNormWith = maximizeGenesWithInPoNZScores();

					//Too few genes?
					if (genesToNormWith.size()==0) {
						IO.pl("\tERROR: too few genes found to normalize! SKIPPING SAMPLE");
						numNotProc++;
						continue;
					}
					//All on one chrom
					else if (chromosomeCheck(genesToNormWith) == false) {
						IO.pl("\tERROR: all genes to normalize were found on one chromosome! SKIPPING SAMPLE");
						numNotProc++;
						continue;
					}
					numProc++;

					//Do they want to use the best scalar method?
					if (useBestScalar) {
						//normalize the test with the best scalar
						IO.pl("\tNormalizing the test sample with the best scalar "+bestScalar);
						meanScaleTestReadCoverageDatasetsWithDefinedCaptureRegions(bestScalar , captureRegions);
					}
					//Nope, want to normalize using just the select genes
					else {
						IO.pl("\tGenes for normalization: "+genesToNormWith);
						// Collect the Capture Regions to use in normalizing
						ArrayList<CaptureRegion> crsToNorm = new ArrayList<CaptureRegion>();
						for (String gn: genesToNormWith) {
							KohliGene kg = genes.get(gn);
							crsToNorm.addAll(kg.getCaptureRegions());
						}

						//Normalize the PoN
						IO.pl("\tNormalizing the PoN with the select genes...");
						meanScaleGermlineReadCoverageDatasetsAllGenes(1000.0, crsToNorm, false);

						//Normalize the test
						IO.pl("\tNormalizing the test with the select genes...");
						meanScaleTestReadCoverageDatasetsWithDefinedCaptureRegions(1000.0, crsToNorm);
					}
				}

				//Generate some plots for each gene
				generateBoxPlots("("+Misc.treeSetToString(genesToNormWith, ",")+")");
				generateZScoreHistograms("("+Misc.treeSetToString(genesToNormWith, ",")+")");

				//Calculate ttest pvalues
				calculateCaptureRegionTTestPValues();

				//Stat copy changes
				generateGeneStatistics();

				//Reset the PoN to the original coverage values, might or might not be needed depending on the norm proceedure above
				IO.pl("\tResetting the PoN to the original values...");
				for (CaptureRegion cr: captureRegions) cr.copyInitialPoNToCurrentPon();
				
			//}
			}
		}
		IO.pl("\tNum samples successfully normalized:\t"+numProc+"\tNum skipped:\t"+numNotProc+"\n");
	}
	
	private void calculateCaptureRegionTTestPValues() throws IOException {
		IO.pl("\tCalculating ttest pvalues for each capture region...");
		double numberCaptureRegions = captureRegions.size();
		//for each gene
		for (KohliGene kg: genes.values()) {
			double[][] crsScores = kg.fetchCaptureRegionTestAndPoNScaledScores();
			double[] pvals = fetchTTestPValues(crsScores, kg.getGeneName());
			ArrayList<CaptureRegion> crs = kg.getCaptureRegions();
			
			double numPassingPVals = 0;
			for (int i=0; i< crs.size(); i++) {
				double adjPval = pvals[i] * numberCaptureRegions;
				if (adjPval>1) adjPval = 1.0;
				crs.get(i).setpValTTest(adjPval);
				pvals[i] = adjPval;
				if (adjPval <= minimumAdjustedCRTTestPval) numPassingPVals++;
			}
			//set combine pval
			double combinePVal = combinePValues.calculateCombinePValues(pvals);
			kg.setCombinePValueTTest(combinePVal);
			kg.setFractionPassingAdjPvals(numPassingPVals/(double)crs.size());
		}
	}

	private boolean chromosomeCheck(TreeSet<String> geneNames) {
		HashSet<String> chrs = new HashSet<String>();
		for (String geneName: geneNames) chrs.add(genes.get(geneName).getChromosome());
		if (chrs.size()>1) return true;
		return false;
	}

	private void generateGeneStatistics() {
		IO.pl("\tGenerating statistic spreadsheet...");
		StringBuilder sb = new StringBuilder();
		sb.append("Patient\tTest Sample\tGermline\tUsed To Norm\tGene\t# Capture Regions\t"
				+ "# CRs Test Outside PoN\tMean CR Z-Score\tPseudoMedian CR Z-Score\tMean CR Scaled Test\tMean CR Scaled PoN\tLog2(Test/PoN)\tTTestPVals\tTTestCombineAdjPVal\tTTestFracPassAdjPVals\n");
		for (KohliGene kg: genes.values()) {
			sb.append(testSample.getPatient().getHciPatientId()); sb.append("\t");
			sb.append(testSample.getSampleId()); sb.append("\t");
			sb.append(testSample.isGermlineSample()); sb.append("\t");
			sb.append(genesToNormWith.contains(kg.getGeneName()));  sb.append("\t");
			sb.append(kg.getGeneName()); sb.append("\t");
			sb.append(kg.getCaptureRegions().size());  sb.append("\t");
			sb.append(kg.getNumberTestValuesOutsideCaptureRegions()); sb.append("\t");
			sb.append(kg.getMeanCaptureRegionZScore()); sb.append("\t");
			sb.append(kg.getPseudoMedianCaptureRegionZScore()); sb.append("\t");
			double meanTest = kg.getMeanCaptureRegionScaledTestScore();
			sb.append(meanTest); sb.append("\t");
			double meanPoN = kg.getMeanCaptureRegionPoNScore();
			sb.append(meanPoN); sb.append("\t");
			sb.append(Num.log2(meanTest/meanPoN)); sb.append("\t");
			
			//ttest pvalues
			ArrayList<CaptureRegion> crs = kg.getCaptureRegions();
			double[] ttPVals = new double[crs.size()];
			for (int i=0; i< ttPVals.length; i++) ttPVals[i] = crs.get(i).getpValTTest();
			sb.append(Num.doubleArrayToString(ttPVals, ","));
			sb.append("\t");
			sb.append(kg.getCombinePValueTTest() );
			sb.append("\t");
			sb.append(kg.getFractionPassingAdjPvals() );
			
			sb.append("\n");
		}
		
		File spreadsheet = new File( workingSaveDir, workingSaveDir.getName()+".geneStats.xls" );
		IO.writeString(sb.toString(), spreadsheet);
	}
	
	private void generateDetailedGeneStatistics() {
		IO.pl("\tGenerating detailed statistic spreadsheet...");
		
		IO.pl("Patient: "+testSample.getPatient().getHciPatientId()+"\tSample: "+testSample.getSampleId());
		for (KohliGene kg: genes.values()) {
			double[][] crsScores = kg.fetchCaptureRegionTestAndPoNScaledScores();
			int numScores = crsScores[0].length;
			int numCrs = crsScores.length;
			
			//header
			IO.pl("\nGene: "+kg.getGeneName());
			for (int i=0; i< numCrs; i++) {
				IO.p("\t");
				IO.p("CapReg_"+i);
			}
			IO.pl();

			//for each score
			for (int j=0; j< numScores; j++) {
				if (j==0) IO.p("tumor");
				else IO.p("PoN_"+j);
				//for each capture region
				for (int i=0; i<numCrs; i++) {
					IO.p("\t");
					IO.p(crsScores[i][j]);
				}
				IO.pl();
			}
		}
		IO.pl();
	
	}

	private void generateBoxPlots(String normGenes) throws IOException {
		IO.pl("\tGenerating boxplots...");
		//write out data files for R
		File rDir = new File(workingSaveDir, "R");
		rDir.mkdir();
		File boxPlotDataGermline = new File(rDir, "germlineBoxPlotData.txt");
		saveScaledGermlineReadCoverageDatasetsForBoxPlot(boxPlotDataGermline);
		File boxPlotDataTest = new File(rDir, "testBoxPlotData.txt");
		saveScaledTestReadCoverageDatasetsForBoxPlot(boxPlotDataTest);
		
		//boxPlot dir
		File boxPlots = new File(workingSaveDir, "BoxPlots");
		boxPlots.mkdir();
		
		ArrayList<String> rScriptAL = new ArrayList<String>();
		rScriptAL.add("library(ggplot2)");
		rScriptAL.add("setwd('"+boxPlots.getCanonicalPath()+"')");
		rScriptAL.add("df = read.table('"+boxPlotDataGermline.getCanonicalPath()+"', header = TRUE)");
		rScriptAL.add("testDf = read.table('"+boxPlotDataTest.getCanonicalPath()+"', header=TRUE)");
		rScriptAL.add("normGenes = '"+normGenes+"'");
		rScriptAL.add("testDataset = '"+testSample.getSampleId()+"'");
		rScriptAL.add("genes = c("+testedGenes+")");
		rScriptAL.add("for (x in genes) {");
		rScriptAL.add("  print(x)");
		rScriptAL.add("  toSub=paste(x,'_',sep='')");
		rScriptAL.add("  dfSub = subset(df, grepl(toSub, Gene))");
		rScriptAL.add("  dfSubTest = subset(testDf, grepl(toSub, Gene))");
		rScriptAL.add("  boxPlotObj = ggplot(dfSub, aes(x=Gene, y=ScaledCount)) + ");
		//rScriptAL.add("  geom_boxplot()+ coord_cartesian(ylim = c(250, 1750)) + ");
		rScriptAL.add("  geom_boxplot()+ ");
		rScriptAL.add("  ggtitle( paste(testDataset, x, normGenes)) + ");
		rScriptAL.add("  ylab('1K Mean Scaled Counts') + ");
		rScriptAL.add("  theme( ");
		rScriptAL.add("    plot.title = element_text(size=12,face='bold'), ");
		rScriptAL.add("    axis.title.x = element_blank(), ");
		rScriptAL.add("    axis.title.y = element_text(face='bold', size=12), ");
		rScriptAL.add("    axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + ");
		rScriptAL.add("  geom_point(data=dfSubTest, aes(x=Gene, y=ScaledCount), color = 'red', shape=8) ");
		rScriptAL.add("  ggsave(file=paste(x,'.box.png', sep = ''), plot=boxPlotObj, width=10, height=8) ");
		rScriptAL.add("}");
		rScriptAL.add("print('COMPLETE')");
		String rScript = Misc.stringArrayListToString(rScriptAL, "\n");
		File rScriptFile = new File(rDir, "boxplotScript.R");
		File rLogFile = new File(rDir, "boxplot.Rout");
		IO.writeString(rScript, rScriptFile );
		
		//make command
		String[] command = new String[] {
				fullPathToR.getCanonicalPath(),
				"CMD",
				"BATCH",
				"--no-save",
				"--no-restore",
				rScriptFile.getCanonicalPath(),
				rLogFile.getCanonicalPath()};			
		//execute
		IO.executeCommandLine(command);
		
		//read in results 
		String[] results = IO.loadFile(rLogFile);
		boolean complete = false;
		for (String s: results) if (s.contains("COMPLETE"))complete = true;
		if (complete == false) Misc.printErrAndExit("\nERROR, failed to find 'COMPLETE' in the boxplot R output, see "+rLogFile);
	}
	
	private double[] fetchTTestPValues(double[][] crsScores, String geneName) throws IOException {
		//write out data files for R
		File rDir = new File(workingSaveDir, "R");
		rDir.mkdir();
		
		int numScores = crsScores[0].length;
		int numCrs = crsScores.length;
		
		//for each score
		File tmpCounts = new File(rDir, "ttestCounts_"+geneName+".txt");
		PrintWriter out = new PrintWriter (new FileWriter(tmpCounts));
		for (int j=0; j< numScores; j++) {
			//for each capture region
			for (int i=0; i<numCrs; i++) {
				if (i!=0) out.print("\t");
				out.print(crsScores[i][j]);
			}
			out.println();
		}
		out.close();
		
		File tmpResults = new File(rDir, "ttestPValues_"+geneName+".txt");
		
		ArrayList<String> rScriptAL = new ArrayList<String>();
		rScriptAL.add("counts = read.table(file = '"+tmpCounts.getCanonicalPath()+"', header = FALSE, sep = '\t')");
		rScriptAL.add("nr = nrow(counts)");
		rScriptAL.add("nc = ncol(counts)");
		rScriptAL.add("pvalues = rep(NA, times=nc)");
		rScriptAL.add("for (x in 1:nc){");
		rScriptAL.add("  A = counts[1,x]");
		rScriptAL.add("  B = counts[2:nr, x]");
		rScriptAL.add("  R = t.test(A, B, var.equal=TRUE)");
		rScriptAL.add("  pvalues[x] = R$p.value");
		rScriptAL.add("}");
		rScriptAL.add("write.csv(pvalues,file='"+tmpResults.getCanonicalPath()+"',row.names=F)");
		rScriptAL.add("print('COMPLETE')");
		String rScript = Misc.stringArrayListToString(rScriptAL, "\n");
		File rScriptFile = new File(rDir, "ttestScript.R");
		File rLogFile = new File(rDir, "ttest.Rout");
		IO.writeString(rScript, rScriptFile );
		
		//make command
		String[] command = new String[] {
				fullPathToR.getCanonicalPath(),
				"CMD",
				"BATCH",
				"--no-save",
				"--no-restore",
				rScriptFile.getCanonicalPath(),
				rLogFile.getCanonicalPath()};			
		//execute
		IO.executeCommandLine(command);
		
		//read in log file 
		String[] results = IO.loadFile(rLogFile);
		boolean complete = false;
		for (String s: results) if (s.contains("[1] \"COMPLETE\""))complete = true;
		if (complete == false) Misc.printErrAndExit("\nERROR, failed to find 'COMPLETE' in the ttest R output, see "+rLogFile);
		
		//read in results, skip first line
		String[] pvalues = IO.loadFile(tmpResults);
		double[] pvals = new double[pvalues.length-1];
		int index = 0;
		for (int i=1; i< pvalues.length; i++) pvals[index++] = Double.parseDouble(pvalues[i]);
		
		IO.deleteDirectory(rDir);
		return pvals;
	}
	
	private TreeSet<String> maximizeGenesWithInPoNZScores() {
		bestScalar = -1;
		NormalizationResult[] res = new NormalizationResult[10000 - 100 +1];
		int counter = 0;
		for (int i=100; i<=10000; i++) {
			//scale all to i	
			meanScaleTestReadCoverageDatasetsWithDefinedCaptureRegions(i , captureRegions);
			res[counter++] = countNumberGenesWithTestInPon(i);
		}
		//sort by highest num genes and then sum fraction copy regions
		Arrays.sort(res);
		IO.pl("\tScalar\t#Passing\tGeneNames");
		TreeSet<String> toNorm = new TreeSet<String>();
		if (res[0].numPassingGenes >= minimumNumberPassingGenes ) {
			IO.pl(res[0]);
			for (KohliGene kg: res[0].genesPassingThresholds) {
				toNorm.add(kg.getGeneName());
			}
			bestScalar = res[0].normTarget;
		}
		return toNorm;
	}
		
	private NormalizationResult countNumberGenesWithTestInPon(int normTarget) {
		//for each gene, count the number of capture regions where the test is within the PoN zscores
		ArrayList<KohliGene> testInPoN = new ArrayList<KohliGene>();
		ArrayList<Double> fracCRs = new ArrayList<Double>();
		ArrayList<Double> fracPosZScoresAL = new ArrayList<Double>();
		for (KohliGene kg : genes.values()) {
			double inPoN = 0;
			double positiveZScores = 0;
			double negativeZScores = 0;
			for (CaptureRegion cr: kg.getCaptureRegions()) {
				double zscore = cr.calculateZScoreFromPoN();
				if (Math.abs(zscore) <= maximumZScore) inPoN++;
				if (zscore > 0) positiveZScores++;
				else if (zscore < 0) negativeZScores++;
			}
			
			double frac = inPoN/ (double)kg.getCaptureRegions().size();
			double fracPosZScores = positiveZScores/(positiveZScores+negativeZScores);
			if (frac >= minimumFractionPassingGeneCaptureRegions && (fracPosZScores >= minPosZScores && fracPosZScores <= maxPosZScores)) {
				testInPoN.add(kg);
				fracCRs.add(frac);
				fracPosZScoresAL.add(fracPosZScores);
			}
		}
		return new NormalizationResult(normTarget, testInPoN, fracCRs, fracPosZScoresAL);
	}
	
	private class NormalizationResult implements Comparable<NormalizationResult>{
		
		int normTarget;
		ArrayList<KohliGene> genesPassingThresholds = null;
		ArrayList<Double> fracCRs = new ArrayList<Double>();
		ArrayList<Double> fracPosZScores = new ArrayList<Double>();
		int numPassingGenes;
		double sumFracCRs;
		int numPassingChromosomes;
		
		public NormalizationResult(int normTarget, ArrayList<KohliGene> genesPassingThresholds, ArrayList<Double> fracCRs, ArrayList<Double> fracPosZScores) {
			this.genesPassingThresholds = genesPassingThresholds;
			this.normTarget = normTarget;
			numPassingGenes = genesPassingThresholds.size();
			this.fracCRs = fracCRs;
			sumFracCRs = Num.sumArray(Num.arrayListOfDoubleToArray(fracCRs));
			this.fracPosZScores = fracPosZScores;
			numPassingChromosomes = countNumberChroms();
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder("\t");
			sb.append(normTarget);
			sb.append("\t");
			sb.append(numPassingGenes);
			if (numPassingGenes > 0) {
				sb.append("\t");
				sb.append(genesPassingThresholds.get(0).getGeneName());
				sb.append("(");
				sb.append(Num.formatNumber(fracCRs.get(0), 2));
				sb.append(")");
				for (int i=1; i< numPassingGenes; i++) {
					sb.append(",");
					sb.append(genesPassingThresholds.get(i).getGeneName());
					sb.append("(");
					sb.append(Num.formatNumber(fracCRs.get(i), 2));
					sb.append(")");
				}
			}
			
			return sb.toString();
		}

		public int compareTo(NormalizationResult o) {
			if (this.numPassingGenes > o.numPassingGenes) return -1;
			else if (this.numPassingGenes < o.numPassingGenes) return 1;
			//same number of passing genes
			else {
				//maximize on chr diversity
				if (this.numPassingChromosomes > o.numPassingChromosomes) return -1;
				else if (this.numPassingChromosomes < o.numPassingChromosomes) return 1;
				//same number chrs
				else {
					if (this.sumFracCRs > o.sumFracCRs) return -1;
					else if (this.sumFracCRs < o.sumFracCRs) return 1;
					return 0;
				}
			}
		}

		private int countNumberChroms() {
			HashSet<String> chrs = new HashSet<String>();
			for (KohliGene kg: genesPassingThresholds) chrs.add(kg.getChromosome());
			return chrs.size();
		}
	}
	
	private void generateZScoreHistograms(String normGenes) throws IOException {
		IO.pl("\tGenerating histograms...");
		
		//write out data files for R
		File rDir = new File(workingSaveDir, "R");
		rDir.mkdir();
		File histogramData = new File(rDir, "histogramData.txt");
		saveTestZScoresForHistograms(histogramData);
		
		//Histo dir
		File histogramPlots = new File(workingSaveDir, "HistogramPlots");
		histogramPlots.mkdir();

		ArrayList<String> rScriptAL = new ArrayList<String>();
		rScriptAL.add("library(ggplot2)");
		rScriptAL.add("setwd('"+histogramPlots.getCanonicalPath()+"')");
		rScriptAL.add("df = read.table('"+histogramData.getCanonicalPath()+"', header = TRUE)");
		rScriptAL.add("normGenes = '"+normGenes+"'");
		rScriptAL.add("testDataset = '"+testSample.getSampleId()+"'");
		rScriptAL.add("genes = c("+testedGenes+")");
		rScriptAL.add("for (x in genes) {");
		rScriptAL.add("  print(x)");
		rScriptAL.add("  toSub=paste(x,'_',sep='')");
		rScriptAL.add("  dfSub = subset(df, grepl(toSub, Gene))");
		rScriptAL.add("  obj = ggplot(dfSub, aes(x=ZScore)) +  ");
		rScriptAL.add("  geom_histogram(color='black',fill='gray',bins=20) +  ");
		rScriptAL.add("  geom_vline(aes(xintercept=mean(ZScore)),color='red', linetype='dashed', size=1)+  ");
		rScriptAL.add("  ggtitle( paste(testDataset, x, normGenes)) + ");
		rScriptAL.add("  scale_x_continuous(name = 'Capture Region Z-Scores') + ");
		rScriptAL.add("  scale_y_continuous(name = 'Count')");
		rScriptAL.add("  ggsave(file=paste(x,'.hist.png', sep = ''), plot=obj, width=10, height=8) ");
		rScriptAL.add("}");
		rScriptAL.add("print('COMPLETE')");
		String rScript = Misc.stringArrayListToString(rScriptAL, "\n");
		File rScriptFile = new File(rDir, "histogramScript.R");
		File rLogFile = new File(rDir, "histogram.Rout");
		IO.writeString(rScript, rScriptFile );
		
		//make command
		String[] command = new String[] {
				fullPathToR.getCanonicalPath(),
				"CMD",
				"BATCH",
				"--no-save",
				"--no-restore",
				rScriptFile.getCanonicalPath(),
				rLogFile.getCanonicalPath()};			
		//execute
		IO.executeCommandLine(command);
		
		//read in results 
		String[] results = IO.loadFile(rLogFile);
		boolean complete = false;
		for (String s: results) if (s.contains("COMPLETE"))complete = true;
		if (complete == false) Misc.printErrAndExit("\nERROR, failed to find 'COMPLETE' in the boxplot R output, see "+rLogFile);

	}
	
	private void saveTestZScoresForHistograms(File toSave) {
		//print out all of the scaled counts for boxplots
		StringBuilder sb = new StringBuilder();
		//header row
		sb.append("Gene\tZScore\n");
		String gene = "";
		int counter = 0;
		for (int i=0; i< captureRegions.size(); i++) {
			CaptureRegion cr = captureRegions.get(i);
			if (cr.getGene().equals(gene) == false) {
				gene = cr.getGene();
				counter = 0;
			}
			counter++;
			String counterString = Integer.toString(counter);
			if (counterString.length()==1) counterString = "00"+counterString;
			else if (counterString.length()==2) counterString = "0"+counterString;
			String name = gene+"_"+counterString;
			sb.append(name);
			sb.append("\t");
			sb.append(Double.toString(cr.calculateZScoreFromPoN()));
			sb.append("\n");
		}
		IO.writeString(sb.toString(), toSave);
	}

	private void saveScaledTestReadCoverageDatasetsForBoxPlot(File toSave) {
		//print out all of the scaled counts for boxplots
		StringBuilder sb = new StringBuilder();
		//header row
		sb.append("Gene\tScaledCount\n");
		String gene = "";
		int counter = 0;
		for (int i=0; i< captureRegions.size(); i++) {
			CaptureRegion cr = captureRegions.get(i);
			if (cr.getGene().equals(gene) == false) {
				gene = cr.getGene();
				counter = 0;
			}
			counter++;
			String counterString = Integer.toString(counter);
			if (counterString.length()==1) counterString = "00"+counterString;
			else if (counterString.length()==2) counterString = "0"+counterString;
			String name = gene+"_"+counterString;
			sb.append(name);
			sb.append("\t");
			sb.append(Double.toString(cr.getScaledTestCount()));
			sb.append("\n");
		}
		IO.writeString(sb.toString(), toSave);
	}

	private void meanScaleGermlineReadCoverageDatasetsAllGenes(double target, ArrayList<CaptureRegion> crsToScaleWith, boolean setInitial) {
		//scale the counts using the select capture regions for ALL germline samples
		for (KohliSample ks: samples.values()) {
			if (ks.isGermlineSample()) ks.createScaledRegionCountMap(crsToScaleWith, captureRegions, target);
		}

		//for all capture region, collect the scaled germline counts and save them to the cr for those good for PoN
		for (CaptureRegion cr: captureRegions) {
			String crId = cr.getOriginalInput();
			ArrayList<Double> germlineScaledCounts = new ArrayList<Double>();
			//for each germline pon sample
			for (KohliSample ks: samples.values()) {
				if (ks.isGermlineSample() && ks.isExcludeFromPon()==false) {
					HashMap<String, Double> sampleScaledCaptureRegionMeans = ks.getScaledReadCoverageCounts();
					germlineScaledCounts.add(sampleScaledCaptureRegionMeans.get(crId));
				}
			}
			//convert to double[] and set in cr
			cr.setPonScaledCountsCalculateStats(Num.arrayListOfDoubleToArray(germlineScaledCounts), setInitial);
		}
	}

	private void meanScaleTestReadCoverageDatasetsWithDefinedCaptureRegions(double target, ArrayList<CaptureRegion> forScaling) {
		testSample.createScaledRegionCountMap(forScaling, captureRegions, target);
		//load the capture regions with the scaled counts
		for (CaptureRegion cr: captureRegions) {
			String crId = cr.getOriginalInput();
			HashMap<String, Double> sampleScaledCaptureRegionMeans = testSample.getScaledReadCoverageCounts();
			cr.setScaledTestCount(sampleScaledCaptureRegionMeans.get(crId));
		}
	}
	
	private void saveScaledGermlineReadCoverageDatasetsForBoxPlot(File toSave) {
		//header row
		StringBuilder sb = new StringBuilder("Gene\tScaledCount\n");
		String gene = "";
		int counter = 0;
		for (int i=0; i< captureRegions.size(); i++) {
			CaptureRegion cr = captureRegions.get(i);
			String crId = cr.getOriginalInput();
			if (cr.getGene().equals(gene) == false) {
				gene = cr.getGene();
				counter = 0;
			}
			counter++;
			String counterString = Integer.toString(counter);
			if (counterString.length()==1) counterString = "00"+counterString;
			else if (counterString.length()==2) counterString = "0"+counterString;
			String name = gene+"_"+counterString;
			//for each sample
			for (KohliSample ks: samples.values()) {
				if (ks.isGermlineSample() == false || ks.isExcludeFromPon() == true) continue;
				sb.append(name);
				sb.append("\t");
				sb.append(ks.getScaledReadCoverageCounts().get(crId));
				sb.append("\n");
			}
		}
		IO.writeString(sb.toString(), toSave);
	}
	
	private void loadReadCoverageData() {
		IO.pl("Loading capture region count data...");
		String[] lines = IO.loadFileIntoStringArray(regionReadCountData);
		// parse the header line, create a hashMap of index and KohliSample
		// Chr_Start_Stop_Info 22712X10 22712X11 22712X12 22712X13 ...
		if (lines[0].startsWith("Chr_Start_Stop_Info") == false) Misc.printErrAndExit("First line in the read count data file doesn't start with 'Chr_Start_Stop_Info'? -> "+lines[0]);
		String[] header = Misc.TAB.split(lines[0]);
		KohliSample[] indexSamples = new KohliSample[header.length];
		for (int i=1; i< header.length; i++ ) {
			KohliSample ks = samples.get(header[i]);
			if (ks == null) Misc.printErrAndExit("Failed to find the sample "+header[i]+ "in "+samples);
			indexSamples[i] = ks;
		}
		
		// parse the data lines for each data line
		for (int i=1; i< lines.length; i++) {
			//chr16_72787084_72787208_ZFHX3	2169.5	5336.5	2498.5	2703 ...
			if (lines[i].startsWith("chr") == false) Misc.printErrAndExit("Data line doesn't start with 'chr'? "+lines[i]);
			String[] f = Misc.TAB.split(lines[i]);
			//use the first field to create a Capture Region
			CaptureRegion p = new CaptureRegion(f[0]);
			captureRegions.add(p);

			//for each data point
			for (int j=1; j< f.length; j++) {
				//get the sample, add the values
				KohliSample ks = indexSamples[j];
				ks.addCaptureRegionCount(p, f[j]);
			}
		}
		
		//make the genes
		genes = new LinkedHashMap<String, KohliGene>();
		
		//add in the capture regions
		for (CaptureRegion cr: captureRegions) {
			KohliGene kg = genes.get(cr.getGene());
			if (kg == null) {
				kg = new KohliGene(cr.getGene());
				genes.put(cr.getGene(), kg);
			}
			kg.addCaptureRegion(cr);
		}
		
		//set testedGeneNames
		//"'AR','BRCA2','CHEK2','MYC','NKX3-1','OPHN1','PIK3CA','PIK3CB','TP53','ZBTB16'"
		StringBuilder sb = new StringBuilder("'");
		for (String gn: genes.keySet()) {
			sb.append(gn);
			sb.append("','");
		}
		testedGenes = sb.toString();
		testedGenes = testedGenes.substring(0, testedGenes.length()-2);
		IO.pl("\tTested Genes: "+testedGenes);
		
	}
	
	private void loadSampleNormalizerGenes() {
		if (sampleNormalizerGenes == null) return;
		IO.pl("Loading genes to use in normalizing each sample data...");

		for (String l: IO.loadFileIntoStringArray(sampleNormalizerGenes)) {
			if (l.startsWith("#")) continue;
			//sampleId multiplier	gene1 gene2 ... space delimited
			String[] f = Misc.WHITESPACE.split(l);
			if (f.length < 3) Misc.printErrAndExit("\nError parsing more than one field from "+l+" from "+sampleNormalizerGenes);
			KohliSample ks = samples.get(f[0]);
			if (ks==null) Misc.printErrAndExit("\nError failed to find the sample id "+f[0]+" from "+sampleNormalizerGenes);
			HashSet<String> normGenes = new HashSet<String>();
			for (int i=2; i< f.length; i++) normGenes.add(f[i]);
			ks.setGenesToUseInNormalization(normGenes);
		}
	}


	private void loadPatients() {
		IO.pl("Loading patient and sample meta data...");
		String[] lines = IO.loadFileIntoStringArray(patientAnnotationFile);
		if (lines[0].equals("#HCIPersonID\tSampleID\tType\tDateDrawn") == false) Misc.printErrAndExit("First line in the anno file isn't '#HCIPersonID SampleID Type	DateDrawn'? -> "+lines[0]);
		for (int i=1; i< lines.length; i++) {
			String[] tokens = Misc.TAB.split(lines[i]);
			KohliPatient kp = patients.get(tokens[0]);
			if (kp == null) {
				kp = new KohliPatient(tokens);
				patients.put(tokens[0], kp);
			}
			else {
				boolean isGermline = false;
				if (tokens[2].toLowerCase().contains("germline")) isGermline = true;
				kp.getSamples().add(new KohliSample(kp, tokens[1], isGermline, tokens[3]));
			}
		}
		//load the sample treemap
		for (KohliPatient kp: patients.values()) {
			for (KohliSample ks: kp.getSamples()) {
				if (samples.containsKey(ks.getSampleId())) Misc.printErrAndExit("Duplicate sample found! "+ks.getSampleId());
				samples.put(ks.getSampleId(), ks);
			}
		}
		//any germline samples to exclude from the Pon?
		if (germlineSamplesToRemoveFromPoN.size()!=0) {
			for (String sampleId: germlineSamplesToRemoveFromPoN) {
				KohliSample ks = samples.get(sampleId);
				ks.setExcludeFromPon(true);
			}
		}
	}

	public static void main(String[] args) throws Exception {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SmallPanelCopyAlterationAnalysis(args);
	}	

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args) throws Exception{
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		String germlineSamplesNoPoN = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': patientAnnotationFile = new File(args[++i]).getCanonicalFile(); break;
					case 'b': useBestScalar = false; break;
					case 'c': regionReadCountData = new File(args[++i]); break;
					case 'r': resultsDirectory = new File(args[++i]).getCanonicalFile(); break;
					case 'g': germlineSamplesNoPoN = args[++i]; break;
					case 'e': fullPathToR = new File(args[++i]); break;
					case 'n': sampleNormalizerGenes = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//missing files?
		
		
		//any germline samples to exclude from PoN?
		germlineSamplesToRemoveFromPoN = new ArrayList<String>();
		if (germlineSamplesNoPoN != null) {
			for (String s: Misc.COMMA.split(germlineSamplesNoPoN)) germlineSamplesToRemoveFromPoN.add(s);
		}
		
		/*
		-p /Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/sampleMatchingPersonGNomIDs_13June2024.txt 
		-c /Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/KohliCopyRatio/ReadCoverageCounts/medianCounts150bpDoOverPostQC.txt 
		-r /Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/KohliCopyRatio/CopyAnalysisWithChrXDelme
		-n /Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/KohliCopyRatio/probePickedNormGenesV2.txt 
		-g 23802X4,23802X11
		*/


	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                     Small Panel Copy Ratio Analysis:   August 2024               **\n" +
				"**************************************************************************************\n" +
				"xxxx.\n"+

			"\nOptions:\n"+
			"-p Path to the patient annotation txt file.\n"+
			"-c Path to the capture region count data file.\n"+
			"-r File path to a directory to write the results files.\n"+
			"-g Germline sample IDs to exclude from the PoN, comma delimited, no spaces.\n"+
			"-n (Optional) Path to a txt file describing which genes to use in normalizing each sample.\n"+
			"-b Use scanned genes for normalization instead of the best scalar used to find them.\n"+

			"\nExample: java -Xmx100G -jar pathTo/USeq/Apps/xxxxx\n\n" +

				"**************************************************************************************\n");

	}
}

