package edu.utah.kohli;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.utah.kohli.KohliGene.CopyTestResult;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class AccuGenomicsSpikeins {
	
	private TreeMap<String, KohliPatient> patients = new TreeMap<String,KohliPatient>();
	private TreeMap<String, KohliSample> samples = new TreeMap<String,KohliSample>();
	private ArrayList<AccuGenProbe> accuGenProbes = new ArrayList<AccuGenProbe>();
	private ArrayList<CaptureRegion> captureRegions = new ArrayList<CaptureRegion>();
	private LinkedHashMap<String, KohliGene> genes = null;
	private HashMap<String, NormalizerGene> geneNameNormalizedGene = new HashMap<String, NormalizerGene>();
	private File resultsDirectory = null;
	private ArrayList<String> germlineSamplesToRemoveFromPoN = null;
	private LinkedHashMap<String, CopyAnalysisTest> copyAnalysisTests = new LinkedHashMap<String, CopyAnalysisTest>();
	private String testedGenes = "'ARID1A','PTEN','CCND1','ATM','ZBTB16','CDKN1B','KMT2D','CDK4','MDM2','BRCA2','RB1','FOXA1','AKT1','ZFHX3','TP53','NCOR1','CDK12','BRCA1','SPOP','RNF43','MSH2','ERG','TMPRSS2','CHEK2','CTNNB1','FOXP1','PIK3CB','PIK3CA','PIK3R1','CHD1','APC','CDK6','BRAF','KMT2C','NKX3-1','NCOA2','MYC','COL22A1','CDKN2A','NOTCH1','AR-Enh','AR','OPHN1'";
	private File fullPathToR = new File("/usr/local/bin/R");
	private KohliSample testSample = null;
	private File saveDir = null;
	private TreeSet<String> genesToNormWith = null;
	
	//params to pick genes to use in normalizing a test sample
	private double maximumZScore = 3.0;
	private double minimumFractionPassingGeneCaptureRegions = 0.75;
	private double minimumNumberPassingGenes = 3;
	
	public static void main(String[] args) throws IOException {
		File anno = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/sampleMatchingPersonGNomIDs_13June2024.txt");
		File vcfCountData = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/BcfToolsForAccGen/justVcfSampleCounts.txt");
		File regionCountData = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/CopyRatio/Split150bp/medianCounts150bpDoOverPostQC.txt");
		File resultsDirectory = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/CopyRatio/Split150bp/CnvZScoreNormJustTest");
		
		//bad germline samples for PoN, still want for the paired analysis
		ArrayList<String> germlineSamplesToRemoveFromPoN = new ArrayList<String>();
		germlineSamplesToRemoveFromPoN.add("23802X11");
		germlineSamplesToRemoveFromPoN.add("23802X4");
		new AccuGenomicsSpikeins(anno, vcfCountData, regionCountData, resultsDirectory, germlineSamplesToRemoveFromPoN);
	}
	
	public AccuGenomicsSpikeins(File anno, File vcfCountData, File regionCountData, File resultsDirectory, ArrayList<String> germlineSamplesToRemoveFromPoN) throws IOException {
		this.resultsDirectory = resultsDirectory;
		this.germlineSamplesToRemoveFromPoN = germlineSamplesToRemoveFromPoN;

		//load patient sample annotations
		loadPatients(anno);
		//for (KohliPatient kp: patients.values()) IO.pl(kp);

		//load Accugenomic spike-in data
		/*loadAccugenomicSpikeinData(vcfCountData);
		for (KohliSample ks: samples.values()) IO.pl(ks.toStringProbes(accuGenProbes));
		GOAL: flag probes to ignore
		for each probe generate a histogram and boxplot of AFs for the ctDNA and germline
		IO.pl("\nAFs!");
		printAFs();
		for each probe generate a histogram and boxplot of DPs for the ctDNA and germline
		IO.pl("\nDPs!");
		printDPs();
		GOAL: flag samples to ignore
		remove flagged probes
		histogram and boxplot each sample for AF and DPs*/

		//load read coverage for each sample, region, and gene
		loadReadCoverageData(regionCountData);

		//mean scale the germline samples to 1000, needed for the TN comp
		IO.pl("Initial all germline all region normalization to 1000...");
		if (germlineSamplesToRemoveFromPoN.size()!=0) IO.pl("\tExcluding from PoN: "+germlineSamplesToRemoveFromPoN);
		meanScaleGermlineReadCoverageDatasetsAllGenes(1000.0, captureRegions, true);
		//File boxPlotDataGermline = new File(boxPlotDirectory, "germlineBoxPlotData.txt");
		//saveScaledGermlineReadCoverageDatasetsForBoxPlot(boxPlotDataGermline);

		//define the paired t/n datasets
		//matchedTest = samples.get("23802X1");
		//matchedTest = samples.get("22712X1");
		//genesToNormWith = maximizeGenesWithInPoNZScores();
		
		//runCfDNASamples();
		
		runAllGermlinesAsTest();


	
		IO.pl("\nCOMPLETE!");

	}
	
	private void runAllGermlinesAsTest() throws IOException {

		for (KohliSample ks: samples.values()) {
			if (ks.isGermlineSample()==false) continue;
			testSample = ks;

			//exclude the test from the pon
			germlineSamplesToRemoveFromPoN.add(testSample.getSampleId());

			saveDir = new File(resultsDirectory, testSample.getPatient().getHciPatientId()+"_"+testSample.getSampleId());
			saveDir.mkdirs();
			IO.pl("Testing Sample "+saveDir.getName()+"...");
			
			//mean scale the pon minus the "test" sample
			meanScaleGermlineReadCoverageDatasetsAllGenes(1000.0, captureRegions, true);

			//ZScore gene selection method
			genesToNormWith = maximizeGenesWithInPoNZScores();

			//Too few genes?
			if (genesToNormWith.size()==0) {
				IO.pl("ERROR: too few genes found to normalize to, adjust zscore requirements and restart! SKIPPING");
				continue;
			}
			//All on one chrom
			else if (chromosomeCheck(genesToNormWith) == false) {
				IO.pl("ERROR: all genes to normalize were found on one chromosome, adjust zscore requirements and restart! SKIPPING");
				continue;
			}

			// Collect the Capture Regions to use in normalizing
			ArrayList<CaptureRegion> crsToNorm = new ArrayList<CaptureRegion>();
			for (String gn: genesToNormWith) {
				KohliGene kg = genes.get(gn);
				crsToNorm.addAll(kg.getCaptureRegions());
			}

			//ReNormalize the PoN
			IO.pl("Renormalizing the PoN with the select genes...");
			meanScaleGermlineReadCoverageDatasetsAllGenes(1000.0, crsToNorm, false);

			//Normalize the test
			IO.pl("Normalizing the test with the select genes...");
			meanScaleTestReadCoverageDatasetsWithDefinedCaptureRegions(1000.0, crsToNorm);

			//Generate some plots for each gene
			generateBoxPlots("("+Misc.treeSetToString(genesToNormWith, ",")+")");
			generateZScoreHistograms("("+Misc.treeSetToString(genesToNormWith, ",")+")");

			//Stat copy changes
			generateGeneStatistics();
			
			//remove the test from the PoN blocker
			germlineSamplesToRemoveFromPoN.remove(germlineSamplesToRemoveFromPoN.size()-1);
		}
		
	}

	private void runCfDNASamples() throws IOException {
		for (KohliPatient kp: patients.values()) {
			IO.pl("\nTesting Patient "+kp.getHciPatientId()+"...");

			for (KohliSample ks: kp.getNonGermlineSamples()) {
				testSample = ks;

				saveDir = new File(resultsDirectory, kp.getHciPatientId()+"_"+testSample.getSampleId());
				saveDir.mkdirs();
				IO.pl("Testing Sample "+saveDir.getName()+"...");

				
				//ZScore gene selection method
				genesToNormWith = maximizeGenesWithInPoNZScores();
				
				//Too few genes?
				if (genesToNormWith.size()==0) {
					IO.pl("ERROR: too few genes found to normalize to, adjust zscore requirements and restart! SKIPPING");
					continue;
				}
				//All on one chrom
				else if (chromosomeCheck(genesToNormWith) == false) {
					IO.pl("ERROR: all genes to normalize were found on one chromosome, adjust zscore requirements and restart! SKIPPING");
					continue;
				}
	
				//IO.pl(genesToNormWith.size()+" genes for normalization: "+genesToNormWith);
				// Collect the Capture Regions to use in normalizing
				ArrayList<CaptureRegion> crsToNorm = new ArrayList<CaptureRegion>();
				for (String gn: genesToNormWith) {
					KohliGene kg = genes.get(gn);
					crsToNorm.addAll(kg.getCaptureRegions());
				}

				//ReNormalize the PoN
				IO.pl("Renormalizing the PoN with the select genes...");
				meanScaleGermlineReadCoverageDatasetsAllGenes(1000.0, crsToNorm, false);

				//Normalize the test
				IO.pl("Normalizing the test with the select genes...");
				meanScaleTestReadCoverageDatasetsWithDefinedCaptureRegions(1000.0, crsToNorm);

				//Generate some plots for each gene
				generateBoxPlots("("+Misc.treeSetToString(genesToNormWith, ",")+")");
				generateZScoreHistograms("("+Misc.treeSetToString(genesToNormWith, ",")+")");

				//Stat copy changes
				generateGeneStatistics();
				
				//reset the initial PoN values in the CaptureRegions
				IO.pl("Resetting the PoN to the original values...");
				for (CaptureRegion cr: captureRegions) cr.copyInitialPoNToCurrentPon();
				
			}
		}

	}
	
	private boolean chromosomeCheck(TreeSet<String> geneNames) {
		HashSet<String> chrs = new HashSet<String>();
		for (String geneName: geneNames) chrs.add(genes.get(geneName).getChromosome());
		if (chrs.size()>1) return true;
		return false;
	}

	private void generateGeneStatistics() {
		IO.pl("Generating statistic spreadsheet...");
		StringBuilder sb = new StringBuilder();
		sb.append("Patient\tTest Sample\tGermline\tUsed To Norm\tGene\t# Capture Regions\t"
				+ "# CRs Test Outside PoN\tMean CR Z-Score\tMean CR Scaled Test\tMean CR Scaled PoN\tLog2(Test/PoN)\n");
		for (KohliGene kg: genes.values()) {
			sb.append(testSample.getPatient().getHciPatientId()); sb.append("\t");
			sb.append(testSample.getSampleId()); sb.append("\t");
			sb.append(testSample.isGermlineSample()); sb.append("\t");
			sb.append(genesToNormWith.contains(kg.getGeneName()));  sb.append("\t");
			sb.append(kg.getGeneName()); sb.append("\t");
			sb.append(kg.getCaptureRegions().size());  sb.append("\t");
			sb.append(kg.getNumberTestValuesOutsideCaptureRegions()); sb.append("\t");
			sb.append(kg.getMeanCaptureRegionZScore()); sb.append("\t");
			double meanTest = kg.getMeanCaptureRegionScaledTestScore();
			sb.append(meanTest); sb.append("\t");
			double meanPoN = kg.getMeanCaptureRegionPoNScore();
			sb.append(meanPoN); sb.append("\t");
			sb.append(Num.log2(meanTest/meanPoN)); sb.append("\t");
			sb.append("\n");
		}
		
		File spreadsheet = new File( saveDir, saveDir.getName()+".geneStats.xls" );
		IO.writeString(sb.toString(), spreadsheet);
	}

	private void generateBoxPlots(String normGenes) throws IOException {
		IO.pl("Generating boxplots...");
		//write out data files for R
		File rDir = new File(saveDir, "R");
		rDir.mkdir();
		File boxPlotDataGermline = new File(rDir, "germlineBoxPlotData.txt");
		saveScaledGermlineReadCoverageDatasetsForBoxPlot(boxPlotDataGermline);
		File boxPlotDataTest = new File(rDir, "testBoxPlotData.txt");
		saveScaledTestReadCoverageDatasetsForBoxPlot(boxPlotDataTest);
		
		//boxPlot dir
		File boxPlots = new File(saveDir, "BoxPlots");
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

	private TreeSet<String> maximizeGenesWithCenteredCoverage(KohliSample testSample, int[][] upDownSameGermline) {
		IO.pl("Scanning scalars to identify genes in the test that resemble the matched normal...");
		UpDownResult[] res = new UpDownResult[10000 - 100 +1];
		int counter = 0;
		for (int i=100; i<=10000; i++) {
			int[][] upDownSameTest = countGenesWithCenteredCoverage(i, testSample);
			ArrayList<String>[] numMatchsOneOffsTwoOffs = compare(upDownSameGermline, upDownSameTest, false);
			res[counter++] = new UpDownResult(numMatchsOneOffsTwoOffs, i);
		}
		//sort by matches, oneOffs, twoOffs, threeOffs
		Arrays.sort(res);
		IO.pl("\tScalar\tMatches\tOneOffs\tTwoOffs\tThreeOffs");
		IO.pl(res[0]);
		TreeSet<String> toNorm = new TreeSet<String>();
		for (String geneName: res[0].matchesOneTwoThreeOffs[0]) toNorm.add(geneName);
		for (String geneName: res[0].matchesOneTwoThreeOffs[1]) toNorm.add(geneName);
		if (toNorm.size() ==0) {
			IO.pl("WARNING: no genes found to normalize to, using all!");
			toNorm.addAll(genes.keySet());
		}
		else {
			for (int i=1; i< res.length; i++) {
				if (res[0].compareTo(res[i]) != 0) break;
				IO.pl(res[i]);
				//add in matches and one offs? how about two offs?
				for (String geneName: res[i].matchesOneTwoThreeOffs[0]) toNorm.add(geneName);
				for (String geneName: res[i].matchesOneTwoThreeOffs[1]) toNorm.add(geneName);
			}
		}
		return toNorm;
	}
	
	private TreeSet<String> maximizeGenesWithInPoNZScores() {
		IO.pl("Scanning scalars to identify genes in the test that resemble the matched normal...");
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
		}
		return toNorm;
	}
	
	private NormalizationResult countNumberGenesWithTestInPon(int normTarget) {
		//for each gene, count the number of capture regions where the test is within the PoN zscores
		ArrayList<KohliGene> testInPoN = new ArrayList<KohliGene>();
		ArrayList<Double> fracCRs = new ArrayList<Double>();
		for (KohliGene kg : genes.values()) {
			double inPoN = 0;
			for (CaptureRegion cr: kg.getCaptureRegions()) {
				if (Math.abs(cr.calculateZScoreFromPoN()) <= maximumZScore) inPoN++;
			}
			double frac = inPoN/ (double)kg.getCaptureRegions().size();
			if (frac >= minimumFractionPassingGeneCaptureRegions) {
				testInPoN.add(kg);
				fracCRs.add(frac);
			}
		}
		return new NormalizationResult(normTarget, testInPoN, fracCRs);
	}
	
	private class UpDownResult implements Comparable<UpDownResult>{
		
		int normTarget;
		ArrayList<String>[] matchesOneTwoThreeOffs = null;
		int[] numMatchesOneTwoThreeOffs = null;
		
		public UpDownResult(ArrayList<String>[] matchesOneTwoThreeOffs, int normTarget) {
			this.normTarget = normTarget;
			this.matchesOneTwoThreeOffs = matchesOneTwoThreeOffs;
			numMatchesOneTwoThreeOffs = new int[] {matchesOneTwoThreeOffs[0].size(), matchesOneTwoThreeOffs[1].size(), matchesOneTwoThreeOffs[2].size(), matchesOneTwoThreeOffs[3].size()};
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder("\t");
			sb.append(normTarget);
			sb.append("\t");
			sb.append(matchesOneTwoThreeOffs[0]);
			for (int i=1; i< matchesOneTwoThreeOffs.length; i++) {
				sb.append("\t");
				sb.append(matchesOneTwoThreeOffs[i]);
			}
			return sb.toString();
		}

		public int compareTo(UpDownResult o) {
			//matches
			if (this.numMatchesOneTwoThreeOffs[0] > o.numMatchesOneTwoThreeOffs[0]) return -1;
			else if (this.numMatchesOneTwoThreeOffs[0] < o.numMatchesOneTwoThreeOffs[0]) return 1;
			
			//matches same, thus one offs
			if (this.numMatchesOneTwoThreeOffs[1] > o.numMatchesOneTwoThreeOffs[1]) return -1;
			else if (this.numMatchesOneTwoThreeOffs[1] < o.numMatchesOneTwoThreeOffs[1]) return 1;
			
			//one offs same, thus two offs
			if (this.numMatchesOneTwoThreeOffs[2] > o.numMatchesOneTwoThreeOffs[2]) return -1;
			else if (this.numMatchesOneTwoThreeOffs[2] < o.numMatchesOneTwoThreeOffs[2]) return 1;
			
			//two offs same, thus three offs
			if (this.numMatchesOneTwoThreeOffs[3] > o.numMatchesOneTwoThreeOffs[3]) return -1;
			else if (this.numMatchesOneTwoThreeOffs[3] < o.numMatchesOneTwoThreeOffs[3]) return 1;
			
			return 0;
		}
	}
	
	private class NormalizationResult implements Comparable<NormalizationResult>{
		
		int normTarget;
		ArrayList<KohliGene> genesPassingThresholds = null;
		ArrayList<Double> fracCRs = new ArrayList<Double>();
		int numPassingGenes;
		double sumFracCRs;
		
		public NormalizationResult(int normTarget, ArrayList<KohliGene> genesPassingThresholds, ArrayList<Double> fracCRs) {
			this.genesPassingThresholds = genesPassingThresholds;
			this.normTarget = normTarget;
			numPassingGenes = genesPassingThresholds.size();
			this.fracCRs = fracCRs;
			sumFracCRs = Num.sumArray(Num.arrayListOfDoubleToArray(fracCRs));
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
			else {
				if (this.sumFracCRs > o.sumFracCRs) return -1;
				else if (this.sumFracCRs < o.sumFracCRs) return 1;
				return 0;
			}
		}
	}
	
	private ArrayList<String>[] compare(int[][] upDownSameGermline, int[][] upDownSameTest, boolean printStats) {
		ArrayList<String> matches = new ArrayList<String>();
		ArrayList<String> oneOffs = new ArrayList<String>();
		ArrayList<String> twoOffs = new ArrayList<String>();
		ArrayList<String> threeOffs = new ArrayList<String>();
		
		int i=-1;
		for (String geneName: genes.keySet()) {
			i++;
			if (printStats) IO.p(geneName+"\t"+upDownSameGermline[i][0]+":"+upDownSameGermline[i][1]+":"+upDownSameGermline[i][2]+"\t"+upDownSameTest[i][0]+":"+upDownSameTest[i][1]+":"+upDownSameTest[i][2]);
			
			//No matches for zero values in germline
			if (upDownSameGermline[i][0]==0 || upDownSameGermline[i][1]==0) {
				if (printStats) IO.pl("\tZeroInGermlineGeneSkipping");
			}
			else {
				int diffInUp = Math.abs(upDownSameGermline[i][0] - upDownSameTest[i][0]);
				if (diffInUp == 0) {
					if (printStats) IO.pl("\tMatch");
					matches.add(geneName);
				}
				else if (diffInUp == 1) {
					if (printStats) IO.pl("\tOneOff\t"+(upDownSameGermline[i][0] - upDownSameTest[i][0]));
					oneOffs.add(geneName);
				}
				else if (diffInUp == 2) {
					if (printStats) IO.pl("\tTwoOff\t"+(upDownSameGermline[i][0] - upDownSameTest[i][0]));
					twoOffs.add(geneName);
				}
				else if (diffInUp == 3) {
					if (printStats) IO.pl("\tThreeOff\t"+(upDownSameGermline[i][0] - upDownSameTest[i][0]));
					threeOffs.add(geneName);
				}
				else if (printStats) IO.pl("\tNoMatch");
			}
		}
		if (printStats)IO.pl(matches+ "\t"+ oneOffs+"\t"+twoOffs+"\t"+threeOffs);
		return new ArrayList[] {matches, oneOffs, twoOffs, threeOffs};
	}

	private int[][] countGenesWithCenteredCoverage(double target, KohliSample ks) {
		//scale the counts using all of the capture regions
		ks.createScaledRegionCountMap(captureRegions, captureRegions, target);
		HashMap<String, Double> sampleScaledCaptureRegionMeans = ks.getScaledReadCoverageCounts();
		
		//for each gene count the above and below means from the PoN
		int[][] upDown = new int[genes.size()][];
		int counter = 0;
		for (KohliGene kg : genes.values()) {
			int numAboveMean = 0;
			int numBelowMean = 0;
			int numSame = 0; //never happens!
			for (CaptureRegion cr: kg.getCaptureRegions()) {
				String crId = cr.getOriginalInput();
				double testScaledCount = sampleScaledCaptureRegionMeans.get(crId);
				double ponScaledCount = cr.getPonScaledMean();
				if (testScaledCount > ponScaledCount) numAboveMean++;
				else if (testScaledCount < ponScaledCount) numBelowMean++;
				else numSame++;
			}
			upDown[counter++] = new int[]{numAboveMean, numBelowMean, numSame};
		}
		return upDown;
	}

	private void testOneGeneMethod(String sampleIdToTest) {
		
		CopyAnalysisTest cat = new CopyAnalysisTest(sampleIdToTest);
		copyAnalysisTests.put(sampleIdToTest, cat);
		
		KohliSample ks = samples.get(sampleIdToTest);
		
		//@TODO: print out min gene norm box plots and max gene norm box plots
		//vis inspect to see oddity
		
		//use each gene, use it as a normalizer
		for (KohliGene kg: genes.values()) {
			NormalizerGene ng = geneNameNormalizedGene.get(kg.getGeneName());
			HashMap<String, double[]> geneMinMaxZScores = ng.getGeneNameZScoreMinMax();

			//scale the germline data and load it into the capture regions
			//@TODO really should cache this!
			meanScaleGermlineReadCoverageDatasetsOneGene(1000, kg.getGeneName());
			
			//scale the test
			meanScaleTestReadCoverageDatasetsOneGene(1000, kg.getGeneName(), sampleIdToTest);
File boxPlotDataTest = new File(resultsDirectory, "testing_"+sampleIdToTest+"_"+kg.getGeneName()+"_BoxPlotData.txt");
saveScaledTestReadCoverageDatasetsForBoxPlot(boxPlotDataTest);
			
			//for each gene compare
			int numOutOfRange = 0;
			
			//@TODO: skip TMPRSS2-ERG
			
			for (KohliGene kgs: genes.values()) {
				CopyTestResult testResult = kgs.compareTestvsPoN();
				//is the meanZScore within the germline test zscores?
				double[] minMax = geneMinMaxZScores.get(kgs.getGeneName());
				double testMeanZScore = testResult.getMeanZScore();
				boolean inGermZScoreRange = false;
				if (testMeanZScore < minMax[0] || testMeanZScore > minMax[1]) {
					inGermZScoreRange = true;
					numOutOfRange++;
				}
				//IO.pl(sampleIdToTest+"\t"+kg.getGeneName()+"\t"+kgs.getGeneName()+"\t"+testResult.getMeanZScore()+"\t"+inGermZScoreRange);
			}
			IO.pl(sampleIdToTest+"\t"+ kg.getGeneName()+ "\t"+ numOutOfRange);
			
			//@TODO, convert to show sorted list of normalizer genes
			/*if (numOutOfRange < minNumOutOfRange) {
				minNumOutOfRange = numOutOfRange;
				bestNormalizer = sampleIdToTest+"\t"+ kg.getGeneName();
			}*/
		}
		//IO.pl(bestNormalizer);
			
		
	}

	private void generatePonBackground() {
		IO.pl("Generating all gene normalized PoN...");
		//for each gene
		for (KohliGene kg: genes.values()) {
			String normalizerGeneName = kg.getGeneName(); 
			NormalizerGene ng = new NormalizerGene(normalizerGeneName);
			geneNameNormalizedGene.put(normalizerGeneName, ng);

			//scale the PoN, not necessary, just useful for boxplots, not used in zscore calcs
			meanScaleGermlineReadCoverageDatasetsOneGene(1000, normalizerGeneName);
			File boxPlotDataGermline = new File(resultsDirectory, normalizerGeneName+"_All_Germline_BoxPlotData.txt");
			saveScaledGermlineReadCoverageDatasetsForBoxPlot(boxPlotDataGermline);
			
			//for each germline sample that isn't flagged, use it as a mock test sample
			ArrayList<String> testSamples = new ArrayList<String>();
			for (KohliSample ks: samples.values()) {
				
				//skip any non germline samples and any that have failed a prior determination to be bad
				if (ks.isGermlineSample() == false || ks.isExcludeFromPon() == true) continue;
				
				//set exclude the current test sample from pon generation
				ks.setExcludeFromPon(true);
				//scale the PoN without the mock test sample
				meanScaleGermlineReadCoverageDatasetsOneGene(1000, normalizerGeneName);
				//File boxPlotDataGermline = new File(boxPlotDirectory, normalizerGeneName+"_Minus_"+ks.getSampleId()+"_Germline_BoxPlotData.txt");
				//saveScaledGermlineReadCoverageDatasetsForBoxPlot(boxPlotDataGermline);
				//reset to include for the next round
				ks.setExcludeFromPon(false);

				//scale the test
				String testSampleId = ks.getSampleId(); 
				testSamples.add(testSampleId);
				meanScaleTestReadCoverageDatasetsOneGene(1000, normalizerGeneName, testSampleId);
				//File boxPlotDataTest = new File(boxPlotDirectory, normalizerGeneName+"_"+testSampleId+"_BoxPlotData.txt");
				//saveScaledTestReadCoverageDatasetsForBoxPlot(boxPlotDataTest);
				
				//compare the test to the pon
				for (KohliGene kgs: genes.values()) {
					double compScore = kgs.compareTestvsPoN().getMeanZScore();
					//IO.pl(kgs.getGeneName()+"/t"+compScore);
					ng.addZScore(kgs.getGeneName(), compScore);
				}
			}
			ng.setGermlineSampleNames(testSamples);
			//IO.pl(ng.toString());
		}
		//IO.pl("Calculate the mean and stdev of these zscores for mock germline sample against the PoN to identify and remove outliers.");
		if (germlineSamplesToRemoveFromPoN.size()!=0) IO.pl("\tExcluding user flagged "+Misc.stringArrayListToString(germlineSamplesToRemoveFromPoN, ", "));
	}
	
	private void generateZScoreHistograms(String normGenes) throws IOException {
		IO.pl("Generating histograms...");
		
		//write out data files for R
		File rDir = new File(saveDir, "R");
		rDir.mkdir();
		File histogramData = new File(rDir, "histogramData.txt");
		saveTestZScoresForHistograms(histogramData);
		
		//Histo dir
		File histogramPlots = new File(saveDir, "HistogramPlots");
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
		for (KohliSample ks: samples.values()) if (ks.isGermlineSample()) {
			ks.createScaledRegionCountMap(crsToScaleWith, captureRegions, target);
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
	
	private void meanScaleGermlineReadCoverageDatasetsOneGene(double target, String gene) {

		//scale the counts using one gene to the target value
		KohliGene kg = genes.get(gene);
		ArrayList<CaptureRegion> forScalarRegions = kg.getCaptureRegions();
		for (KohliSample ks: samples.values()) {
			if (ks.isGermlineSample() && ks.isExcludeFromPon()==false) ks.createScaledRegionCountMap(forScalarRegions, captureRegions, target);
		}
	
		//load the capture regions with the scaled counts
		for (CaptureRegion cr: captureRegions) {
			String crId = cr.getOriginalInput();
			ArrayList<Double> germlineScaledCounts = new ArrayList<Double>();
			
			//for each germline pon sample, exclude those that have been flagged
			for (KohliSample ks: samples.values()) {
				if (ks.isGermlineSample()  && ks.isExcludeFromPon()==false) {
					HashMap<String, Double> sampleScaledCaptureRegionMeans = ks.getScaledReadCoverageCounts();
					germlineScaledCounts.add(sampleScaledCaptureRegionMeans.get(crId));
				}
			}
			//convert to double[] and set in cr
			cr.setPonScaledCountsCalculateStats(Num.arrayListOfDoubleToArray(germlineScaledCounts), false);
		}
	}
	
	private void meanScaleTestReadCoverageDatasetsOneGene(double target, String gene, String testSampleId) {

		//scale the counts using one gene to the target value
//IO.pl("Scaling germline sample counts using gene "+gene+" to target "+target);
		KohliGene kg = genes.get(gene);
		KohliSample ks = samples.get(testSampleId);
		ks.createScaledRegionCountMap(kg.getCaptureRegions(), captureRegions, target);
		
		//load the capture regions with the scaled counts
//IO.pl("Pulling counts over each capture region for the test sample....");
		for (CaptureRegion cr: captureRegions) {
			String crId = cr.getOriginalInput();
			HashMap<String, Double> sampleScaledCaptureRegionMeans = ks.getScaledReadCoverageCounts();
			cr.setScaledTestCount(sampleScaledCaptureRegionMeans.get(crId));
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
	
	private void meanScaleAndPrintReadCoverageDatasets(double target) {
		//print header row
		StringBuilder sb = new StringBuilder("CaptureRegions");
		//for each sample
		for (KohliSample ks: samples.values()) {
			sb.append("\t");
			sb.append(ks.getSampleId());
			//if (ks.getSampleId().equals("22712X10")) IO.pl("\n"+Num.doubleArrayToString(ks.getScaledRegionCountArray(captureRegions, target), ","));
		}
		sb.append("\n");

		//add the rows
		//for each region
		String gene = "";
		int counter = 0;
		for (CaptureRegion region: captureRegions) {
			//new gene name?
			if (region.getGene().equals(gene)==false) {
				gene = region.getGene();
				counter = 0;
			}
			else counter++;
			//gene id
			sb.append(gene); sb.append("_"); sb.append(counter);
			//for each sample
			for (KohliSample ks: samples.values()) {
				sb.append("\t");
				ks.createScaledRegionCountMap(captureRegions, captureRegions, target);
				sb.append(ks.getScaledReadCoverageCounts().get(region.getOriginalInput()));
			}
			sb.append("\n");
		}
		IO.pl(sb);
	}
	private void meanScaleReadCoverageDatasetsBoxPlot(double target, boolean germline) {
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
IO.pl(cr.getOriginalInput()+"\t"+name);
			//for each sample 
			for (KohliSample ks: samples.values()) {
				//if (ks.getType().equals(type) == false) continue;
				if ((ks.isGermlineSample() == germline) == false) continue;
				sb.append(name);
				sb.append("\t");
				ks.createScaledRegionCountMap(captureRegions, captureRegions, target);
				sb.append(ks.getScaledReadCoverageCounts().get(crId));
				sb.append("\n");
			}
		}
IO.pl(sb);
	}

	private void printDPs() {
		//for each probe in rows, for each sample in columns, print the AFs
		//header row
		StringBuilder sb = new StringBuilder();
		sb.append("ProbId");
		for (String sampleId: samples.keySet()) {
			KohliSample ks = samples.get(sampleId);
			String type = "germline";
			if (ks.isGermlineSample()==false) type = "tumor";
			sb.append("\t");
			sb.append(sampleId+"_"+type);
		}
		sb.append("\n");
		//data lines
		for (AccuGenProbe p: accuGenProbes) {
			sb.append(p.getOriginalInput());
			//for each sample
			for (KohliSample ks: samples.values()) {
				sb.append("\t");
				sb.append(ks.getAccuGenProbeCount(p.getOriginalInput()).getReadDepth());
			}
			sb.append("\n");
		}
		IO.p(sb);
	}

	private void printAFs() {
		//for each probe in rows, for each sample in columns, print the AFs
		//header row
		StringBuilder sb = new StringBuilder();
		sb.append("ProbId");
		for (String sampleId: samples.keySet()) {
			KohliSample ks = samples.get(sampleId);
			String type = "germline";
			if (ks.isGermlineSample()==false) type = "tumor";
			sb.append("\t");
			sb.append(sampleId+"_"+type);
		}
		sb.append("\n");
		//data lines
		for (AccuGenProbe p: accuGenProbes) {
			sb.append(p.getOriginalInput());
			//for each sample
			for (KohliSample ks: samples.values()) {
				sb.append("\t");
				sb.append(ks.getAccuGenProbeCount(p.getOriginalInput()).getAlleleFraction());
			}
			sb.append("\n");
		}
		IO.p(sb);
	}
	
	private void loadReadCoverageData(File countData) {
		IO.pl("Loading capture region count data...");
		String[] lines = IO.loadFileIntoStringArray(countData);
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
		
		
	}

	private void loadAccugenomicSpikeinData(File vcfCountData) {
		String[] lines = IO.loadFileIntoStringArray(vcfCountData);

		// parse the header line, create a hashMap of index and KohliSample
		// Chr_POS_REF_ALT_GENE 22712X10 22712X11 22712X12 22712X13 ...
		if (lines[0].startsWith("Chr_POS_REF_ALT_GENE") == false) Misc.printErrAndExit("First line in the vcf count data file doesn't start with 'Chr_POS_REF_ALT_GENE'? -> "+lines[0]);
		String[] header = Misc.TAB.split(lines[0]);
		KohliSample[] indexSamples = new KohliSample[header.length];
		for (int i=1; i< header.length; i++ ) {
			KohliSample ks = samples.get(header[i]);
			if (ks == null) Misc.printErrAndExit("Failed to find the sample "+header[i]+ "in "+samples);
			indexSamples[i] = ks;
		}
		
		// parse the data lines for each data line
		for (int i=1; i< lines.length; i++) {
			//chr11_114064031_A_T_ZBTB16	0,255,255:3150,152	0,99,255:8440,685	0,136,255:3177,188
			if (lines[i].startsWith("chr") == false) Misc.printErrAndExit("Data line doesn't start with 'chr'? "+lines[i]);
			String[] f = Misc.TAB.split(lines[i]);
			//use the first field to create a probe
			AccuGenProbe p = new AccuGenProbe(f[0]);
			accuGenProbes.add(p);

			//for each data point
			for (int j=1; j< f.length; j++) {
				//get the sample, add the values
				KohliSample ks = indexSamples[j];
				ks.addAccuGenProbeCounts(p, f[j]);
			}
		}
	}

	private void loadPatients(File anno) {
		IO.pl("Loading patient and sample meta data...");
		String[] lines = IO.loadFileIntoStringArray(anno);
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

}
