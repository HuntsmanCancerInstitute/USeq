package edu.utah.kohli;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.TreeMap;

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
	private File boxPlotDirectory = null;
	private ArrayList<String> germlineSamplesToRemoveFromPoN = null;
	private LinkedHashMap<String, CopyAnalysisTest> copyAnalysisTests = new LinkedHashMap<String, CopyAnalysisTest>();

	public static void main(String[] args) {
		File anno = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/sampleMatchingPersonGNomIDs_13June2024.txt");
		File vcfCountData = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/BcfToolsForAccGen/justVcfSampleCounts.txt");
		File regionCountData = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/CopyRatio/Split150bp/medianMergedSplit150ProbeCountsQCed.txt");
		File boxPlotDirectory = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/CopyRatio/Split150bp/BoxPlotDataFiles");
		ArrayList<String> germlineSamplesToRemoveFromPoN = new ArrayList<String>();
		germlineSamplesToRemoveFromPoN.add("23802X11");
		new AccuGenomicsSpikeins(anno, vcfCountData, regionCountData, boxPlotDirectory, germlineSamplesToRemoveFromPoN);
	}
	
	public AccuGenomicsSpikeins(File anno, File vcfCountData, File regionCountData, File boxPlotDirectory, ArrayList<String> germlineSamplesToRemoveFromPoN) {
		this.boxPlotDirectory = boxPlotDirectory;
		this.germlineSamplesToRemoveFromPoN = germlineSamplesToRemoveFromPoN;
		
		//load patient sample annotations
		loadPatients(anno);
		//for (KohliPatient kp: patients.values()) IO.pl(kp);
		
		//load Accugenomic spike-in data
		//loadAccugenomicSpikeinData(vcfCountData);
		//for (KohliSample ks: samples.values()) IO.pl(ks.toStringProbes(accuGenProbes));
		
		//GOAL: flag probes to ignore
		//for each probe generate a histogram and boxplot of AFs for the ctDNA and germline
		//IO.pl("\nAFs!");
		//printAFs();
		//for each probe generate a histogram and boxplot of DPs for the ctDNA and germline
		//IO.pl("\nDPs!");
		//printDPs();
		
		//GOAL: flag samples to ignore
		//remove flagged probes
		//histogram and boxplot each sample for AF and DPs
		
		//load read coverage for each sample, region, and gene
		loadReadCoverageData(regionCountData);

		//for (KohliSample ks: samples.values()) IO.pl(ks.meanRegionCounts(captureRegions));
		//meanScaleAndPrintReadCoverageDatasets(1000.0);
		//meanScaleReadCoverageDatasetsBoxPlot(1000.0,"germline");
		
		//mean scale the germline samples
		//meanScaleGermlineReadCoverageDatasetsAllGenes(1000.0);
		
		//print out the genes
		//IO.pl("Genes, their capture regions, the mean of all of the germline samples and the individual values");
		//for (KohliGene kg: genes.values()) {
		//	IO.pl(kg.toString());
		//}
		
		
		//scan the pon to define the accepted range for each normalizer gene
		generatePonBackground();
		
		//test a sample against the PoN
		test("22712X1");
		//test("23802X1");
		//for (String sampleId : samples.keySet()) {
			//test(sampleId);
		//}
		
	}
	
	private void test(String sampleIdToTest) {
		
//IO.pl("Testing, "+sampleIdToTest+"...");
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
File boxPlotDataTest = new File(boxPlotDirectory, "testing_"+sampleIdToTest+"_"+kg.getGeneName()+"_BoxPlotData.txt");
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
			File boxPlotDataGermline = new File(boxPlotDirectory, normalizerGeneName+"_All_Germline_BoxPlotData.txt");
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
			if (counterString.length()==1) counterString = "0"+counterString;
			String name = gene+"_"+counterString;
			sb.append(name);
			sb.append("\t");
			sb.append(Double.toString(cr.getScaledTestCount()));
			sb.append("\n");
		}
		IO.writeString(sb.toString(), toSave);
	}

	private void meanScaleGermlineReadCoverageDatasetsAllGenes(double target) {

		//scale the counts using all of the capture regions
//IO.pl("Scaling germline sample counts....");
		for (KohliSample ks: samples.values()) {
			if (ks.isGermlineSample()) ks.createScaledRegionCountMap(captureRegions, captureRegions, target);
		}
		
		//for each capture region
//IO.pl("Pulling counts over each capture region for each germline sample....");
		for (CaptureRegion cr: captureRegions) {
			String crId = cr.getOriginalInput();
			ArrayList<Double> germlineScaledCounts = new ArrayList<Double>();
			//for each germline pon sample
			for (KohliSample ks: samples.values()) {
				if (ks.isGermlineSample()) {
					HashMap<String, Double> sampleScaledCaptureRegionMeans = ks.getScaledReadCoverageCounts();
					germlineScaledCounts.add(sampleScaledCaptureRegionMeans.get(crId));
				}
			}
			//convert to double[] and set in cr
//IO.pl("Setting counts in "+cr.getOriginalInput());
			cr.setPonScaledCountsCalculateStats(Num.arrayListOfDoubleToArray(germlineScaledCounts));
//IO.pl("\t"+cr.getPonScaledMean());
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
			cr.setPonScaledCountsCalculateStats(Num.arrayListOfDoubleToArray(germlineScaledCounts));
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
			if (counterString.length()==1) counterString = "0"+counterString;
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
			if (counterString.length()==1) counterString = "0"+counterString;
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
				kp.getSamples().add(new KohliSample(tokens[1], isGermline, tokens[3]));
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
