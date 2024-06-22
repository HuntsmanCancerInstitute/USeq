package edu.utah.kohli;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.TreeMap;

import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class AccuGenomicsSpikeins {
	
	private TreeMap<String, KohliPatient> patients = new TreeMap<String,KohliPatient>();
	private TreeMap<String, KohliSample> samples = new TreeMap<String,KohliSample>();
	private ArrayList<AccuGenProbe> accuGenProbes = new ArrayList<AccuGenProbe>();
	private ArrayList<CaptureRegion> captureRegions = new ArrayList<CaptureRegion>();
	private LinkedHashMap<String, KohliGene> genes = null;

	public static void main(String[] args) {
		File anno = new File ("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/sampleMatchingPersonGNomIDs_13June2024.txt");
		File vcfCountData = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/BcfToolsForAccGen/justVcfSampleCounts.txt");
		File regionCountData = new File("/Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/CopyRatio/Split150bp/medianMergedSplit150ProbeCountsQCed.txt");
		
		new AccuGenomicsSpikeins(anno, vcfCountData, regionCountData);
	}
	
	public AccuGenomicsSpikeins(File anno, File vcfCountData, File regionCountData) {
		
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
		meanScaleGermlineReadCoverageDatasets(1000.0);
		
		//print out the genes
		//IO.pl("Genes, their capture regions, the mean of all of the germline samples and the individual values");
		//for (KohliGene kg: genes.values()) {
		//	IO.pl(kg.toString());
		//}
		
		//pick a sample to test and a gene to normalize to
		String testSampleId = "23802X10";  //one of the germline samples used in the boxplot
		String normalizerGeneName = "TP53";
		String boxplotData = scaleTestCountsUsingGene(testSampleId, normalizerGeneName);
		IO.pl(boxplotData);
		
		//compare
		for (KohliGene kg: genes.values()) {
			double compScore = kg.compareTestvsPoN();
			IO.pl(kg.getGeneName()+"\t"+ kg.compareTestvsPoN());
			
IO.pl("Here");
			
		}
		
		
		//scan the pon to define the accepted range for each gene
	}
	
	private String scaleTestCountsUsingGene(String testSampleId, String normalizerGeneName) {
		//get the testSample
		KohliSample testSample = samples.get(testSampleId);
		//add in raw test values to all of the CaptureRegions
		for (CaptureRegion cr: captureRegions) {
			String crId = cr.getOriginalInput();
			double count = testSample.getReadCoverageCounts().get(crId);
			cr.setRawTestCount(count);
		}
		//use the normalizer gene to calculate a scalar to normalize the test raw counts
		KohliGene normalizerGene = genes.get(normalizerGeneName);
		double[] ponGeneCounts = normalizerGene.getScaledPonCounts();
		//get the testCounts for that gene
		double[] testGeneCounts = normalizerGene.getTestCounts();
		//calculate the scalar
		double scalar = Num.calculateScalar(ponGeneCounts, testGeneCounts);
		//scale all of the test counts
		for (CaptureRegion cr: captureRegions) {
			double rawCount = cr.getRawTestCount();
			cr.setScaledTestCount(rawCount*scalar);
		}
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
		return sb.toString();
	}

	private void meanScaleGermlineReadCoverageDatasets(double target) {

		//scale the counts
//IO.pl("Scaling germline sample counts....");
		for (KohliSample ks: samples.values()) {
			if (ks.getType().equals("germline")) ks.getScaledRegionCountMap(captureRegions, target);
		}
		
		//for each capture region
//IO.pl("Pulling counts over each capture region for each germline sample....");
		for (CaptureRegion cr: captureRegions) {
			String crId = cr.getOriginalInput();
			ArrayList<Double> germlineScaledCounts = new ArrayList<Double>();
			//for each germline pon sample
			for (KohliSample ks: samples.values()) {
				if (ks.getType().equals("germline")) {
					HashMap<String, Double> sampleScaledCaptureRegionMeans = ks.getScaledRegionCountMap(captureRegions, target);
					germlineScaledCounts.add(sampleScaledCaptureRegionMeans.get(crId));
				}
			}
			//convert to double[] and set in cr
//IO.pl("Setting counts in "+cr.getOriginalInput());
			cr.setPonScaledCounts(Num.arrayListOfDoubleToArray(germlineScaledCounts));
//IO.pl("\t"+cr.getPonScaledMean());
		}
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
				sb.append(ks.getScaledRegionCountMap(captureRegions, target).get(region.getOriginalInput()));
			}
			sb.append("\n");
		}
		IO.pl(sb);
	}
	
	private void meanScaleReadCoverageDatasetsBoxPlot(double target, String type) {
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
				if (ks.getType().equals(type) == false) continue;
				sb.append(name);
				sb.append("\t");
				sb.append(ks.getScaledRegionCountMap(captureRegions, target).get(crId));
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
			sb.append("\t");
			sb.append(sampleId+"_"+ks.getType());
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
			sb.append("\t");
			sb.append(sampleId+"_"+ks.getType());
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
	
	private void loadReadCoverageData(File vcfCountData) {
		String[] lines = IO.loadFileIntoStringArray(vcfCountData);
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
		String[] lines = IO.loadFileIntoStringArray(anno);
		if (lines[0].equals("#HCIPersonID\tSampleID\tType\tDateDrawn") == false) Misc.printErrAndExit("First line in the anno file isn't '#HCIPersonID SampleID Type	DateDrawn'? -> "+lines[0]);
		for (int i=1; i< lines.length; i++) {
			String[] tokens = Misc.TAB.split(lines[i]);
			KohliPatient kp = patients.get(tokens[0]);
			if (kp == null) {
				kp = new KohliPatient(tokens);
				patients.put(tokens[0], kp);
			}
			else kp.getSamples().add(new KohliSample(tokens[1], tokens[2], tokens[3]));
		}
		//load the sample treemap
		for (KohliPatient kp: patients.values()) {
			for (KohliSample ks: kp.getSamples()) {
				if (samples.containsKey(ks.getSampleId())) Misc.printErrAndExit("Duplicate sample found! "+ks.getSampleId());
				samples.put(ks.getSampleId(), ks);
			}
		}
	}



}
