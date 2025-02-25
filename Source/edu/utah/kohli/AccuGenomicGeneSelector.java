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
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class AccuGenomicGeneSelector {
	
	//user input
	private File patientAnnotationFile = null;
	private File vcfCountData = null;
	private File resultsDirectory = null;
	private File fullPathToR = new File("/usr/local/bin/R");
	private String probesToSkip= "";
	int numGenesToSumForBestScalar = 2;
	
	
	//internal 
	private double maxAFForBoxplot = 0.25;
	private String testedGenes = null;
	private TreeMap<String, KohliPatient> patients = new TreeMap<String,KohliPatient>();
	private TreeMap<String, KohliSample> samples = new TreeMap<String,KohliSample>();
	private ArrayList<AccuGenProbe> accuGenProbes = new ArrayList<AccuGenProbe>();
	private TreeMap<String, AccuGenGene> genes = null;
	private ArrayList<String> germlineSamplesToRemoveFromPoN = null;
	private KohliSample testSample = null;
	private File workingSaveDir = null;
	private int numberGermlineSamples = 0;
	private ArrayList<String> bestScaleResultLines = new ArrayList<String>();
	
	public AccuGenomicGeneSelector(String[] args) throws Exception {
		processArgs(args);


		//load patient sample annotations
		loadPatients();

		//load Accugenomic spike-in data
		loadAccugenomicSpikeinData(vcfCountData);
		
		// calculates mean AFs for each sample, then the mean of the means, then scales all of the sample afs to the mean
		normalizeGermlineAFs();
		loadGenesWithProbes();
		
		//for (KohliSample ks: samples.values()) IO.pl(ks.toStringProbes(accuGenProbes));
		//IO.pl("\nAFs!");
		//printAFs();
		//IO.pl("\nDPs!");
		//printDPs();
		
		runCfDNASamples();  
		
		//testScalars = "1.0";
		runAllGermlinesAsTest();

		IO.pl("\nCOMPLETE!");

	}
	
	private void runAllGermlinesAsTest() throws IOException {
		IO.pl("\nRunning germline samples as tests...");
		for (KohliSample ks: samples.values()) {
			if (ks.isGermlineSample()==false) continue;
			testSample = ks;

			//exclude the test from the pon
			testSample.setExcludeFromPon(true);
			germlineSamplesToRemoveFromPoN.add(testSample.getSampleId());
			
			//write out AFs from just those not excluded for box plotting
			File boxPlotDataGermline = new File(resultsDirectory, "germlineAF_PoN_BoxPlotData.txt");
			saveGermlineAFsForBoxPlot(boxPlotDataGermline);

			workingSaveDir = new File(resultsDirectory, testSample.getPatient().getHciPatientId()+"_"+testSample.getSampleId());
			workingSaveDir.mkdirs();
			IO.pl("Testing Sample "+workingSaveDir.getName()+"...");
			
			//Generate some plots for each gene
			generateBoxPlots(boxPlotDataGermline);
			
			//Stat copy changes
			//generateGeneStatistics();
			
			//remove the test from the PoN blocker
			testSample.setExcludeFromPon(false);
			germlineSamplesToRemoveFromPoN.remove(germlineSamplesToRemoveFromPoN.size()-1);
		}
	}

	private void runCfDNASamples() throws IOException {
		IO.pl("\nRunning cfDNA samples...");

		//write out scaled AFs from just those not excluded for box plotting
		File boxPlotDataGermline = new File(resultsDirectory, "germlineAF_PoN_BoxPlotData.txt");
		saveGermlineAFsForBoxPlot(boxPlotDataGermline);
		
		//for each patient
		for (KohliPatient kp: patients.values()) {
			IO.pl("\nTesting Patient "+kp.getHciPatientId()+"...");

			//for each sample
			for (KohliSample ks: kp.getNonGermlineSamples()) {
				testSample = ks;
				String id = ks.getSampleId();
				
				scaleSampleAfs(false);
				
				//"22712X61") || id.equals("22712X62") || 
				//if (id.equals("22712X63") ) scaleSampleAfs(true);

				///if (id.equals("22712X28") || id.equals("22712X9") || id.equals("22712X16") ) {
					workingSaveDir = new File(resultsDirectory, kp.getHciPatientId()+"_"+testSample.getSampleId());
					workingSaveDir.mkdirs();
					IO.pl("Testing Sample "+workingSaveDir.getName()+"...");

					//Generate some plots for each gene
					generateBoxPlots(boxPlotDataGermline);

					//Stat copy changes
					//generateGeneStatistics();
				//}

				
			}
		}
		
		//print the best scalar data
		IO.pl("\nBest scalar info, inspect for outliers (cfDNAId SumAbsZScores BestScalar Gene:MeanAbsZScore etc)...");
		IO.pl(Misc.stringArrayListToString(bestScaleResultLines, "\n"));
	}
	
	
	private void scaleSampleAfs(boolean verbose) {
		
		ScalarResult bestResult = null;
		double sumBest = 0;
		
		for (double scalar=0.1; scalar< 3; scalar+=0.001) {
			ScalarResult sr = new ScalarResult(scalar);

			//for each gene
			for (String geneName: genes.keySet()) {
				AccuGenGene g = genes.get(geneName);
				//for each probe in the gene
				double totalAbsZScore = 0;
				for (AccuGenProbe p: g.getProbes()) {
					double testAf = scalar * testSample.getAccuGenProbeCount(p.getOriginalInput()).getAlleleFraction();
					double zscore = p.calculateZScore(testAf);
					totalAbsZScore+= Math.abs(zscore);
				}
				double meanGeneAbsZscore = totalAbsZScore/ (double) g.getProbes().size();
				sr.geneAbsZScore.add(new GeneScore(geneName, meanGeneAbsZscore));
			}

			if (verbose) IO.pl(testSample.getSampleId()+"\t"+sr.toString());
			if (bestResult == null) {
				bestResult = sr;
				sumBest = bestResult.getSumLowestZScores(numGenesToSumForBestScalar);
			}
			else {
				double sumLowest = sr.getSumLowestZScores(numGenesToSumForBestScalar);
				if (sumLowest< sumBest) {
					bestResult = sr;
					sumBest = sumLowest;
				}
			}
		}
		bestScaleResultLines.add(testSample.getSampleId()+"\t"+sumBest+"\t"+bestResult.toString());
		
		//scale the scores
		for (AccuGenProbe p: accuGenProbes) {
			AccuGenProbeCounts c = testSample.getAccuGenProbeCount(p.getOriginalInput());
			double af = c.getAlleleFraction();
			c.setScaledAlleleFraction(af * bestResult.scalar);
		}

	}
	
	private class ScalarResult {
		
		private double scalar;
		private ArrayList<GeneScore> geneAbsZScore = new ArrayList<GeneScore>();
		GeneScore[] sortedGeneAbsZScore = null;
		
		private ScalarResult (double scalar) {
			this.scalar = scalar;
		}
		
		public double getSumLowestZScores(int numToSum) {
			if (sortedGeneAbsZScore == null) {
				sortedGeneAbsZScore = new GeneScore[geneAbsZScore.size()];
				geneAbsZScore.toArray(sortedGeneAbsZScore);
				Arrays.sort(sortedGeneAbsZScore);
			}
			double sum = 0;
			for (int i=0; i< numToSum; i++) sum+=sortedGeneAbsZScore[i].score;
			return sum;
		}

		public String toString() {
			sortedGeneAbsZScore = new GeneScore[geneAbsZScore.size()];
			geneAbsZScore.toArray(sortedGeneAbsZScore);
			Arrays.sort(sortedGeneAbsZScore);
			StringBuilder sb = new StringBuilder();
			sb.append(scalar); 
			for (GeneScore pairs: sortedGeneAbsZScore) {
				sb.append("\t");
				sb.append(pairs.toString());
			}
			return sb.toString();
		}
	}
	
	private class GeneScore implements Comparable<GeneScore>{
		private String geneName = null;
		private double score = 0;
		
		public GeneScore(String geneName, double score) {
			this.geneName = geneName;
			this.score = score;
		}

		public int compareTo(GeneScore o) {
			if (o.score > this.score) return -1;
			if (o.score < this.score) return 1;
			return 0;
		}
		
		public String toString() {
			return geneName+":"+score;
		}
	}

	private void generateBoxPlots(File germline) throws IOException {
		IO.pl("\tGenerating boxplots...");
		
		//write out test data files for R
		File rDir = new File(workingSaveDir, "R");
		rDir.mkdir();
		File boxPlotDataTest = new File(rDir, testSample.getSampleId()+"_PoN_BoxPlotData.txt");
		saveTestAFsForBoxPlot(boxPlotDataTest);
		
		//boxPlot dir 
		File boxPlots = new File(workingSaveDir, "BoxPlots");
		boxPlots.mkdir();
		
		//
		
		ArrayList<String> rScriptAL = new ArrayList<String>();
		rScriptAL.add("library(ggplot2)");
		rScriptAL.add("setwd('"+boxPlots.getCanonicalPath()+"')");
		rScriptAL.add("df = read.table('"+germline.getCanonicalPath()+"', header=TRUE)");
		rScriptAL.add("testDf = read.table('"+boxPlotDataTest.getCanonicalPath()+"', header=TRUE)");
		rScriptAL.add("testDataset = '"+testSample.getSampleId()+"'");
		rScriptAL.add("genes = c("+testedGenes+")");
		rScriptAL.add("for (x in genes) {");
		rScriptAL.add("  print(x)");
		rScriptAL.add("  toSub=paste(x,'_',sep='')");
		rScriptAL.add("  dfSub = subset(df, grepl(toSub, Gene))");
		rScriptAL.add("    dfSubTest = subset(testDf, grepl(toSub, Gene))");
		rScriptAL.add("    boxPlotObj = ggplot(dfSub, aes(x=Gene, y=AlleleFraction)) + ");
		rScriptAL.add("    geom_boxplot()+ coord_cartesian(ylim = c(0, "+maxAFForBoxplot+")) + ");
		rScriptAL.add("    geom_boxplot()+ ");
		rScriptAL.add("    ggtitle( paste(testDataset, x)) + ");
		rScriptAL.add("    ylab('Allele Fraction') + ");
		rScriptAL.add("    theme( ");
		rScriptAL.add("      plot.title = element_text(size=12,face='bold'), ");
		rScriptAL.add("      axis.title.x = element_blank(), ");
		rScriptAL.add("      axis.title.y = element_text(face='bold', size=12), ");
		rScriptAL.add("      axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + ");
		rScriptAL.add("    geom_point(data=dfSubTest, aes(x=Gene, y=AlleleFraction), color = 'red', shape=8) ");
		rScriptAL.add("    ggsave(file=paste(x,'.box.png', sep = ''), plot=boxPlotObj, width=10, height=8) ");

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

	private void saveTestAFsForBoxPlot(File toSave) {
		//header row
		StringBuilder sb = new StringBuilder("Gene\tAlleleFraction\n");
		String gene = "";
		int counter = 0;
		for (int i=0; i< accuGenProbes.size(); i++) {
			AccuGenProbe probe = accuGenProbes.get(i);
			if (probe.getGene().equals(gene) == false) {
				gene = probe.getGene();
				counter = 0;
			}
			counter++;
			String counterString = Integer.toString(counter);
			if (counterString.length()==1) counterString = "00"+counterString;
			else if (counterString.length()==2) counterString = "0"+counterString;
			String name = gene+"_"+counterString;
			//for the test sample
			sb.append(name);
			sb.append("\t");
			//sb.append(testSample.getAccuGenProbeCount(probe.getOriginalInput()).getAlleleFraction());
			sb.append(testSample.getAccuGenProbeCount(probe.getOriginalInput()).getScaledAlleleFraction());
			sb.append("\n");
		}
		IO.writeString(sb.toString(), toSave);
	}

	private void saveGermlineAFsForBoxPlot(File toSave) {
		if (germlineSamplesToRemoveFromPoN.size()!=0) IO.pl("\tExcluding from PoN: "+germlineSamplesToRemoveFromPoN);
		//header row
		StringBuilder sb = new StringBuilder("Gene\tAlleleFraction\n");
		String gene = "";
		int counter = 0;
		for (int i=0; i< accuGenProbes.size(); i++) {
			AccuGenProbe probe = accuGenProbes.get(i);
			if (probe.getGene().equals(gene) == false) {
				gene = probe.getGene();
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
				//sb.append(ks.getAccuGenProbeCount(probe.getOriginalInput()).getAlleleFraction());
				sb.append(ks.getAccuGenProbeCount(probe.getOriginalInput()).getScaledAlleleFraction());
				sb.append("\n");
			}
		}
		IO.writeString(sb.toString(), toSave);
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
	
	private void normalizeGermlineAFs() {
		IO.pl("\nNormalizing Germline Allele Fractions....");
		IO.pl("\nOriginal Germline Allele Fractions....");
		//print header
		IO.p("SampleId");
		for (AccuGenProbe p: accuGenProbes) {
			IO.p("\t");
			IO.p(p.getOriginalInput());
		}
		IO.pl("\tMean");

		//for each sample
		ArrayList<Double> means = new ArrayList<Double>();
		for (KohliSample ks: samples.values()) {
			if (ks.isGermlineSample() == false) continue;
			numberGermlineSamples++;
			IO.p(ks.getSampleId());
			double total = 0;
			for (AccuGenProbe p: accuGenProbes) {
				double af = ks.getAccuGenProbeCount(p.getOriginalInput()).getAlleleFraction();
				IO.p("\t");
				IO.p(Num.formatNumber(af, 3));
				total +=af;
			}
			double mean = total/(double)accuGenProbes.size() ;
			IO.pl("\t"+mean);
			means.add(mean);
		}

		//calculate the target
		double targetMean = Num.meanDouble(means);
		IO.pl("\nTargetMean: "+targetMean);

		IO.pl("\nScaling the individual germline samples (sampleID meanAF scalar)....");
		//Scale all of the scores to the means
		int meanIndex = 0;
		for (KohliSample ks: samples.values()) {
			if (ks.isGermlineSample() == false) continue;

			//scale all of the probe afs
			double mean = means.get(meanIndex++);
			double scalar = targetMean/mean;
			IO.pl(ks.getSampleId()+"\t"+mean+"\t"+scalar);
			for (AccuGenProbe p: accuGenProbes) {
				AccuGenProbeCounts c = ks.getAccuGenProbeCount(p.getOriginalInput());
				double af = c.getAlleleFraction();
				c.setScaledAlleleFraction(af * scalar);
			}
		}

		IO.pl("\nScaled Germline Allele Fractions....");
		//print header
		IO.p("SampleId");
		for (AccuGenProbe p: accuGenProbes) {
			IO.p("\t");
			IO.p(p.getOriginalInput());
		}
		IO.pl("\tMean");

		//for each sample
		for (KohliSample ks: samples.values()) {
			if (ks.isGermlineSample() == false) continue;
			IO.p(ks.getSampleId());
			double total = 0;
			for (AccuGenProbe p: accuGenProbes) {
				double af = ks.getAccuGenProbeCount(p.getOriginalInput()).getScaledAlleleFraction();
				IO.p("\t");
				IO.p(Num.formatNumber(af, 3));
				total +=af;
			}
			double mean = total/(double)accuGenProbes.size() ;
			IO.pl("\t"+mean);
			means.add(mean);
		}
	}
	
	private void loadAccugenomicSpikeinData(File vcfCountData) {
		String[] lines = IO.loadFileIntoStringArray(vcfCountData);

		//create a hash of probes to skip, might be empty
		HashSet<String> toSkip = new HashSet<String>();
		for (String s: Misc.COMMA.split(probesToSkip)) toSkip.add(s);
		
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
		TreeSet<String> geneNames = new TreeSet<String>();
		// parse the data lines for each data line
		for (int i=1; i< lines.length; i++) {
			//chr11_114064031_A_T_ZBTB16	0,255,255:3150,152	0,99,255:8440,685	0,136,255:3177,188
			if (lines[i].startsWith("chr") == false) Misc.printErrAndExit("Data line doesn't start with 'chr'? "+lines[i]);
			String[] f = Misc.TAB.split(lines[i]);
			
			//skip it?
			if (toSkip.contains(f[0])) continue;
			
			//use the first field to create a probe
			AccuGenProbe p = new AccuGenProbe(f[0]);
			accuGenProbes.add(p);
			geneNames.add(p.getGene());

			//for each data point
			for (int j=1; j< f.length; j++) {
				//get the sample, add the values
				KohliSample ks = indexSamples[j];
				ks.addAccuGenProbeCounts(p, f[j]);
			}
		}
		//set testedGeneNames
		//"'AR','BRCA2','CHEK2','MYC','NKX3-1','OPHN1','PIK3CA','PIK3CB','TP53','ZBTB16'"
		StringBuilder sb = new StringBuilder("'");
		for (String gn: geneNames) {
			sb.append(gn);
			sb.append("','");
		}
		testedGenes = sb.toString();
		testedGenes = testedGenes.substring(0, testedGenes.length()-2);
		IO.pl("\tTested Genes: "+testedGenes);
		
	}

	private void loadGenesWithProbes() {
		genes = new TreeMap<String, AccuGenGene>(); 
		for (AccuGenProbe p: accuGenProbes) {
			String geneName = p.getGene();
			AccuGenGene gene = genes.get(geneName);
			if (gene==null) {
				gene = new AccuGenGene(geneName);
				genes.put(geneName, gene);
			}
			gene.getProbes().add(p);
		}
		
		//for each sample
		for (KohliSample ks: samples.values()) {
			if (ks.isGermlineSample() == false) continue;
			//for each probe
			for (AccuGenProbe p: accuGenProbes) {
				double scaledAf = ks.getAccuGenProbeCount(p.getOriginalInput()).getScaledAlleleFraction();
				p.getPonAfs().add(scaledAf);
			}
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
		new AccuGenomicGeneSelector(args);
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
					case 'v': vcfCountData = new File(args[++i]).getCanonicalFile(); break;
					case 'r': resultsDirectory = new File(args[++i]).getCanonicalFile(); break;
					case 'g': germlineSamplesNoPoN = args[++i]; break;
					case 's': probesToSkip = args[++i]; break;
					case 'e': fullPathToR = new File(args[++i]); break;
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
		if (patientAnnotationFile==null || patientAnnotationFile.exists()==false) Misc.printErrAndExit("\nFailed to find your -p Patient Annotation File, correct and restart.\n");
		if (vcfCountData==null || vcfCountData.exists()==false) Misc.printErrAndExit("\nFailed to find your -v Vcf Alt Ref Count File, correct and restart.\n");
		if (fullPathToR==null || fullPathToR.exists()==false) Misc.printErrAndExit("\nFailed to find your -e R executible path, correct and restart.\n");
		if (resultsDirectory==null) Misc.printErrAndExit("\nFailed to find your -r results directory path, correct and restart.\n");
		resultsDirectory.mkdirs();

		//any germline samples to exclude from PoN?
		germlineSamplesToRemoveFromPoN = new ArrayList<String>();
		if (germlineSamplesNoPoN != null) {
			for (String s: Misc.COMMA.split(germlineSamplesNoPoN)) germlineSamplesToRemoveFromPoN.add(s);
		}
		
		/*
		-p /Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/sampleMatchingPersonGNomIDs_13June2024.txt 
		-v /Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/BcfToolsForAccGen/justVcfSampleCounts.txt 
		-r /Users/u0028003/HCI/Labs/Kohli_Manish/ProstateCfDNAAnalysis/FinalAggregateAnalysis/KohliCopyRatio/AccuGenomicGeneSelector
		-g 23802X4,23802X11
		*/
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                     AccuGenomic Gene Selector:   August 2024               **\n" +
				"**************************************************************************************\n" +
				"xxxx.\n"+

			"\nOptions:\n"+
			"-p Path to the patient annotation txt file.\n"+
			"-v Path to the vcf count data txt file.\n"+
			"-g Germline sample IDs to exclude from the PoN, comma delimited, no spaces.\n"+
			"-r File path to a directory to write the results files.\n"+

			"\nExample: java -Xmx100G -jar pathTo/USeq/Apps/xxxxx\n\n" +

				"**************************************************************************************\n");

	}
}

