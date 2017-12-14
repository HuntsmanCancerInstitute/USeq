package edu.utah.seq.parsers.mpileup;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;

/**@author davidnix*/
public class MpileupVariantComparator {

	//user defined fields
	private File pileupFile;
	private int minSnvDP = 25;
	private double minAFForHom = 0.95;
	private double minAFForMatch = 0.90;
	private double minAFForHis = 0.05;
	private int minBaseQuality = 20;

	//internal fields
	private double maxIndel = 0.01;
	private Histogram[] afHist;
	private Histogram[] chrXAfHist;
	private Similarity[] similarities;
	private String[] sampleNames = null;

	//constructor
	public MpileupVariantComparator(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);

		printThresholds();

		parseFile();

		printStats();

		printHist();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}

	public void printThresholds(){
		IO.p("Settings:");
		IO.p(minSnvDP+"\tMinSnvDP");
		IO.p(minAFForHom+"\tMinAFForHom");
		IO.p(minAFForMatch+"\tminAFForMatch");
		IO.p(minAFForHis+"\tMinAFForHis");
		IO.p(minBaseQuality+"\tMinBaseQuality");
		IO.p(maxIndel+"\tMaxIndel");
	}

	public void printStats(){
		//calculate max match and sort
		for (int i=0; i< similarities.length; i++) similarities[i].calculateMaxMatch();
		Arrays.sort(similarities);
		System.out.println("\nStats:");
		for (int i=0; i< similarities.length; i++) System.out.println(similarities[i].toString(sampleNames));
	}

	public void printHist(){
		IO.p("\nHistograms of AFs observed - gender check:");
		for (int i=0; i< afHist.length; i++){
			System.out.println("Sample\t"+sampleNames[i]);
			System.out.println("All NonRef Allele Frequencies:");
			afHist[i].printScaledHistogram();
			System.out.println("chrX NonRef Allele Frequencies:");
			chrXAfHist[i].printScaledHistogram();
			System.out.println();
		}
		
		IO.p("\nHistograms of the AFs for mis matched homozygous variants - contamination check:");
		for (int i=0; i< similarities.length; i++) {
			similarities[i].toStringMismatch(sampleNames);
			System.out.println();
		}
		
	}

	public void parseFile(){
		if (pileupFile != null) System.out.println("\nParsing "+pileupFile.getName());
		else System.out.println("\nParsing standard in");
		String line = null;
		BufferedReader in = null;
		
		try {
			
			if (pileupFile != null) in = IO.fetchBufferedReader(pileupFile);
			else in = new BufferedReader (new InputStreamReader(System.in));
			ParsedSample[] parsedSamples = null;

			int counter = 0;
			while ((line=in.readLine())!= null){
				if (line.startsWith("#")){
					System.out.println("\t"+line);
					continue;
				}

				//parse line
				MpileupLine ml = new MpileupLine(line, minBaseQuality);
				boolean chrX = ml.getChr().equals("X");
				if (ml.getChr() == null || ml.getRef().equals("N")) continue;

				//get samples
				MpileupSample[] samples = ml.getSamples();
				if (samples == null) continue;

				//first set?
				if (afHist == null) {
					makeCounters(samples.length);
					parsedSamples = new ParsedSample[samples.length];
				}

				//check for proper number of samples
				if (samples.length != afHist.length) continue;

				//parse and check samples, this increments the histograms
				boolean fail = false;
				for (int i=0; i<samples.length; i++){
					//watch out for failing samples
					if (samples[i].isPass() == false) {
						fail = true;
						break;
					}
					parsedSamples[i] = new ParsedSample(samples[i], chrX, i);
				}
				if (fail) continue;

				//contrast samples
				for (int i=0; i< similarities.length; i++) similarities[i].contrast(parsedSamples);

				//if (counter++ > 50000) return;
			}


		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing "+ line);
		}
		finally {
			if (in != null)
				try {
					in.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
		}
	}

	//a helper classes
	class ParsedSample {
		boolean passDpIndel;
		MpileupSample ms;
		double[] maxAFIndex = null;

		ParsedSample (MpileupSample ms, boolean chrX, int index){
			this.ms = ms;

			//check snv DP
			int snvDp = ms.getReadCoverageForwardBases() + ms.getReadCoverageReverseBases();
			if (snvDp < minSnvDP) passDpIndel = false;
			else {
				//check indel freq
				double indelRto = 1.0 - ((double)snvDp/(double)ms.getReadCoverageAll());
				if (indelRto > maxIndel) passDpIndel = false;
				else {
					passDpIndel = true;
					maxAFIndex = ms.findMaxSnvAFAndIndex();
					if (maxAFIndex[0] >= minAFForHis){
						afHist[index].count(maxAFIndex[0]);
						if (chrX) chrXAfHist[index].count(maxAFIndex[0]);
					}
				}
			}
		}
	}


	public void makeCounters(int numSamples) throws Exception{
		if (sampleNames!= null && sampleNames.length != numSamples) throw new Exception("\nERROR: the number of samples in the mpileup file ("+numSamples+") don't equal the number of provided sample names.");
		//make names?
		if (sampleNames == null){
			sampleNames = new String[numSamples];
			for (int i=0; i<numSamples; i++) sampleNames[i] = ""+i;
		}
		//make histograms
		afHist = new Histogram[numSamples];
		chrXAfHist = new Histogram[numSamples];
		for (int i=0; i< numSamples; i++){
			afHist[i] = new Histogram(minAFForHis,1,25);
			chrXAfHist[i] = new Histogram(minAFForHis,1,25);
		}

		//make Sims and PCs
		ArrayList<Similarity> al = new ArrayList<Similarity>();
		for (int i=0; i< numSamples; i++){
			for (int j=i+1; j< numSamples; j++){
				al.add(new Similarity(i,j,minAFForHom, minAFForMatch));
			}
		}
		similarities = new Similarity[al.size()];
		al.toArray(similarities);
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MpileupVariantComparator(args);
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
					case 'p': pileupFile = new File(args[++i]); break;
					case 'n': sampleNames = Misc.COMMA.split(args[++i]); break;
					case 'd': minSnvDP = Integer.parseInt(args[++i]); break; 
					case 'a': minAFForHom = Double.parseDouble(args[++i]); break;
					case 'm': minAFForMatch = Double.parseDouble(args[++i]); break;
					case 's': minAFForHis = Double.parseDouble(args[++i]); break;
					case 'b': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Mpileup Variant Comparator: Dec 2017                   **\n" +
				"**************************************************************************************\n" +
				"Identifies homozygous snvs in each sample then looks for them in the other samples to\n"+
				"calculate a similarity score for variants that meet the read depth threshold. MVC\n"+
				"generates a variety of AF histograms for checking gender and sample contamination.\n"+
				"This is slow, best to run it on just the related sample bams (e.g Tumor & Normal\n"+
				"Exome, Tumor RNASeq) plus one unrelated sample bam for comparison. Use the USeq\n"+
				"ClusterMultiSampleVCF app on large batches of vcfs to check for sample misnaming.\n\n"+
				
				"WARNING! Mpileup does not check that the chr order is the same across samples and the\n"+
				"fasta reference, be certain they are or mpileup will produce incorrect counts.\n\n"+

				"Options:\n"+
				"-p Path to a multi-sample mpileup file. e.g.:\n"+
				"      'echo \"#\"$(ls *bam) | gzip > mpileup.gz; samtools mpileup -B\n"+
				"      -q 13 -d 1000000 -f $fastaIndex -l $bedFile *.bam | gzip >> mpileup.gz'\n"+
				"            Alternatively, skip and pipe in the mpileup output to this app.\n"+
				"-d Minimum read depth, defaults to 25\n"+
				"-a Minimum allele frequency to count as a homozygous variant, defaults to 0.95\n"+
				"-m Minimum allele frequency to count a homozygous match, defaults to 0.9\n"+
				"-s Minimim allele frequency for adding a variant to the histograms, defaults to 0.05\n"+
				"-b Minimum base quality, defaults to 20\n"+
				"-n Comma delimited list of sample names as ordered in the mpileup file\n"+


				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/MpileupVariantComparator -p mpileup.gz -v\n"+
				"      -d 20 -b 13 -n 7GF,3GM,7TF,3TM\n\n" +

				"**************************************************************************************\n");

	}

}
