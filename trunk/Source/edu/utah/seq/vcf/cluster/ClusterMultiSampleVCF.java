package edu.utah.seq.vcf.cluster;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.vcf.*;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Clusters vcf samples based on genotype.  Use to confirm relatedness, or the lack there of.*/
public class ClusterMultiSampleVCF {

	//fields
	private File vcfFile;
	private File saveDirectory;
	private float minimumRecordQuality = 20;
	private float minimumSampleGenotypeQuality = 20;
	
	VCFClusterSample[] vcfClusterSample = null;
	private int numberRecords = 0;
	private int numberPassingRecords = 0;
	private HashMap<String, VCFClusterPair> clusterPairs = new HashMap<String, VCFClusterPair>();
	double[][] similarityMatrix = null;

	public ClusterMultiSampleVCF(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		
		processArgs(args);
		
		System.out.println("Loading vcf samples...");
		buildClusterSamples();
		System.out.println("\t"+ numberRecords+ "\t# VCF Records");
		System.out.println("\t"+ numberPassingRecords+ "\t# Passing VCF Records (minQUAL "+Num.formatNumber(minimumRecordQuality, 1) +", minGT "+Num.formatNumber(minimumSampleGenotypeQuality, 1)+")");
		System.out.println("\t"+ vcfClusterSample.length+ "\t# Samples");

		statClusterSamples();
		
		System.out.println("\nMaking pairwise comparisons...");
		makePairwiseComparisons();
		
		System.out.println("\nClustering samples...");
		clusterSamples();
		
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}
	
	private void clusterSamples() {
		//make initial Clusters, one per sample 
		VCFCluster[] clusters = new VCFCluster[vcfClusterSample.length];
		for (int i=0; i< clusters.length; i++){
			clusters[i] = new VCFCluster();
			clusters[i].getSamples().add(vcfClusterSample[i]);
		}
		
		
		while (true){
			//find best pair
			if (findBestPair(clusters) == false) break;
		}
		
	}

	private boolean findBestPair(VCFCluster[] clusters) {
		//make all pair comparison
		int bestA = -1;
		int bestB = -1;
		double maxSim = -1;
		for (int i=0; i< clusters.length; i++){
			VCFCluster first = clusters[i];
			if (first == null) continue;
			for (int j=i+1; j< clusters.length; j++){
				VCFCluster second = clusters[j];
				if (second == null) continue;
				double sim = compareClusters (first, second);
				if (sim > maxSim) {
					maxSim = sim;
					bestA = i;
					bestB = j;
				}
			}
		}
		//anything found?
		if (bestA == -1) return false;
		//output
		StringBuilder sb = new StringBuilder("\t");
		sb.append(maxSim);
		sb.append("\t");
		clusters[bestA].fetchIndexNames(sb);
		sb.append("\t");
		clusters[bestB].fetchIndexNames(sb);
		System.out.println (sb);
		//nuke em
		//add b's to a
		clusters[bestA].getSamples().addAll(clusters[bestB].getSamples());
		clusters[bestB] = null;
		return true;
	}
	
	private double compareClusters(VCFCluster first, VCFCluster second){
		//make array of pairwise sim scores
		ArrayList<VCFClusterSample> firstSamples = first.getSamples();
		ArrayList<VCFClusterSample> secondSamples = second.getSamples();
		int numFirst = firstSamples.size();
		int numSecond = secondSamples.size();
		double[] sims = new double[numFirst*numSecond];
		int index = 0;
		for (int i=0; i< numFirst; i++){
			int f = firstSamples.get(i).getSampleIndex();
			for (int j=0; j< numSecond; j++){
				sims[index++] = similarityMatrix[f][secondSamples.get(j).getSampleIndex()];
			}
		}
		//mean  
		return Num.mean(sims);
	}

	private void makePairwiseComparisons() {
		//build hash
		for (int i=0; i< vcfClusterSample.length; i++){
			for (int j=i+1; j< vcfClusterSample.length; j++){
				String name = vcfClusterSample[i].getSampleIndex()+"_"+vcfClusterSample[j].getSampleIndex();
				VCFClusterPair pair = new VCFClusterPair(vcfClusterSample[i], vcfClusterSample[j]);
				clusterPairs.put(name, pair);
			}
		}
		//build out matrix
		similarityMatrix = new double[clusterPairs.size()][clusterPairs.size()];
		for (VCFClusterPair p: clusterPairs.values()){
			int f = p.getFirst().getSampleIndex();
			int s = p.getSecond().getSampleIndex();
			double sim = p.getSimilarity();
			similarityMatrix[f][s] = sim;
			similarityMatrix[s][f] = sim;
		}
	}

	private void statClusterSamples() {
		System.out.println("\nSampleIndex\tName\t#NoCall\tFracNoCall");
		double numPassRec = numberPassingRecords;
		for (int i=0; i< vcfClusterSample.length; i++){
			StringBuilder sb = new StringBuilder();
			sb.append(vcfClusterSample[i].getSampleIndex());
			sb.append("\t");
			sb.append(vcfClusterSample[i].getSampleName());
			sb.append("\t");
			int numNoCall = vcfClusterSample[i].makeCalls();
			sb.append(numNoCall);
			sb.append("\t");
			double noCallFrac = ((double)numNoCall) / numPassRec;
			sb.append(noCallFrac);
			System.out.println(sb);
		}
		
	}

	private void buildClusterSamples() {
		VCFParser parser = new VCFParser(vcfFile);
		
		//make vcf samples
		String[] sampleName = parser.getSampleNames();
		vcfClusterSample = new VCFClusterSample[sampleName.length];
		for (int i=0; i< vcfClusterSample.length; i++) vcfClusterSample[i] = new VCFClusterSample(new String(sampleName[i]), i);
		 
		//walk through each record, check, if passes add genotype to samples
		VCFRecord[] records = parser.getVcfRecords();
		numberRecords = records.length;
		HashSet<String> genotypes = new HashSet<String>();
		for (VCFRecord vr : records){
			//pass record QUAL?
			if (vr.getQuality() < minimumRecordQuality) continue;
			
			//different genotypes?
			VCFSample[] sample = vr.getSample();
			genotypes.clear();
			for (int i=0; i< sample.length; i++){
				//no call? pass qc?
				if (sample[i].isNoCall() || sample[i].getGenotypeQualityGQ() < minimumSampleGenotypeQuality) continue;
				//add genotype to hash
				genotypes.add(sample[i].getGenotypeGT());
			}
			//more than one genotype?
			if (genotypes.size() < 2) continue;
			numberPassingRecords++;
			
			//add info to VCFClusterSamples
			///0=homozygous ref, 1=heterozygous, 2=homozygous alt, 3=no call
			for (int i=0; i< sample.length; i++){
				byte gt = -1;
				//no call? pass qc?
				if (sample[i].isNoCall() || sample[i].getGenotypeQualityGQ() < minimumSampleGenotypeQuality) gt = 3;
				//figure out genotype
				else {
					String genotype = sample[i].getGenotypeGT();
					if (genotype.equals("0/0")) gt = 0;
					else if (genotype.equals("0/1")) gt = 1;
					else if (genotype.equals("1/1")) gt = 2;
					else {
						//System.err.println("WARNING: Odd genotype for sample "+i+" "+sample[i].getUnmodifiedSampleString()+" skipping");
						gt = 3;
					}
				}
				vcfClusterSample[i].getCallsAL().add(new Byte(gt));
			}
			vr = null;
		}
		//clean up, this app is a memory pig!
		parser = null;
		System.gc();
	}

	public static void main(String[] args) {
		if (args.length<1){
			printDocs();
			System.exit(0);	
		}
		new ClusterMultiSampleVCF(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': vcfFile = new File(args[++i]); break;
					//case 's': saveDirectory = new File(args[++i]); break;
					case 'r': minimumRecordQuality = Float.parseFloat(args[++i]); break;
					case 'g': minimumSampleGenotypeQuality = Float.parseFloat(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		//check args
		if (vcfFile == null || vcfFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your vcf file -> "+vcfFile+"\n");
		//if (saveDirectory == null) Misc.printErrAndExit("\nError: please provide a path to a directory to save the results.\n");
		//saveDirectory.mkdirs();
		//if (saveDirectory.isDirectory() == false || saveDirectory.canWrite() == false) Misc.printErrAndExit("\nError: cannot access results directory?\n");
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Cluster Multi Sample VCF: Nov 2014                      **\n" +
				"**************************************************************************************\n" +
				"Clusters samples based on the genotypes of each that differ in one or more samples.\n\n" +

				"Options:\n"+
				//"-s Save directory, full path.\n"+
				"-v Full path to a multi sample vcf file (xxx.vcf/xxx.vcf.gz)). Note, Java often fails\n"+
				"       to parse tabix compressed vcf files.  Best to uncompress.\n\n"+
				"-r Minimum record QUAL score, defaults to 20.\n" +
				"-g Minimum sample genotype GT score, defaults to 20.\n" +
				
				"\n"+
				"Example: java -Xmx2G -jar pathTo/USeq/Apps/ClusterMultiSampleVCF -v ~/UGP/suicide.vcf\n" +
				"\n" +

		"**************************************************************************************\n");		
	}
}
