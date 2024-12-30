package edu.utah.seq.cnv.wgs;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.useq.data.RegionScoreText;
import util.bio.annotation.Bed;

/**
 * @author Nix
 * */
public class CopyAnalysisZScoreMinimizer {

	//user defined fields
	private File tumorPointDir;
	private File[] ponPointDirs;
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private int minimumNumberObservationsInRegion = 100;
	private File bedFile = null;
	private double targetMedian = 1000;

	//internal fields
	private HashMap<String, RegionScoreText[]> chrRegions = null;
	private HashMap<String, PointData[]> tumorPointData = null;
	private HashMap<String, PointData[]>[] ponPointData = null;
	private ArrayList<CopyRegion> copyRegions = new ArrayList<CopyRegion>();
	private int numRegionsFailingReadDepth = 0;
	

	
	//by chromosome
	private String chromosome;

	//constructor	
	/**Stand alone.*/
	public CopyAnalysisZScoreMinimizer(String[] args){
		long startTime = System.currentTimeMillis();
		
		//set fields
		processArgs(args);
		
		//launch
		run();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		IO.pl("\nDone! "+Math.round(diffTime)+" min\n");
	}

	public void run(){	

		registerBedAndPointData();

		//For each chromosome in the bed file
		scanAllChromosomes();
		
		//Normalize pon to the targetMedian
		normalizePoN();
		
		
		
		//find best target and scalar to minimize median z-score
		findBestTarget();
		
		//write out per region stats
		writeOutStats();
		
		//write out bed file of regions with zscores outside pon
		writeOutBed();
		
	}

	private void writeOutBed() {
		IO.pl("\nWriting out bed file of passing regions...");
		Gzipper out = null;
		File results = new File(saveDirectory, "results.bed.gz");
		try {
			out = new Gzipper(results);
			
			for (int j=0; j< copyRegions.size(); j++) {
				CopyRegion cr = copyRegions.get(j);
				if (cr.isTumorZScoreOutsidePonZScores() == false) continue;
				StringBuilder sb = new StringBuilder(cr.getChr());
				sb.append("\t");
				sb.append(cr.getRegion().getStart());
				sb.append("\t");
				sb.append(cr.getRegion().getStop());
				sb.append("\t.\t");
				sb.append((float)cr.calculateZScore());
				sb.append("\t.");
				out.println(sb.toString());
			}
		} catch (Exception e) {
			e.printStackTrace();
			results.delete();
			System.exit(1);
		} finally {
			out.closeNoException();
		}
		
		
	}


	private void writeOutStats() {
		IO.pl("\nWriting out matrix table...");
		Gzipper out = null;
		File results = new File(saveDirectory, "resultsMatrix.txt");
		try {
			out = new Gzipper(results);
			out.println("#Chr:Start-Stop\tZScore\tOutsidePonZScores\tNormTumor\tNormPon1\tNormPon2\tNormPon3...");
			for (int j=0; j< copyRegions.size(); j++) {
				CopyRegion cr = copyRegions.get(j);
				StringBuilder sb = new StringBuilder(cr.getChr());
				sb.append(":");
				sb.append(cr.getRegion().getStart());
				sb.append("-");
				sb.append(cr.getRegion().getStop());
				sb.append("\t");
				sb.append((float)cr.calculateZScore());
				sb.append("\t");
				sb.append(cr.isTumorZScoreOutsidePonZScores());
				sb.append("\t");
				sb.append((float)cr.getNormalizedTumorObs());
				for (float pon: cr.getNormalizedPonObs()) {
					sb.append("\t");
					sb.append(pon);
				}
				out.println(sb.toString());
			}
		} catch (Exception e) {
			e.printStackTrace();
			results.delete();
			System.exit(1);
		} finally {
			out.closeNoException();
		}
	}

	private void findBestTarget() {
		IO.pl("\nScanning tumor for the smallest median z-score...");
		int minTarget = 100;
		int maxTarget = 10000;
		
		//fetch all of the tumor observations, just do once
		float[] tumorObs = new float[copyRegions.size()];
		for (int j=0; j< tumorObs.length; j++) {
			CopyRegion cr = copyRegions.get(j);
			tumorObs[j] = cr.getTumorObs();
		}
		Arrays.sort(tumorObs);
		double medianTumorObs = Num.median(tumorObs);
		IO.pl("\tName\t"+tumorPointDir.getName()+"\n\tMedianRegionObs\t"+medianTumorObs);
		
		
		ZScoreScanResult[] scanResults = new ZScoreScanResult[1+maxTarget-minTarget];
		
		//crude scan, could speed up as threaded
		int index = 0;
		for (int i=minTarget; i<= maxTarget; i++) {
			
			//calculate the scalar
			double regionCountScalar = (double)i / medianTumorObs;
			
			//scale the tumor counts
			double[] zscores = new double[tumorObs.length];
			for (int j=0; j< tumorObs.length; j++) {
				CopyRegion cr = copyRegions.get(j);
				double scaledTumorCount = regionCountScalar * (double)cr.getTumorObs();
				cr.setNormalizedTumorObs((float)scaledTumorCount);
				zscores[j] = cr.calculateZScore();
			}
			
			Arrays.sort(zscores);
			double medianZScores = Num.median(zscores);
			
			scanResults[index++] = new ZScoreScanResult(i,regionCountScalar, medianZScores);
			
		}
		
		//print best 6 results
		Arrays.sort(scanResults);
		int[] topTargetMedians = new int[6];
		IO.pl("\n\tTop results from initial scan. (targetMedian scalar medianZScore)...");
		for (int i=0; i< 6; i++) {
			IO.pl("\t\t"+scanResults[i].targetMedian+"\t"+scanResults[i].scalar+"\t"+scanResults[i].medianZScore);
			topTargetMedians[i] = (int)scanResults[i].targetMedian;
		}
		
		//check that top medians are 1 apart
		Arrays.sort(topTargetMedians);
		int prior = topTargetMedians[0];
		for (int i=1; i< topTargetMedians.length; i++) {
			int diff = topTargetMedians[i] - prior;
			if (diff != 1) Misc.printErrAndExit("\nERROR! The top target medians are not sequential, multiple minimums? Contact developer! Aborting!\n");
			prior = topTargetMedians[i];
		}
		
		//fine scan
		IO.pl("\n\tTop results from the fine scan (targetMedian scalar medianZScore)...");
		double bestTargetMinOne = scanResults[0].targetMedian- 1;
		double bestTargetPlusOne = scanResults[0].targetMedian+ 1;
		ZScoreScanResult[] scanResultsFine = new ZScoreScanResult[201];
		index = 0;
		for (double i=bestTargetMinOne; i<= bestTargetPlusOne; i+=0.01) {

			//calculate the scalar
			double regionCountScalar = i / medianTumorObs;

			//scale the tumor counts
			double[] zscores = new double[tumorObs.length];
			for (int j=0; j< tumorObs.length; j++) {
				CopyRegion cr = copyRegions.get(j);
				double scaledTumorCount = regionCountScalar * (double)cr.getTumorObs();
				cr.setNormalizedTumorObs((float)scaledTumorCount);
				zscores[j] = cr.calculateZScore();
			}
			//calc median
			Arrays.sort(zscores);
			double medianZScores = Num.median(zscores);
			scanResultsFine[index++] = new ZScoreScanResult(i,regionCountScalar, medianZScores);
		}
		//print best 6 results
		Arrays.sort(scanResultsFine);
		
		for (int i=0; i< 6; i++) {
			IO.p("\t\t"+scanResultsFine[i].targetMedian+"\t"+scanResultsFine[i].scalar+"\t"+scanResultsFine[i].medianZScore);
			if (i==0)IO.pl("\t<- Best Result!");
			else IO.pl();
		}
		
		//set final scaled tumor obs
		double bestRegionCountScalar = scanResultsFine[0].scalar;
		for (int j=0; j< tumorObs.length; j++) {
			CopyRegion cr = copyRegions.get(j);
			double scaledTumorCount = bestRegionCountScalar * (double)cr.getTumorObs();
			cr.setNormalizedTumorObs((float)scaledTumorCount);
		}
	}
	
	private class ZScoreScanResult implements Comparable<ZScoreScanResult>{
		private double targetMedian;
		private double scalar;
		private double medianZScore;
		private double absMedianZScore;

		private ZScoreScanResult (double targetMedian, double scalar, double medianZScore) {
			this.targetMedian = targetMedian;
			this.scalar = scalar;
			this.medianZScore = medianZScore;
			absMedianZScore = Math.abs(medianZScore);
		}
		public int compareTo(ZScoreScanResult o) {
			if (o.absMedianZScore< this.absMedianZScore) return 1;
			if (o.absMedianZScore> this.absMedianZScore) return -1;
			return 0;
		}
	}

	private void scanAllChromosomes() {
		IO.p("\nFetching observations over each chromosome...\n\t");
		for (String chr: chrRegions.keySet()) {
			IO.p(chr+" ");
			chromosome = chr;
			scanChromosome();
		}
		IO.pl();
		int totalRegions = copyRegions.size()+numRegionsFailingReadDepth;
		double fractionPassing = (double)copyRegions.size() / (double)totalRegions;
		IO.pl("\t"+copyRegions.size()+ "("+Num.formatNumber(fractionPassing, 2)+")\tRegions passing minimum observations, "+ minimumNumberObservationsInRegion);
		
		
	}

	private void normalizePoN() {
		IO.pl("\nNormalizing the PoN (name medianRegionObs targetMedian scalar)...");
		
		//for each PoN replica
		for (int i=0; i<ponPointDirs.length; i++) {
			IO.p("\t"+ ponPointDirs[i].getName());
			
			//fetch all of the counts
			float[] regionCounts = new float[copyRegions.size()];
			for (int j=0; j< regionCounts.length; j++) {
				CopyRegion cr = copyRegions.get(j);
				regionCounts[j] = cr.getPonObs()[i];
			}
			
			//calculate the scalar
			Arrays.sort(regionCounts);
			double median = Num.median(regionCounts);
			double regionCountScalar = targetMedian / median;
			
			//scale the counts
			for (int j=0; j< regionCounts.length; j++) {
				CopyRegion cr = copyRegions.get(j);
				regionCounts[j] = cr.getPonObs()[i];
				double scaledCount = regionCountScalar * (double)regionCounts[j];
				float[] normCounts = cr.getNormalizedPonObs();
				normCounts[i] = (float) scaledCount;
			}
			IO.pl("\t"+median+"\t1000\t"+regionCountScalar);
		}
		
	}

	/**Window scans a chromosome collecting read count data and calculating binomial p-values.*/
	public void scanChromosome(){
		//load and merge strands
		PointData[] tumorStranded = tumorPointData.get(chromosome);
		PointData tumor = PointData.mergePairedPointData(tumorStranded[0], tumorStranded[1], true);
		
		PointData[] pon = new PointData[ponPointData.length];
		for (int i=0; i< ponPointData.length; i++) {
			PointData[] ponStranded = ponPointData[i].get(chromosome);
			pon[i] = PointData.mergePairedPointData(ponStranded[0], ponStranded[1], true);
		}
		
		scanRegions(tumor, pon);
	}
	

	private void scanRegions(PointData tumor, PointData[] pon){
		
		//for each bed region in the chromosome 
		for (RegionScoreText region: chrRegions.get(chromosome)){
			//tumor
			float tumorObs = tumor.sumScoreBP(region.getStart(), region.getStop());
			if (tumorObs < minimumNumberObservationsInRegion) {
				numRegionsFailingReadDepth++;
				continue;
			}
			
			//pon
			float[] ponObs = new float[pon.length];
			for (int i=0; i< pon.length; i++) {
				ponObs[i] = pon[i].sumScoreBP(region.getStart(), region.getStop());
				if (ponObs[i] < minimumNumberObservationsInRegion) {
					numRegionsFailingReadDepth++;
					continue;
				}
			}
			//save it
			copyRegions.add(new CopyRegion(chromosome, region, tumorObs, ponObs));
		}
	}


	
	/**Collects and calculates a bunch of stats re the PointData.*/
	private void registerBedAndPointData(){
		IO.pl("Loading bed, tumor, and pon files...");
		//parse bed file 
		chrRegions = Bed.parseBedFile(bedFile, true, false);
		
		//fetch tumor PointData, does not load the actual scores
		tumorPointData = PointData.fetchStrandedPointData(tumorPointDir);
		
		//PoN
		ponPointData = PointData.fetchStrandedPointData(ponPointDirs);

		//check chrs
		for (String chr: chrRegions.keySet()) {
			if (tumorPointData.containsKey(chr) == false || tumorPointData.get(chr).length!=2) Misc.printErrAndExit("\tFailed to find chr '"+chr+"' in "+tumorPointDir+" or two stranded datasets.");
			for (int i=0; i< ponPointData.length; i++) {
				if (ponPointData[i].containsKey(chr) == false || ponPointData[i].get(chr).length!=2) Misc.printErrAndExit("\tFailed to find chr '"+chr+"' in "+ponPointDirs[i]+" or two stranded datasets.");
			}
		}
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CopyAnalysisZScoreMinimizer(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': tumorPointDir = new File(args[++i]); break;
					case 'p': ponPointDirs = IO.extractOnlyDirectories(new File(args[++i])); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'm': minimumNumberObservationsInRegion = Integer.parseInt(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//look bed file
		if (bedFile == null ) Misc.printExit("\nError: cannot find your bed file of regions to scan!\n");

		//look for point directories
		if (tumorPointDir == null || tumorPointDir.isDirectory() == false) Misc.printExit("\nError: cannot find your tumor PointData directory!\n");
		
		//pon data
		if (ponPointDirs == null || ponPointDirs.length < 3) Misc.printExit("\nError: cannot find three or more normal PointData directories!\n");

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory to save results.\n");
		if (saveDirectory.exists() == false && saveDirectory.mkdir() == false)  Misc.printExit("\nError: could not find nor make the save directory.\n");
		
		//check for R and required libraries
		/*
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}
		else {
			String errors = IO.runRCommandLookForError("library(DESeq2); library(gplots)", fullPathToR, saveDirectory);
			if (errors == null || errors.length() !=0){
				Misc.printExit("\nError: Cannot find the required R libraries?  Did you install DESeq2? gplots? Once installed, " +
						"launch an R terminal and type 'library(DESeq2); library(gplots)' to see if present. R error message:\n\t\t"+errors+"\n\n");
			}
		}*/
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                       Copy Analysis ZScore Minimizer: May 2024                   **\n" +
				"**************************************************************************************\n" +
				"Beta\n"+
				"Takes a bed file of regions, counts the number of observations under each for the single tumor and panel of normals."+
				"Scales each normal in the PoN to 1000, then calculates a scalar for the tumor that minimizes the median z-score by scanning."+
				"Use this scalar in downstream apps to normalize the tumor when calling copy alterations."+

				"\nOptions:\n"+
				"-s Save directory, full path.\n"+
				"-t Treatment replica PointData directories, full path, comma delimited, no spaces,\n" +
				"       one per biological replica. Use the PointDataManipulator app to merge same\n" +
				"       replica and technical replica datasets. Each directory should contain stranded\n" +
				"       chromosome specific xxx_-/+_.bar.zip files. Alternatively, provide one\n" +
				"       directory that contains multiple biological replical PointData directories.\n" +
				"-c Control replica PointData directories, ditto. \n" +
				"-r Full path to 64bit R loaded with DESeq library, defaults to '/usr/bin/R' file, see\n" +
				"       http://www-huber.embl.de/users/anders/DESeq/ . Type 'library(DESeq2)' in\n" +
				"       an R terminal to see if it is installed.\n"+
				"-p Peak shift, average distance between + and - strand peaks for chIP-Seq data, see\n" +
				"       PeakShiftFinder or set it to 100bp. For RNA-Seq set it to 0. It will be used\n" +
				"       to shift the PointData by 1/2 the peak shift.\n"+
				"-w Window size, defaults to the peak shift. For chIP-Seq data, a good alternative \n" +
				"       is the peak shift plus the standard deviation, see the PeakShiftFinder app.\n" +
				"       For RNA-Seq data, set this to 100-250.\n"+
				
				"\nAdvanced Options:\n"+
				"-m Minimum number of reads in a window, defaults to 15\n"+
				"-d Don't delete temp files\n" +

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/MultipleReplicaScanSeqs -t\n" +
				"      /Data/PolIIRep1/,/Data/PolIIRep2/ -c /Data/Input1/,Data/Input2/ -s\n" +
				"      /Data/PolIIResults/ -p 150 -w 250 -b \n\n" +

		"**************************************************************************************\n");

	}
}
