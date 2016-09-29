package edu.utah.seq.analysis.ase;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import org.apache.commons.math3.stat.inference.ChiSquareTest;

import htsjdk.samtools.*;
import util.bio.annotation.Bed;
import util.gen.*;
import edu.utah.seq.analysis.BisSeq;
import edu.utah.seq.analysis.multi.PairedCondition;
import edu.utah.seq.data.SmoothingWindow;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.data.sam.SamLayoutForMutation;

/** Application for identifying allelic expression based on a table of snps and bam alignments filtered for alignment bias (see ReferenceMutator and SamComparator).
 * @author Nix
 * */
public class AllelicExpressionComparator {

	//user defined fields
	private File[] nameSnpFiles;
	private File fullPathToR = new File ("/usr/bin/R");
	private File saveDirectory;
	private int minimumSamples = 2;
	private float minFDR = 13;
	private float minLog2Ratio = 1;
	private LinkedHashMap<String, AllelicExpressionSnp> nameSnpA = null;
	private LinkedHashMap<String, AllelicExpressionSnp> nameSnpB = null;
	private boolean deleteTempFiles = true;
	
	private double cacheNumber = 5000;
	private FisherExact fischerExact = new FisherExact((int)cacheNumber);
	private ChiSquareTest chiSquare = new ChiSquareTest();
	

	//constructor
	/**Stand alone.*/
	@SuppressWarnings("unchecked")
	public AllelicExpressionComparator(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//for each pairing
		System.out.println("Contrasting...");
		for (int i=0; i< nameSnpFiles.length; i++){
			File a = nameSnpFiles[i];
			String nameA = Misc.removeExtension(a.getName());
			nameSnpA = (LinkedHashMap<String, AllelicExpressionSnp>) IO.fetchObject(a);
			for (int j=i+1; j< nameSnpFiles.length; j++){
				File b = nameSnpFiles[j];
				String nameB = Misc.removeExtension(b.getName());
				nameSnpB = (LinkedHashMap<String, AllelicExpressionSnp>) IO.fetchObject(b);
				compare(nameSnpA, nameSnpB, nameA, nameB);
			}
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Num.formatNumber(diffTime, 2)+" minutes\n");
	}



	
	private void compare(LinkedHashMap<String, AllelicExpressionSnp> nameSnpA, LinkedHashMap<String, AllelicExpressionSnp> nameSnpB, String nameA, String nameB) {
		System.out.print("\t"+nameA+" vs "+nameB+" ");
		//identify allelic snps with enough data and pass the thresholds
		ArrayList<String> toScore = new ArrayList<String>();
		for (String snpName : nameSnpA.keySet()){
			//good in A?
			AllelicExpressionSnp aesA = nameSnpA.get(snpName);
			if (aesA.getSampleName().size() < minimumSamples) continue;
			//good in B?
			AllelicExpressionSnp aesB = nameSnpB.get(snpName);
			if (aesB == null || aesB.getSampleName().size() < minimumSamples) continue;
			
			//is either ae?
			boolean aFails = (aesA.getFdr() < minFDR || Math.abs(aesA.getLog2Ratio()) < minLog2Ratio); 
			boolean bFails = (aesB.getFdr() < minFDR || Math.abs(aesB.getLog2Ratio()) < minLog2Ratio);
			if (aFails == true && bFails == true) continue;
			
			toScore.add(snpName);
		}
		System.out.println(toScore.size()+" snps");
		
		//calc pvals and log2Ratios (a/b)
		double[][] pvalsRatios = calculatePValues(toScore, nameSnpA, nameSnpB);
		float[] fdrs = Num.benjaminiHochbergCorrectUnsorted(Num.doubleArrayToFloatArray(pvalsRatios[0]));
		
		//print!
		printReport(toScore, nameSnpA, nameSnpB, nameA, nameB, fdrs, pvalsRatios[1]);
		
	}

	private void printReport(ArrayList<String> toScore, LinkedHashMap<String, AllelicExpressionSnp> nameSnpA2, LinkedHashMap<String, AllelicExpressionSnp> nameSnpB2,
			String nameA, String nameB, float[] fdrs, double[] log2Ratios) {
		
		File f = new File (saveDirectory, nameA+"_"+nameB+"_AEC.xls.gz");
		try {
			Gzipper out = new Gzipper(f);
			//header
			out.println("Chr\tPos\tRef\tAlt\tName\tSampleNames A\tRefAltCounts A\tTotalRef A\tTotalAlt A\tSampleGCScores A\tFDR A\tLog2Rto A\tBadHetCall A\tSampleNames B\tRefAltCounts B\tTotalRef B\tTotalAlt B\tSampleGCScores B\tFDR B\tLog2Rto B\tBadHetCall B\tFDR\tLog2Rto");
			
			for (int i = 0; i< toScore.size(); i++){
				String snpName = toScore.get(i);
				AllelicExpressionSnp aesA = nameSnpA.get(snpName);
				AllelicExpressionSnp aesB = nameSnpB.get(snpName);
				out.print(aesA.toString(false));
				out.print("\t");
				out.print(aesB.toString(true));
				out.print("\t");
				out.print(fdrs[i]);
				out.print("\t");
				out.println(log2Ratios[i]);
			}
			
			out.close();
		} catch (Exception e) {
			
			e.printStackTrace();
		} 
	}




	public double[][] calculatePValues(ArrayList<String> toScore, LinkedHashMap<String, AllelicExpressionSnp> nameSnpA, LinkedHashMap<String, AllelicExpressionSnp> nameSnpB){
		double[] pvals = new double[toScore.size()];
		double[] log2Ratios = new double[toScore.size()];
		ArrayList<long[][]> forR = new ArrayList<long[][]>();
		for (int i=0; i< pvals.length; i++){
			String snpName = toScore.get(i);
			AllelicExpressionSnp aesA = nameSnpA.get(snpName);
			AllelicExpressionSnp aesB = nameSnpB.get(snpName);
			
			int[] refAltA = aesA.getTotalRefAltCounts();
			int[] refAltB = aesB.getTotalRefAltCounts();
			
			long[][] counts = new long[][]{
					{refAltA[0], refAltA[1]},
					{refAltB[0],refAltB[1]}
				};
			//calc log2Ratios
			double a = (double) (refAltA[1]+1)/ (double)(2+refAltA[0]+refAltA[1]);
			double b = (double) (refAltB[1]+1)/ (double)(2+refAltB[0]+refAltB[1]);
			log2Ratios[i] = Num.log2(a/b);

			//calculate p-value for differences in unstranded methylation
			boolean sendToR = false;
			double totalObservations = refAltA[0] + refAltA[1] + refAltB[0] + refAltB[1];
			
			//use fischer's?
			if (totalObservations < cacheNumber) {
				double pNoLog = fischerExact.getTwoTailedP(refAltA[0], refAltA[1], refAltB[0], refAltB[1]);
				pvals[i] = Num.minus10log10Float(pNoLog);
			}

			//use chi-square?
			else {	
				//too many for java chi-square?
				double chiSquareStat = Num.chiSquareTestStatistic(counts);
				if (chiSquareStat > 1400) sendToR = true;
				//use apache chi-square
				else{
					double pNoLog = -1;
					try {
						pNoLog = chiSquare.chiSquareTest(counts);
					} catch (Exception e) {}
					if (BisSeq.pValCheck(pNoLog) == false) sendToR = true;
					else pvals[i] = Num.minus10log10Float(pNoLog);
				}
			}
			//assign
			if (sendToR) {
				pvals[i] = Double.MIN_VALUE;
				forR.add(counts);
			}
		}
		
		//any for R
		if (forR.size() !=0){
			int[][] a = new int[forR.size()][2];
			int[][] b = new int[forR.size()][2];
			for (int i=0; i<a.length; i++){
				long[][] counts = forR.get(i);
				a[i] = new int[]{(int)counts[0][0], (int)counts[0][1]};
				b[i] = new int[]{(int)counts[1][0], (int)counts[1][1]};
			}
			double[] pvalues = Num.chiSquareIndependenceTest(a, b, saveDirectory, fullPathToR, true);
			//add back
			int index = 0;
			for (int i=0; i< pvals.length; i++){
				if (pvals[i] == Double.MIN_VALUE) pvals[i] = pvalues[index++];
			}
		}
		
		return new double[][]{pvals, log2Ratios};
	}



	public static void main(String[] args) {

		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new AllelicExpressionComparator(args);
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
					case 'd': nameSnpFiles = IO.extractFiles(new File(args[++i]), ".obj"); break;
					case 'm': minimumSamples = Integer.parseInt(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'f': minFDR = Float.parseFloat(args[++i]); break;
					case 'l': minLog2Ratio = Float.parseFloat(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

			
		if (saveDirectory == null) Misc.printErrAndExit("\nError: please enter a directory in which to save the results.\n");
		saveDirectory.mkdirs();
		if (saveDirectory.canWrite() == false || saveDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: please enter a directory in which to save the results.\n");

		//check for R and required libraries
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printErrAndExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                      Allelic Expression Comparator:  Oct 2014                    **\n" +
				"**************************************************************************************\n" +
				"Looks for changes in allelic expression between two conditions. First run the\n"+
				"AllelicExpressionDetector on each condition. Only snps with minimum # samples and\n"+
				"where in one of the conditions it also passes FDR and log2Rto thresholds.\n\n"+

				"Required Options:\n"+
				"-d Directory containing nameSnp.obj files to compare from the AllelicExpressionDetector.\n"+
				"-s Save directory.\n"+
				
				"\nDefault Options:\n"+
				"-f Minimum -10Log10(FDR) for individual condition allelic expression, default 13\n"+
				"-l Minimum abs(log2Ratio) for individual condition allelic expression, default 1\n"+
				"-m Minimum samples in each condition to compare Snp count data, defaults to 2.\n"+
				"-r Full path to R. Defaults to '/usr/bin/R'\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/AllelicExpressionComparator -s EyeAEC/\n"+
				"       -d EyeAED\n\n" +

				"**************************************************************************************\n");

	}
}
