package edu.utah.seq.mes;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Java implementation of the Max Ent Scan score3 algorithm. 
 * See Yeo and Burge 2004, http://www.ncbi.nlm.nih.gov/pubmed/15285897
 * Scores 3' splice site 20bpIntron:3bpExon
 * @author David.Nix@hci.utah.edu.*/
public class MaxEntScanScore3 {
	//user defined fields
	private File spliceModelDirectory;
	private File user23merFile;

	//internal fields
	private double[][] seqScores;
	private HashMap<String, double[]> baseScores;
	private static final double[] hashValues = {1,4,16,64,256,1024,4096,16384};
	private double[] scoredScores = null;
	private String[] scoredSubSequences = null;

	//constructor
	/**Stand alone*/
	public MaxEntScanScore3(String[] args) {
		processArgs(args);
		loadBaseScores();
		loadSeqScores();
		scoreSequences();
	}

	/**For incorp into other apps.*/
	public MaxEntScanScore3(File spliceModelDirectory) {
		this.spliceModelDirectory = spliceModelDirectory;
		loadBaseScores();
		loadSeqScores();
	}

	/**Checks for 23mers and non GATCgatc bases.*/
	private void scoreSequences() {
		try {
			BufferedReader in = IO.fetchBufferedReader(user23merFile);
			String seq;
			Matcher mat;
			while ((seq = in.readLine()) != null){
				if (seq.startsWith(">")) System.out.println(seq);
				else {
					System.out.print(seq +"\t");
					seq = seq.trim();
					if (seq.length() != 23) {
						System.out.println("skipping, non 23Mer");
						continue;
					}
					String ucSeq = seq.toUpperCase();
					mat = MaxEntScanScore5.NonGATC.matcher(ucSeq);
					if (mat.find()) {
						System.out.println("skipping, non GATCgatc base ");
						continue;
					}
					double score = scoreSequenceNoChecks(ucSeq);
					System.out.println(score);
				}
			}
			in.close();
		} catch (IOException e) {
			System.err.println("\nError loading user test sequences from "+user23merFile+"\n");
			e.printStackTrace();
			System.exit(1);
		}

	}

	/**Scores an upper case 23Mer for 5' splicing potential using MaxEntScan algorithm. 
	 * Only upper case GATC bases, nothing else. Does not check.*/
	public double scoreSequenceNoChecks(String upperCase23Mer) {
		double consensus = scoreConsensus(upperCase23Mer);
		//reduce to 21mer
		String sub = upperCase23Mer.substring(0,18) + upperCase23Mer.substring(20,23);
		double mes = maxEntScore(sub);
		//calc final 
		return consensus + mes;
	}
	
	/**Looks to see if 23mer and no nonGATC bases before scoring, if fail, returns Double.MIN_VALUE*/
	public double scoreSequenceWithChecks(String seq){
		if (seq.length() != 23) return Double.MIN_VALUE;
		Matcher mat = MaxEntScanScore5.NonGATC.matcher(seq);
		if (mat.find()) return Double.MAX_VALUE;
		return scoreSequenceNoChecks(seq);
	}
	
	
	/**Scores each 23mer in the sequence moving 5' to 3', only GATC, upper case sensitive. Does not check.
	 * Assumes you've upper cased and removed non GATC bases.*/
	public double[] scanSequenceNoChecks(String seq){
		int num = seq.length() - 22;
		double[] scores = new double[num];
		for (int i=0; i< num; i++){
			String subSeq = seq.substring(i, i+23);
			scores[i] = scoreSequenceNoChecks(subSeq);
		}
		return scores;
	}
	
	/**Scores each 23mer in the sequence moving 5' to 3' skipping 23mers with non GATCgatc bases.
	 * Case insensitive.*/
	public double[] scanSequence(String seq){
		int num = seq.length() - 22;
		Matcher mat;
		ArrayList<Double> al = new ArrayList<Double>();
		String ucSeq = seq.toUpperCase();
		for (int i=0; i< num; i++){
			String subSeq = ucSeq.substring(i, i+23);
			mat = MaxEntScanScore5.NonGATC.matcher(subSeq);
			if (mat.find() == false) {
				double score = scoreSequenceNoChecks(subSeq);
				al.add(score);
			}
		}
		return Num.arrayListOfDoubleToArray(al);
	}
	
	/**Scores each 23mer in the sequence moving 5' to 3' 23mers with a non GATCgatc base are not scored and are assigned the defaultScore.
	 * Use the getScoredSubSequences, and getScoredScores methods to fetch.
	 * Case insensitive.*/
	public double[] scanSequence(String seq, double defaultScore){
		int num = seq.length() - 22;
		Matcher mat;
		scoredScores = new double[num];
		scoredSubSequences = new String[num];
		String ucSeq = seq.toUpperCase();
		for (int i=0; i< num; i++){
			String subSeq = ucSeq.substring(i, i+23);
			scoredSubSequences[i] = subSeq;
			mat = MaxEntScanScore5.NonGATC.matcher(subSeq);
			if (mat.find() == false) {
				double score = scoreSequenceNoChecks(subSeq);
				scoredScores[i] = score;
			}
			else scoredScores[i] = defaultScore;
		}
		return scoredScores;
	}
	
	/**Scores each 23mer in the sequence moving 5' to 3' skipping seqs with non GATCgatc bases.
	 * Prints results to screen. Case insensitive.*/
	public void scanPrintSequence(String seq){
		int num = seq.length() - 22;
		Matcher mat;
		ArrayList<Double> al = new ArrayList<Double>();
		String ucSeq = seq.toUpperCase();
		System.out.println(seq);
		StringBuilder spaces = new StringBuilder();
		for (int i=0; i< num; i++){
			String subSeq = ucSeq.substring(i, i+23);
			mat = MaxEntScanScore5.NonGATC.matcher(subSeq);
			if (mat.find() == false) {
				double score = scoreSequenceNoChecks(subSeq);
				al.add(score);
				System.out.println(spaces+subSeq+"\t"+score);
			}
			else System.out.println(spaces+subSeq+"\tNonGATCgatc base, skipping");
			spaces.append(" ");
		}
	}
	
	
	

	private double maxEntScore(String seq) {
		double[] sc = new double[9];

		sc[0] = seqScores[0][seqIndex(seq.substring(0,7))];
		sc[1] = seqScores[1][seqIndex(seq.substring(7,14))];
		sc[2] = seqScores[2][seqIndex(seq.substring(14,21))];
		sc[3] = seqScores[3][seqIndex(seq.substring(4,11))];
		sc[4] = seqScores[4][seqIndex(seq.substring(11,18))];
		sc[5] = seqScores[5][seqIndex(seq.substring(4,7))];
		sc[6] = seqScores[6][seqIndex(seq.substring(7,11))];
		sc[7] = seqScores[7][seqIndex(seq.substring(11,14))];
		sc[8] = seqScores[8][seqIndex(seq.substring(14,18))];
		
		return (sc[0] + sc[1] + sc[2] + sc[3] + sc[4]) - (sc[5] + sc[6] + sc[7] + sc[8]);
	}

	private int seqIndex(String seq) {
		int sum = 0;
		int len = seq.length();
		for (int i=0; i< len; i++){
			char b = seq.charAt(i);
			if (b == 'A') continue;
			double baseValue = 1;
			if (b == 'G') baseValue = 2;
			else if (b == 'T') baseValue = 3;
			sum += baseValue * hashValues[len- i- 1];
		}
		return sum;
	}

	private double scoreConsensus(String seq) {
		//fetch base scores bdg, cons1, cons2
		double[] posA = baseScores.get(seq.substring(18, 19));
		double[] posB = baseScores.get(seq.substring(19, 20));
		return (posA[1] + posB[2]) - (posA[0] + posB[0]);
	}

	//methods
	/**Base and bgd, cons1, cons2 scores in log2 space*/
	public void loadBaseScores(){
		baseScores = new HashMap<String, double[]>();
		double[] s= new double[]{0.27,0.9903,0.0027};
		Num.log2(s);
		baseScores.put("A", s);
		s= new double[]{0.27,0.0030,0.0030};
		Num.log2(s);
		baseScores.put("T", s);
		s= new double[]{0.23,0.0034,0.9905};
		Num.log2(s);
		baseScores.put("G", s);
		s= new double[]{0.23,0.0032,0.0037};
		Num.log2(s);
		baseScores.put("C", s);
	}


	/**Precomputed  log2 sequence scores*/
	public void loadSeqScores(){
		String[] tables = {"me2x3acc1","me2x3acc2","me2x3acc3","me2x3acc4","me2x3acc5","me2x3acc6","me2x3acc7","me2x3acc8","me2x3acc9"};
		seqScores = new double[tables.length][];
		for (int i=0; i< tables.length; i++){
			File f = new File(spliceModelDirectory, tables[i]);
			if (f.exists() == false) Misc.printErrAndExit("\nError: cannot find precomputed seq scores for -> "+f);
			seqScores[i] = Num.loadDoubles(f);
			Num.log2(seqScores[i]);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MaxEntScanScore3(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-zA-Z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i];
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 's': spliceModelDirectory = new File(args[++i]); break;
					case 't': user23merFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for required args
		if (spliceModelDirectory == null || spliceModelDirectory.isDirectory()== false) Misc.printExit("\nError: please enter a directory containing the me2x3acc1-9 files.\n");
		if (user23merFile == null) Misc.printExit("\nError: please enter your 23mer sequences to test, one per line, fasta format OK too.\n");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               MaxEntScanScore3: Nov 2013                         **\n" +
				"**************************************************************************************\n" +
				"Implementation of Max Ent Scan's score3 algorithm for human splice site detection. See\n"+
				"Yeo and Burge 2004, http://www.ncbi.nlm.nih.gov/pubmed/15285897 \n\n" +

				"Options:\n"+
				"-s Full path directory name containing the me2x3acc1-9 splice model files. See\n"+
				"     USeq/Documentation/ or http://genes.mit.edu/burgelab/maxent/download/ \n"+
				"-t Full path file name for 23mer test sequences, GATCgatc only, one per line. Fasta OK.\n"+

				"\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/MaxEntScanScore3 -s ~/MES/splicemodels -t\n"+
				"     ~/MES/seqsToTest.fasta \n\n" +

				"**************************************************************************************\n");

	}
	public double[] getScoredScores() {
		return scoredScores;
	}

	public String[] getScoredSubSequences() {
		return scoredSubSequences;
	}

}
