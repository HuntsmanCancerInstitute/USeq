package edu.utah.seq.mes;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Java implementation of the Max Ent Scan score5 algorithm. 
 * See Yeo and Burge 2004, http://www.ncbi.nlm.nih.gov/pubmed/15285897
 * Scores 5' splice site 2bpExon:7bpIntron
 * @author David.Nix@hci.utah.edu.*/
public class MaxEntScanScore5 {
	//user defined fields
	private File spliceModelDirectory;
	private File user9merFile;

	//internal fields
	private HashMap<String, Double> seqScores;
	private HashMap<String, double[]> baseScores;
	public static final Pattern TAB = Pattern.compile("\t");
	public static final Pattern NonGATC = Pattern.compile("[^GATC]");

	//constructor
	/**Stand alone*/
	public MaxEntScanScore5(String[] args) {
		processArgs(args);
		loadBaseScores();
		loadSeqScores();
		scoreSequences();
	}
	
	/**For incorp into other apps.*/
	public MaxEntScanScore5(File spliceModelDirectory) {
		this.spliceModelDirectory = spliceModelDirectory;
		loadBaseScores();
		loadSeqScores();
	}

	/**Checks for 9mers and non GATCgatc bases.*/
	private void scoreSequences() {
		try {
			BufferedReader in = IO.fetchBufferedReader(user9merFile);
			String seq;
			Matcher mat;
			while ((seq = in.readLine()) != null){
				if (seq.startsWith(">")) System.out.println(seq);
				else {
					System.out.print(seq +"\t");
					seq = seq.trim();
					if (seq.length() != 9) {
						System.out.println("skipping, non 9Mer");
						continue;
					}
					String ucSeq = seq.toUpperCase();
					mat = NonGATC.matcher(ucSeq);
					if (mat.find()) {
						System.out.println("skipping, non GATCgatc base ");
						continue;
					}
					double score = scoreSequence(ucSeq);
					System.out.println(score);
				}
			}
			in.close();
		} catch (IOException e) {
			System.err.println("\nError loading user test sequences from "+user9merFile+"\n");
			e.printStackTrace();
			System.exit(1);
		}
		
	}

	/**Scores an upper case 9Mer for 5' splicing potential using MaxEntScan algorithm. 
	 * Only upper case GATC bases, nothing else. Does not check.*/
	public double scoreSequence(String upperCase9Mer) {
		double consensus = scoreConsensus(upperCase9Mer);
		//reduce to 7mer
		String sub = upperCase9Mer.substring(0,3) + upperCase9Mer.substring(5);
		double preComputedScore = seqScores.get(sub);
		//calc final 
		double score = consensus + preComputedScore;
		return score;
	}
	
	/**Scores each 9mer in the sequence moving 5' to 3', only GATC, upper case sensitive. Does not check.
	 * Assumes you've upper cased and removed non GATC bases.*/
	public double[] scanSequenceNoChecks(String seq){
		int num = seq.length() - 8;
		double[] scores = new double[num];
		for (int i=0; i< num; i++){
			String subSeq = seq.substring(i, i+9);
			scores[i] = scoreSequence(subSeq);
		}
		return scores;
	}
	
	/**Scores each 9mer in the sequence moving 5' to 3' skipping 9mers with non GATCgatc bases.
	 * Case insensitive.*/
	public double[] scanSequence(String seq){
		int num = seq.length() - 8;
		Matcher mat;
		ArrayList<Double> al = new ArrayList<Double>();
		String ucSeq = seq.toUpperCase();
		for (int i=0; i< num; i++){
			String subSeq = ucSeq.substring(i, i+9);
			mat = NonGATC.matcher(subSeq);
			if (mat.find() == false) {
				double score = scoreSequence(subSeq);
				al.add(score);
			}
		}
		return Num.arrayListOfDoubleToArray(al);
	}

	private double scoreConsensus(String seq) {
		//fetch base scores bdg, cons1, cons2
		double[] pos3 = baseScores.get(seq.substring(3, 4));
		double[] pos4 = baseScores.get(seq.substring(4, 5));
		double score = (pos3[1] + pos4[2]) - (pos3[0] + pos4[0]);
		return score;
	}

	//methods
	/**Base and bgd, cons1, cons2 scores in log2 space*/
	public void loadBaseScores(){
		baseScores = new HashMap<String, double[]>();
		double[] s= new double[]{0.27,0.004,0.0034};
		Num.log2(s);
		baseScores.put("A", s);
		s= new double[]{0.27,0.0032,0.9884};
		Num.log2(s);
		baseScores.put("T", s);
		s= new double[]{0.23,0.9896,0.0042};
		Num.log2(s);
		baseScores.put("G", s);
		s= new double[]{0.23,0.0032,0.0039};
		Num.log2(s);
		baseScores.put("C", s);
	}

	/**7mer sequence and precomputed score*/
	public void loadSeqScores(){
		try {
			seqScores = new HashMap<String, Double>();
			File seqs = new File (spliceModelDirectory, "splice5sequences");
			File scores = new File (spliceModelDirectory, "me2x5");
			if (seqs.exists() == false || scores.exists() == false) throw new IOException("Failed to find splice5sequences and/ or me2x5 files!");
			BufferedReader inSeqs = IO.fetchBufferedReader(seqs);
			BufferedReader inScores = IO.fetchBufferedReader(scores);
			String seq;
			String scoreString;
			while ((seq = inSeqs.readLine()) != null){
				seq = seq.trim();
				scoreString = inScores.readLine();
				if (scoreString == null) throw new IOException("Score file is shorter than Seq file?!");
				seqScores.put(seq, Num.log2(Double.parseDouble(scoreString.trim())));
			}
			inSeqs.close();
			inScores.close();
		} catch (IOException e) {
			System.err.println("\nError loading splice5sequences or me2x5 from "+spliceModelDirectory+"\n");
			e.printStackTrace();
			System.exit(1);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MaxEntScanScore5(args);
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
					case 't': user9merFile = new File(args[++i]); break;
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
		if (spliceModelDirectory == null) Misc.printExit("\nError: please enter the full  path to the directory containing the 5Seq2Me2x5Score.txt matrix file.\n");
		if (user9merFile == null) Misc.printExit("\nError: please enter your 9mer sequences to test, one per line, fasta format OK too.\n");
	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               MaxEntScanScore5: Nov 2013                         **\n" +
				"**************************************************************************************\n" +
				"Implementation of Max Ent Scan's score5 algorithm for human splice site detection. See\n"+
				"Yeo and Burge 2004, http://www.ncbi.nlm.nih.gov/pubmed/15285897 \n\n" +

				"Options:\n"+
				"-s Full path directory containing the splice5sequences and me2x5 splice model files.\n"+
				"     See USeq/Documentation/ or http://genes.mit.edu/burgelab/maxent/download/ \n"+
				"-t Full path file name for 9mer test sequences, GATCgatc only, one per line. Fasta OK.\n"+

				"\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/MaxEntScanScore5 -s ~/MES/splicemodels -t\n"+
				"     ~/MES/seqsToTest.fasta \n\n" +

				"**************************************************************************************\n");

	}
}
