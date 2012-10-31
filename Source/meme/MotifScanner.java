package meme;
import java.util.*;

import util.bio.seq.*;

/**MotifScanner objects are used to scan sequences with a particular weight matrix and generate 
 an array of MotifHit objects.*/
public class MotifScanner {
	private double[][] matrix;  //log likelihood position specific probability matrix
	private int lenMotif;
	private MemeResults results;
	private ArrayList motifHits = new ArrayList(); //arraylist containing an array of MotifHits for each seq found to have hits
	double[] baseScores;		//centered motif hit score for every base in sequence scanned
	
	public MotifScanner(double[][] logLikelihoodMatrix, MemeResults res){
		matrix = logLikelihoodMatrix;
		lenMotif = matrix.length;
		results = res;
	}
	public MotifScanner(double[][] logLikelihoodMatrix){
		matrix = logLikelihoodMatrix;
		lenMotif = matrix.length;
	}    
	
	public ArrayList getMotifHits(){
		return motifHits;
	}
	/**Calculates the maximum number of hits that occur within a given window size.
	 * Assumes that the array hasn't been sorted and is ordered by position.*/
	public static int calculateMaxCluster(MotifHit[] hits, int window){
		int len = hits.length;
		if (len < 2) return len;
		int counter = 1;
		int oldCounter=0;
		int index = 0;
		MotifHit start = hits[index];
		for (int i=1; i<len; i++){
			//if within window, increment counter
			if ( (hits[i].getStop() - start.getStart()) <= window){
				counter++;
			}
			//not within window
			else{
				if (counter> oldCounter) oldCounter = counter;
				start = hits[++index];
				i = index;
				counter = 1;
			}
		}
		if (counter> oldCounter) oldCounter = counter;
		return oldCounter;
	}
	
	/**Returns the lowest score from scoring all the seqs.  Useful in getting a
	 cutOff score when the seqs are the seqs used to generate the matrix.
	 All the seqs must be the same size as the matrix!*/
	public double findLowestScoringSeq(String[] seqs){
		int len= seqs.length;
		double lowestScore = 10000000;
		for (int i=0; i<len; i++){
			double test = scoreSubSeq(seqs[i]);
			//System.out.println(test + " "+ seqs[i]);
			if (test<lowestScore) lowestScore = test;
		}
		return lowestScore;
	}
	
	/**Scans a set of sequences for the motif, saves results in a MemeResults object.*/
	public int scanEm(String[] seqs, String[] names, double cutOff){
		int numSeqsWithMotifs = 0;
		for (int i=0; i<seqs.length; i++){
			MotifHit[] hits = scoreSequence(cutOff, seqs[i]);
			Arrays.sort(hits);
			int len = hits.length;
			if (len != 0){
				numSeqsWithMotifs++;
				results.appendMotifSearchResults(
						"Name: "+names[i]+     //printing motif hits
						"  Seq: "+seqs[i]+"\n");
				for (int j=0; j<len; j++){
					results.appendMotifSearchResults(hits[j]+"\n");
				}
				
				
				motifHits.add(hits);
			}
		}
		return numSeqsWithMotifs;
	}
	
	/**Scans a set of sequences for the motif, prints as it goes.*/
	public void scanPrintSequences(String[] seqs, String[] names, double cutOff){
		for (int i=0; i<seqs.length; i++){
			MotifHit[] hits = scoreSequence(cutOff, seqs[i].toUpperCase());
			Arrays.sort(hits);
			int len = hits.length;
			if (len != 0){
				System.out.print(
						"Name: "+names[i]+     //printing motif hits
						"  Seq: "+seqs[i]+" ");
				for (int j=0; j<len; j++){
					System.out.print(hits[j]+" ");
				}
				System.out.println();
			}
		}
		
	}
	/**Scans a set of sequences for the motif, prints as it goes.*/
	public void scanPrintAllSequences(String[] seqs, String[] names){
		for (int i=0; i<seqs.length; i++){
			MotifHit[] hits = scoreSequence(0, seqs[i].toUpperCase());
			Arrays.sort(hits);
			int len = hits.length;
				System.out.println(
						"Name: "+names[i]+     //printing motif hits
						"  Seq: "+seqs[i]+" ");
				for (int j=0; j<len; j++){
					System.out.println("\t"+hits[j]+" ");
				}
				System.out.println();
		}
		
	}
	
	/**Scans a set of sequences for the motif, prints as it goes.*/
	public void scanPrintNumberHits(double scoreCutOff, String[] seqs, String[] names){
		int numWithHits = 0;
		long totalNumberHits = 0;
		long totalLength = 0;
		for (int i=0; i<seqs.length; i++){
			MotifHit[] hits = scoreSequence(scoreCutOff, seqs[i].toUpperCase());
			System.out.println("Seq: "+names[i]);
			System.out.println("\tLength: "+seqs[i].length());
			System.out.println("\t# Hits: "+hits.length);
			totalLength+= seqs[i].length();
			if (hits.length !=0) {
				//System.out.println(names[i]+"\t"+len+"\t");
				numWithHits++;
				totalNumberHits += hits.length;
			}
		}
		System.out.println("\nNumber sequences with hits: "+numWithHits);
		System.out.println("Total number hits: "+totalNumberHits);
		System.out.println("Total length scanned: "+totalLength);
		
	}
	
	
	/**Returns an array of MotifHits that are > or = to the loglikelihood cutOff.
	 * Set cutOff to a very negative number (ie -20*lengthOfMotif) if you want all the MotifHits.
	 * Will not score any non GATC containing sequences. Searches both strands but returns coordinates 
	 * for either in forward direction. 
	 * Also assigning the higher score to the center base for every window (baseWindow double[]) where 
	 * index 0 in the double[] is base 1 in the seq. 0s will be found at 
	 * the beginning and stop. Scores < 0 are assigned 0.
	 * */
	public MotifHit[] scoreSequence(double cutOff, String seq){
		double score;
		double revCompScore;
		String test;
		String revCompTest;
		int ori;
		ArrayList hits = new ArrayList();
		int seqLen = seq.length();
		if (seqLen < lenMotif){
			System.err.println("Warning: cannot score sequence, sequence is smaller than motif!");
			return new MotifHit[0];
		}
		baseScores = new double[seqLen];
		int halfLenMotif = Math.round((float)lenMotif/2f);
		//fragment text
		int len = seqLen-lenMotif+1;
		for (int i=0; i<len; i++){
			test = seq.substring(i, i+lenMotif);
			revCompTest = Seq.reverseComplementDNA(test);
			//score the subfrags
			score = scoreSubSeq(test);
			revCompScore = scoreSubSeq(revCompTest);
			//find max score
			ori = 0;
			if (revCompScore > score) {
				score = revCompScore;
				ori = 1;
				test = revCompTest;
			}
			//add to array of MotifHits, interbase coordinates
			if (score>= cutOff ) {
				hits.add(new MotifHit(i, i+lenMotif+1, score, test, ori));
			}
			//record score in double[]
			if (score<0) score = 0;
			baseScores[i+halfLenMotif] = score;
		}
		MotifHit[] mh= new MotifHit[hits.size()];
		hits.toArray(mh);
		return mh;
	}
	
	/**Scans a sequence for the motif, scores every window, forward and reverse complement, assigning the 
	 * higher score to the center base. Index 0 in the double[] is base 1 in the seq. 0s will be found at 
	 * the beginning and stop.  Scores < 0 are assigned 0*/
	public void scanSequence(String seq){
		double score;
		double revCompScore;
		String test;
		String revCompTest;
		int ori;
		int seqLen = seq.length();
		//fragment text
		baseScores = new double[seqLen];
		int len = seqLen-lenMotif+1;
		int halfLenMotif = Math.round((float)lenMotif/2f);
		for (int i=0; i<len; i++){
			test = seq.substring(i, i+lenMotif);
			revCompTest = Seq.reverseComplementDNA(test);
			//score the subfrags
			score = scoreSubSeq(test);
			revCompScore = scoreSubSeq(revCompTest);
			//find max score
			if (revCompScore > score) score = revCompScore;
			if (score<0) score = 0;
			baseScores[i+halfLenMotif] = score;
		}
	}
	
	
	
	/**ignores non GATC chars*/
	public double scoreSubSeq(String seq){
		double score = 0;
		for (int j=0; j<lenMotif; j++){
			char test = seq.charAt(j);
			switch (test){
			case 'A': score+= matrix[j][0]; break;
			case 'C': score+= matrix[j][1]; break;
			case 'G': score+= matrix[j][2]; break;
			case 'T': score+= matrix[j][3]; break;
			//default :return 0;
			}
		}
		return score;
	}
	
	
	public double[] getBaseScores() {
		return baseScores;
	}
}
