package edu.utah.seq.vcf.cluster;

import java.util.Random;

import util.gen.Misc;
import util.gen.Num;

public class VCFClusterPair {
	//fields
	private VCFClusterSample first;
	private VCFClusterSample second;
	private double observedScore = 0;
	private double maximumScore = 0;
	private double similarity = 0;
	
	//constructor
	public VCFClusterPair(VCFClusterSample first, VCFClusterSample second){
		this.first = first;
		this.second = second;
		score();
	}

	private void score() {
		byte[] fCalls = first.getCalls();
		byte[] sCalls = second.getCalls();
		/**0=homozygous ref, 1=heterozygous, 2=homozygous alt, 3=no call*/
		for (int i=0; i< fCalls.length; i++){
			//either no call? skip
			if (fCalls[i] == 3 || sCalls[i] == 3) continue;
			//add two for max
			maximumScore += 2.0;
			//the same?
			if (fCalls[i] == sCalls[i]) observedScore += 2.0;
			//either a het?
			else if (fCalls[i] == 1 || sCalls[i] == 1) observedScore += 1.0;
			//must be a mismatch
			else observedScore -= 2.0;
		}
		similarity = observedScore/maximumScore;
	}
	
	/*Not used....*/
	private double averageRandom(int numberPermutations){
		//copy first calls
		byte[] rCalls = new byte[first.getCalls().length];
		System.arraycopy(first.getCalls(), 0, rCalls, 0, rCalls.length);
		byte[] sCalls = second.getCalls();
		double[] rScores = new double[numberPermutations];
		Random r = new Random();
		for (int j=0; j< numberPermutations; j++){
			Num.randomize(rCalls, r);
			//score em
			for (int i=0; i< rCalls.length; i++){
				//either no call? skip
				if (rCalls[i] == 3 || sCalls[i] == 3) continue;
				//the same?
				if (rCalls[i] == sCalls[i]) rScores[j] += 2.0;
				//either a het?
				else if (rCalls[i] == 1 || sCalls[i] == 1) rScores[j] += 1.0;
				//must be a mismatch
				else rScores[j] -= 2.0;
			}
		}
		rCalls = null;
		System.out.println(similarity +" vs mean Random "+(observedScore/Num.mean(rScores)));
		return Num.mean(rScores);
	}

	public VCFClusterSample getFirst() {
		return first;
	}

	public VCFClusterSample getSecond() {
		return second;
	}

	public double getObservedScore() {
		return observedScore;
	}

	public double getMaximumScore() {
		return maximumScore;
	}

	public double getSimilarity() {
		return similarity;
	}
}
