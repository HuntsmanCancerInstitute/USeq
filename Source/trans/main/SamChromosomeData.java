package trans.main;
import java.util.*;
import util.gen.*;

public class SamChromosomeData {
	
	//fields
	private float[] realTStats;
	private float[] realMeanDiff;
	private float[] realVariance;
	
	private float[][] permTStats;
	private float[][] permMeanDiff;
	private float[][] permVariance;
	
	//constructor
	public SamChromosomeData(int numberOfWindows, int numberOfPermutations){
		//instantiate arrays
		realTStats = new float[numberOfWindows];
		realMeanDiff = new float[numberOfWindows];
		realVariance = new float[numberOfWindows];
		permTStats = new float[numberOfPermutations][numberOfWindows];
		permMeanDiff = new float[numberOfPermutations][numberOfWindows];
		permVariance = new float[numberOfPermutations][numberOfWindows];
	}
	
	//methods
	public void setRealSamScores(float[] samScores, int windowIndex){
		realTStats[windowIndex] = samScores[0];
		realMeanDiff[windowIndex] = samScores[1];
		realVariance[windowIndex] = samScores[1];
	}
	
	public void setPermSamScores(float[] samScores, int permutationIndex, int windowIndex){
		permTStats[permutationIndex][windowIndex] = samScores[0];
		permMeanDiff[permutationIndex][windowIndex] = samScores[1];
		permVariance[permutationIndex][windowIndex] = samScores[1];
	}
	
	public void calculateFudgeFactor(int windowIndex){
		//fetch perm error terms
		float[] permWindowVar = new float[permVariance.length];
		for (int x=0; x< permVariance.length; x++ ){
			permWindowVar[x] = permVariance[x][windowIndex];
		}
		//sort
		Arrays.sort(permWindowVar);
		//create a vector containing d statistics using every 5th percentile of the perm error terms and the real scores
		float [] dScores = new float[19];
		double percentile =0;
		for (int i=0; i<dScores.length; i++){
			percentile+=5;
			float fudge = new Double(Num.percentile(permWindowVar, percentile)).floatValue();
			float dScore = realMeanDiff[windowIndex]/ (realVariance[windowIndex] + fudge);
			dScores[i] = dScore;
		}
	}
	
	
	//getters setters
	public float[][] getPermMeanDiff() {
		return permMeanDiff;
	}
	public void setPermMeanDiff(float[][] permMeanDiff) {
		this.permMeanDiff = permMeanDiff;
	}
	public float[][] getPermTStats() {
		return permTStats;
	}
	public void setPermTStats(float[][] permTStats) {
		this.permTStats = permTStats;
	}
	public float[][] getPermVariance() {
		return permVariance;
	}
	public void setPermVariance(float[][] permVariance) {
		this.permVariance = permVariance;
	}
	public float[] getRealMeanDiff() {
		return realMeanDiff;
	}
	public void setRealMeanDiff(float[] realMeanDiff) {
		this.realMeanDiff = realMeanDiff;
	}
	public float[] getRealTStats() {
		return realTStats;
	}
	public void setRealTStats(float[] realTStats) {
		this.realTStats = realTStats;
	}
	public float[] getRealVariance() {
		return realVariance;
	}
	public void setRealVariance(float[] realVariance) {
		this.realVariance = realVariance;
	}

}
