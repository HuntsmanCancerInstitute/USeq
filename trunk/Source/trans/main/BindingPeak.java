package trans.main;

import java.io.Serializable;

import util.gen.*;
/**
 * Holds information about a binding peak within an interval. 
 */
public class BindingPeak implements Serializable{
	//fields
	private int firstPeakIndex; //there may be additional oligos to the right with the same peak score
	private int leftFlankingIndex;
	private int rightFlankingIndex;
	private int peakBP;
	private double score;
	
	//constructor
	public BindingPeak(int firstPeakIndex, int leftFlankingIndex, int rightFlankingIndex, 
			int peakBP, double score){
		this.firstPeakIndex = firstPeakIndex;
		this.leftFlankingIndex = leftFlankingIndex;
		this.rightFlankingIndex = rightFlankingIndex;
		this.peakBP = peakBP;
		this.score = score;
	}
	
	/**Score, Left, Peak, Right*/
	public String toString(Oligo[] oligos, int oligoLength){
		StringBuffer s = new StringBuffer();
		s.append(Num.formatNumberOneFraction(score)); s.append("\t");
		s.append(oligos[leftFlankingIndex].getStart());s.append("\t");
		s.append(peakBP);s.append("\t");
		s.append(oligos[rightFlankingIndex].getStart()+oligoLength);
		return s.toString();
	}
	
	
	public int getFirstPeakIndex() {
		return firstPeakIndex;
	}
	public int getLeftFlankingIndex() {
		return leftFlankingIndex;
	}
	public int getPeakBP() {
		return peakBP;
	}
	public int getRightFlankingIndex() {
		return rightFlankingIndex;
	}
	public double getScore() {
		return score;
	}
}
