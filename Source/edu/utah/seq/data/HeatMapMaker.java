package edu.utah.seq.data;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.*;
import trans.misc.*;
import java.io.*;
import java.text.*;

/**Converts a SmoothingWindow[] into a heat map specific bar file for import into IGB. 
 * Assumes windows are sorted by chromosome and start position
 * Score index is the index to use when fetching a score for a window from it's double[] of scores.*/
public class HeatMapMaker {

	//fields
	private int scoreIndex = 0;	
	private int numWindows; 
	private float baseScore = 0;
	private SmoothingWindow[] windows;

	//constructor
	public HeatMapMaker(int scoreIndex, float baseScore){
		this.scoreIndex = scoreIndex;
		this.baseScore = baseScore;
	}

	/**Returns a PointData object with loaded position and values representing a top score heatmap/ stairstep.
	 * Be sure to expand the SmoothingWindows so their start and stops represent the true read size and NOT the
	 * middle position. Interbase coordinates.*/
	public PointData makeHeatMapPositionValues(SmoothingWindow[] windows, boolean oneDecimal){
		PointData pd = new PointData();
		
		//clone and reduce sig figs to one decimal?
		if (oneDecimal){
			DecimalFormat df = new DecimalFormat("#.#");
			SmoothingWindow[] clnWins = new SmoothingWindow[windows.length];
			//for each window, round score and place in new sw
			for (int i=0; i< windows.length; i++){
				float[] scores = windows[i].getScores();
				float[] rndScrs = new float[scores.length];
				rndScrs[scoreIndex] = Float.parseFloat(df.format(scores[scoreIndex]));
				clnWins[i] = new SmoothingWindow(windows[i].getStart(), windows[i].getStop(), rndScrs);
			}
			this.windows = clnWins;
		}
		//nope just reference
		else this.windows = windows;
		
		//any windows?
		numWindows = windows.length;	
		if (numWindows < 1) return pd;
		
		//make AL to hold data
		ArrayList<Integer> positions = new ArrayList(numWindows/5);
		ArrayList<Float> values = new ArrayList(numWindows/5);
		
		//set base score
		positions.add(new Integer(0));
		values.add(new Float(baseScore));
		
		//assign 1st score and 1st position
		int base = windows[0].getStart();
		if (base == 0) base = 1;
		positions.add(new Integer((base -1)));
		values.add(new Float(baseScore));
		double score = windows[0].getScores()[scoreIndex];
		positions.add(new Integer((base)));
		values.add(new Float(score));

		//for each window  
		for (int y = 0; y< numWindows; y++){
			//for each base in the window
			int endBase = windows[y].getStop();
			for (int z=base; z<endBase; z++){
				double highestScore = highestScoringWindow(z, y);
				if (highestScore != score){
					//close old block start new block
					positions.add(new Integer((z-1)));
					values.add(new Float(score));
					positions.add(new Integer((z)));
					values.add(new Float(highestScore));
					score = highestScore;
				}
			}
			//is this the last window in the chromosome?
			if (numWindows == (y+1)){
				//close old block
				positions.add(new Integer((endBase-1)));
				values.add(new Float(score));
				positions.add(new Integer((endBase)));
				values.add(new Float(baseScore));
			}
			//Does the next base overlap next window?
			else if (baseOverlapsWindow(endBase, windows[y+1]) == false){
				//close old block
				positions.add(new Integer((endBase-1)));
				values.add(new Float(score));
				positions.add(new Integer((endBase)));
				values.add(new Float(baseScore));
				//reset score and base
				score = windows[y+1].getScores()[scoreIndex];
				base = windows[y+1].getStart();
				//start new block
				positions.add(new Integer((base-1)));
				values.add(new Float(baseScore));
				positions.add(new Integer((base)));
				values.add(new Float(score));
			}
			else{
				base = endBase;
			}
		}
		//set position and values in PointData object
		pd.setPositions(Num.arrayListOfIntegerToInts(positions));
		pd.setScores(Num.arrayListOfFloatToArray(values));
		return pd;
	}	



	/**Looks forward at overlapping windows, returns highest score for the given bp position.*/
	public double highestScoringWindow (int base, int currentWindowIndex){
		double score = windows[currentWindowIndex].getScores()[scoreIndex];
		for (int i=currentWindowIndex; i<numWindows; i++){
			if (baseOverlapsWindow(base,windows[i])) {
				//check score
				double testScore = windows[i].getScores()[scoreIndex];
				if (testScore> score) score = testScore;				
			}
			else return score;
		}
		return score;
	}

	/**Is a given base contained within the window?
	 * Interbase coordinates so last window position is not included.*/
	public boolean baseOverlapsWindow(int bp, SmoothingWindow win){
		if (bp>= win.getStart() && bp < win.getStop()) return true;
		return false;
	}

	public int getScoreIndex() {
		return scoreIndex;
	}
	public void setScoreIndex(int scoreIndex) {
		this.scoreIndex = scoreIndex;
	}
	public float getBaseScore() {
		return baseScore;
	}
	public void setBaseScore(float baseScore) {
		this.baseScore = baseScore;
	}
}
