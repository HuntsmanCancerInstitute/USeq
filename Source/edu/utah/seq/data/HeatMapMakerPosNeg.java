package edu.utah.seq.data;

import java.util.*;
import java.util.regex.*;

import util.gen.*;
import trans.misc.*;
import java.io.*;

/**Converts a SmoothingWindow[] into a heat map specific bar file for import into IGB. 
 * Assumes windows are sorted by start position
 * Score index is the index to use when fetching a score for a window from it's float[] of scores.
 * Note, this will place at the 0 and 1 base two values to maintain symmetry in the values.*/
public class HeatMapMakerPosNeg {
	//fields
	private int scoreIndex = -1;	
	private SmoothingWindow[] windows; 
	private float maxValue;
	private float minValue;
	private float minWindowFilter = 0;
	private float maxWindowFilter = 0;
	private ArrayList<Integer> bases = new ArrayList<Integer>();
	private ArrayList<Float> values = new ArrayList<Float>();

	//constructors
	public HeatMapMakerPosNeg(int scoreIndex, float minWindowFilter, float maxWindowFilter){
		this.scoreIndex = scoreIndex;
		this.minWindowFilter = minWindowFilter;
		this.maxWindowFilter = maxWindowFilter;
	}

	/**Tosses any windows within the min max filter range.*/
	public SmoothingWindow[] filterWindows(SmoothingWindow[] windows){
		if (maxWindowFilter == 0 && minWindowFilter == 0) return windows;
		ArrayList<SmoothingWindow> winAL = new ArrayList<SmoothingWindow> (windows.length/2);
		for (int i=0; i< windows.length; i++){
			double score = windows[i].getScores()[scoreIndex];
			if (score <= 0 && score <= minWindowFilter) winAL.add(windows[i]);
			else if (score >= 0 && score >= maxWindowFilter) winAL.add(windows[i]);
		}
		SmoothingWindow[] filtered = new SmoothingWindow[winAL.size()];
		winAL.toArray(filtered);
		return filtered;
	}

	/**Makes a stair step heat map from an array of windows in bar format.
	 * Don't forget to set the barDirectory and score Index!!!!!!!*/
	public PointData makeHeatMapPositionValues(SmoothingWindow[] win){
		//filter?
		windows = filterWindows(win);
		//add blocks
		assembleHeatMapBlocks();
		//balance by adding max or min at zero base provided this isn't covered by a window!
		if (win[0].start > 1) balanceValues();
		//make PD
		PointData pd = new PointData();
		pd.setPositions(Num.arrayListOfIntegerToInts(bases));
		pd.setScores(Num.arrayListOfFloatToArray(values));
		//clear ArrayLists
		bases.clear();
		values.clear();
		return pd;
	}

	public void balanceValues(){
		float value2Balance = 0;
		if ((-1*minValue) > maxValue) value2Balance = -1* minValue;
		else value2Balance = -1 * maxValue;
		bases.add(0, new Integer(1)); values.add(0, new Float(0));
		bases.add(0, new Integer(0)); values.add(0, new Float(value2Balance));
	}


	public void assembleHeatMapBlocks(){
		//reset max and min
		maxValue = 0;
		minValue = 0;
		//find windowed block
		int startIndex = 0;
		int endIndex = 0;
		SmoothingWindow leftWindow = windows[startIndex];
		for (int i=1; i< windows.length; i++){
			SmoothingWindow rightWindow = windows[i];
			//do they overlap?
			//no block found!
			if (overlapOrAbut(leftWindow, rightWindow) == false){			
				endIndex = i;
				//note endIndex is not included when making heatMap blocks
				addHeatMapBlocks(startIndex, endIndex);
				//start new block
				startIndex = i;

			}
			leftWindow = rightWindow;
		}
		//find last
		addHeatMapBlocks(startIndex, windows.length);
	}

	/**Make a heatmap blocks from given window indexes, stop not included.*/
	public void addHeatMapBlocks(int startIndex, int stopIndex){

		//create arrays of float, one per base to hold max positive and max negative scores
		int startBase = windows[startIndex].getStart();
		int stopIndexMin1 = stopIndex -1;
		int stopBase = windows[stopIndexMin1].getStop();
		int numBases = stopBase- startBase;
		float[] maxPosValues = new float[numBases];
		float[] maxNegValues = new float[numBases];

		//load max arrays with max scores
		//for each window
		for (int i=startIndex; i< stopIndex; i++){
			float score = windows[i].getScores()[scoreIndex];
			int baseIndex = windows[i].getStart() - startBase;
			int length = windows[i].getStop() -  windows[i].getStart() + baseIndex;
			//is the score positive 
			if (score > 0.0f ){
				//It's positive, attempt to increase max
				for (int j=baseIndex; j< length; j++){
					if (score > maxPosValues[j]){ 
						maxPosValues[j] = score;
						if (score > maxValue) maxValue = score;
					}
				}
			}
			//is it negative
			else if (score < 0.0f){
				for (int j=baseIndex; j< length; j++){
					if (score < maxNegValues[j]) {
						maxNegValues[j] = score;
						if (score < minValue) minValue = score;
					}
				}
			}
		}

		//average overlapping bases setting in maxPosValues[]
		for (int i=0; i<numBases; i++){
			if (maxPosValues[i] == 0 && maxNegValues[i] != 0) maxPosValues[i] = maxNegValues[i];
			else if (maxPosValues[i] != 0 && maxNegValues[i] != 0 ) maxPosValues[i] = (maxPosValues[i]+maxNegValues[i])/2.0f;
		}		

		//build blocks		
		//open first block
		//set zero mark
		int previousBase = startBase -1;
		if (previousBase < 0) previousBase = 0;
		add(previousBase, 0);
		//set block value
		float blockValue = maxPosValues[0];
		add(startBase, blockValue);

		//advance each base opening and closing blocks
		for (int i=1; i<numBases; i++){
			float testValue = maxPosValues[i];
			if (testValue != blockValue){
				//close old
				add(i-1+startBase,blockValue);
				//open new
				blockValue = testValue;
				add(i+startBase,blockValue);
			}
		}
		//close last block
		add(numBases-1+startBase,blockValue);
		add(numBases+startBase,0);


	}

	/**Adds a line to the global ArrayLists.*/
	public void add(int base, float value){
		bases.add(new Integer(base));
		values.add(new Float(value));
	}

	/**Assumes right window is to right of left or abuts or doesn't overlap.*/
	public boolean overlapOrAbut (SmoothingWindow left, SmoothingWindow right){
		//overlap or abut, note not using sizeOfOligoMinusOne to include abuts
		if (left.getStop() >= right.getStart()) return true;
		return false;
	}

	public float getMaxWindowFilter() {
		return maxWindowFilter;
	}

	public void setMaxWindowFilter(float maxWindowFilter) {
		this.maxWindowFilter = maxWindowFilter;
	}

	public double getMinWindowFilter() {
		return minWindowFilter;
	}

	public void setMinWindowFilter(float minWindowFilter) {
		this.minWindowFilter = minWindowFilter;
	}
}
