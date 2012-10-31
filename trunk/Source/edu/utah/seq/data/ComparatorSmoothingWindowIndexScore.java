package edu.utah.seq.data;

import java.util.Comparator;
/**Sorts by score, largest to smallest*/
public class ComparatorSmoothingWindowIndexScore implements Comparator<SmoothingWindowIndex> {
	//fields
	private int index;
	
	public ComparatorSmoothingWindowIndexScore (int scoreIndex){
		index = scoreIndex;
	}
	
		/**Sorts by score, largest to smallest*/
		public int compare(SmoothingWindowIndex first, SmoothingWindowIndex second) {
			float firstScore = first.getSmoothingWindow().getScores()[index];
			float secondScore = second.getSmoothingWindow().getScores()[index];
			if (firstScore > secondScore) return -1;
			if (firstScore < secondScore) return 1;
			return 0;
		}

	}


