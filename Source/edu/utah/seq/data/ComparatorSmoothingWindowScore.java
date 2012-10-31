package edu.utah.seq.data;

import java.util.Comparator;

public class ComparatorSmoothingWindowScore implements Comparator {
	//fields
	private int index;
	
	public ComparatorSmoothingWindowScore (int scoreIndex){
		index = scoreIndex;
	}
	
		/**Sorts by score, smallest to largest*/
		public int compare(Object arg0, Object arg1) {
			SmoothingWindow first = (SmoothingWindow)arg0;
			SmoothingWindow second = (SmoothingWindow)arg1;
			if (first.scores[index] < second.scores[index]) return -1;
			if (first.scores[index] > second.scores[index]) return 1;
			return 0;
		}

	}


