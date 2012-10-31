package edu.utah.seq.data;

import java.util.Comparator;

public class ComparatorEnrichedRegionScore implements Comparator {
	
		private int index;
		
		public ComparatorEnrichedRegionScore(int scoreIndex){
			this.index = scoreIndex;
		}
		
		/**Sorts by score, largest to smallest*/
		public int compare(Object arg0, Object arg1) {
			SmoothingWindow first = ((EnrichedRegion)arg0).getBestWindow();
			SmoothingWindow second = ((EnrichedRegion)arg1).getBestWindow();
			if (first.scores[index] > second.scores[index]) return -1;
			if (first.scores[index] < second.scores[index]) return 1;
			return 0;
		}

	}
