package edu.utah.seq.data;

import java.util.Comparator;

public class ComparatorEnrichedRegionAbsoluteScore implements Comparator<EnrichedRegion> {
	
		private int index;
		
		public ComparatorEnrichedRegionAbsoluteScore(int scoreIndex){
			this.index = scoreIndex;
		}
		
		/**Sorts by score, largest to smallest*/
		public int compare(EnrichedRegion arg0, EnrichedRegion arg1) {
			float first = Math.abs(arg0.getBestWindow().getScores()[index]);
			float second = Math.abs(arg1.getBestWindow().getScores()[index]);
			if (first > second) return -1;
			if (first < second) return 1;
			return 0;
		}

	}
