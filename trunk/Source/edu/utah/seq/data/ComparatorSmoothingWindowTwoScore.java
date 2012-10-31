package edu.utah.seq.data;
import java.util.Comparator;

/**Sorts by scores[0], then scores[1]*/
public class ComparatorSmoothingWindowTwoScore implements Comparator {
	
			public int compare(Object arg0, Object arg1) {
				float[] first = ((SmoothingWindow)arg0).getScores();
				float[] second = ((SmoothingWindow)arg1).getScores();
				//by first score
				if (first[0] < second[0]) return -1;
				if (first[0] > second[0]) return 1;
				//by second score
				if (first[1] < second[1]) return -1;
				if (first[1] > second[1]) return 1;
				return 0;
			}

		



}
