package edu.utah.seq.data;

import java.util.Comparator;

public class ComparatorSmoothingWindowPosition implements Comparator {

		/**Sorts by start position, smallest to largest*/
		public int compare(Object arg0, Object arg1) {
			SmoothingWindow first = (SmoothingWindow)arg0;
			SmoothingWindow second = (SmoothingWindow)arg1;
			if (first.start < second.start) return -1;
			if (first.start > second.start) return 1;
			return 0;
		}

	}


