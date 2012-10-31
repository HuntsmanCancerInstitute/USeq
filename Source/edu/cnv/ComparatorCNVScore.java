package edu.cnv;
import java.util.*;

public class ComparatorCNVScore implements Comparator {
		
			/**Sorts by score, largest to smallest*/
			public int compare(Object arg0, Object arg1) {
				CNV first = (CNV) arg0;
				CNV second = (CNV) arg1;
				if (first.medianScore > second.medianScore) return -1;
				if (first.medianScore < second.medianScore) return 1;
				return 0;
			}

		}
