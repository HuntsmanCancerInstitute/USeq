package edu.cnv;
import java.util.*;


public class ComparatorCNVChromPosition implements Comparator {
		
			/**Sorts by score, largest to smallest*/
			public int compare(Object arg0, Object arg1) {
				CNV first = (CNV) arg0;
				CNV second = (CNV) arg1;
				//sort by chromosome
				int compare = first.getGrGraph().getChromosome().compareTo(second.getGrGraph().getChromosome());
				if (compare !=0) return compare * -1;
				//sort by start position
				if (first.getStart() > second.getStart()) return -1;
				if (first.getStart() < second.getStart()) return 1;
				return 0;
			}

		}
