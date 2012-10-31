package edu.utah.seq.data;

import java.util.Comparator;

public class ComparatorEnrichedRegionPosition implements Comparator {
		
		/**Sorts by chromosome and position*/
		public int compare(Object arg0, Object arg1) {
			EnrichedRegion first = (EnrichedRegion)arg0;
			EnrichedRegion second = (EnrichedRegion)arg1;
			String chrom1 = first.getChromosome();
			String chrom2 = second.getChromosome();
			if (chrom1.equals(chrom2) == false) return chrom1.compareTo(chrom2);
			int start1 = first.getStart();
			int start2 = second.getStart();
			if (start1 > start2) return 1;
			if (start1 < start2) return -1;
			return 0;
		}

	}
