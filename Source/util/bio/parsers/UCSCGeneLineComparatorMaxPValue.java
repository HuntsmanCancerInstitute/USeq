package util.bio.parsers;

import java.util.Comparator;

public class UCSCGeneLineComparatorMaxPValue implements Comparator<UCSCGeneLine> {
		
		public int compare(UCSCGeneLine first, UCSCGeneLine second) {
			if (first.getMaxPValue() > second.getMaxPValue()) return -1;
			if (first.getMaxPValue() < second.getMaxPValue()) return 1;
			return 0;
		}

	}
