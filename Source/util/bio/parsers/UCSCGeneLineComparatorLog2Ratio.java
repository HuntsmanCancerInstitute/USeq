package util.bio.parsers;

import java.util.Comparator;

public class UCSCGeneLineComparatorLog2Ratio implements Comparator<UCSCGeneLine> {
		
		public int compare(UCSCGeneLine first, UCSCGeneLine second) {
			if (first.getLog2Ratio() > second.getLog2Ratio()) return -1;
			if (first.getLog2Ratio() < second.getLog2Ratio()) return 1;
			return 0;
		}

	}
