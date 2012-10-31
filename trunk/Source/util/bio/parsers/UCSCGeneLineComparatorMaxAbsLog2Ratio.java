package util.bio.parsers;

import java.util.Comparator;

public class UCSCGeneLineComparatorMaxAbsLog2Ratio implements Comparator<UCSCGeneLine> {
		
		public int compare(UCSCGeneLine first, UCSCGeneLine second) {
			if (first.getMaxAbsLog2Ratio() > second.getMaxAbsLog2Ratio()) return -1;
			if (first.getMaxAbsLog2Ratio() < second.getMaxAbsLog2Ratio()) return 1;
			return 0;
		}

	}
