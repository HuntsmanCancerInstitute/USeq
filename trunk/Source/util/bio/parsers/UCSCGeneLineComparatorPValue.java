package util.bio.parsers;

import java.util.Comparator;

public class UCSCGeneLineComparatorPValue implements Comparator<UCSCGeneLine> {
		
		public int compare(UCSCGeneLine first, UCSCGeneLine second) {
			if (first.getpValue() > second.getpValue()) return -1;
			if (first.getpValue() < second.getpValue()) return 1;
			return 0;
		}

	}
