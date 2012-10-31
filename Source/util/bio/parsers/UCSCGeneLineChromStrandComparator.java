package util.bio.parsers;

import java.util.Comparator;

public class UCSCGeneLineChromStrandComparator implements Comparator {

			/**Sorts by chromosomeStrand -> start position -> length, shortest first*/
			public int compare(Object arg0, Object arg1) {
				UCSCGeneLine first = (UCSCGeneLine)arg0;
				UCSCGeneLine second = (UCSCGeneLine)arg1;
				String firstChromStrand = first.getChrom()+first.getStrand();
				String secondChromStrand = second.getChrom()+second.getStrand();
				int i = firstChromStrand.compareTo(secondChromStrand);
				if (i !=0 ) return i;
				if (first.getTxStart()< second.getTxStart()) return -1;
				if (first.getTxStart()> second.getTxStart()) return 1;
				if ( first.getTxEnd() < second.getTxEnd() ) return -1;
				if ( first.getTxEnd() > second.getTxEnd()  ) return 1;
				return 0;
			}

		}

