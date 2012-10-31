package util.bio.parsers;
import java.util.*;

public class UCSCGeneLineChromComparator  implements Comparator {

		/**Sorts by chromosome-> start position -> length, shortest first*/
		public int compare(Object arg0, Object arg1) {
			UCSCGeneLine first = (UCSCGeneLine)arg0;
			UCSCGeneLine second = (UCSCGeneLine)arg1;
			int i = first.getChrom().compareTo(second.getChrom());
			if (i !=0 ) return i;
			if (first.getTxStart()< second.getTxStart()) return -1;
			if (first.getTxStart()> second.getTxStart()) return 1;
			if ( first.getTxEnd() < second.getTxEnd() ) return -1;
			if ( first.getTxEnd() > second.getTxEnd()  ) return 1;
			return 0;
		}

	}
