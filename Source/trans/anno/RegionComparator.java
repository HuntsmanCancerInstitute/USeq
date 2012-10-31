package trans.anno;
import java.util.Comparator;

public class RegionComparator  implements Comparator {
		/**Sorts by chromosome-> start oligo position -> length, shortest first*/
		public int compare(Object arg0, Object arg1) {
			GenomicRegion first = (GenomicRegion)arg0;
			GenomicRegion second = (GenomicRegion)arg1;
			int i = first.getChromosome().compareTo(second.getChromosome());
			if (i !=0 ) return i;
			if (first.getStart()< second.getStart()) return -1;
			if (first.getEnd()> second.getEnd()) return 1;
			if ( (first.getEnd() - first.getStart()) < 
					(second.getEnd() - second.getStart()) ) return -1;
			if ( (first.getEnd() - first.getStart()) < 
					(second.getEnd() - second.getStart()) ) return 1;
			return 0;
		}

	}



