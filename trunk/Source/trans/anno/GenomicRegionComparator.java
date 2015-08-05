package trans.anno;
import java.util.Comparator;

public class GenomicRegionComparator  implements Comparator<GenomicRegion> {
		/**Sorts by chromosome-> start oligo position -> length, shortest first*/
		public int compare(GenomicRegion first, GenomicRegion second) {
			int i = first.getChromosome().compareTo(second.getChromosome());
			if (i !=0 ) return i;
			if (first.getStart()< second.getStart()) return -1;
			if (first.getStart()> second.getStart()) return 1;
			//ok same start now check length
			if ( (first.getEnd() - first.getStart()) < 
					(second.getEnd() - second.getStart()) ) return -1;
			if ( (first.getEnd() - first.getStart()) > 
					(second.getEnd() - second.getStart()) ) return 1;
			return 0;
		}

	}



