package trans.main;
import java.util.Comparator;

/**Class to help sort an array of Interval[] by chromosome, start position, length shortest first.*/
public class IntervalComparator implements Comparator {

	/**Sorts by chromosome-> start oligo position -> length, shortest first*/
	public int compare(Object arg0, Object arg1) {
		Interval first = (Interval)arg0;
		Interval second = (Interval)arg1;
		int i = first.getChromosome().compareTo(second.getChromosome());
		if (i !=0 ) return i;
		if (first.getStart1stOligo()< second.getStart1stOligo()) return -1;
		if (first.getStart1stOligo()> second.getStart1stOligo()) return 1;
		if ( (first.getStartLastOligo() - first.getStart1stOligo()) < 
				(second.getStartLastOligo() - second.getStart1stOligo()) ) return -1;
		if ( (first.getStartLastOligo() - first.getStart1stOligo()) < 
				(second.getStartLastOligo() - second.getStart1stOligo()) ) return 1;
		return 0;
	}

}
