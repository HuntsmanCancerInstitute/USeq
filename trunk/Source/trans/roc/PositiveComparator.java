package trans.roc;
import java.util.Comparator;

public class PositiveComparator   implements Comparator {
	/**Sorts by chromosome-> start oligo position -> length, shortest first*/
	public int compare(Object arg0, Object arg1) {
		Positive first = (Positive)arg0;
		Positive second = (Positive)arg1;
		int i = first.getChromosome().compareTo(second.getChromosome());
		if (i !=0 ) return i;
		if (first.getStart()< second.getStart()) return -1;
		if (first.getStop()> second.getStop()) return 1;
		if ( (first.getStop() - first.getStart()) < 
				(second.getStop() - second.getStart()) ) return -1;
		if ( (first.getStop() - first.getStart()) < 
				(second.getStop() - second.getStart()) ) return 1;
		return 0;
	}
	
}