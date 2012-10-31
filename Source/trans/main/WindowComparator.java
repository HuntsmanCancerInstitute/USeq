package trans.main;
import java.util.Comparator;

/**Class to help sort an array of Window[] by chromosome, start position, length shortest first.*/
public class WindowComparator implements Comparator {

	/**Sorts by chromosome-> start oligo position -> length, shortest first*/
	public int compare(Object arg0, Object arg1) {
		Window first = (Window)arg0;
		Window second = (Window)arg1;
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
