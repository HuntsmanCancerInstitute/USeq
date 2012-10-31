package trans.misc;
import java.util.*;
import trans.roc.*;

public class SgrComparator implements Comparator {
	
	/**Sorts by chromosome then position*/
	public int compare(Object arg0, Object arg1) {
		Sgr first = (Sgr)arg0;
		Sgr second = (Sgr)arg1;
		int i = first.getChromosome().compareTo(second.getChromosome());
		if (i !=0 ) return i;
		if (first.getPosition()< second.getPosition()) return -1;
		if (first.getPosition()> second.getPosition()) return 1;
		return 0;
	}
	
}


