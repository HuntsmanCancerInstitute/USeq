package trans.roc;
import java.util.Comparator;

public class GrComparator   implements Comparator {
	
	/**Sorts by start position */
	public int compare(Object arg0, Object arg1) {
		Gr first = (Gr)arg0;
		Gr second = (Gr)arg1;
		if (first.getPosition() < second.getPosition()) return -1;
		if (first.getPosition() > second.getPosition()) return 1;
		return 0;
	}
	
}