package util.gen;
import java.util.Comparator;

/**Sorts by index, small to large.*/
public class RankSampleIndexComparator implements Comparator {
	public int compare(Object arg0, Object arg1) {
		RankSample first = (RankSample)arg0;
		RankSample second = (RankSample)arg1;
		if (second.index < first.index) return 1;
		if (second.index > first.index) return -1;
		return 0;
	}
	
}
