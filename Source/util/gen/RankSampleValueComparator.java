package  util.gen;
import java.util.Comparator;

/**Sorts by value, small to large.*/
public class RankSampleValueComparator implements Comparator {
	public int compare(Object arg0, Object arg1) {
		RankSample first = (RankSample)arg0;
		RankSample second = (RankSample)arg1;
		if (second.value < first.value) return 1;
		if (second.value > first.value) return -1;
		return 0;
	}
	
}




