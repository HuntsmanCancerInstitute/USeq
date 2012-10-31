package trans.anno;
import java.util.Comparator;

public class BindingRegionComparator  implements Comparator {
		/**Sorts by rank*/
		public int compare(Object arg0, Object arg1) {
			BindingRegion first = (BindingRegion)arg0;
			BindingRegion second = (BindingRegion)arg1;
			if (first.getRank() < second.getRank()) return -1;
			if (first.getRank() > second.getRank()) return 1;
			return 0;
		}

	}