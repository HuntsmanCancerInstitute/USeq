package trans.cel;
import java.util.*;

/**For sorting ComparatorQuantile objects by value or chromosome number and position*/
public class ComparatorQuantile implements Comparator {
	
	//fields
	private boolean sortByValue = true;
	
	/**Set boolean sortByValue = true to sort by value, false to sort by chromosome number and position.*/
	public ComparatorQuantile(boolean sortByValue){
		this.sortByValue = sortByValue;
	}
	
	/**Sorts by value (default) or position based on sortByValue boolean, all smallest to biggest.*/
	public int compare(Object arg0, Object arg1) {
		Quantile first = (Quantile)arg0;
		Quantile second = (Quantile)arg1;
		if (sortByValue){
			//sort by value, smallest to biggest
			if (first.value > second.value) return 1;
			if (first.value < second.value) return -1;
		}
		else{
			//sort by position
			if (first.position > second.position ) return 1;
			if (first.position < second.position ) return -1;			
		}
		//they are the same
		return 0;
	}

	public void setSortByValue(boolean sortByValue) {
		this.sortByValue = sortByValue;
	}
	
	
	
}
