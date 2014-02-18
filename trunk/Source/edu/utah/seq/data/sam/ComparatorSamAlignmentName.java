package edu.utah.seq.data.sam;
import java.util.Comparator;


/**Class to help sort an array of Point[] by score, largest to smallest.*/
public class ComparatorSamAlignmentName implements Comparator<SamAlignment> {

	/**Sorts by position*/
	public int compare(SamAlignment first, SamAlignment second) {
		return first.getName().compareTo(second.getName());
	}
}


