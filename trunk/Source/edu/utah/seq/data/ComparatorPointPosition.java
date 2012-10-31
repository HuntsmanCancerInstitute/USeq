package edu.utah.seq.data;
import java.util.Comparator;

/**Class to help sort an array of Point[] by position, smallest to largest.*/
public class ComparatorPointPosition implements Comparator<Point> {

	/**Sorts by position*/
	public int compare(Point first, Point second) {

		if (first.position < second.position) return -1;
		if (first.position > second.position) return 1;
		return 0;
	}

}


