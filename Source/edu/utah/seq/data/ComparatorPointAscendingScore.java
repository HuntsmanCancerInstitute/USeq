package edu.utah.seq.data;
import java.util.Comparator;

/**Class to help sort an array of Point[] by score, smallest to largest.*/
public class ComparatorPointAscendingScore implements Comparator<Point> {

	/**Sorts by position*/
	public int compare(Point first, Point second) {
		if (first.score < second.score) return -1;
		if (first.score > second.score) return 1;
		return 0;
	}
}


