package edu.utah.seq.vcf;
import java.util.Comparator;

import edu.utah.seq.data.Point;

/**Class to help sort an array of VCFRecord[] by score, smallest to largest.*/
public class ComparatorVCFRecordScore implements Comparator<VCFRecord> {

	/**Sorts by score*/
	public int compare(VCFRecord first, VCFRecord second) {
		if (first.getScore() < second.getScore()) return -1;
		if (first.getScore() > second.getScore()) return 1;
		return 0;
	}
}


