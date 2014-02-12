package edu.utah.seq.useq.data;

import java.util.Comparator;

public class RegionScoreTextComparator implements Comparator<RegionScoreText> {
	
	/**Sorts by longest length then name*/
	public int compare(RegionScoreText a, RegionScoreText b) {
		int lenA = a.getLength();
		int lenB = b.getLength();
		if (lenA < lenB) return 1;
		if (lenA > lenB) return -1;
		return a.getText().compareTo(b.getText());
	}

}
