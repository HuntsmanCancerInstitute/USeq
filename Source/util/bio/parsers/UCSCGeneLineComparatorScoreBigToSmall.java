package util.bio.parsers;

import java.util.Comparator;

/**Sorts by score, largest to smallest.*/
public class UCSCGeneLineComparatorScoreBigToSmall implements Comparator {
	
	private int scoreIndex;
	
	public UCSCGeneLineComparatorScoreBigToSmall (int scoreIndex){
		this.scoreIndex = scoreIndex;
	}
	
	public int compare(Object arg0, Object arg1) {
		float first = ((UCSCGeneLine)arg0).getScores()[scoreIndex];
		float second = ((UCSCGeneLine)arg1).getScores()[scoreIndex];
		if (first > second) return -1;
		if (first < second) return 1;
		return 0;
	}

}
