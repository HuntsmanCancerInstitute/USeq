package gata.plotter;

/**
 * @author Nix
 Simple object used for sorting an arrays of start stops
 */
public class StartStop implements Comparable {
	//field
	int start;
	int stop;

	public StartStop(int start, int stop) {
		this.start = start;
		this.stop = stop;
	}

	public int compareTo(Object otherObject) {
		StartStop other = (StartStop) otherObject;
		if (start < other.start)
			return -1;
		if (start > other.start)
			return 1;
		return 0;
	}
	public int[] getStartStop(){
		return new int[]{start,stop};
	}
}
