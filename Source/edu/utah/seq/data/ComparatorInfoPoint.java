package edu.utah.seq.data;
import java.util.Comparator;

/**Class to sort Info by chromosome.*/
public class ComparatorInfoPoint implements Comparator {

	/**Sorts by text*/
	public int compare(Object arg0, Object arg1) {
		InfoPoint first = (InfoPoint)arg0;
		InfoPoint second = (InfoPoint)arg1;
		return first.getInfo().getChromosome().compareTo(second.getInfo().getChromosome());
	}

}


