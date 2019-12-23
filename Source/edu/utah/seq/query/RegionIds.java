package edu.utah.seq.query;

import java.util.HashSet;
import util.gen.Misc;

/**For loading prior bed indexes.*/
public class RegionIds {
	int start;
	int stop;
	HashSet<Integer> ids = new HashSet<Integer>();
	public RegionIds(String line) {
		String[] fields = Misc.TAB.split(line);
		start = Integer.parseInt(fields[1]);
		stop = Integer.parseInt(fields[2]);
		for (String s: Misc.COMMA.split(fields[3])) ids.add(Integer.parseInt(s));
	}
}