package edu.utah.seq.useq.data;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import util.gen.Num;

/**Class to rapidly scan regions for intersection. Assumes interbase coordinates throughout.*/
public class IntersectingRegions{
	public int startBp;
	public int lastBp;
	public Region[] regions;
	public int[] indexedBps;
	
	public IntersectingRegions (Region[] regions){
		this.regions = regions;
		//find first and last bp
		startBp = regions[0].getStart();
		lastBp = regions[0].getStop();
		for (int i=1; i< regions.length; i++){
			if (regions[i].getStart() < startBp) startBp = regions[i].getStart();
			if (regions[i].getStop() > lastBp) lastBp = regions[i].getStop();
		}
		//flip booleans to true if in a region
		int len = lastBp - startBp;
		indexedBps = new int[len];
		Arrays.fill(indexedBps, -1);
		//for each region
		for (int i=0; i< regions.length; i++){
			int s = regions[i].getStart()- startBp;
			int e = regions[i].getStop()- startBp;
			for (int j= s; j < e; j++) indexedBps[j] = i;
		}
	}
	
	/**Returns null if no intersection or a hash of region indexes intersected by these coordinates.*/
	public HashSet<Integer> fetchIntersectingIndexes (int start, int end){
		//watch for out of range
		if (end <= startBp || start >= lastBp) return null;
		//watch for overlaps
		if (start < startBp) start = startBp;
		if (end > lastBp) end = lastBp;
		start = start - startBp;
		end = end - startBp;
		//scan
		HashSet<Integer> hash = new HashSet<Integer>();
		for (int i=start; i< end; i++) if (indexedBps[i] != -1) hash.add(indexedBps[i]);
		if (hash.size() == 0) return null;
		return hash;
		
	}
	
	public static HashMap<String, IntersectingRegions> parseRegions(File bedFormat){
		HashMap<String, Region[]> regions = Region.parseStartStops(bedFormat, 0, 0, 1);
		HashMap<String, IntersectingRegions> br = new HashMap<String, IntersectingRegions>();
		for (String chr: regions.keySet()){
			br.put(chr, new IntersectingRegions(regions.get(chr)));
		}
		return br;
	}

	public Region[] getRegions() {
		return regions;
	}
}