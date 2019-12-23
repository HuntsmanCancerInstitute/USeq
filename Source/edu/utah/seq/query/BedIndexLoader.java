package edu.utah.seq.query;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import util.gen.IO;

public class BedIndexLoader {
		
		private BufferedReader in;
		private RegionIds lastRegionIds = null;
		private ArrayList<IndexRegion>[] workingIndex  = null;
		private HashSet<Integer> oldDataSourceIDsToDelete = null;
		private int numOldIds = 0;
		private ArrayList<RegionIds> riBlock =  new ArrayList<RegionIds>();
		private long numRegionsLoaded = 0;
		
		//constructor
		public BedIndexLoader(File bed, HashSet<Integer> oldDataSourceIDsToDelete, ArrayList<IndexRegion>[] workingIndex) throws IOException {
			this.oldDataSourceIDsToDelete = oldDataSourceIDsToDelete;
			this.workingIndex = workingIndex;
			//oldDataSourceIDsToDelete.add(4);
			//oldDataSourceIDsToDelete.add(8);
			numOldIds = oldDataSourceIDsToDelete.size();
			in = IO.fetchBufferedReader(bed);
			 
			while (true) {
				fetchBlock(riBlock);
				if (riBlock.size()==0) {
					in.close();
					break;
				}
				else processBlock();
			} 
		}
		
		public void fetchBlock(ArrayList<RegionIds> regionIds) throws IOException {
			regionIds.clear();
			String line;

			while ((line = in.readLine())!= null) {
				//remove old ids
				RegionIds ris = new RegionIds(line);
				if (numOldIds !=0) {
					ris.ids.removeAll(oldDataSourceIDsToDelete);
					if (ris.ids.size() == 0) continue;
				}
				//first one?
				if (lastRegionIds  == null) lastRegionIds = ris;
				else {
					regionIds.add(lastRegionIds);
					//continuous block?
					if (lastRegionIds.stop == ris.start) lastRegionIds = ris;
					//nope break
					else {
						lastRegionIds = ris;
						return;
					}
				} 
			}
			//last one?
			if (lastRegionIds != null) {
				regionIds.add(lastRegionIds);
				lastRegionIds =  null;
			}
		}
	
	private void processBlock() {

		/*
		IO.pl("\nBlocks");
		for (RegionIds ris: riBlock) {
			IO.pl("\t"+ris.start+"-"+ris.stop+":"+ris.ids);
		}*/

		
		int numInBlock = riBlock.size();
		//just one line?
		if (numInBlock == 1) {
			addIndexRegions(riBlock.get(0));
		}
		else {
			//multi line, block here!
			for (int x=0; x< numInBlock; x++) {
				//start new
				RegionIds ri = riBlock.get(x);
//IO.pl("Starting at top "+x+" "+ri.start+"-"+ri.stop+"\t"+ri.ids);
				
				//any ids left?
				int numIds = ri.ids.size();
				if (numIds == 0) {
//IO.pl("\tNo Ids, continuing");
					continue;
				}
				
				int start = ri.start;
				int stop = ri.stop;
				//for each id
				Iterator<Integer> it = ri.ids.iterator();
				while (it.hasNext()) {


					int id = it.next();
					//ri.ids.remove(id);
//IO.pl("\tProc Ids "+id);

					//for each region downstream, go till break
					boolean saveIt = true;
					for (int i=x+1; i< numInBlock; i++) {
						RegionIds downstreamRI = riBlock.get(i);
//IO.pl("\tInner "+x+" "+downstreamRI.start+"-"+downstreamRI.stop+"\t"+downstreamRI.ids);
						//present in downstream?
						if (downstreamRI.ids.contains(id)) {
//IO.pl("\t\tContains");
							downstreamRI.ids.remove(id);
							stop = downstreamRI.stop;
						}
						//not present or no id's
						else {
							//save region
							IndexRegion ir = new IndexRegion(start, stop, id);
//IO.pl("\t\t\tSaving "+ir);
							numRegionsLoaded++;
							//add start and stop irs
							if (workingIndex[start] == null) workingIndex[start] = new ArrayList<IndexRegion>();
							workingIndex[start].add(ir);
							if (workingIndex[stop] == null) workingIndex[stop] = new ArrayList<IndexRegion>();
							workingIndex[stop].add(ir);
							saveIt = false;
							break;
						}
					}
					//save it?
					if (saveIt) {
						//save region
						IndexRegion ir = new IndexRegion(start, stop, id);
//IO.pl("\t\t\tSaving last "+ir);
						numRegionsLoaded++;
						//add start and stop irs
						if (workingIndex[start] == null) workingIndex[start] = new ArrayList<IndexRegion>();
						workingIndex[start].add(ir);
						if (workingIndex[stop] == null) workingIndex[stop] = new ArrayList<IndexRegion>();
						workingIndex[stop].add(ir);
					}
				}
			}

		}
		
	}
	private void addIndexRegions(RegionIds regionIds) {
		//for each id, add a IndexRegion
		Iterator<Integer> iter = regionIds.ids.iterator();
		int start  = regionIds.start;
		int stop = regionIds.stop;
		while (iter.hasNext()) {
			IndexRegion ir = new IndexRegion(start, stop, iter.next());
//IO.pl("\t\t\tSaving "+ir);
			numRegionsLoaded++;
			//add start and stop irs
			if (workingIndex[start] == null) workingIndex[start] = new ArrayList<IndexRegion>();
			workingIndex[start].add(ir);
			if (workingIndex[stop] == null) workingIndex[stop] = new ArrayList<IndexRegion>();
			workingIndex[stop].add(ir);
		}
	}

	public long getNumRegionsLoaded() {
		return numRegionsLoaded;
	}

	
	public static void main(String[] args)  {

		try {
			ArrayList<IndexRegion>[] workingIndex  = new ArrayList[100000000];
			HashSet<Integer> oldDataSourceIDsToDelete = new HashSet<Integer>();
			//File regions = new File("/Users/u0028003/Downloads/QueryIndexDev/QueryIndex/20.qi.bed");
			//File  regions = new File("/Users/u0028003/Downloads/QueryIndexDev/testDelme.txt");
			File regions = new File("/Users/u0028003/Downloads/QueryIndexDev/testBlockParsing.txt");
			
			BedIndexLoader dm = new BedIndexLoader(regions, oldDataSourceIDsToDelete, workingIndex);
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}






}

