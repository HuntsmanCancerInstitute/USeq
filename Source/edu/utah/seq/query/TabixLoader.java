package edu.utah.seq.query;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import htsjdk.tribble.readers.TabixReader;
import util.gen.Misc;

public class TabixLoader implements Runnable{

	//fields
	private File tabixFile = null;
	private TabixReader reader = null;
	private ArrayList<TabixQuery> toQuery;
	private boolean failed = false;
	private boolean printWarnings;
	private TQuery tQuery;
	private long numberRetrievedResults = 0;
	private int numberQueriesWithResults = 0;
	private long numberSavedResults = 0;
	private boolean matchVcf = false;

	public TabixLoader (TQuery tQuery) throws IOException{
		this.tQuery = tQuery;
		this.printWarnings = tQuery.isPrintWarnings();
		this.matchVcf = tQuery.getTabixFilter().isMatchVcf();
	}

	public void run() {	
		TabixQuery tq = null;
		try {
			//get next chunk of queries and their associated file
			while (getChunk()){ 

				int size = toQuery.size();
				for (int i=0; i< size; i++){
					tq = toQuery.get(i);			
					TabixReader.Iterator it = reader.query(tq.getTabixCoordinates());
					String hit;

					ArrayList<String> al = new ArrayList<String>();
					int numRes = 0;
					//check that the vcf pos, ref, alt is the same, only one of the alts need match
					if (matchVcf && tabixFile.getName().endsWith(".vcf.gz")){
						while ((hit = it.next()) != null) {
							String[] t = Misc.TAB.split(hit);
							if (tq.compareVcf(t[1], t[3], Misc.COMMA.split(t[4]))) {
								al.add(hit);
								numberSavedResults++;
							}
							numRes++;
						}
					}
					//nope add everything
					else {
						while ((hit = it.next()) != null) al.add(hit);
						numRes = al.size();
						numberSavedResults+= numRes;
					}
					tq.addResults(tabixFile, al);
					
					//stats
					if (numRes != 0) {
						numberQueriesWithResults++;
						numberRetrievedResults += numRes;
					}
					//good to know if a query was created yet no results came back, this should be rare, comes form how tabix defines the foot print and how TQuery does it.
					else if (printWarnings) System.err.println("\tWARNING: failed to return any data from "+tabixFile+" for region "+tq.getInterbaseCoordinates()+" -> converted tbx query-> "+tq.getTabixCoordinates() );
				}
			}		
		} catch (IOException e) {
			failed = true;
			System.err.println("\nError: searching "+tabixFile+" for "+tq.getTabixCoordinates() );
			e.printStackTrace();
		}
	}


	/**This gets a chunk of data to process as well as a Tabix reader.*/
	private boolean getChunk() throws IOException {
		TabixChunk c = tQuery.getChunk();

		//anything there?
		if (c == null) {
			//return the reader, this is important since TQuery handles the close()
			tQuery.getSetTabixReader(tabixFile, reader);
			//signal a loader shutdown
			return false;
		}

		toQuery = c.queries;

		//first launch?
		if (tabixFile == null){
			tabixFile = c.tabixFile;
			reader = tQuery.getSetTabixReader(tabixFile, null);
			if (reader == null) reader = new TabixReader(tabixFile.toString());
		}

		//need new reader?
		else if (tabixFile.equals(c.tabixFile) == false){
			//return old reader, this is important since TQuery handles the close()
			tQuery.getSetTabixReader(tabixFile, reader);
			//get new reader
			tabixFile = c.tabixFile;
			reader = tQuery.getSetTabixReader(tabixFile, null);
			if (reader == null) reader = new TabixReader(tabixFile.toString());
		}
		return true;
	}

	public String toString(){
		return "TabixLoader:\n\t"+tabixFile.toString()+"\n\t"+toQuery.size();
	}
	public boolean isFailed() {
		return failed;
	}

	public long getNumberRetrievedResults() {
		return numberRetrievedResults;
	}

	public long getNumberQueriesWithResults() {
		return numberQueriesWithResults;
	}

	public long getNumberSavedResults() {
		return numberSavedResults;
	}

}
