package edu.utah.seq.query;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import htsjdk.tribble.readers.TabixReader;
import util.gen.Misc;

public class TabixLoader implements Runnable{

	//fields
	private boolean printWarnings;
	private boolean failed = false;
	private QueryLoader queryLoader;
	
	//working fields
	private TabixReader reader = null;
	private File tabixFile = null;
	private TabixQuery tq = null;
	private TabixChunk tc = null;
	private ArrayList<TabixQuery> toQuery;
	private long numberRetrievedResults = 0;
	private int numberQueriesWithResults = 0;
	private long numberSavedResults = 0;
	private boolean matchVcf = false;
	

	public TabixLoader (QueryLoader queryLoader, boolean printWarnings) throws IOException{
		this.queryLoader = queryLoader;
		this.printWarnings = printWarnings;
	}
	
	public void run() {	
		try {
			//get next chunk of work
			while ((tc = queryLoader.getChunk()) != null){ 
				loadWorkChunk();
				
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
				
				//update the QueryRequest
				updateQueryStats();
			}	
			//return the reader
			if (tabixFile != null) queryLoader.getSetTabixReader(tabixFile, reader);
		} catch (IOException e) {
			failed = true;
			System.err.println("\nError: searching "+tabixFile+" for "+tq.getTabixCoordinates() );
			e.printStackTrace();
		}
	}
	
	/*This loads a chunk of work to do and deals with the readers, returning the old or fetching a new one.*/
	private void loadWorkChunk() throws IOException {
		//check reader
		if (reader != null){
			//fetch new?
			if (tabixFile.equals(tc.getTabixFile()) == false){
				//return old
				queryLoader.getSetTabixReader(tabixFile, reader);
				//get new
				tabixFile = tc.getTabixFile();
				reader = queryLoader.getSetTabixReader(tabixFile, null);
			}
		}
		else {
			tabixFile = tc.getTabixFile();
			reader = queryLoader.getSetTabixReader(tabixFile, null);
		}
		//set working objects
		toQuery = tc.getQueries();
		matchVcf = tc.getQueryRequest().getFilter().isMatchVcf();
		
		//reset counters
		numberRetrievedResults = 0;
		numberQueriesWithResults = 0;
		numberSavedResults = 0;
	}
	
	private void updateQueryStats() {
		QueryRequest qr = tc.getQueryRequest();
		qr.incrementNumQueriesWithResults(numberQueriesWithResults);
		qr.incrementNumRetrievedResults(numberRetrievedResults);
		qr.incrementNumSavedResults(numberSavedResults);
	}

	public String toString(){
		return "TabixLoader:\n\t"+tabixFile.toString()+"\n\t"+toQuery.size();
	}
	public boolean isFailed() {
		return failed;
	}
}
