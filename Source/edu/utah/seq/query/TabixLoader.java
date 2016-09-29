package edu.utah.seq.query;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import htsjdk.tribble.readers.TabixReader;

public class TabixLoader implements Runnable{

	//fields
	private File tabixFile = null;
	private TabixReader reader = null;
	private ArrayList<TabixQuery> toQuery;
	private boolean complete = false;
	private boolean failed = false;
	private boolean printWarnings;
	private TQuery tQuery;
	private long numberRetrievedResults = 0;
	private int numberQueriesWithResults = 0;


	public TabixLoader (TQuery tQuery) throws IOException{
		this.tQuery = tQuery;
		this.printWarnings = tQuery.isPrintWarnings();
		
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
					//TODO: these shouldn't be saved but streamed out?
					//while ((hit = it.next()) != null){};
					ArrayList<String> al = new ArrayList<String>();
					while ((hit = it.next()) != null) al.add(hit);
					tq.addResults(tabixFile, al);
					//stats
					int numRes = al.size();
					numberRetrievedResults += numRes;
					if (numRes != 0) numberQueriesWithResults++;
					else if (printWarnings) System.err.println("\tWARNING: failed to return any data from "+tabixFile+" for region "+tq.getInterbaseCoordinates()+" -> converted tbx query-> "+tq.getTabixCoordinates() );
					
				}
			}
			complete = true;
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
	public boolean isComplete() {
		return complete;
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

}
