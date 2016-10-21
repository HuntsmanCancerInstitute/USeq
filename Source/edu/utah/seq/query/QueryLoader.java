package edu.utah.seq.query;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import htsjdk.tribble.readers.TabixReader;

public class QueryLoader {
	
	private TQuery tQuery;
	
	/*Container for active TabixReaders*/
	private HashMap<File, TabixReader> tabixReaders = new HashMap<File, TabixReader>();

	private QueryRequest currentQueryRequest;
	
	
	public QueryLoader(TQuery tQuery) {
		this.tQuery = tQuery;
	}

	/**Synchronized method for getting or setting a TabixReader.
	 * This is a primitive cache of readers, only one kept per file. Would be better to time these out and permit many readers for big files, e.g. dbSNP
	 * @throws IOException */
	public synchronized TabixReader getSetTabixReader(File tabixFile, TabixReader reader) throws IOException {
		//do they want an active reader?
		if (reader == null) {
			//does it exist?
			if (tabixReaders.containsKey(tabixFile)) return tabixReaders.remove(tabixFile);
			else return new TabixReader(tabixFile.toString());
		}
		//nope, want to put one back, do so if it doesn't exist
		else {
			//save it?
			if (tabixReaders.containsKey(tabixFile) == false) {
				tabixReaders.put(tabixFile, reader);
			}
			//nope, close it
			else reader.close();
			return null;
		}
	}

	/**Method for freeing up the file handles.*/
	public void closeTabixReaders(){
		for (TabixReader tr: tabixReaders.values()) tr.close();
	}
	
	/**This is the method the TabixLoaders use to pull a chunk of data to process.
	 * @return a TabixChunk or null when nothing is left to do. This shuts down the TabixLoader.*/
	public synchronized TabixChunk getChunk() throws IOException {
		//at some point will want to queue up many QueryRequest and keep the loaders busy
		//for now all just processing one QueryRequest at a time
		return currentQueryRequest.getChunk();
	}

	/**This takes QueryRequests and loads their TabixQueries with data. */
	public void loadTabixQueriesWithData(QueryRequest qr) throws IOException {
		long startTime = System.currentTimeMillis();
		currentQueryRequest = qr;
		
		//try to make a loader for each chunk
		int numToMake= qr.getTabixChunks().size();
		if (numToMake > tQuery.getNumberThreads()) numToMake = tQuery.getNumberThreads();
		TabixLoader[] loader = new TabixLoader[numToMake];
		ExecutorService executor = Executors.newFixedThreadPool(numToMake);
		for (int i=0; i< loader.length; i++){
			loader[i] = new TabixLoader(this, tQuery.isPrintWarnings());
			executor.execute(loader[i]);
		}
		executor.shutdown();

		//spins here until the executer is terminated, e.g. all threads complete
		while (!executor.isTerminated()) {}

		//check loaders 
		for (TabixLoader c: loader) {
			if (c.isFailed()) throw new IOException("\nERROR: TabixLoader issue! \n"+c);
		}
		qr.setTime2Load(System.currentTimeMillis() -startTime);
	}

}
