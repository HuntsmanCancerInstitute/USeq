package edu.utah.seq.query;

import java.io.File;
import java.util.ArrayList;

public class TabixChunk {
	
	private File tabixFile;
	private ArrayList<TabixQuery> queries;
	private QueryRequest queryRequest;
	
	public TabixChunk (File tabixFile, ArrayList<TabixQuery> queries, QueryRequest queryRequest){
		this.tabixFile = tabixFile;
		this.queries = queries;
		this.queryRequest = queryRequest;
	}

	public File getTabixFile() {
		return tabixFile;
	}

	public void setTabixFile(File tabixFile) {
		this.tabixFile = tabixFile;
	}

	public ArrayList<TabixQuery> getQueries() {
		return queries;
	}

	public void setQueries(ArrayList<TabixQuery> queries) {
		this.queries = queries;
	}

	public QueryRequest getQueryRequest() {
		return queryRequest;
	}

	public void setQueryRequest(QueryRequest queryRequest) {
		this.queryRequest = queryRequest;
	}
	
}
