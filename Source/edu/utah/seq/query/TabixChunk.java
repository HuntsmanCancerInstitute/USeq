package edu.utah.seq.query;

import java.io.File;
import java.util.ArrayList;

public class TabixChunk {
	
	public File tabixFile;
	public ArrayList<TabixQuery> queries;
	
	public TabixChunk (File tabixFile, ArrayList<TabixQuery> queries){
		this.tabixFile = tabixFile;
		this.queries = queries;
	}
	
}
