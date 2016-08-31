package edu.utah.seq.query;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import htsjdk.tribble.readers.TabixReader;

public class TabixLoader implements Runnable{

	//fields
	private File tabixFile;
	private TabixReader reader = null;
	private ArrayList<TabixQuery> toQuery;
	private boolean complete = false;
	private boolean failed = false;


	public TabixLoader (File tabixFile, ArrayList<TabixQuery> toQuery) throws IOException{
		this.tabixFile = tabixFile;
		reader = new TabixReader( this.tabixFile.toString() );
		this.toQuery = toQuery;
	}

	public void run() {
		TabixQuery tq = null;
		try {
			int size = toQuery.size();
			for (int i=0; i< size; i++){
				tq = toQuery.get(i);
				TabixReader.Iterator it = reader.query(tq.getCoordinates());
				String hit;
				ArrayList<String> al = new ArrayList<String>();
				while ((hit = it.next()) != null) al.add(hit);
			}
			reader.close();
			complete = true;
		} catch (IOException e) {
			failed = true;
			System.err.println("\n:Error searching "+tabixFile+" for "+tq.getCoordinates() );
			e.printStackTrace();
		}
	}

	public boolean isComplete() {
		return complete;
	}

	public boolean isFailed() {
		return failed;
	}

}
