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
	private boolean printWarnings;


	public TabixLoader (File tabixFile, ArrayList<TabixQuery> toQuery, boolean printWarnings) throws IOException{
		this.tabixFile = tabixFile;
		this.printWarnings = printWarnings;
		reader = new TabixReader( this.tabixFile.toString() );
		this.toQuery = toQuery;
	}

	public void run() {
		TabixQuery tq = null;
		try {
			int size = toQuery.size();
			for (int i=0; i< size; i++){
				tq = toQuery.get(i);
				TabixReader.Iterator it = reader.query(tq.getTabixCoordinates());
				String hit;
				ArrayList<String> al = new ArrayList<String>();
				while ((hit = it.next()) != null) al.add(hit);
				tq.addResults(tabixFile, al);
				//TODO: debug issue with tabix not returning some results
				if (al.size() == 0 && printWarnings) System.err.println("\nError: failed to return any data from "+tabixFile+" for "+tq.getInterbaseCoordinates()+" -> tbx query-> "+tq.getTabixCoordinates() );
			}
			complete = true;
		} catch (IOException e) {
			failed = true;
			System.err.println("\nError: searching "+tabixFile+" for "+tq.getTabixCoordinates() );
			e.printStackTrace();
		} finally {
			reader.close();
		}
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

}
