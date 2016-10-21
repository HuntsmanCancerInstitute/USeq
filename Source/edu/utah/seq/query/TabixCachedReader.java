package edu.utah.seq.query;

import htsjdk.tribble.readers.TabixReader;

public class TabixCachedReader {
	
	//fields
	TabixReader reader;
	long lastCalled;
	
	public TabixCachedReader (TabixReader reader){
		this.reader = reader;
		lastCalled = System.currentTimeMillis();
	}
	
}
