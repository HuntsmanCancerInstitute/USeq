package edu.utah.seq.parsers.jpileup;

import java.util.HashSet;

public class RCBaseCount {
	int bpPosition;
	int baseCount = 0;
	HashSet<String> readNames = new HashSet<String>();

	public RCBaseCount(int bpPosition){
		this.bpPosition = bpPosition;
	}

	
}
