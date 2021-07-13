package edu.utah.seq.parsers.jpileup;

import util.gen.Misc;

public class BpileupLine {
	
	//fields
	private String line; 
	private String chr;
	private int zeroPos;
	private char ref;
	private BaseCount[] samples;
	
	public BpileupLine(String line) throws Exception{
		this.line = line;
		
		//chr 1BasePos Ref A,C,G,T,N,Del,Ins,FailBQ A,C,G,T,N,Del,Ins,FailBQ ...
		String[] fields = Misc.TAB.split(line);
		int numFields = fields.length;
		
		//parse standard fields
		chr = fields[0];
		zeroPos = Integer.parseInt(fields[1]) -1;
		ref = fields[2].charAt(0);
	
		//parse samples
		samples = new BaseCount[fields.length-3];
		int index = 0;
		for (int i=3; i< numFields; i++) samples[index++] = new BaseCount(zeroPos, ref, fields[i]);	
	}
	
	public String getBed(){
		StringBuilder sb = new StringBuilder(chr);
		sb.append("\t");
		sb.append(zeroPos);
		sb.append("\t");
		sb.append(zeroPos+1);
		return sb.toString();
	}
	
	//getters and setters
	public String getChr() {
		return chr;
	}
	public int getZeroPos() {
		return zeroPos;
	}
	public char getRef() {
		return ref;
	}
	public BaseCount[] getSamples() {
		return samples;
	}
	public String getLine() {
		return line;
	}
}
