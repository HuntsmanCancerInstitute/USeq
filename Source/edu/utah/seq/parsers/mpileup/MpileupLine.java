package edu.utah.seq.parsers.mpileup;

import java.io.BufferedReader;
import java.io.File;

import util.gen.IO;
import util.gen.Misc;

public class MpileupLine {
	
	//fields
	private String line; 
	private String chr;
	private int zeroPos;
	private String ref;

	private MpileupSample[] samples;
	private int minBaseQuality;
	
	public MpileupLine(String line, int minBaseQuality) throws Exception{
		this.line = line;
		this.minBaseQuality = minBaseQuality;
		String[] fields = Misc.TAB.split(line);
		int numFields = fields.length;
		
		//watch out for lines where last is 0 obs and no bases or quals so skip it
		int fieldsToParse = numFields-3;
		if (fieldsToParse % 3 !=0) {
			numFields--;
			fieldsToParse--;
			//other junkers
			if (fieldsToParse % 3 !=0) {
				System.err.println("Malformed, skipping: "+line);
				return;
			}
		}
		
		//parse standard fields
		chr = fields[0];
		zeroPos = Integer.parseInt(fields[1]) -1;
		ref = fields[2].toUpperCase();
	
		//parse samples
		int numSamples = fieldsToParse/3;
		samples = new MpileupSample[numSamples];
		numSamples = 0;
		for (int i=3; i< numFields; i+=3) {
			samples[numSamples++] = new MpileupSample(fields[i], fields[i+1], fields[i+2], this);
		}
		
	}

	/*
	public static void main(String[] args) throws Exception{
		BufferedReader in = IO.fetchBufferedReader(new File ("/Users/u0028003/Downloads/All8_143351842.mpileup"));
		
		//String line = "chrX	156	A	11	.$......+2AG.+2AG.+2AGGG	<975;:<<<<<";
		//line= "seq3	200	A	20	,,,,,..,.-4CACC.-4CACC....,.,,.^~.	==<<<<<<<<<<<::<;2<<";
		//line= "seq1	276	G	22	...T,,.,.,...,,,.,....	33;+<<7=7<<7<&<<1;<<6<";
		String line;
		long counter = 0;
		while ((line=in.readLine())!= null){
			//System.out.print("\n"+line);
			MpileupLine ml = new MpileupLine(line, 20);
			//ml.getSamples()[0].debug();
			//if ((ml.getZeroPos()+1) == 143351882) break;
			counter++;
		}
		in.close();
		System.out.println(counter);
	}*/
	
	//getters and setters
	public String getChr() {
		return chr;
	}
	public int getZeroPos() {
		return zeroPos;
	}
	public String getRef() {
		return ref;
	}
	public MpileupSample[] getSamples() {
		return samples;
	}
	public String getLine() {
		return line;
	}

	public int getMinBaseQuality() {
		return minBaseQuality;
	}
}
