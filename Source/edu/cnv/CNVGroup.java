package edu.cnv;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import util.gen.IO;
import util.gen.Misc;

/**Represents CNVs from a particular dataset coming from the CNVScanner.*/
public class CNVGroup{
	private String name;
	private CNV[] upCNVs = null;
	private CNV[] downCNVs = null;
	private static final ComparatorCNVChromPosition comp = new ComparatorCNVChromPosition();

	public String toString(){
		StringBuilder sb = new StringBuilder(name);
		if (upCNVs!=null){
			sb.append("\n\tUp:\n\t\t");
			for (int i=0; i< upCNVs.length; i++){
				sb.append(upCNVs[i].toStringBed(""));
				sb.append(", ");
			}
		}
		if (downCNVs!=null){
			sb.append("\n\tDn:\n\t\t");
			for (int i=0; i< downCNVs.length; i++){
				sb.append(downCNVs[i].toStringBed(""));
				sb.append(", ");
			}
		}
		return sb.toString();
	}

	public CNVGroup (File cnvFile, double scoreDivider){
		name = Misc.removeExtension(cnvFile.getName());
		CNV[] all = (CNV[])IO.fetchObject(cnvFile);
		//split
		ArrayList<CNV> up = new ArrayList<CNV>();
		ArrayList<CNV> down = new ArrayList<CNV>();
		for (int i=0; i< all.length; i++){
			if ( all[i].getMedianScore() > scoreDivider) up.add(all[i]);
			else down.add(all[i]);
		}
		if (up.size() !=0) {
			upCNVs = new CNV[up.size()];
			up.toArray(upCNVs);
			Arrays.sort(upCNVs,comp);
		}
		if (down.size() !=0) {
			downCNVs = new CNV[down.size()];
			down.toArray(downCNVs);
			Arrays.sort(downCNVs, comp);
		}
	}

	public String getName() {
		return name;
	}

	public CNV[] getUpCNVs() {
		return upCNVs;
	}

	public CNV[] getDownCNVs() {
		return downCNVs;
	}
}

