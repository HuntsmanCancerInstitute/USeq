package edu.utah.seq.analysis.multi;

import java.io.*;

import util.gen.IO;


public class Condition implements Serializable{

	//fields
	private String name;
	private Replica[] replicas;
	private static final long serialVersionUID = 1L;
	
	//constructor
	public Condition (File conditionDirectory){
		name = conditionDirectory.getName().replaceAll(" ", "_");
		File[] bamFiles = IO.extractFiles(conditionDirectory, ".bam");
		replicas = new Replica[bamFiles.length];
		if (bamFiles.length == 1) replicas[0] = new Replica(name, bamFiles[0]);
		else for (int i=0; i< bamFiles.length; i++) replicas[i] = new Replica(name+i, bamFiles[i]);
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder("\t");
		sb.append(name);
		sb.append("\n");
		for (Replica tcr: replicas){
			sb.append("\t\t");
			sb.append(tcr.getNameNumber());
			sb.append("\t");
			sb.append(tcr.getBamFile().toString());
			sb.append("\t");
			sb.append(tcr.getTotalCounts());
			sb.append("\n");
		}
		return sb.toString();
	}

	public String getName() {
		return name;
	}

	public Replica[] getReplicas() {
		return replicas;
	}

}


