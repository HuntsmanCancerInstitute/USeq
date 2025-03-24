package edu.utah.kegg.pathway;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;

import edu.utah.kegg.api.KeggApiNetwork;
import edu.utah.kegg.api.KeggApiPathway;
import util.gen.Misc;
import util.gen.Num;

public class CombinePathway {
	
	//Common 
	private KeggApiPathway keggPathway = null;
	private double combinePValue = -1;
	private double combineFdr = -1;
	private TreeMap<String, AnalyzedNetwork> networkIdsAnalyzedNetworks = new TreeMap<String, AnalyzedNetwork>();

	//Constructor
	public CombinePathway (KeggApiPathway keggPathway) {
		this.keggPathway = keggPathway;
		
	}
	
	public double getCombinePValue() {
		return combinePValue;
	}
	public void setCombinePValue(double combinePValue) {
		this.combinePValue = combinePValue;
	}
	public double getCombineFdr() {
		return combineFdr;
	}
	public void setCombineFdr(double combineFdr) {
		this.combineFdr = combineFdr;
	}
	public KeggApiPathway getKeggPathway() {
		return keggPathway;
	}

	public TreeMap<String, AnalyzedNetwork> getNetworkIdsAnalyzedNetworks() {
		return networkIdsAnalyzedNetworks;
	}
	
}
