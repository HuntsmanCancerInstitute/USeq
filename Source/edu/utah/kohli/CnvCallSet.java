package edu.utah.kohli;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;

import util.gen.IO;
import util.gen.Misc;

public class CnvCallSet {
	
	//fields
	private String methodName;
	private String sampleId;
	private int numberCalls;
	private TreeMap<String, GeneCallResult> geneCallResults;
	private ArrayList<String> copyAlteredGenesWithDirection = null;
	
	public CnvCallSet (String methodName, String sampleId, int numberCalls, TreeMap<String, GeneCallResult> geneCallResults) {
		this.methodName = methodName;
		this.sampleId = sampleId;
		this.numberCalls = numberCalls;
		this.geneCallResults = geneCallResults;
	}

	public String getMethodName() {
		return methodName;
	}

	public String getSampleId() {
		return sampleId;
	}

	public TreeMap<String, GeneCallResult> getGeneCallResults() {
		return geneCallResults;
	}

	public String getCopyAlteredGeneString(String[] testedGenes) {
		getCopyAlteredGenes(testedGenes);
		if (copyAlteredGenesWithDirection.size() == 0) return "None";
		return Misc.stringArrayListToString(copyAlteredGenesWithDirection, ",");
	}

	public ArrayList<String> getCopyAlteredGenes(String[] testedGenes) {
		if (copyAlteredGenesWithDirection == null) {
			copyAlteredGenesWithDirection = new ArrayList<String>();

			for (String geneName: testedGenes) {
				GeneCallResult gcr = geneCallResults.get(geneName);
				if (gcr!=null && gcr.isCopyAltered()) {
					if (gcr.getIchorCode()!=null) copyAlteredGenesWithDirection.add(geneName+ gcr.getIchorCode());
					else {
						String dir = "-";
						if (gcr.isAmplified()) dir = "+";
						copyAlteredGenesWithDirection.add(geneName+dir);
					}
				}
			}
		}
		return copyAlteredGenesWithDirection;
		
	}

	public int getNumberCalls() {
		return numberCalls;
	}

	public void setNumberCalls(int numberCalls) {
		this.numberCalls = numberCalls;
	}

}
