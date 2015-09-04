package edu.utah.seq.barcodes;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class BarcodeCluster {
	
	//fields
	private HashMap<String, char[]> familyMembers; 
	private BarcodeClusterEngine csa;
	
	//constructors
	public BarcodeCluster(String bcsToCluster, BarcodeClusterEngine csa){
		familyMembers = new HashMap<String, char[]>();
		familyMembers.put(bcsToCluster, bcsToCluster.toCharArray());
		this.csa = csa;
	}
	
	//methods
	public String toString(){
		StringBuilder sb = new StringBuilder();
		Iterator<String> keys = familyMembers.keySet().iterator();
		sb.append(keys.next());
		while (keys.hasNext()){
			sb.append(", ");
			sb.append(keys.next());
		}
		return sb.toString();
	}
	
	public boolean isMember(String testString, char[] testChar){
		//look for it in the hash
		if (familyMembers.containsKey(testString)) return true;
		
		//ok, an exact match isn't present so need to calculate the fraction similarity to all members
		//if any is less than min identity return false
		Iterator<char[]> it = familyMembers.values().iterator();
		while (it.hasNext()){
			char[] famChar = it.next();
			double fractionMatches = scoreMatches(testChar, famChar);
			if (fractionMatches < csa.minFractionIdentity) return false;
		}
		return true;
	}
	
	private double scoreMatches(char[] test, char[] member) {
		int num = member.length;
		double matches = 0;
		double nonMatches = 0;
		for (int i=0; i< num; i++){
			if (test[i] != 'N' && member[i] != 'N'){
				if (test[i] == member[i]) matches++;
				else nonMatches++;
			}
			
		}
		double total = matches+nonMatches;
		if (total < csa.minNumBases) return 0;
		return matches/(matches+nonMatches);
	}

	/**Iterates through each barcode in the cluster comparing them to this, if any are not members returns false.*/
	public boolean isMemberOf(BarcodeCluster test) {
		HashMap<String, char[]> testBarcodes = test.familyMembers;
		//for each entry
		for (Map.Entry<String, char[]> entry: testBarcodes.entrySet()){
			if (this.isMember(entry.getKey(), entry.getValue()) == false) return false;
		}
		return true;
	}

	public HashMap<String, char[]> getFamilyMembers() {
		return familyMembers;
	}

	public void setFamilyMembers(HashMap<String, char[]> familyMembers) {
		this.familyMembers = familyMembers;
	}
	
}
