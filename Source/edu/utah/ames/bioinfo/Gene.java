package edu.utah.ames.bioinfo;

import java.util.ArrayList;

/**
 * 
 * @author darren.ames@hci.utah.edu
 *
 */

public class Gene {

	//fields
	ArrayList<Double>[] vals;
	
	//constructor
	public Gene() {
		
		vals = new ArrayList[2];
		vals[0] = new ArrayList<Double>();
		vals[1] = new ArrayList<Double>();
	}
	
	public void addIndex(int index, double value) {
		vals[index].add(value);
	}
	
	public double[] getArray(int index) {
		Double[] retVal = new Double[vals[index].size()];
		double[] retVal2 = new double[vals[index].size()];
		vals[index].toArray(retVal);
		for (int i = 0; i < retVal.length; i++) {
			retVal2[i] = (double)retVal[i];
		}
		return retVal2;
	}
	
	public double getArray(int index, int index2) {
		return vals[index].get(index2);
	}
	
	public boolean isValid(int index) {
		
		boolean valid = false;
		for (double v : vals[index]) {
			if (v != 0.0) {
				valid = true;
			}
		}
		return valid;
	}
}
