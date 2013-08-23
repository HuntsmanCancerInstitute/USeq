package edu.utah.ames.bioinfo;

import java.text.DecimalFormat;
import java.util.Comparator;

import util.gen.FisherExact;


public class Cartridge {
	
	public String seq;
	public int count[] = {0,0}; //count per seq
	static int tCount[] = {0,0}; //total count
	public double pValue;
	public double FDR;
	
	//constructor
	public Cartridge(String seq, int fileNum) {
		this.seq = seq;
		this.increment(fileNum);
	}
	
	public void increment(int sampleNum) {
		count[sampleNum]++;
		tCount[sampleNum]++;
	}
	
	public void calculateFisher(FisherExact fe) {
		//setup the 2x2 table
		//System.out.println(count[0] + " " + count[1] + " " + (tCount[0]-count[0]) + " " + (tCount[1]-count[1]));
		double pVal = fe.getTwoTailedP(count[0], count[1], tCount[0]-count[0], tCount[1]-count[1]);
		pValue = pVal;
	}
	
	//perform a two-tailed Fisher Exact test
	public String getFisher(double cutOff) {
		
		String out = null;
		
		if (this.FDR <= cutOff) {
			//convert p-value to scientific notation
			DecimalFormat df = new DecimalFormat("0.000E00");
			
			//return results
			out = String.format("%s\t%d\t%d\t%s\t%s\n", this.seq, count[0], count[1], df.format(FDR), df.format(pValue)); 
		}
		return out;
	}
	

	public String getSeq() {
		return seq;
	}

	public int getCount(int sampleNum) {
		return count[sampleNum];
	}
	
	//sets max matrix size
	public double getMax() {
		return count[0] + count[1] + (tCount[0]-count[0]) + (tCount[1]-count[1]);
	}

	public double getpValue() {
		return pValue;
	}

	public void setFDR(double fDR) {
		FDR = fDR;
	}
	
	static class CartridgeComparator implements Comparator<Cartridge> {

		@Override
		public int compare(Cartridge arg0, Cartridge arg1) {
			
			if (arg1.getpValue() > arg0.getpValue()) {
				return 1;
			} 
			else if (arg1.getpValue() < arg0.getpValue()) {
				return -1;
			} 
			else {
				return 0;
			}
		
			//int value = (int)(arg1.getpValue() - arg0.getpValue());
			//return value;
		}
		
	}
}
