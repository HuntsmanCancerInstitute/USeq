package trans.main;
import java.util.*;

import util.bio.calc.*;
import util.gen.*;

/**
 * Implementation of a  Wilcoxon Rank Sum Test.
 * Note, if the number of samples for either the treatment or control group is <10, 
 * no test statistics are calculated and the pValue is set to 1, logTPvalue is set to 0, 
 * and z to -10.
 */

public class WilcoxonRankSumTest {
	private WilcoxonSample[] samples;
	private double treatmentSum;
	private double controlSum;
	private int numTreatments;
	private int numControls;
	private double z;
	private double pValue;
	private double logTPvalue;	
	private boolean twoTailed = false; //set to true if you want to consider cases where control> treatment, by default negative z scores are ignored
	
	/**Returns a -10LogBase10(pValue) if n1 and n2 >9.
	 * if twoTailed is true and z is neg, will return a - logTPValue.*/
	public double test(float[] treatment, float[] control){
		//initialize params
		treatmentSum =0;
		controlSum =0;
		z=-10;
		pValue = 1; 
		logTPvalue =0;
		
		//make samples
		makeSamples(treatment, control);
		
		//rank samples
		Arrays.sort(samples);
		rankSamples();
		
		//sum the ranks for the treatment and control samples
		sumRanks();
		
		//calc -10*Log10(p-value)	
		if (numTreatments>9 && numControls>9){
			z = calculateZ();
			if (z>0 || twoTailed) {
				pValue = Alignment.estimatePValue(z);
				//System.out.println("z: "+z+ " p-val "+pValue);			
				if (pValue==0) logTPvalue = 90; 
				else logTPvalue = Alignment.transform10Log10(pValue);
				if (z<0) logTPvalue *= -1;
			}
			else {
				pValue = 1;
				logTPvalue = 0;
			}
		}		
		return logTPvalue;
	}
	
	/**Returns a -10LogBase10(pValue) if n1 and n2 >9  */
	public double test(float[][] treatment, float[][] control, int startIndex, int stopIndex){
		//initialize params
		treatmentSum =0;
		controlSum =0;
		z=-10;
		pValue = 1;
		logTPvalue =0;
		
		//make samples
		makeSamples(treatment, control, startIndex, stopIndex);
		
		//rank samples
		Arrays.sort(samples);
		rankSamples();
		
		//sum the ranks for the treatment and control samples
		sumRanks();
		
		//calc -10*Log10(p-value)	
		if (numTreatments>9 && numControls>9){
			z = calculateZ();
			if (z>0 || twoTailed) {
				pValue = Alignment.estimatePValue(z);
				//System.out.println("z: "+z+ " p-val "+pValue);			
				if (pValue==0) logTPvalue = 90; 
				else logTPvalue = Alignment.transform10Log10(pValue);
			}
			else {
				pValue = 1;
				logTPvalue = 0;
			}
		}		
		return logTPvalue;
	}

	/**Only valid when # treatment samples and # control samples >9 each.*/
	public double  calculateZ(){
		//for some odd reason must sum seperately before multiplying?!!
		double sum = numTreatments+numControls+1;
		double stndiv = Math.sqrt(numTreatments* numControls* sum /12);
		double top = treatmentSum - (numTreatments*(numTreatments+ numControls+1))/2;
		return top/stndiv;
	}

	public void sumRanks(){
		int numSamples = samples.length;
		for (int i= 0; i<numSamples; i++){
			if (samples[i].treatment) treatmentSum+= samples[i].rank;
			else controlSum+= samples[i].rank;
		}
	}

	/**Ranks a sorted array of WilcoxonSamples based on value.
	 * If no ties are found then this is simply their array index number+1.
	 * (ie 1,2,3,4...)
	 * If ties are encountered, ties are assigned the average of their index
	 * positions+1. (ie if index+1's: 2,3,4 have the same absolute difference, all are assigned a
	 * rank of 3).*/
	public void rankSamples(){
		int num = samples.length;
		//assign ranks as index+1
		for (int i=0; i<num; i++)samples[i].rank=i+1;
		//check for ties
		int start;
		int end;
		for (int i=0; i<num; i++){
			start = i;
			end = i;
			//advance stop until the former and latter don't have the same value and the stop
			//	of the array hasn't been reached
			while (++end < num && samples[start].value==samples[end].value){}
			//check if i was advanced
			if (end-start!=1){// ties found
				//get average of ranks		
				float ave = Num.getAverageInts((int)samples[start].rank, (int)samples[end-1].rank);
				//assign averages
				for (int x=start; x<end; x++) samples[x].rank = ave;
				//reset i
				i=end-1;
			}		
		}
	}	
	
	public void makeSamples(float[] treatment, float[] control){
		numTreatments = treatment.length;
		numControls = control.length;
		int counter = numTreatments;
		samples = new WilcoxonSample[numTreatments+numControls];
		for (int i=0; i<numTreatments; i++) samples[i] = new WilcoxonSample(treatment[i],true);
		for (int i=0; i<numControls; i++) samples[counter++] = new WilcoxonSample(control[i], false);
	}
	
	/** Builds samples from part of a float[replica][intensities] for treat and control.
	 * @param float[replicas][intensities]*/
	public void makeSamples(float[][] treatment, float[][] control, int startIndex, int stopIndex){
		int numIntensities = 1+stopIndex-startIndex;
		numTreatments = treatment.length * numIntensities;
		numControls = control.length * numIntensities;
		samples = new WilcoxonSample[numTreatments + numControls];
		//for each intensity/ oligo
		int counter = 0;
		for (int j=startIndex; j<=stopIndex; j++){
			//make treatments
			for (int i=0; i<treatment.length; i++) { 
				samples[counter++] = new WilcoxonSample(treatment[i][j],true);
			}
			for (int i=0; i<control.length; i++) {
				samples[counter++] = new WilcoxonSample(control[i][j], false);
			}
		}
	}
	
	public static void main(String[] args) {
		WilcoxonRankSumTest w = new WilcoxonRankSumTest();
		w.setTwoTailed(true);
		float[] c = {45,	60,	52,	48,	101,	47,	70,	28,	65,	57};
		float[] t = {34,	23,	36,	25,	69,	34,	26,	44,	32,	45};
		System.out.println(w.test(t,c));
	}
	/*
		pValue 0.0036423256759730016
		 	z 2.6835477583655134
		 	24.386212246711295
	*/	
	
	public double getLogTPvalue() {
		return logTPvalue;
	}
	public int getNumberTreatmentsControls(){
		return numTreatments+numControls;
	}
	public void setTwoTailed(boolean twoTailed) {
		this.twoTailed = twoTailed;
	}
}
