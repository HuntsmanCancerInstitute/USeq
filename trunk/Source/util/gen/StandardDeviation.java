package util.gen;

import java.io.Serializable;

/**Uses running totals to calculate the variance and standard deviation of an array of numbers.*/
public class StandardDeviation implements Serializable {
	
	//running fields
	double numberObservations = 0;
	double sum = 0;
	double sumSquares = 0;
	double mean = Double.MIN_NORMAL;
	double standardDeviation = Double.MIN_NORMAL;
	private static final long serialVersionUID = 1L;
	
	/**Count each number as it becomes available.*/
	public void count(double number){
		numberObservations++;
		sum+= number;
		sumSquares+= (number*number);
		//watch out for negative values indicating that they have exceeded the size 
		if (sumSquares < 0) System.err.println("You have exceeded the size of the StandardDeviation counters with value "+number);
	}
	
	public void count(StandardDeviation other){
		numberObservations+= other.numberObservations;
		sum+= other.sum;
		sumSquares+= other.sumSquares;
		//watch out for negative values indicating that they have exceeded the size 
		if (sumSquares < 0) System.err.println("You have exceeded the size of the StandardDeviation counter!");
	}

	public double getMean() {
		if (mean == Double.MIN_NORMAL) mean = sum/numberObservations;
		return mean;
	}

	public double getVariance() {
		return (sumSquares - sum * getMean()) / (numberObservations - 1.0);
	}

	public double getStandardDeviation() {
		if (standardDeviation == Double.MIN_NORMAL) standardDeviation = Math.sqrt(getVariance());
		return standardDeviation;
	}
	
	/**@return mean tab standard deviation*/
	public String toString(){
		return getMean()+"\t"+getStandardDeviation();
	}

	public double getNumberObservations() {
		return numberObservations;
	}
	
	public double getCoefficientOfVariation(){
		return getStandardDeviation()/ getMean();
	}
	
	public double getSignal2Noise(){
		return getMean()/ getStandardDeviation();
	}
	/**returns (Score - Mean) / Standard Deviation*/
	public double getZScore(double score){
		if (mean == Double.MIN_NORMAL){
			getMean();
			getStandardDeviation();
		}
		return (score - mean) / standardDeviation;
	}
	
}
