package util.gen;

/**Uses running totals to calculate the variance and standard deviation of an array of numbers.*/
public class StandardDeviation {
	
	//running fields
	long numberObservations = 0;
	long sum = 0;
	long sumSquares = 0;
	
	/**Count each number as it becomes available.*/
	public void count(double number) throws Exception{
		numberObservations++;
		sum+= number;
		sumSquares+= (number*number);
		//watch out for negative values indicating that they have exceeded the size 
		if (sumSquares < 0) throw new Exception("You have exceeded the size of the counters.");
	}

	public double getMean() {
		return ((double)sum)/((double)numberObservations);
	}

	public double getVariance() {
		return (((double)sumSquares) - ((double)sum) * getMean()) / (((double)numberObservations) - 1.0);
	}

	public double getStandardDeviation() {
		return Math.sqrt(getVariance());
	}
	
	/**@return mean tab standard deviation*/
	public String toString(){
		return getMean()+"\t"+getStandardDeviation();
	}

	public long getNumberObservations() {
		return numberObservations;
	}
	
	public double getCoefficientOfVariation(){
		return getStandardDeviation()/ getMean();
	}
	
	public double getSignal2Noise(){
		return getMean()/ getStandardDeviation();
	}
	
}
