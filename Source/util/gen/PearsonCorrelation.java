package util.gen;

/**Calculates Pearson Correlations (R) in two ways, using a straight static method on two arrays or
 * by 1st adding paired arrays sequentially and 2nd calculating the correlation on the cumulative 
 * sum.  This later method is good if you want to keep the memory requirements down.  Otherwise 
 * just use the static methods.*/
public class PearsonCorrelation {

	//fields
	private double xTot=0;
	private double yTot=0;
	private double sqrXTot =0;
	private double sqrYTot =0;
	private double xYTot=0;
	private double totalN = 0;

	/**Use this method to sequentially add multiple pairs for correlation.*/
	public void addPairsToCorrelate (float[] x, float[] y){
		double N = x.length;
		for (int i=0; i<N; i++){
			xTot += x[i];
			yTot += y[i];
			sqrXTot += Math.pow(x[i],2);
			sqrYTot += Math.pow(y[i],2);
			xYTot += (x[i] * y[i]);
		} 
		totalN += N;
	}

	/**Use this method to sequentially add pairs for correlation.*/
	public void addPairsToCorrelate (float x, float y){
		xTot += x;
		yTot += y;
		sqrXTot += Math.pow(x,2);
		sqrYTot += Math.pow(y,2);
		xYTot += (x * y);
		totalN++;
	}

	/**Following loading of this class using the addPairsToCorrelate(),
	 * call this method to calculate to correlation coefficient (R) on the total.
	 * Returns null if an error is encountered.*/
	public Double calculateAdditivePairCorrelation(){
		double top = (totalN * xYTot) - (xTot * yTot);
		double botLeft = Math.sqrt( (totalN * sqrXTot) - Math.pow(xTot,2) );
		double botRight = Math.sqrt( (totalN * sqrYTot) - Math.pow(yTot,2) );
		double bot = botLeft*botRight;
		if (bot == 0) {
			System.err.println("\nERROR calculating correlation:");
			System.err.println("\tTop "+top);
			System.err.println("\tBotL "+botLeft);
			System.err.println("\tBotR "+botRight);
			System.err.println("\tBot "+bot);
			return null;
		}
		return new Double(top/bot);
	}

	/**Calculates Pearson correlation coefficient, r, from two float[]s. 
	 * Cannot have one float[] be uniform values, returns -2 if error.*/
	public static double correlationCoefficient (float[] x, float[] y){
		double N = x.length;
		double xTot=0;
		double yTot=0;
		double sqrXTot =0;
		double sqrYTot =0;
		double xYTot=0;
		for (int i=0; i<N; i++){
			xTot += x[i];
			yTot += y[i];
			sqrXTot += Math.pow(x[i],2);
			sqrYTot += Math.pow(y[i],2);
			xYTot += (x[i] * y[i]);
		}
		double top = (N * xYTot) - (xTot * yTot);
		double botLeft = Math.sqrt( (N * sqrXTot) - Math.pow(xTot,2) );
		double botRight = Math.sqrt( (N * sqrYTot) - Math.pow(yTot,2) );
		double bot = botLeft*botRight;
		if (bot == 0) {
			String message = "Warning: calc of corr coef error, Num.java";
			System.out.println(message);
			return -2;
		}
		return top/bot;
	}

	/**Calculates Pearson correlation coefficient, r, from two double[]s*/
	public static double correlationCoefficient (double[] x, double[] y){
		double N = x.length;
		double xTot=0;
		double yTot=0;
		double sqrXTot =0;
		double sqrYTot =0;
		double xYTot=0;
		for (int i=0; i<N; i++){
			xTot += x[i];
			yTot += y[i];
			sqrXTot += Math.pow(x[i],2);
			sqrYTot += Math.pow(y[i],2);
			xYTot += (x[i] * y[i]);
		}
		double top = (N * xYTot) - (xTot * yTot);
		double botLeft = Math.sqrt( (N * sqrXTot) - Math.pow(xTot,2) );
		double botRight = Math.sqrt( (N * sqrYTot) - Math.pow(yTot,2) );
		return top/(botLeft*botRight);
	}
}
