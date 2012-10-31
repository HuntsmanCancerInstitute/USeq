package trans.qc;

import java.text.NumberFormat;

/**Info related to a particular statistic used to measure a cel file.*/
public class StatFlag {
	
	//fields, not all necessarily set
	private String name;	//text of the statistic
	private double minimum;	//minimum accepted value
	private double maximum;	//maximum accepted value
	private double median;	//median value of the all of the statistics for a set of arrays
	private double stndDev;
	private static final NumberFormat rnd = NumberFormat.getNumberInstance();
	private boolean checkMin = true;
	private boolean checkMax = true;
	
	//constructors
	public StatFlag(String name, double minimum, double maximum){
		this.name = name;
		this.minimum = minimum;
		this.maximum = maximum;
		//make output formatter for pretty printing
		rnd.setMaximumFractionDigits(3);
	}

	//methods
	
	/**Checks min  max to see if value is contained withing, adds a 'Warning' line to the StringBuffer if so.
	 * Can disable checking minimum by setting boolean.*/
	public void rangeCheck (double value, StringBuffer sb){
		if ((checkMax && value > maximum) || (checkMin && value < minimum)) sb.append(name+", ");
	}
	
	/**@return text median stndDev min max adjustor*/
	public String toString(){
		StringBuffer sb = new StringBuffer(name);
		sb.append("\t");
		sb.append(rnd.format(median));
		sb.append("\t");
		sb.append(rnd.format(stndDev));
		sb.append("\t");
		if (checkMin) sb.append(rnd.format(minimum));
		else sb.append("-");
		sb.append("\t");
		if (checkMax)sb.append(rnd.format(maximum));
		return sb.toString();
	}
	
	//getters setters

	public double getMaximum() {
		return maximum;
	}

	public void setMaximum(double maximum) {
		this.maximum = maximum;
	}

	public double getMedian() {
		return median;
	}

	public void setMedian(double median) {
		this.median = median;
	}

	public double getMinimum() {
		return minimum;
	}

	public void setMinimum(double minimum) {
		this.minimum = minimum;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public void setCheckMin(boolean checkMin) {
		this.checkMin = checkMin;
	}	
	
	public void setCheckMax(boolean checkMax) {
		this.checkMax = checkMax;
	}

	public double getStndDev() {
		return stndDev;
	}

	public void setStndDev(double stndDev) {
		this.stndDev = stndDev;
	}
}
