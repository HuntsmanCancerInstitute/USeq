package util.gen;
import java.util.*;

/**Info and methods about a particular bin, used by the Histogram class.*/
public class Bin {
	
	private double start;
	private double stop;
	private double middle;
	private int hits =0;
	private ArrayList scores = new ArrayList();
	private Histogram histogram;
	
	public Bin (double start, double stop, Histogram histogram){
		this.start = start;
		this.stop = stop;
		middle = ((stop-start)/2.0) + start;
		this.histogram = histogram;
	}
	
	/**Returns true if the value is >= bin.start and < bin.stop, thus stop is exclusive, start is inclusive.*/
	public boolean contains(double value){
		if (value >= start && value < stop) return true;
		return false;
	}
	
	/**If value is contains() returns true and increments the bin hit counter, otherwise returns false.*/
	public boolean count( double value){
		if (contains(value)) {
			hits ++;
			return true;
		}
		return false;
	}
	
	/**If value is contains() returns true and increments the bin hit counter and total otherwise returns false.*/
	public boolean count( double value, double score){
		if (contains(value)) {
			hits ++;
			scores.add(new Double (score));
			return true;
		}
		return false;
	}
	
	public void addCountToHits(int counts){
		hits += counts;
	}
	
	public String toString(){
		String tab = "\t";
		StringBuffer sb = new StringBuffer(getLabel());
		sb.append(tab);
		sb.append(start);
		sb.append(tab);
		sb.append(middle);
		sb.append(tab);
		sb.append(stop);
		sb.append(tab);
		sb.append(hits);
		if (histogram.isMeanScoresPresent()){
			sb.append(tab);
			sb.append(getMean());
			sb.append(tab);
			sb.append(getMedian());
		}
		return sb.toString();
	}
	
	public String getLabel(){
		return ">="+Num.formatNumber(start,3)+ "<"+Num.formatNumber(stop,3);
	}
	
	public String getIntLabel(){
		return ""+(int)start;
	}
	

	public int getHits() {
		return hits;
	}
	public double getStart() {
		return start;
	}
	public double getStop() {
		return stop;
	}
	public double getMiddle() {
		return middle;
	}
	public Double getMean(){
		if (hits == 0) return null;
		double[] s = Num.arrayListOfDoubleToArray(scores);
		return new Double (Num.mean(s));
	}
	public Double getMedian(){
		if (hits == 0) return null;
		double[] s = Num.arrayListOfDoubleToArray(scores);
		Arrays.sort(s);
		return new Double (Num.median(s));
	}
	
}
