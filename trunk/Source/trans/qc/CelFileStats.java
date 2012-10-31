package trans.qc;

import util.gen.*;
import java.io.*;
import java.text.NumberFormat;
import java.util.*;
import trans.graphics.*;

/**Container for cel file and it's associated statistics.*/
public class CelFileStats {
	//fields
	private double mean = 0;
	private double rawMedian = 0;
	private double medianScalar = 1;
	private double coefficientOfVariation = 0; //coefficient of the variation
	private double[] quartiles; //25th, 50th, 75th
	private String groupName;
	private File serializedCelFile;
	private float[] intensities;
	float[][] virtualCel;
	private static NumberFormat rnd;
	private boolean ok = true;
	private StringBuffer notes = new StringBuffer();
	private double medianNoSynthControls = 0;
	private double medianDimControls = 0;
	private double medianBrightControls = 0;
	private int numNoSynthOutliers = 0;
	private int numDimOutliers = 0;
	private int numBrightOutliers = 0;
	
	//constructors
	/**Loads virtualCel float[][] and calculates a rawMedian.*/
	public CelFileStats(File serializedCelFile, String groupName){
		//make output formatter for pretty printing
		rnd = NumberFormat.getNumberInstance();
		rnd.setMaximumFractionDigits(3);
		//load intensities
		this.groupName = groupName;
		this.serializedCelFile = serializedCelFile;
		virtualCel = (float[][])IO.fetchObject(serializedCelFile);
		intensities = Num.collapseFloatArray(virtualCel);
		//calculate raw median
		Arrays.sort(intensities);
		rawMedian = Num.median(intensities);
	}
	
	/**Reloads intensites and calculates mean, coef var, percentiles.*/
	public void calculateStats(){
		intensities = Num.collapseFloatArray(virtualCel);
		mean = Num.mean(intensities);
		coefficientOfVariation = Num.standardDeviation(intensities, mean)/mean;
		
		//Copy array, sort and pitch
		float[] copy = new float[intensities.length];
		System.arraycopy(intensities,0,copy,0,intensities.length);
		Arrays.sort(copy);
		quartiles = Num.quartiles(copy);
		copy = null;
	}

	/**Calculates the median value for each of the different classes of controls. 
	 * NoSynt, Dim, Bright, using the virtualCel and the supplied xy coordinates.
	 * Also counts the number of feature outliers based on the number of standard deviations*/
	public void calculateControlStats(int[][][] controls, double noSynthMultiplier, 
			double dimMultiplierHigh, double dimMultiplierLow, double brightMultiplier){
		//noSynth
		float[] floats = Num.fetchMatrixValues(virtualCel, controls[0]);
		//watch out for no values
		if (floats!=null && floats.length!=0){
			Arrays.sort(floats);
			medianNoSynthControls = Num.median(floats);
			double threshold = medianNoSynthControls * noSynthMultiplier;	
			numNoSynthOutliers = Num.countOutliers(floats, threshold, true);
		}
		//dim
		floats = Num.fetchMatrixValues(virtualCel, controls[1]);
		if (floats!=null && floats.length!=0){
			Arrays.sort(floats);
			medianDimControls = Num.median(floats);
			double threshold = medianDimControls * dimMultiplierHigh;	
			int high = Num.countOutliers(floats, threshold, true);
			threshold = medianDimControls * dimMultiplierLow;
			int low = Num.countOutliers(floats, threshold, false);
			//set highest num outliers
			if (high > low ) numDimOutliers = high;
			else numDimOutliers = low;
		}
		//bright
		if (floats!=null && floats.length!=0){
			floats = Num.fetchMatrixValues(virtualCel, controls[2]);
			Arrays.sort(floats);
			medianBrightControls = Num.median(floats);		
			double threshold = medianBrightControls * brightMultiplier;	
			numBrightOutliers = Num.countOutliers(floats, threshold, false);
		}
	}
	
	//methods
	/**Nulls the virtualCel and collapsed intensities array.*/
	public void nullIntensityArrays(){
		intensities = null;
		virtualCel = null;
	}
	
	/**Loads the serialized virtualCel (float[][]) from disk and then 
	 * collapses it to generate the intensities (float[]).*/
	public void loadFloatArrays(boolean makeCollapsedFloatArray){
		virtualCel = (float[][])IO.fetchObject(serializedCelFile);
		if (makeCollapsedFloatArray) intensities = Num.collapseFloatArray(virtualCel);
	}
	
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append(serializedCelFile.getName());
		sb.append("\t");
		sb.append(rnd.format(medianScalar));
		sb.append("\t");
		sb.append(rnd.format(mean));
		sb.append("\t");
		sb.append(rnd.format(coefficientOfVariation));
		sb.append("\t");
		sb.append(rnd.format(quartiles[0]));
		sb.append("\t");
		sb.append(rnd.format(quartiles[1]));
		sb.append("\t");
		sb.append(rnd.format(quartiles[2]));
		sb.append("\t");
		sb.append(numNoSynthOutliers);
		sb.append("\t");
		sb.append(numDimOutliers);
		sb.append("\t");
		sb.append(numBrightOutliers);
		return sb.toString();
	}

	public double getCoefficientOfVariation() {
		return coefficientOfVariation;
	}
	public float[] getIntensities() {
		return intensities;
	}
	public double getMean() {
		return mean;
	}
	public double[] getQuartiles() {
		return quartiles;
	}
	public String getGroupName() {
		return groupName;
	}
	public File getSerializedCelFile() {
		return serializedCelFile;
	}
	public boolean isOk() {
		return ok;
	}
	public void setOk(boolean ok) {
		this.ok = ok;
	}
	public StringBuffer getNotes() {
		return notes;
	}
	public void setNotes(StringBuffer notes) {
		this.notes = notes;
	}
	public void appendNotes(String comment){
		notes.append(comment);
	}
	public float[][] getVirtualCel() {
		return virtualCel;
	}
	public double getMedianBrightControls() {
		return medianBrightControls;
	}
	public double getMedianDimControls() {
		return medianDimControls;
	}
	public double getMedianNoSynthControls() {
		return medianNoSynthControls;
	}
	public int getNumBrightOutliers() {
		return numBrightOutliers;
	}
	public void setNumBrightOutliers(int numBrightOutliers) {
		this.numBrightOutliers = numBrightOutliers;
	}
	public int getNumDimOutliers() {
		return numDimOutliers;
	}
	public void setNumDimOutliers(int numDimOutliers) {
		this.numDimOutliers = numDimOutliers;
	}
	public int getNumNoSynthOutliers() {
		return numNoSynthOutliers;
	}
	public void setNumNoSynthOutliers(int numNoSynthOutliers) {
		this.numNoSynthOutliers = numNoSynthOutliers;
	}

	public double getMedianScalar() {
		return medianScalar;
	}

	public void setMedianScalar(double medianScalar) {
		this.medianScalar = medianScalar;
	}

	public double getRawMedian() {
		return rawMedian;
	}

	public void setRawMedian(double rawMedian) {
		this.rawMedian = rawMedian;
	}
}
