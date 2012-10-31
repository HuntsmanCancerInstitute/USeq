package trans.cel;
import java.util.*;
import java.io.*;

import trans.tpmap.MummerMapper;
import util.gen.*;

/**
 *Performs quantile normalization on arrays of double values.  Arrays must be of
 *equal length.
 *
 *Note: median scaling can introduce an annoying periodicity into the intensity measurements that is apparent in histograms.
 *This is not a bug.
 *
 *An implementation of:
 *	Bolstad BM, Irizarry RA, Astrand M, Speed TP.
 *	A comparison of normalization methods for high density oligonucleotide array data based on variance and bias.
 *	Bioinformatics. 2003 Jan 22;19(2):185-93.
 *
 */

public class QuantileNormalization {
	
	//fields
	private Quantile[][] quantiles;
	private File[] celaMapFiles; 		//full path text for files being normalized
	private File[] celpFiles;
	private float[][] intensities;		//raw then normalized
	private int[][] controlIndexes;
	private float targetMedian =250; 
	private boolean useMMData = false;
	private ComparatorQuantile quantileComparator= null;
	private boolean scaleByControlStats = false;	
	private ControlStats[] controlStats;
	private  ArrayList tpmapInfo;
	private String[] controlChromosomeNames = {"chrGapdh", "chrActb"};
	private double[] modelValuesForScaling = {30, 900, 1200};

	public QuantileNormalization(){}
	
	public QuantileNormalization(File[] files){
		System.out.println("\nLoading intensity arrays...");
		this.celaMapFiles = files;
		intensities = new float[files.length][];
		for (int i=0; i< files.length; i++){
			intensities[i] = Num.loadFloats(files[i],0);
		}
		printIntensityStats();
		System.out.println("\nQuantile normalizing...");
		quantileNormalize();
		intensities = extractNormalizedValues();
		printIntensityStats();
		System.out.println("\nSaving...");
		for (int i=0; i< files.length; i++){
			File norm = new File (files[i].getParentFile(), Misc.removeExtension(files[i].getName()) +"_Norm.xls");
			Num.writeToFile(intensities[i], norm);
			File cela = new File (files[i].getParentFile(), Misc.removeExtension(files[i].getName()) +"_Norm.cela");
			IO.saveObject(cela, intensities[i]);
		}
		
	}
	
	public static void main(String[] args){
		File[] files = IO.extractFiles(new File (args[0]));
		new QuantileNormalization(files);
	}
	
	//methods
	/**Performs a median scaling then quantile normalization of the loaded float[replicas][intensities].*/
	public void quantileNormalize(){
		//convert arrays to Quantile arrays
		System.out.println("\tMaking Quantiles...");
		quantiles = makeQuantiles(intensities);
		intensities = null;
		
		//sort Quantile arrays by intensity value
		quantileComparator = new ComparatorQuantile(true);
		System.out.println("\tSorting Quantiles by intensity...");
		sortQuantileArrays(quantiles);
		
		//no control indexes therefore median scale by all mapped intensities  
		if (controlIndexes == null){
			System.out.println("\tScaling intensities to a median of "+targetMedian+"...");
			scaleQuantileArraysByMedian(quantiles, targetMedian);
		}
 
		//normalize the scores
		System.out.println("\tAveraging Quantiles...");
		normalizeQuantileArrayScores(quantiles); 
		
		//sort Quantile Arrays by original position
		System.out.println("\tResorting Quantiles by original position...");
		quantileComparator.setSortByValue(false);
		sortQuantileArrays(quantiles);

		//median scale by control intensities?
		if (controlIndexes != null){
			System.out.println("\tScaling intensities based on the median of the control oligos. CompositCNV-> "+targetMedian+"....");
			scaleQuantileArraysByControlIntensities();
		}
	}
	
	/**Uses a cubic spline and user defined regions to non-linearly scale the data.*/
	public void splineNormalize(){
		//convert arrays to Quantile arrays
		System.out.println("\tMaking Quantiles...");
		quantiles = makeQuantiles(intensities);

		//median scale by control intensities?
		if (controlIndexes != null){
			System.out.println("\tScaling intensities based on the median of the control oligos. CompositCNV-> "+targetMedian+"....");
			scaleQuantileArraysByControlIntensities();
		}
		
		//calculate median of control intensities from POSITION SORTED Quantile[]s and scale using Splines
		makeControlStats();
		calculateChromControlSpecificMedian(controlChromosomeNames);
		prependTargetMedianInControlStats();
		//scale based on controlStats		
		SplineScalar model = new SplineScalar(modelValuesForScaling);
		System.out.println("Control Stats\t"+Misc.stringArrayToString(controlStats[0].getNames(),"\t"));
		for (int i=0; i<controlStats.length; i++){
			System.out.println(celaMapFiles[i].getName()+"\t"+controlStats[i]);
			//scale based on model using splines and linear regression
			SplineScalar raw = new SplineScalar(controlStats[i].getStats());
			for (int x=0; x< quantiles[i].length; x++){
				double scalar = model.scalar(quantiles[i][x].value, raw);
				quantiles[i][x].value = new Double (quantiles[i][x].value * scalar).floatValue();
			}
		}
		
		//check?
		if (false){
			int columns = quantiles.length;
			double actualMedian;
			//for each column
			for (int i=0; i<columns; i++){
				//fetch control intensities
				float[] ctrlInt = fetchControlIntensities(quantiles[i]);
				//too few control values?
				if (ctrlInt.length<25) Misc.printExit("\nError: too few control intensities! (<25) Did you " +
						"include a 'chromosome' fasta file named and headed '"+
						MummerMapper.controlChromosomeName+" loaded with control genomic sequences " +
				"when making your tpmap?\n");
				Arrays.sort(ctrlInt);
				actualMedian = Num.median(ctrlInt);
				System.out.println("\t\tFinal Control Median "+actualMedian);
			}
			
			calculateChromControlSpecificMedian(controlChromosomeNames);
			System.out.println("Control Stats\t"+Misc.stringArrayToString(controlStats[0].getNames(),"\t"));
			for (int i=0; i<controlStats.length; i++){
				System.out.println(celaMapFiles[i].getName()+"\t"+controlStats[i]);
			}
		}
		
		
	}
	
	public void prependTargetMedianInControlStats(){
		for (int i=0; i<controlStats.length; i++){
			double[] stats = Num.prependDouble(controlStats[i].getStats(), targetMedian);
			controlStats[i].setStats(stats);
			String[] names = Misc.prependString(controlStats[i].getNames(), "targetMedian");
			controlStats[i].setNames(names);
		}
	}
	
	/**Given chromosome names (e.g. chrGapdh, chrActb) will calculate the median intensities.
	 * Sets the values and names in the ControlStats[]. Don't forget to call makeControlStats().*/
	public void calculateChromControlSpecificMedian(String[] chrNameControls){
		//for each Quantile[]
		for (int i=0; i< controlStats.length; i++){
			ArrayList medians = new ArrayList();
			//for each chrNameControls
			for (int j=0; j<chrNameControls.length; j++){
				double median = calculateMedianOfChromosome(quantiles[i], chrNameControls[j]);
				if (median == -1) Misc.printExit("\nError calculating control median.\n");
				medians.add(new Double(median));
			}
			controlStats[i].setNames(chrNameControls);
			controlStats[i].setStats(Num.arrayListOfDoubleToArray(medians));
		}
	}
	
	/**Instantiates a ControlStats object for each Quantile[] object.*/
	public void makeControlStats(){
		controlStats = new ControlStats[quantiles.length];
		for (int i=0; i< quantiles.length; i++) controlStats[i] = new ControlStats();
	}
	
	/**Looks in the tpmapInfo for the appropriate chromosome, extracts the intensities from the Quantile[], calculates a median.  
	 * Returns -1 if no chromosome found.*/
	public double calculateMedianOfChromosome(Quantile[] intensities, String chromosomeName){
		//start at 1 since 1st num in the number of data lines
		for (int i=1; i< tpmapInfo.size(); i+=4){
			//fetch fields
			String chromName = (String)tpmapInfo.get(i);
			//correct chromosome?
			if (chromName.equals(chromosomeName) == false) continue;
			int startIndex = ((Integer)tpmapInfo.get(i+1)).intValue();
			int length = ((Integer)tpmapInfo.get(i+3)).intValue();
			int end = startIndex + length;
			
			//make float[]
			float[] ifs = new float[length];
			int counter = 0;
			for (int x=startIndex; x<end; x++){
				ifs[counter++] = intensities[x].value;
			}	
			Arrays.sort(ifs);
			return Num.median(ifs);
		}
		System.out.println("\nError: no chromosome found for "+chromosomeName+". Cannot calculate median, returning -1.\n");
		return -1;
	}
	
	/**Median scale current intensity values.*/
	public void medianScaleIntensities(){
		int numFiles = celaMapFiles.length;
		System.out.println("\tScaling intensities to a median of "+targetMedian+"...");
		for (int i=0; i<numFiles; i++){
			intensities[i] = Num.medianNormalize(intensities[i], targetMedian);
		}
		
	}
	
	/**Prints statistics about the current intensity float[][] arrays*/
	public void printIntensityStats(){
		int numFiles = celaMapFiles.length;
		for (int i=0; i<numFiles; i++){
			System.out.println("\nStatistics for intensity values from -> "+celaMapFiles[i]);
			Num.statFloatArray(intensities[i], false);
		}
	}

	/**Fetch intensity arrays from .celaMap files.*/
	public void loadIntensityArrays (File[] celaMapFiles){
		this.celaMapFiles = celaMapFiles;
		intensities = null;
		try {		
			int num = celaMapFiles.length;
			//initialize directory field		
			System.out.println ("\tFetching intensity arrays from "+num+" files...");	
			intensities = new float[num][];
			for (int i = 0; i < num; i++) {
				intensities[i]= (float[])IO.fetchObject(celaMapFiles[i]);
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		}			
	}		
	
	/**Saves intensities, one float[] for each replica.  */
	public void saveIntensities(){
		int numColumns = intensities.length;
		celpFiles = new File[numColumns];
		for (int i=0; i<numColumns; i++){
			String name = Misc.replaceEnd(celaMapFiles[i].getName(), "celaMap", "celp");
			celpFiles[i] = new File(celaMapFiles[i].getParentFile(),name);
			if (useMMData) IO.saveObject(celpFiles[i], transformPMMM(intensities[i]));
			else IO.saveObject(celpFiles[i], intensities[i]);
		}
	}
	
	/**Performs a max(PM-MM,1) on the float array where MM is just after each PM.*/
	public static float[] transformPMMM(float[] f){
		int num = f.length;
		float[] t = new float[num/2];
		float[] pmF = new float[num/2];
		float[] mmF = new float[num/2];
		int counter = 0;
		for (int i=0; i<num; i++){
			float pm = f[i];
			i++;
			float mm = f[i];
			pmF[counter] = pm;
			mmF[counter] = mm;
			float value = pm - mm;
			if (value < 1) value = 1;
			t[counter++] = value;
		}
		double cc = PearsonCorrelation.correlationCoefficient(pmF, mmF);
		System.out.println("\tCorrelation Coefficient for PM vs MM: "+cc);
		return t;
	}
	
	/**Takes the Quantile[][], fetches the float[] values for each, nulls the Quantile[], saves the float[] to disk.
	 * Transforms the data if mm is present.*/
	public void saveAndNullQuantiles(){
		int numColumns = quantiles.length;
		int numRows = quantiles[0].length;
		intensities = new float[numColumns][];
		celpFiles = new File[numColumns];
		for (int i=0; i<numColumns; i++){
			intensities[i] = new float[numRows];
			for (int j=0; j<numRows; j++){
				intensities[i][j]= quantiles[i][j].value;
			}
			//kill Quantile[]
			quantiles[i]=null;
			//save its normalized values
			String name = Misc.replaceEnd(celaMapFiles[i].getName(), "celaMap", "celp");
			celpFiles[i] = new File(celaMapFiles[i].getParentFile(),name);
			if (useMMData) IO.saveObject(celpFiles[i], transformPMMM(intensities[i]));
			else IO.saveObject(celpFiles[i], intensities[i]);
		}
	}	
	
	/**Makes a Quantile[][] given float[][] where the quantile.position is assigned its
	 * array index position and the quantile.value is assigned the float value.*/
	public static Quantile[][] makeQuantiles(float[][] values){
		int numColumns = values.length;
		int numValues = values[0].length;
		
		Quantile[][] q = new Quantile[numColumns][numValues];
		for (int i=numColumns-1; i>=0; i--){
			for (int j=numValues-1; j>=0; j--){
				q[i][j]= new Quantile(j,values[i][j]);
			}
		}
		return q;		
	}
	
	/**Takes arrays of Quantile, and averages the values for a given array index.  */
	public static void normalizeQuantileArrayScores(Quantile[][] quantiles){
		float columns = quantiles.length;
		int values = quantiles[0].length;
		float total;
		float average;
		//run thru all values/ rows
		for (int i=0; i<values; i++){
			total = 0;
			//run thru each dataset/ column
			for (int j=0; j<columns; j++){
				total +=quantiles[j][i].value;
			}
			average= total/columns;
			//assign average to each Quantile
			for (int j=(int)columns-1; j>=0; j--){			
				quantiles[j][i].value = average;
			}	
		}
	}
	
	/**Takes arrays of Quantile that have been sorted by value, and scales their values to the target median.
	 * This can introduce an annoying periodicity into the data.*/
	public static void scaleQuantileArraysByMedian(Quantile[][] quantiles, float targetMedian){
		int columns = quantiles.length;
		int values = quantiles[0].length;
		double actualMedian;
		double scalar;
		//for each column
		for (int i=0; i<columns; i++){
			actualMedian = fetchMedianQuantileValue(quantiles[i]);
			System.out.println("\t\tRaw Median "+actualMedian);
			if (actualMedian == 0) Misc.printExit("\nError: median of a file is zero! Check your xxx.cela conversions. Use the -s flag to find the file.\n");
			scalar = ((double)targetMedian)/actualMedian;
			//run thru and scale values
			for (int j=values-1; j>=0; j--) {
				quantiles[i][j].value *=scalar;
			} 
			
		}
	}
	
	/**Scales a quantile array that is sorted by original position in the MapFeature[] 
	 * using the median value of the control intensities.
	 * Be sure the qn array is sorted by position!*/
	public void scaleQuantileArraysByControlIntensities(){
		int columns = quantiles.length;
		int values = quantiles[0].length;
		double actualMedian;
		double scalar;
		//for each column
		for (int i=0; i<columns; i++){
			//fetch control intensities
			float[] ctrlInt = fetchControlIntensities(quantiles[i]);
			//too few control values?
			if (ctrlInt.length<25) Misc.printExit("\nError: too few control intensities! (<25) Did you " +
					"include a 'chromosome' fasta file named and headed '"+
					MummerMapper.controlChromosomeName+" loaded with control genomic sequences " +
							"when making your tpmap?\n");
			Arrays.sort(ctrlInt);
			actualMedian = Num.median(ctrlInt);
			

			// actualMedian = calculateMedianOfChromosome(quantiles[i], "chrGapdh");
			System.out.println("\t\tRaw Control Median "+actualMedian);
			if (actualMedian == 0) Misc.printExit("\nError: median of control sequences is zero! " +
					"Check your xxx.cela conversions. Use the -s flag to find the file.\n");
			scalar = ((double)targetMedian)/actualMedian;
			//run thru and scale values
			for (int j=values-1; j>=0; j--) {
				quantiles[i][j].value *=scalar;
			} 
			
		}
	}
	
	/**Uses the control indexes to fetch the control intensities from the Quantile[].
	 * Assumes the Quantile[] is sorted by position, not value!
	 * Takes the median of duplicates.*/
	public float[] fetchControlIntensities(Quantile[] q){
		float[] intensities = new float[controlIndexes.length];
		for (int i=0; i< controlIndexes.length; i++){
			int[] dupIndexes = controlIndexes[i];
			float[] values = new float[dupIndexes.length];
			for (int j=0; j< dupIndexes.length; j++){
				if (useMMData) {
					int index = dupIndexes[j]*2;
					values[j]=q[index].value;
				}
				else values[j]=q[dupIndexes[j]].value;
			}
			Arrays.sort(values);
			intensities[i] = new Double(Num.median(values)).floatValue();
		}
		return intensities;
	}
	
	/**Calculates the median intensity of a sorted Quantile[].*/
	public static double fetchMedianQuantileValue(Quantile[] m) {
		int length = m.length;
		int middle = length/2;  // subscript of middle element
		//Odd number of elements -- return the middle one.
		if (m.length%2 == 1) return m[middle].value;
		// Even number -- return average of middle two
		return ((double)m[middle-1].value + (double)m[middle].value) / 2.0;
	}	
	
	/**Returns the value of a given percentile from a SORTED Quantile array.
	 * Percentile is from 0-1, ie 0.95 and is according to Lentner, 1982.*/
	public static float percentile(Quantile[] sortedQs, double percentile){
		//calculate index
		double index = ((percentile * (double)sortedQs.length)) - 0.5;
		int trunkIndex = (int)index;
		//is it a whole number?
		double rnd = index%1;
		if (rnd<0.00001 || rnd> 0.99999){
			return sortedQs[trunkIndex].value;
		}
		//otherwise average trunk and trunk +1
		else {
			return (sortedQs[trunkIndex].value + sortedQs[trunkIndex+1].value)/2;
		}
	}
	
	
	/**Sorts quantile arrays based on comparator, either position or intensity value.*/
	public void sortQuantileArrays(Quantile[][] q){
		for (int i= q.length-1; i>=0; i--) Arrays.sort(q[i], quantileComparator);
	}
	
	/**Extracts out the normalized values regenerating the original float[][] except
	 * the scores have been quantile normalized.*/
	public float[][] extractNormalizedValues(){
		int numColumns = quantiles.length;
		int numValues = quantiles[0].length;
		
		float[][] q = new float[numColumns][numValues];
		for (int i=numColumns-1; i>=0; i--){
			for (int j=numValues-1; j>=0; j--){
				q[i][j]= quantiles[i][j].value;
			}
		}
		return q;		
	}
	
	/** Fetches an array of sorted normalized Quantiles*/
	public Quantile[][] getQuantiles() {
		return quantiles;
	}
	public float[][] getIntensities() {
		return intensities;
	}
	public void setIntensities(float[][] intensities) {
		this.intensities = intensities;
	}
	public float getTargetMedian() {
		return targetMedian;
	}
	public void setTargetMedian(float targetMedian) {
		this.targetMedian = targetMedian;
	}
	public void setQuantiles(Quantile[][] quantiles) {
		this.quantiles = quantiles;
	}

	public void setControlIndexes(int[][] controlIndexes) {
		this.controlIndexes = controlIndexes;
	}

	public boolean isUseMMData() {
		return useMMData;
	}

	public void setUseMMData(boolean useMMData) {
		this.useMMData = useMMData;
	}

	public ControlStats[] getControlStats() {
		return controlStats;
	}

	public void setTpmapInfo(ArrayList tpmapInfo) {
		this.tpmapInfo = tpmapInfo;
	}

	public void setScaleByControlStats(boolean scaleByControlStats) {
		this.scaleByControlStats = scaleByControlStats;
	}

	public void setControlChromosomeNames(String[] controlChromosomeNames) {
		this.controlChromosomeNames = controlChromosomeNames;
	}

	public void setModelValuesForScaling(double[] modelValuesForScaling) {
		this.modelValuesForScaling = modelValuesForScaling;
	}

	public File[] getCelpFiles() {
		return celpFiles;
	}
}
