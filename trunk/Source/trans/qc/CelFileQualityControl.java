package trans.qc;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.regex.*;
import java.util.*;
import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.data.category.*;
import util.bio.cluster.*;
import trans.graphics.VirtualCelPanel;
import util.gen.*;

/**
 * Provides automated qc checking for cel files, flags those of poor quality 
 * based on a variety of statistics. 
 * 
 * For every cel file, the median intensity is calculated, then a median determined for all the cel files' medians,
 * then a scaler is calculated for each cel file required to get it to the overall median.  This is the cel file specific 
 * median scalar and indicates the degree of brighness relative to the other chips.  After this calculation the intensities
 * on all the chips are median scaled to the targetMedian defined in the application, ie 100.  Lastly, all the stats are 
 * calculated on these scaled values.
 * 
 * For the dim bright no synth calculations, these are made using the particular chip intensities, not any global values.
 * 
 * @author David Nix
 */
public class CelFileQualityControl {
	
	//fields
	private File[][] celFileGroups;
	private String[] groupNames;
	private CelFileStats[][] celFileStats;
	private StatFlag[] statFlags;
	private HashMap statFlagsHash;
	private NumberFormat rnd;
	private double minWithinGroupCorCoeff = 0.75;		//R, not R squared
	private int maxClusterSizeToFlag = 1;
	private int[][][] controlCoordinates; 				//noSynth, dim, bright, pm x,y coordinates
	private File controlCoordinateFile;
	private boolean writeOutlierPNGs = false;
	private boolean writeAllPNGs = false;
	private boolean displayCharts = false;
	private double targetMedian = 100;
	private boolean switchDimBright = false;
	private File saveDirectory = null;
	private StringBuffer output = new StringBuffer();
	
	//hard coded range values
	private double medianScalerLow = 0.5;
	private double medianScalerHigh = 1.75;
	private double meanLow = 120;
	private double meanHigh = 190;
	private double coeffVarHigh = 5.5;
	private double p25thLow = 70;
	private double p25thHigh = 110;
	private double p75thLow = 110;
	private double p75thHigh = 130;
	
	//for determining what is an outlier, multipliers for median control value
	private double noSynthMultiplierHigh = 3;
	private double dimMultiplierHigh = 3;
	private double dimMultiplierLow = 0.333;
	private double brightMultiplierLow = 0.333;
	
	//percent of total number of pixes that represents too many outliers
	private double noSynthHigh = 0.02041858;
	private double dimHigh = 0.39246468;
	private double brightHigh = 0.07861635;
	
	//primary methods
	/**Main app controler.*/
	public void testCelFiles(){
		System.out.print("\nLaunched ");
		//make output formatter for pretty printing
		rnd = NumberFormat.getNumberInstance();
		rnd.setMaximumFractionDigits(3);
		
		//load control xy coordinates, derived from CoordinateExtractor1lq.java
		output.append("\nReading coordinates...\n");
		System.out.print(".");
		controlCoordinates = (int[][][])IO.fetchObject(controlCoordinateFile);
		int[] numCoor = CoordinateExtractor1lq.countNumberIntensities(controlCoordinates);
		
		// for some damn reason, the dim and bright coordinates switch occasionally
		if (switchDimBright){
			output.append("\tSwitching dim and bright coordinates\n");
			int[][] dim = controlCoordinates[1];
			int[][] bright = controlCoordinates[2];
			controlCoordinates[1] = bright;
			controlCoordinates[2] = dim;
		}
		
		output.append("\t# No Synth (0)\t"+numCoor[0]);
		output.append("\n\t# Dim (-3,-6)\t"+numCoor[1]);
		output.append("\n\t# Bright (3,6)\t"+numCoor[2]);
		output.append("\n\t# PM (-111,111)\t"+numCoor[3]+"\n");
		controlCoordinates[3] = null;
		
		//create CelFileStats[][] from celFileGroups
		output.append("\nLoading cel files...\n");
		System.out.print(".");
		createCelFileStatsArray();
		
		//check median No Synth, Dim, Bright; Dim and Bright are often inverted!
		output.append("\nChecking control intensities...\n");
		System.out.print(".");
		checkForFlippedControls();
		
		//create array of StatFlags to test the cel files
		instantiateStatFlags();
		
		//calculate median stat flag values and set range
		output.append("\nCalculating global median values...\n");
		System.out.print(".");
		calculateMedianStatFlagValues();
		
		//check singles
		output.append("\nChecking single cel file stats...\n");
		System.out.print(".");
		checkSingleCelFileStats();
		
		//hierarchical cluster groups of cel files using median scaled PM intensities
		int[][] pmCoordinates = ((int[][][])IO.fetchObject(controlCoordinateFile))[3];
		System.out.print(".");
		output.append("\nHierarchically clustering cel file groups based on "+pmCoordinates.length+" PM intensities...\n");
		cluster(pmCoordinates);
		pmCoordinates = null;
		
		//print results on single cel file analysis
		output.append("\nGlobal median (intensity and count) values...\n");
		printStatFlags();
		output.append("\nSingle cel file statistics...\n");
		printCelFileStats();
		
		//save output to text file
		File results = new File(saveDirectory, "qcResults.txt");
		IO.writeString(output.toString(), results);
		System.out.println(output);
		
		//make charts for diff values
		makeCharts();
		
		//write outlier PNGs?
		if(writeOutlierPNGs || writeAllPNGs) {
			output.append("\nWriting PNGs...\n");
			System.out.print(".");
			writePNGs(writeOutlierPNGs);
		}
	}
	
	/**Clusters cel files based on pm intensities by group.*/
	public void cluster(int[][] pmCoordinates){
		//for each group
		for (int i=0; i<celFileGroups.length; i++ ){
			File[] vcs = celFileGroups[i];
			output.append("\n"+celFileStats[i][0].getGroupName()+" group...\n");
			File pngClusterPlot =  new File (saveDirectory,celFileStats[i][0].getGroupName()+"ClusterPlot.png");
			HierarchicalClustering vhc = new HierarchicalClustering (minWithinGroupCorCoeff, maxClusterSizeToFlag, output, pngClusterPlot, displayCharts);
			vhc.setTitle(celFileStats[i][0].getGroupName()+"Group (RSqr x 100)");
			File[] badFiles = vhc.clusterVirtualCelFiles(vcs, pmCoordinates);
			if (badFiles != null){
				output.append("\n***** Correlation coefficient warning! One or more clusters has an R value = 1 or < "+minWithinGroupCorCoeff +". *****\n\n");
				//flag appropriate celFileStats
				for (int x=0; x<badFiles.length; x++){
					File badFile = badFiles[x];
					for (int z=0; z<celFileStats[i].length; z++){
						if (badFile.equals(celFileStats[i][z].getSerializedCelFile())){
							celFileStats[i][z].appendNotes("HC");
							celFileStats[i][z].setOk(false);
						}
					}
				}
			}
			
		}
	}
	
	/**Writes charts to the save directory and if so indicated displays them.*/
	public void makeCharts(){
		//kill the need to display so it can run on the cluster without a DISPLAY
		if (displayCharts == false )System.setProperty("java.awt.headless","true");
		//make datasets
		DefaultCategoryDataset medianScalarDS = new DefaultCategoryDataset();
		DefaultCategoryDataset meanDS = new DefaultCategoryDataset();
		DefaultCategoryDataset coeffVarDS = new DefaultCategoryDataset();
		DefaultCategoryDataset p25thDS = new DefaultCategoryDataset();
		DefaultCategoryDataset p75thDS = new DefaultCategoryDataset();
		DefaultCategoryDataset noSynthDS = new DefaultCategoryDataset();
		DefaultCategoryDataset dimDS = new DefaultCategoryDataset();
		DefaultCategoryDataset brightDS = new DefaultCategoryDataset();
		
		//make eight graphs each
		int numFiles = 0;
		
		//make trunkated file names
		String[][] fileNames = IO.getTruncatedNames(celFileGroups);
		
		for (int i=0; i< celFileStats.length; i++){
			
			//for each file load different datasets
			for (int j=0; j< celFileStats[i].length; j++){
				medianScalarDS.setValue(celFileStats[i][j].getMedianScalar(), "Median Scalar", fileNames[i][j]);
				meanDS.setValue(celFileStats[i][j].getMean(), "Mean", fileNames[i][j]);
				coeffVarDS.setValue(celFileStats[i][j].getCoefficientOfVariation(), "Coeff Var", fileNames[i][j]);
				double[] quartiles = celFileStats[i][j].getQuartiles();
				p25thDS.setValue(quartiles[0], "25th Percentile", fileNames[i][j]);
				p75thDS.setValue(quartiles[2], "75th Percentile", fileNames[i][j]);
				noSynthDS.setValue(celFileStats[i][j].getNumNoSynthOutliers(), "Number No Synth Outliers", fileNames[i][j]);
				dimDS.setValue(celFileStats[i][j].getNumDimOutliers(), "Number Dim Outliers", fileNames[i][j]);
				brightDS.setValue(celFileStats[i][j].getNumBrightOutliers(), "Number Bright Outliers", fileNames[i][j]);
				numFiles++;
			}
		}
		//make charts
		JFreeChart median = ChartFactory.createBarChart("Median Scalar",null, null, medianScalarDS, PlotOrientation.HORIZONTAL,false, true, false);
		JFreeChart mean = ChartFactory.createBarChart("Mean",null, null, meanDS, PlotOrientation.HORIZONTAL,false, true, false);
		JFreeChart coeffVar = ChartFactory.createBarChart("Coefficient of Variation",null, null, coeffVarDS, PlotOrientation.HORIZONTAL,false, true, false);
		JFreeChart p25th = ChartFactory.createBarChart("25th Percentile",null, null, p25thDS, PlotOrientation.HORIZONTAL,false, true, false);
		JFreeChart p75th = ChartFactory.createBarChart("75th Percentile",null, null, p75thDS, PlotOrientation.HORIZONTAL,false, true, false);
		JFreeChart noSynth = ChartFactory.createBarChart("No Synthesis Outliers",null, null, noSynthDS, PlotOrientation.HORIZONTAL,false, true, false);
		JFreeChart dim = ChartFactory.createBarChart("Dim Outliers",null, null, dimDS, PlotOrientation.HORIZONTAL,false, true, false);
		JFreeChart bright = ChartFactory.createBarChart("Bright Outliers",null, null, brightDS, PlotOrientation.HORIZONTAL,false, true, false);
		//set boundary lines
		
		//Display?
		if (displayCharts){
			JFreeChart[] charts = {median, mean, coeffVar, p25th, p75th, noSynth, dim, bright};
			new ChartFrame(charts, "Individual Cel File Statistics", 2,4); 
		}
		
		//save
		try {
			//set height at 15 pix per file
			int height = numFiles * 17;
			if (height < 200 ) height = 200;
			ChartUtilities.saveChartAsPNG(new File(saveDirectory, "MedianScalar.png"), median, 800, height);
			ChartUtilities.saveChartAsPNG(new File(saveDirectory, "Mean.png"), mean, 800, height);
			ChartUtilities.saveChartAsPNG(new File(saveDirectory, "CoeffVar.png"), coeffVar, 800, height);
			ChartUtilities.saveChartAsPNG(new File(saveDirectory, "25th.png"), p25th, 800, height);
			ChartUtilities.saveChartAsPNG(new File(saveDirectory, "75th.png"), p75th, 800, height);
			ChartUtilities.saveChartAsPNG(new File(saveDirectory, "NoSynth.png"), noSynth, 800, height);
			ChartUtilities.saveChartAsPNG(new File(saveDirectory, "Dim.png"), dim, 800, height);
			ChartUtilities.saveChartAsPNG(new File(saveDirectory, "Bright.png"), bright, 800, height);
		} catch (IOException e) {
			System.out.println("Problem occurred creating chart!");
			e.printStackTrace();
		}
		
		
	}
	
	
	
	/**Runs through groups of cel files, calculates the median no Synth, Dim, and Bright intensities,
	 * checks to see if they are inverted, for some reason they sometimes are!  Note these medians 
	 * are calculated on the raw, pre scaled chip intensities.*/
	public void checkForFlippedControls(){
		//calculate median intensities for each class of controls
		double[][] noSynth = new double[celFileStats.length][];
		double[][] dim = new double[celFileStats.length][];
		double[][] bright = new double[celFileStats.length][];
		for (int i=0; i< celFileStats.length; i++){
			noSynth[i] = new double[celFileStats[i].length];
			dim[i] = new double[celFileStats[i].length];
			bright[i] = new double[celFileStats[i].length];
			
			for (int j=0; j< celFileStats[i].length; j++){
				noSynth[i][j] = celFileStats[i][j].getMedianNoSynthControls();
				dim[i][j] = celFileStats[i][j].getMedianDimControls();
				bright[i][j] = celFileStats[i][j].getMedianBrightControls();
			}
		}
		//collapse and sort
		double[] colNoSynth = Num.collapseDoubleArray(noSynth);
		double[] colDim = Num.collapseDoubleArray(dim);
		double[] colBright = Num.collapseDoubleArray(bright);
		
		Arrays.sort(colNoSynth);
		Arrays.sort(colDim);
		Arrays.sort(colBright);
		//calculate median median
		double medianNoSynth = Num.median(colNoSynth);
		double medianDim = Num.median(colDim);
		double medianBright = Num.median(colBright);
		
		output.append("\t"+rnd.format(medianNoSynth)+"\t Median no synth intensity\n");
		output.append("\t"+rnd.format(medianDim)+"\t Median dim intensity\n");
		output.append("\t"+rnd.format(medianBright)+"\t Median bright intensity\n");
		if (medianNoSynth > medianDim) output.append("***** WARNING! No synth intensities are greater than dim intensities! *****\n");
		if (medianNoSynth > medianBright) output.append("***** WARNING! No synth intensities are greater than bright intensities! *****\n");
		if (medianDim > medianBright) output.append("***** WARNING! Dim intensities are greater than bright intensities! Use the -s flag and restart. *****\n");
	}
	
	/**Runs through groups of cel files, calling statSingleCelFile() on each
	 * and prints the results.*/
	public void printCelFileStats(){
		output.append("\tCelFile\tMedian Scalar\tMean\tCoefVar\t25th\t50th\t75th\t#NoSynthOutlrs\t#DimOutlrs\t#BrightOutlrs\tStatus\n\n");
		for (int i=0; i< celFileStats.length; i++){
			output.append("\tGroup: "+celFileStats[i][0].getGroupName()+"\n");
			for (int j=0; j< celFileStats[i].length; j++){
				output.append("\t"+celFileStats[i][j]);	
				if (celFileStats[i][j].isOk()) {
					output.append("\tPass\n");
				}
				else {
					output.append("\tFail\t"+celFileStats[i][j].getNotes().toString()+"\n");
				}
			}
			output.append("\n");
		}
	}
	
	/**Runs through groups of cel files, calling statSingleCelFile() on each
	 * and prints the results.*/
	public void clusterCelFiles(int[][] pmCoordinates){
		
		for (int i=0; i< celFileStats.length; i++){
			//for each group
			output.append("\tGroup: "+celFileStats[i][0].getGroupName()+"\n");
			//build array of 
			for (int j=0; j< celFileStats[i].length; j++){
				output.append("\t"+celFileStats[i][j]);	
				if (celFileStats[i][j].isOk()) {
					output.append("\tPass\n");
				}
				else {
					output.append("\tFail\t"+celFileStats[i][j].getNotes().toString()+"\n");
				}
			}
			output.append("\n");
		}
	}
	
	/**Runs through groups of cel files, calling statSingleCelFile() on each
	 * and prints the results.*/
	public void checkSingleCelFileStats(){
		for (int i=0; i< celFileStats.length; i++){
			for (int j=0; j< celFileStats[i].length; j++){
				String res = statSingleCelFile(celFileStats[i][j]).toString();
				if (Misc.isEmpty(res)) {
					celFileStats[i][j].setOk(true);
				}
				else {
					celFileStats[i][j].setOk(false);
					celFileStats[i][j].appendNotes(res.toString());
				}
			}
		}
	}
	
	public File[] fetchCelFileStatsFiles(){
		ArrayList al = new ArrayList();
		for (int i=0; i< celFileStats.length; i++){
			for (int j=0; j< celFileStats[i].length; j++){
				al.add(celFileStats[i][j].getSerializedCelFile());
			}
		}
		File[] files = new File[al.size()];
		al.toArray(files);
		return files;
	}
	
	/**Writes PNG files for  CelFileStats.*/
	public void writePNGs(boolean restrictToOutlier){
		//kill the need to display
		System.setProperty("java.awt.headless","true");
		VirtualCelPanel panel; 
		for (int i=0; i< celFileStats.length; i++){
			for (int j=0; j< celFileStats[i].length; j++){
				if (restrictToOutlier == false || celFileStats[i][j].isOk() == false) {
					String name = celFileStats[i][j].getSerializedCelFile().getName();
					celFileStats[i][j].loadFloatArrays(false);
					float[][] intensities = celFileStats[i][j].getVirtualCel();
					intensities = Num.medianNormalize(intensities, 100);
					name = Misc.replaceEnd(name, "cela","png");
					File png = new File (saveDirectory, name);
					panel = new VirtualCelPanel(intensities, 2000, png);
					celFileStats[i][j].nullIntensityArrays();
					intensities = null;
				}
			}
		}
	}
	
	
	/**Prints text, median min max for each stat flag*/
	public void printStatFlags(){
		output.append("\t\tMedian\tStndDev\tMin\tMax\n");
		for (int i=0; i<statFlags.length; i++){
			output.append("\t"+statFlags[i]+"\n");
		}
	}
	
	/**Calculates the median value for each stat flag*/
	public void calculateMedianStatFlagValues(){
		//initialize arrays
		int len = Misc.totalLength(celFileStats);
		double[] scalars = new double[len];
		double[] means = new double[len];
		double[] coefVars = new double[len];
		double[] p25ths = new double[len];
		double[] p50ths = new double[len];
		double[] p75ths = new double[len];
		int[] numNoSynths = new int[len];
		int[] numDims = new int[len];
		int[] numBrights = new int[len];
		int num = celFileStats.length;
		
		//run through all of the CelFileStats
		int c = 0;
		for (int i=0; i<num; i++){
			int num2 = celFileStats[i].length;
			for (int j=0; j<num2; j++){
				scalars[c] = celFileStats[i][j].getMedianScalar();
				means[c] = celFileStats[i][j].getMean();
				coefVars[c] = celFileStats[i][j].getCoefficientOfVariation();
				p25ths[c] = celFileStats[i][j].getQuartiles()[0];
				p50ths[c] = celFileStats[i][j].getQuartiles()[1];
				p75ths[c] = celFileStats[i][j].getQuartiles()[2];
				numNoSynths[c] = celFileStats[i][j].getNumNoSynthOutliers();
				numDims[c] = celFileStats[i][j].getNumDimOutliers();
				numBrights[c] = celFileStats[i][j].getNumBrightOutliers();
				c++;
			}
		}
		if (c!=1){
			//sort arrays
			Arrays.sort(scalars);
			Arrays.sort(means);
			Arrays.sort(coefVars);
			Arrays.sort(p25ths);
			Arrays.sort(p50ths);
			Arrays.sort(p75ths);
			Arrays.sort(numNoSynths);
			Arrays.sort(numDims);
			Arrays.sort(numBrights);
			
			//calculate and set medians and standard deviations
			StatFlag sf = (StatFlag)statFlagsHash.get("Mean");
			sf.setMedian(Num.median(means));
			sf.setStndDev(Num.standardDeviation(means));
			
			sf = (StatFlag)statFlagsHash.get("Median Scalar");
			sf.setMedian(Num.median(scalars));
			sf.setStndDev(Num.standardDeviation(scalars));
			
			sf = (StatFlag)statFlagsHash.get("Coefficient of Variation");
			sf.setMedian(Num.median(coefVars));
			sf.setStndDev(Num.standardDeviation(coefVars));
			
			sf =(StatFlag)statFlagsHash.get("25th Percentile");
			sf.setMedian(Num.median(p25ths));
			sf.setStndDev(Num.standardDeviation(p25ths));
			
			sf= (StatFlag)statFlagsHash.get("75th Percentile");
			sf.setMedian(Num.median(p75ths));
			sf.setStndDev(Num.standardDeviation(p75ths));
			
			sf= (StatFlag)statFlagsHash.get("No Synthesis");
			sf.setMedian(Num.median(numNoSynths));
			sf.setStndDev(Num.standardDeviation(numNoSynths));
			
			sf= (StatFlag)statFlagsHash.get("Dim");
			sf.setMedian(Num.median(numDims));
			sf.setStndDev(Num.standardDeviation(numDims));
			
			sf= (StatFlag)statFlagsHash.get("Bright");
			sf.setMedian(Num.median(numBrights));
			sf.setStndDev(Num.standardDeviation(numBrights));
			
		}
		//only one cel file
		else {
			StatFlag sf = (StatFlag)statFlagsHash.get("Mean");
			sf.setMedian(means[0]);
			
			sf = (StatFlag)statFlagsHash.get("Coefficient of Variation");
			sf.setMedian(coefVars[0]);
			
			sf =(StatFlag)statFlagsHash.get("25th Percentile");
			sf.setMedian(p25ths[0]);
			
			sf= (StatFlag)statFlagsHash.get("Median");
			sf.setMedian(p50ths[0]);
			
			sf= (StatFlag)statFlagsHash.get("75th Percentile");
			sf.setMedian(p75ths[0]);
			
			sf= (StatFlag)statFlagsHash.get("No Synthesis");
			sf.setMedian(numNoSynths[0]);
			
			sf= (StatFlag)statFlagsHash.get("Dim");
			sf.setMedian(numDims[0]);
			
			sf= (StatFlag)statFlagsHash.get("Bright");
			sf.setMedian(numBrights[0]);
			
		}
	}
	
	/**Makes the Stat Flags, add new stats here. Set */
	public void instantiateStatFlags(){
		statFlags = new StatFlag[8];
		statFlagsHash = new HashMap(8);
		
		statFlags[0] = new StatFlag("Median Scalar",medianScalerLow, medianScalerHigh);
		statFlagsHash.put(statFlags[0].getName(),statFlags[0]);
		
		statFlags[1] = new StatFlag("Mean",meanLow, meanHigh); 
		statFlagsHash.put(statFlags[1].getName(),statFlags[1]);
		
		statFlags[2] = new StatFlag("Coefficient of Variation",0,coeffVarHigh);
		statFlags[2].setCheckMin(false);
		statFlagsHash.put(statFlags[2].getName(),statFlags[2]);
		
		statFlags[3] = new StatFlag("25th Percentile",p25thLow, p25thHigh);
		statFlagsHash.put(statFlags[3].getName(),statFlags[3]);
		
		statFlags[4] = new StatFlag("75th Percentile",p75thLow, p75thHigh);
		statFlagsHash.put(statFlags[4].getName(),statFlags[4]);
		
		//what's given are percents so need to calc nums //noSynth, dim, bright x,y coordinates
		double num = controlCoordinates[0].length * noSynthHigh;
		statFlags[5] = new StatFlag("No Synthesis", 0, num);
		statFlags[5].setCheckMin(false);
		statFlagsHash.put(statFlags[5].getName(),statFlags[5]);
		
		num = controlCoordinates[1].length * dimHigh; 
		statFlags[6] = new StatFlag("Dim",0, num);
		statFlags[6].setCheckMin(false);
		statFlagsHash.put(statFlags[6].getName(),statFlags[6]);
		
		num = controlCoordinates[2].length * brightHigh;
		statFlags[7] = new StatFlag("Bright",0, num);
		statFlags[7].setCheckMin(false);
		statFlagsHash.put(statFlags[7].getName(),statFlags[7]);	
	}
	
	/**Creates a cel file stats array, then saves and nulls the float arrays.*/
	public void createCelFileStatsArray(){
		celFileStats = new CelFileStats[celFileGroups.length][];
		int totalNumberCelFiles = Num.countObjects(celFileGroups);
		
		//first load and calculate rawMedians
		double[] rawMedians = new double[totalNumberCelFiles];
		int counter = 0;
		//for each group
		for (int i=0; i<celFileGroups.length; i++){
			celFileStats[i] = new CelFileStats[celFileGroups[i].length];
			String groupName = groupNames[i];
			//check number of cels
			if (celFileGroups[i].length<1 )Misc.printExit("\nError! No 'xxx.cela' files found in group -> "+groupName+". Check your directories.\n");
			//for each file
			for (int j=0; j< celFileGroups[i].length; j++){
				//make celFileStats and calculate median
				celFileStats[i][j] = new CelFileStats(celFileGroups[i][j], groupName);
				rawMedians[counter] = celFileStats[i][j].getRawMedian();
				//check raw median
				if (rawMedians[counter] == 0){
					Misc.printExit("\nError! The raw median for "+celFileGroups[i][j].getName()+ " is zero!  Did you " +
					"set the correct number of rows/ columns when converting this 'xxx.cel' file to a 'xxx.cela' file?\n");
				}
				counter++;
				celFileStats[i][j].nullIntensityArrays();
			}
		}
		
		//calculate overall median
		Arrays.sort(rawMedians);
		double overallMedian = Num.median(rawMedians);
		//median scale each to overallMedian
		
		//for each group
		for (int i=0; i<celFileGroups.length; i++){
			//for each file
			for (int j=0; j< celFileGroups[i].length; j++){
				//calculate and set scalar for average
				double scalar = overallMedian/ celFileStats[i][j].getRawMedian();
				celFileStats[i][j].setMedianScalar(scalar);
				//rescale to target median so statistics look similar
				//set scalar to targetMedian
				scalar = targetMedian/ celFileStats[i][j].getRawMedian();
				//scale float[][]
				celFileStats[i][j].loadFloatArrays(false);
				float[][] f = celFileStats[i][j].getVirtualCel();
				int num = f.length;
				for (int x=0; x<num; x++){
					for (int y=0; y<num; y++){
						Double scaled = new Double(((double)f[x][y]) * scalar);
						f[x][y] = scaled.floatValue();
					}
				}
				//calculate stats
				celFileStats[i][j].calculateStats();
				//calculate stats on the noSynth, dim, bright controls
				celFileStats[i][j].calculateControlStats(controlCoordinates, noSynthMultiplierHigh, 
						dimMultiplierHigh, dimMultiplierLow, brightMultiplierLow);
				celFileStats[i][j].nullIntensityArrays();
			}
		}
	}
	
	
	/**Checks a CelStats object for deviations from the min max settings.
	 * appends nothing if no errors.*/
	public StringBuffer statSingleCelFile(CelFileStats cs){
		StringBuffer sb = new StringBuffer();
		StatFlag statFlag;
		double value;
		
		//check median scalar
		statFlag = (StatFlag)(statFlagsHash.get("Median Scalar"));
		value = cs.getMedianScalar();
		statFlag.rangeCheck(value, sb);
		
		//check mean
		statFlag = (StatFlag)(statFlagsHash.get("Mean"));
		value = cs.getMean();
		statFlag.rangeCheck(value, sb);
		
		//check Coefficient of the variation (stndDev/mean)
		statFlag = (StatFlag)(statFlagsHash.get("Coefficient of Variation"));
		value = cs.getCoefficientOfVariation();
		statFlag.rangeCheck(value, sb);
		
		//check percentiles
		statFlag = (StatFlag)(statFlagsHash.get("25th Percentile"));
		value = cs.getQuartiles()[0];
		statFlag.rangeCheck(value, sb);
		
		statFlag = (StatFlag)(statFlagsHash.get("75th Percentile"));
		value = cs.getQuartiles()[2];
		statFlag.rangeCheck(value, sb);
		
		//check noSynth
		statFlag = (StatFlag)(statFlagsHash.get("No Synthesis"));
		value = cs.getNumNoSynthOutliers();
		statFlag.rangeCheck(value, sb);
		
		//check dim
		statFlag = (StatFlag)(statFlagsHash.get("Dim"));
		value = cs.getNumDimOutliers();
		statFlag.rangeCheck(value, sb);
		
		//check bright
		statFlag = (StatFlag)(statFlagsHash.get("Bright"));
		value = cs.getNumBrightOutliers();
		statFlag.rangeCheck(value, sb);
		
		return sb;
	}
	
	
	
	//main
	public static void main(String[] args) {
		if (args.length!=0) {
			CelFileQualityControl qc = new CelFileQualityControl();
			qc.processArgs(args);
			qc.testCelFiles();
		}
		else {
			printDocs();
		}
		
	}
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		String fileGrouping = null;
		File parameterFile = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': fileGrouping = args[i+1]; i++; break;
					case 'p': writeOutlierPNGs = true;  break;
					case 'a': writeAllPNGs = true;  break;
					case 'd': displayCharts = true;  break;
					case 's': switchDimBright = true; break;
					case 'x': controlCoordinateFile = new File(args[i+1]); i++; break;
					case 'r': saveDirectory = new File(args[i+1]); i++; break;
					case 'l': parameterFile = new File(args[i+1]); i++; break;
					case 'e': Misc.printExit(fetchQCParams());
					case 'c': minWithinGroupCorCoeff= Double.parseDouble(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: {
						Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter request: -"+test+"\n");
				}
			}
		}
		//check for controlCoordinateFile
		if (controlCoordinateFile == null || controlCoordinateFile.canRead() == false){
			Misc.printExit("\nCannot find your serialized 1lq control coordinate file ('xxx.llq.SerCont')?! -> "+controlCoordinateFile);			
		}
		//check for results directory
		if (saveDirectory == null || saveDirectory.isDirectory()==false) Misc.printExit("\nPlease provide a full path directory text for saving the results.\n");
		if (saveDirectory.exists() == false) saveDirectory.mkdir();
		
		//build groupings of files
		ArrayList results = IO.buildFileGroups(fileGrouping,"cela");	
		if (results == null) {
			Misc.printExit("\nSorry could not find any 'xxx.cela' files?!\n");
		}
		groupNames = (String[])results.get(0);
		celFileGroups = (File[][])results.get(1);
		
		//parse params
		if (parameterFile!=null){
			LinkedHashMap hash = IO.loadKeyValueFile(parameterFile);
			if (hash==null)Misc.printExit("\nNo parameters found in file?!\n");
			if (hash.containsKey("medianScalerLow")) medianScalerLow = Double.parseDouble((String)hash.get("medianScalerLow"));
			if (hash.containsKey("medianScalerHigh")) medianScalerHigh = Double.parseDouble((String)hash.get("medianScalerHigh"));
			if (hash.containsKey("meanLow")) meanLow = Double.parseDouble((String)hash.get("meanLow"));
			
			if (hash.containsKey("meanHigh")) meanHigh = Double.parseDouble((String)hash.get("meanHigh"));
			if (hash.containsKey("coeffVarHigh")) coeffVarHigh = Double.parseDouble((String)hash.get("coeffVarHigh"));
			if (hash.containsKey("p25thLow")) p25thLow = Double.parseDouble((String)hash.get("p25thLow"));
			if (hash.containsKey("p25thHigh")) p25thHigh = Double.parseDouble((String)hash.get("p25thHigh"));
			if (hash.containsKey("p75thLow")) p75thLow = Double.parseDouble((String)hash.get("p75thLow"));
			if (hash.containsKey("p75thHigh")) p75thHigh = Double.parseDouble((String)hash.get("p75thHigh"));
			
			if (hash.containsKey("noSynthMultiplierHigh")) noSynthMultiplierHigh = Double.parseDouble((String)hash.get("noSynthMultiplierHigh"));
			if (hash.containsKey("dimMultiplierHigh")) dimMultiplierHigh = Double.parseDouble((String)hash.get("dimMultiplierHigh"));
			if (hash.containsKey("dimMultiplierLow")) dimMultiplierLow = Double.parseDouble((String)hash.get("dimMultiplierLow"));
			if (hash.containsKey("brightMultiplierLow")) brightMultiplierLow = Double.parseDouble((String)hash.get("brightMultiplierLow"));
			
			if (hash.containsKey("noSynthHigh")) noSynthHigh = Double.parseDouble((String)hash.get("noSynthHigh"));
			if (hash.containsKey("dimHigh")) dimHigh = Double.parseDouble((String)hash.get("dimHigh"));
			if (hash.containsKey("brightHigh")) brightHigh = Double.parseDouble((String)hash.get("brightHigh"));
			
			//print out params
			output.append("\nParsed qc params...\n\n");
			output.append(fetchQCParams());
		}
	}
	
	public String fetchQCParams(){
		StringBuffer sb = new StringBuffer();
		sb.append("//hard coded range values, any chip with statistics outside these ranges will be failed");
		sb.append("\nmedianScalerLow = "); sb.append(medianScalerLow);
		sb.append("\nmedianScalerHigh = "); sb.append(medianScalerHigh);
		sb.append("\nmeanLow = "); sb.append(meanLow);
		sb.append("\nmeanHigh = "); sb.append(meanHigh);
		sb.append("\ncoeffVarHighHigh = "); sb.append(coeffVarHigh);
		sb.append("\np25thLow = "); sb.append(p25thLow);
		sb.append("\np25thHigh = "); sb.append(p25thHigh);
		sb.append("\np75thLow = "); sb.append(p75thLow);
		sb.append("\np75thHigh = "); sb.append(p75thHigh);
		
		sb.append("\n\n//multipliers for control value outlier thresholding");
		sb.append("\n//intensity threshold = multiplier X median control value (noSynth, dim, or bright)");
		sb.append("\nnoSynthMultiplierHigh = "); sb.append(noSynthMultiplierHigh);
		sb.append("\ndimMultiplierHigh = "); sb.append(dimMultiplierHigh);
		sb.append("\ndimMultiplierLow = "); sb.append(dimMultiplierLow);
		sb.append("\nbrightMultiplierLow = "); sb.append(brightMultiplierLow);
		
		sb.append("\n\n//fraction of total control counts at which a flag should be thrown");
		sb.append("\n//number of outlier intensities = fraction X total number of controls (noSynth, dim, or bright)");
		sb.append("\nhcNoSynthHigh = "); sb.append(noSynthHigh);
		sb.append("\nhcDimHigh = "); sb.append(dimHigh);
		sb.append("\nhcBrightHigh = "); sb.append(brightHigh);
		sb.append("\n\n");	
		
		return sb.toString();
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                     Cel File Quality Control: May  2006                          **\n" +
				"**************************************************************************************\n" +
				"Calculates a variety of qc statistics on serialized float[][] version cel files.\n" +
				"Uses hierarchical clustering of Pearson correlation coefficients based on cel files'\n" +
				"PM intensities to group like cel files and flag outliers.\n\n" +
				
				"Required Parameters:\n\n" +
				
				"-f Cel files (serialized float[][], 'xxx.cela', output from the CelFileConverter app)\n" +
				"     to process, no spaces. Use key=value grouping to associate\n" +
				"     files for comparison (e.g. grp1=/data/file1.cela,grp1=/data/file2.cela,\n" +
				"     grp2=/data/file3.cela,grp2=/data/file4.cela,/data/file5.cela). Grouped/ ungrouped\n" +
				"     mixing permitted. Directories can also be specified instead of individual files,\n" +
				"     be sure each file you want examined ends in 'cela' (ie grp5=/data/treatDir/)\n" +
				"-x Full path to the serialized 1lq control coordinate file ('xxx.1lq.SerCont')\n" +
				"     generated by the trans/qc/CoordinateExtractor1lq app.\n" +
				"-r Full path directory text where results should be saved.\n\n"+
				
				"Optional Parameters:\n\n" +
				
				"-l Load parameters file, full path file text, see -e option.\n"+
				"-s Switch Dim and Bright coordinates.\n"+
				"-e Print example parameters file with current defaults.\n"+
				"-c Minimum acceptable within group correlation coefficient (R) (defaults to 0.75).\n"+
				"-p Write outlier image PNG files to disk, defaults to false.\n"+
				"-a Write image PNG files for all files, defaults to no.\n"+
				"-d Display charts, default to no.\n\n"+
				
				"Example: java -Xmx512M -jar pathTo/T2/Apps/CelFileQualityControl -f tr=/data/file1.cel,\n" +
				"     tr=/data/file2.cel,cont=/data/file3.cel,cont=/data/file4.cel,/data/file5.cel\n" +
				"     -s -c 0.9 -p -d -x /maps/dmel.1lq.SerCont -r /data/QCDmelRNAResults/\n" +
				
				"\n" +
		        "**************************************************************************************\n");
	}	
	
}
