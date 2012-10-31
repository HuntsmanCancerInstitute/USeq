package trans.cel;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import trans.tpmap.*;
import util.gen.*;
import util.apps.*;

/**
 * Application to process '.cela' files, wrapper for {@link CelMapper}, 
 * {@link QuantileNormalization}.
 */
public class CelProcessor {
	//fields 
	private float targetMedian = 100;
	private String mapDirectory = null;
	private boolean useMMData = false;
	private boolean printStats = false;
	private boolean controlProbeMedianScale = false;
	private boolean splineScale = false;
	private boolean quantileNormalize = true;
	private File[] celaFiles = null;
	private File[] celaMapFiles = null;
	private boolean breakSaveCelpsByChrom = false;
	private boolean deleteCelps = true;
	private File clusterPlotName = null; 
	
	public CelProcessor (String[] args) {
		processArgs(args);
		
		//get tpmap MapFeature[] file 
		System.out.println("\nLoading tpmap MapFeature[]...");
		File tpmapFeatureFile = new File(mapDirectory, "tpmap.fa");
		MapFeature[] features =  (MapFeature[])IO.fetchObject(tpmapFeatureFile);
		
		//Map cel files to tpmap with CelMapper on each directory, will max((PM-MM),1) if useMMData == true, otherwise mm is tossed
		System.out.println("\nLaunching CelMapper...");
		CelMapper cm = new CelMapper(features, useMMData, printStats);
		Arrays.sort(celaFiles);
		//map em	
		celaMapFiles = cm.mapCelFiles(celaFiles);
		features = null;
		
		//Normalize each directory of celMapped files.
		if (useMMData) System.out.println("\nLaunching normalization and PMMM transformation...");
		else System.out.println("\nLaunching normalization...");
		
		QuantileNormalization qn = new QuantileNormalization();
		qn.setUseMMData(useMMData);
		qn.setTargetMedian(targetMedian);
		File infoFile = new File(mapDirectory,"tpmap.faInfo");
		ArrayList info = (ArrayList)IO.fetchObject(infoFile);
		qn.setTpmapInfo(info);
		
		//Use control indexes for median scaling?
		if (controlProbeMedianScale){
			File file = new File(mapDirectory, "tpmap.controlIndexes");
			if (file.exists() == false ){
				Misc.printExit("\nCannot seem to find your tmpap.controlIndexes in the TPMapFiles folder?  " +
				"Did you include a control fasta file when making a tpmap in MummerMapper? If not then you cannot use the -c option.\n");
			}
			int[][] indexes =  (int[][])IO.fetchObject(file);
			qn.setControlIndexes(indexes);
		}
		
		//load float[replicas][intensities] into qn object
		qn.loadIntensityArrays (celaMapFiles);
		
		//quantileNormalize?
		if (quantileNormalize){
			//median scale, quant norm
			qn.quantileNormalize();
			//save null quantiles
			qn.saveAndNullQuantiles();
		}
		//scale non linearly based on control chromosome sequences?
		else if (splineScale) {
			//median scale, quant norm
			qn.splineNormalize();
			//save null quantiles
			qn.saveAndNullQuantiles();
		}
		//just median scale
		else {
			qn.medianScaleIntensities();
			qn.saveIntensities();
		}
		//print stats?
		if (printStats)qn.printIntensityStats();
		
		//cluster and correlate celaMapFiles?
		if (clusterPlotName != null) {
			System.out.println("\nCorrelating and Hierarchical clustering by cel file raw PM intensities...");
			new Correlate (celaMapFiles, clusterPlotName);
		}
		
		//delete old xxx.celaMap files
		IO.deleteFiles(celaMapFiles);
		
		//break and save
		if (breakSaveCelpsByChrom){
			System.out.println("\nSplitting by chromosome...");
			//load features, slow
			features =  (MapFeature[])IO.fetchObject(tpmapFeatureFile);
			//for each celp file
			File[] celps = qn.getCelpFiles();
			for (int i=0; i<celps.length; i++){
				breakSaveFeatureArray(info, features, celps[i]);
				if (deleteCelps) celps[i].delete();
			}
		}
		System.out.println("\nDone!\n");
	}
	/**Use this method to break and save a celp float[] file into chromosomal components in a 
	 * directory named after the celpFile minus the .celp extension.*/
	public static void breakSaveFeatureArray(ArrayList bpmapInfo, MapFeature[] mapFeatures, File celpFile){
		//load intensities
		float[] intensities = (float[])IO.fetchObject(celpFile);
		//start at 1 since 1st num in the number of data lines
		int size = bpmapInfo.size();
		//make folder to hold results
		File dir = new File(celpFile.getParent(), Misc.removeExtension(celpFile.getName()));		
		dir.mkdir();
		for (int i=1; i< size; i+=4){
			//fetch fields
			String chromName = (String)bpmapInfo.get(i);
			int startIndex = ((Integer)bpmapInfo.get(i+1)).intValue();
			int length = ((Integer)bpmapInfo.get(i+3)).intValue();
			int end = startIndex + length;
			
			//make new IntensityFeature[]
			IntensityFeature[] ifs = new IntensityFeature[length];
			int counter = 0;
			for (int x=startIndex; x<end; x++){
				ifs[counter++] = new IntensityFeature(mapFeatures[x], intensities[x]);
			}
			//save IntensityFeature[]
			if (anyNaNs(ifs)) Misc.printExit("\nNaN in intensity features "+celpFile);			
			IO.saveObject(new File(dir,chromName), ifs);				
		}
	}
	
	
	
	/**Looks for NaN in IntensityFeatures. Returns true if NaNs found, false if clean.*/
	public static boolean anyNaNs(IntensityFeature[] ints){
		for (int i=0; i< ints.length; i++){
			if(Float.isNaN(ints[i].intensity)) return true;
		}
		return false;
	}
	
	/**Looks for NaN in IntensityFeatures. Returns true if NaNs found, false if clean.*/
	public static boolean anyNaNs(float[] ints){
		for (int i=0; i< ints.length; i++){
			if(Float.isNaN(ints[i])) return true;
		}
		return false;
	}
	
	/**Given a float[] celp file and the corresponding MapFeature[] will make a 
	 * minimal */
	public static void breakAndSaveCelp(File celpFile, MapFeature[] features){
		float[] ints = (float[])IO.fetchObject(celpFile);
		int num = ints.length;
		IntensityFeature[] f = new IntensityFeature[num];
		for (int i=0; i< num; i++){
			f[i] = new IntensityFeature(features[i], ints[i]);
		}
	}
	
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		File directory = null;
		String filesString = null;
		
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': directory = new File(args[i+1]); i++; break;
					case 'd': filesString = args[i+1]; i++; break;
					case 'u': useMMData = true; break;
					case 's': printStats = true; break;
					case 'c': controlProbeMedianScale = true; break;
					case 'r': breakSaveCelpsByChrom = true; break;
					case 'q': quantileNormalize = false; break;
					case 'p': splineScale = true; controlProbeMedianScale = true; quantileNormalize = false; break;
					case 'm': targetMedian=Float.parseFloat(args[i+1]); i++; break;
					case 'i': clusterPlotName = new File(args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		
		try {
			//check to see if they entered required params
			if (directory==null || directory.isDirectory() == false){
				System.out.println("\nCannot find your tpmap directory!\n");
				System.exit(0);
			}
			else mapDirectory = directory.getCanonicalPath();
		}
		catch(IOException e){
			e.printStackTrace();
		}
		
		//get .cela files
		//any commas
		if (filesString.indexOf(",") == -1) celaFiles = IO.extractFiles(new File(filesString), ".cela");
		else celaFiles = IO.extractFiles(filesString);
		
	}
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Cel Processor: Sept 2008                              **\n" +
				"**************************************************************************************\n" +
				"CP tpmaps, quantile normalizes, and median scales raw intensity values from float[][]\n" +
				"'xxx.cela' files.  Group files using different directories.\n" +
				"\n" +
				"Use the following options when running CP:\n\n" +
				"-t Full path directory text for the 'TPMapFiles' directory generated by the\n" +
				"      TPMapProcessor.\n" +
				"-d A directory containing xxx.cela files for normalization or a comma delimited list.\n" +
				"-u Use MisMatch data, perform max((PM-MM),1) transformation, default is no.\n" +
				"-q Skip quantile normalization, default is to perform QN.\n"+
				"-c Median scale using control sequences instead of all the mapped intensities.\n"+
				"-p Scale using control genes and cubic splines (chip set dependent).\n"+
				"-m The target median value for scaling, defaults to 100. \n" +
				"-r Break and save each normalized cel file as chromosome specific IntensityFeature[]s\n" +
				"      for downstream multi-chip merging.\n"+
				"-s Calculate and print statistics on each .cel file, default is no.\n" +
				"-i Correlate and hierarchical cluster raw cel file PM intensities. Provide a full path\n" +
				"      text to save the chart.\n"+
				"\n" +
				"Example: java -Xmx1500M -jar pathTo/T2/Apps/CelProcessor -t /affy/TPMapFiles/ -d\n" +
				"      /affy/tCels/,/affy/cCels/ -u -s -m 50 -i /affy/QC/hcPlot.png\n" +
				"\n" +
		"**************************************************************************************\n");
	}	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CelProcessor(args);
	}
	
}
