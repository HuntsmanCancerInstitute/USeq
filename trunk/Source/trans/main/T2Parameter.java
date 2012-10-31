package trans.main;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.*;

/**Container with info related to runing multiple T2 programs.*/
public class T2Parameter {
	
	//fields
	private LinkedHashMap map;
	
	//general
	private File java;
	private File appsDirectory;
	private int maxMemory;
	private File resultsDirectory;
	private String resultsName;
	private File logFile;
	private int sizeOfOligo;
	
	//maps and cel files
	private File[] tpmapFiles;
	private File[][] treatmentCelFiles;	//File[chipSet#][replicas]
	private File[][] controlCelFiles;		//File[chipSet#][replicas]
	private File[] normalizedTreatmentReplicas;
	private File[] normalizedControlReplicas;
	
	//CelProcessor
	private double targetMedian;
	private boolean quantileNormalizeAll;
	private boolean medianScaleUsingControlSeqs;
	
	//ScanChromosomesCNV
	private int windowSize;
	private int minimumNumberOfOligos;
	private String genomeVersion;
	private boolean usePseudoMedian;
	private boolean useRandomPermutation;
	private String permutationType;
	private int numberOfRandomPermutations;
	private boolean randomizeData;
	private boolean useSymmetricNull;
	private File symmetricNullApp;
	private File rApp;
	private String strand;
	
	//merge window arrays
	private boolean makeReducedWindows = false;
	
	//SetNumberIntervalMaker
	private String numberOfIntervals = null;
	private int maxGap = 24;
	
	//LoadChipSetIntervalOligoInfo
	private File genomeFastaDirectory = null;
	
	//FindSubBindingRegions
	private int subWindowSize;			
	private int minNumberOligosInSubWin;		
	private int peakPickerWindowSize;		
	private int maxNumberPeaks;	
	
	//IntervalFilter
	private File[] repeatRegionFiles = null;
	
	//Constructor
	public T2Parameter (File paramFile){
		//parse parameter file
		map = IO.loadKeyValueSpaceTabDelimited(paramFile, true);
		if (map == null) Misc.printExit("Error: cannot parse your parameter file. Be sure each line is either a comment or contains two value strings seperated by white space.\n");
		
		//load general info
		loadGeneral();
		
		//load maps and raw cel files
		loadMapsAndRawCelFiles();
		
		//load program params
		loadApps();
		
	}
	
	//Methods
	/**Replaces the File pointers to the converted xxx.cela files preserving the input structure.*/
	public void renameCelFiles(){
		treatmentCelFiles = renameToACels(treatmentCelFiles);
		controlCelFiles = renameToACels(controlCelFiles);
	}
	
	public File[][] renameToACels( File[][] f){
		File dirACels = new File(resultsDirectory, "ACels");
		Pattern p = Pattern.compile("\\.cel$", Pattern.CASE_INSENSITIVE);
		File[][] acels = new File[f.length][];
		for (int i=0; i< f.length; i++){
			File[] renamed = new File[f[i].length];
			for (int j=0; j< f[i].length; j++){
				Matcher mat = p.matcher(f[i][j].getName());
				String name = mat.replaceFirst(".cela");
				File test = new File (dirACels, name);
				if (test.exists() == false) Misc.printExit("\nError: converted xxx.cel file doesn't exist! "+test);
				renamed[j] = test;
			}
			acels[i] = renamed;
		}
		return acels;
	}
	
	public void loadApps(){
		targetMedian = Double.parseDouble( fetchValue ("targetMedian"));
		quantileNormalizeAll = fetchBoolean("quantileNormalizeAll");
		medianScaleUsingControlSeqs = fetchBoolean("medianScaleUsingControlSeqs");
		windowSize = Integer.parseInt( fetchValue ("windowSize"));
		minimumNumberOfOligos = Integer.parseInt( fetchValue ("minimumNumberOfOligos"));
		genomeVersion = fetchValue ("genomeVersion");
		strand = fetchPossibleValue ("strand");
		if (strand == null) strand = ".";
		usePseudoMedian = fetchBoolean("usePseudoMedian");
		useRandomPermutation = fetchBoolean("useRandomPermutation");
		if (useRandomPermutation) {
			numberOfRandomPermutations = Integer.parseInt( fetchValue ("numberOfRandomPermutations"));
			permutationType = fetchValue("permutationType").toLowerCase();
		}
		randomizeData = fetchBoolean("randomizeData");
		useSymmetricNull = fetchBoolean("useSymmetricNull");
		if (useSymmetricNull) {
			symmetricNullApp = fetchFile("symmetricNullApp");
			rApp = fetchFile("rApp");
		}
		makeReducedWindows = fetchBoolean("makeReducedWindows");
		numberOfIntervals = fetchPossibleValue("numberOfIntervals");
		if (numberOfIntervals != null) {
			maxGap = Integer.parseInt( fetchValue ("maxGap"));
			String dir = fetchPossibleValue("genomicFastaSeqDirectory");
			if (dir !=null) genomeFastaDirectory = fetchFile("genomicFastaSeqDirectory");
			if (map.containsKey("subWindowSize")) {
				subWindowSize = Integer.parseInt( fetchValue ("subWindowSize"));
				minNumberOligosInSubWin = Integer.parseInt( fetchValue ("minNumberOligosInSubWin"));		
				peakPickerWindowSize = Integer.parseInt( fetchValue ("peakPickerWindowSize"));		
				maxNumberPeaks = Integer.parseInt( fetchValue ("maxNumberPeaks"));
			}
				
			//IntervalFilter
			String repFiles = fetchPossibleValue("repeatRegionFiles");
			if (repFiles != null) {
				repeatRegionFiles = IO.extractFiles(repFiles);
				if (repeatRegionFiles == null) Misc.printExit("\nError: cannot read/ find one or more of your repeat region files? "+repFiles);
			}
			
		}
	}
	
	public void loadGeneral(){
		java = new File (fetchValue("java"));
		appsDirectory = new File (fetchValue("appsDirectory"));
		maxMemory = Integer.parseInt( fetchValue ("maxMemory"));
		if (appsDirectory.exists() == false || appsDirectory.isDirectory() == false) {
			Misc.printExit("\nError: cannot find or read your TiMAT2 Apps directory -> "+appsDirectory);
		}
		resultsDirectory = new File (fetchValue("resultsDirectory"));
		if (resultsDirectory.exists() == false) {
			if (resultsDirectory.mkdir() == false) Misc.printExit("Error: could not make results directory. Check that the folder in which you want to make the results directory actually exists.\n");
		}
		resultsName = fetchValue("resultsName");
		logFile = new File (resultsDirectory, resultsName+"_log.txt");
		sizeOfOligo = Integer.parseInt( fetchValue ("sizeOfOligo"));
	}
	
	public void loadMapsAndRawCelFiles(){
		//tpmaps
		File mapDirectory = new File (fetchValue("mapDirectory"));
		tpmapFiles = fetchNumberedOneFile("m", mapDirectory);
		//treatments
		File celDirectory = new File (fetchValue("celDirectory"));
		treatmentCelFiles = fetchNumberedFiles("t", celDirectory);
		//controls
		controlCelFiles = fetchNumberedFiles("c", celDirectory);
		
		//check files
		if (tpmapFiles.length != treatmentCelFiles.length || treatmentCelFiles.length != controlCelFiles.length){
			Misc.printExit("\nError: the number chips in the chip set is off:\n\t"+tpmapFiles.length+" map files\n\t"+
					treatmentCelFiles.length+" treatment chip sets\n\t"+controlCelFiles.length+" control chip sets.\n");
		}
		if (checkNumberOfReplicas(treatmentCelFiles) == false){
			Misc.printExit("\nError: the number of treatment replicas differ between the different chip sets.\n");
		}
		if (checkNumberOfReplicas(controlCelFiles) == false){
			Misc.printExit("\nError: the number of control replicas differ between the different chip sets.\n");
		}
	}
	public boolean checkNumberOfReplicas(File[][] files){
		int numberOfReplicas = files[0].length;
		for (int i =1; i< files.length; i++){
			if (files[i].length != numberOfReplicas) return false;
		}
		return true;
	}
	public File[][] fetchNumberedFiles(String prepended, File directory){
		ArrayList files = new ArrayList();
		int counter = 1;
		while (true){
			String key = prepended+counter;			
			if (map.containsKey(key)) {
				String fileName = (String) map.get(key);
				String[] names = fileName.split(",");
				File[] reps = new File[names.length];
				for (int i=0; i< names.length; i++){
					reps[i] = new File (directory,names[i]);
					//does it exist?
					if (reps[i].canRead() == false){
						Misc.printExit("Error: cannot find or read '"+names[i]+"' in the '"+directory +"' directory.\n");
					}
				}
				files.add(reps);
				counter++;
			}
			else break;
		}
		int num = files.size();
		File[][] f = new File[num][];
		for (int i=0; i< num; i++){
			f[i] = (File[]) files.get(i);
		}
		return f;
	}
	
	public File[] fetchNumberedOneFile(String prepended, File directory){
		ArrayList files = new ArrayList();
		int counter = 1;
		while (true){
			String key = prepended+counter;			
			if (map.containsKey(key)) {
				String fileName = (String) map.get(key);
				File file = new File (directory,fileName);
				if (file.canRead() == false){
					Misc.printExit("Error: cannot find or read '"+fileName+"' in the '"+directory +"' directory.\n");
				}
				files.add(file);
				counter++;
			}
			else break;
		}
		int num = files.size();
		File[] f = new File[num];
		for (int i=0; i< num; i++){
			f[i] = (File) files.get(i);
		}
		return f;
	}
	
	
	public String fetchValue (String key) {
		if (key == null || map.containsKey(key) == false) Misc.printExit("\nError: cannot find your "+key+" value in the parameter file.");
		return (String) map.get(key);
	}
	
	public String fetchPossibleValue (String key) {
		if (map.containsKey(key) == false) return null;
		return (String) map.get(key);
	}
	
	public File fetchFile (String key) {
		if (map.containsKey(key) == false) Misc.printExit("\nError: cannot find your file/ app? "+key);
		File f = new File ((String) map.get(key));
		if (f.exists() == false) Misc.printExit("\nError: cannot find your file/ app? "+f);
		return f;
	}
	
	public boolean fetchBoolean (String key) {
		if (map.containsKey(key) == false) return false;
		String x = ((String) map.get(key)).toLowerCase();
		if (x.startsWith("t")) return true;
		return false;
	}

	public File[][] getControlCelFiles() {
		return controlCelFiles;
	}

	public File getLogFile() {
		return logFile;
	}

	public File getResultsDirectory() {
		return resultsDirectory;
	}

	public String getResultsName() {
		return resultsName;
	}

	public int getSizeOfOligo() {
		return sizeOfOligo;
	}

	public File[] getTpmapFiles() {
		return tpmapFiles;
	}

	public File[][] getTreatmentCelFiles() {
		return treatmentCelFiles;
	}


	public double getTargetMedian() {
		return targetMedian;
	}


	public File getAppsDirectory() {
		return appsDirectory;
	}


	public int getMaxMemory() {
		return maxMemory;
	}

	public File[] getNormalizedControlReplicas() {
		return normalizedControlReplicas;
	}

	public void setNormalizedControlReplicas(File[] normalizedControlReplicas) {
		this.normalizedControlReplicas = normalizedControlReplicas;
	}

	public File[] getNormalizedTreatmentReplicas() {
		return normalizedTreatmentReplicas;
	}

	public void setNormalizedTreatmentReplicas(File[] normalizedTreatmentReplicas) {
		this.normalizedTreatmentReplicas = normalizedTreatmentReplicas;
	}

	public int getMinimumNumberOfOligos() {
		return minimumNumberOfOligos;
	}

	public int getNumberOfRandomPermutations() {
		return numberOfRandomPermutations;
	}

	public File getRApp() {
		return rApp;
	}

	public File getSymmetricNullApp() {
		return symmetricNullApp;
	}

	public boolean isUsePseudoMedian() {
		return usePseudoMedian;
	}

	public boolean isUseRandomPermutation() {
		return useRandomPermutation;
	}

	public boolean isUseSymmetricNull() {
		return useSymmetricNull;
	}

	public int getWindowSize() {
		return windowSize;
	}

	public String getGenomeVersion() {
		return genomeVersion;
	}

	public boolean isMakeReducedWindows() {
		return makeReducedWindows;
	}

	public String getNumberOfIntervals() {
		return numberOfIntervals;
	}

	public int getMaxGap() {
		return maxGap;
	}

	public File getGenomeFastaDirectory() {
		return genomeFastaDirectory;
	}

	public int getMaxNumberPeaks() {
		return maxNumberPeaks;
	}

	public int getMinNumberOligosInSubWin() {
		return minNumberOligosInSubWin;
	}

	public int getPeakPickerWindowSize() {
		return peakPickerWindowSize;
	}

	public int getSubWindowSize() {
		return subWindowSize;
	}

	public File[] getRepeatRegionFiles() {
		return repeatRegionFiles;
	}

	public File getJava() {
		return java;
	}

	public boolean isQuantileNormalizeAll() {
		return quantileNormalizeAll;
	}

	public boolean isMedianScaleUsingControlSeqs() {
		return medianScaleUsingControlSeqs;
	}

	public String getPermutationType() {
		return permutationType;
	}

	public boolean isRandomizeData() {
		return randomizeData;
	}

	public String getStrand() {
		return strand;
	}

}
