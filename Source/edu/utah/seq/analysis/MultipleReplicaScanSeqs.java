package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import trans.tpmap.WindowMaker;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.BarParser;
import util.bio.cluster.*;



/**
 * @author Nix
 * */
public class MultipleReplicaScanSeqs {

	//user defined fields
	private File[] treatmentPointDirs;
	private File[] controlPointDirs;
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private int windowSize = -1;
	private int peakShift = -1;
	private int halfPeakShift = 0;
	private int minimumNumberReadsInWindow = 10;
	private float minimumFDR = 0.9f;
	private boolean cluster = true;
	private boolean poolReplicas = false;
	private boolean bypassVarianceOutlierFiltering = false;

	//internal fields
	private WindowMaker windowMaker; 
	private int[][] windows;
	private int totalNumberWindows = 0;
	private int totalNumberWindowsPassingThresholds = 0;
	private SmoothingWindowInfo[] smoothingWindowInfo;
	private int totalInterrogatedTreatmentObservations = 0;
	private int totalInterrogatedControlObservations = 0;
	private double numberTreatmentObservations;
	private double numberControlObservations;
	private int numberTreatmentReplicas;
	private int numberControlReplicas;
	private int totalNumberReplicas;
	private int[] numberTreatmentReads;
	private int[] numberControlReads;
	private HashMap<String,PointData[]> treatmentPlusPointData;
	private HashMap<String,PointData[]> treatmentMinusPointData;
	private HashMap<String,PointData[]> controlPlusPointData;
	private HashMap<String,PointData[]> controlMinusPointData;
	private String genomeVersion;
	private boolean deleteTempFiles = true;
	private String[] scoreNames;
	private String[] scoreDescriptions;
	private String[] scoreUnits;
	private String adapterName = "chrAdapt";
	private File matrixFile;
	private File tempRDirectory;
	private PrintWriter matrixFileOut;
	private File swiFile;
	private boolean verbose = true;
	
	//by chromosome
	private String chromosome;
	private PointData[] treatmentChromPlus = null;
	private PointData[] treatmentChromMinus = null;	
	private PointData[] controlChromPlus = null;
	private PointData[] controlChromMinus = null;

	//constructor	
	/**Stand alone.*/
	public MultipleReplicaScanSeqs(String[] args){
		long startTime = System.currentTimeMillis();
		//set fields
		processArgs(args);
		//launch
		run();
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		if (verbose) System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}
	
	/**For integration with ChIPSeq and RNASeq*/
	public MultipleReplicaScanSeqs(File[] treatmentPointDirs, File[] controlPointDirs, File saveDirectory, File fullPathToR, int windowSize, int peakShift, int minimumNumberReadsInWindow, float minimumFDR, boolean bypassVarianceOutlierFiltering, boolean verbose){
		this.treatmentPointDirs = treatmentPointDirs;
		this.controlPointDirs = controlPointDirs;
		this.saveDirectory = saveDirectory;
		this.fullPathToR = fullPathToR;
		this.windowSize = windowSize;
		this.peakShift = peakShift;
		this.minimumNumberReadsInWindow = minimumNumberReadsInWindow;
		this.minimumFDR = minimumFDR;
		this.bypassVarianceOutlierFiltering = bypassVarianceOutlierFiltering;
		this.verbose = verbose;
		
		//check DESeq?
		DefinedRegionDifferentialSeq.estimateDispersions(fullPathToR, saveDirectory);
		
		//System.out.println("ED? "+useEstimateDispersions);
		
		halfPeakShift = (int)Math.round( ((double)peakShift)/2 );
		if (peakShift == 0 && windowSize == -1) windowSize = peakShift;
		saveDirectory.mkdir();
		tempRDirectory = new File (saveDirectory, "TempRDir_"+Passwords.createRandowWord(7));
		tempRDirectory.mkdir();
		if (deleteTempFiles) tempRDirectory.deleteOnExit();
		matrixFile = new File (tempRDirectory, "matrixFile.txt");
		
		run();
	}
	

	public void run(){	
		//make window maker 
		windowMaker = new WindowMaker(windowSize,minimumNumberReadsInWindow);
		
		//set descriptors
		setScoreStrings();

		//fetch counts
		System.out.println("Stats and parameters...");
		calculateReadCountStatistics();
		System.out.println("\t"+peakShift+"\tPeak shift");
		System.out.println("\t"+windowSize+"\tWindow size");
		System.out.println("\t"+minimumFDR+"\tMinimum Window FDR");
		System.out.println("\t"+minimumNumberReadsInWindow+"\tMinimum number reads in window");
		System.out.println("\t"+bypassVarianceOutlierFiltering+"\tBypass DESeq's variance outlier filtering");

		//make print writer for generating the matrix file
		try{
			matrixFileOut = new PrintWriter ( new FileWriter( matrixFile));
		} catch (Exception e){
			System.err.println("\nProblem making PrintWriter for window count matrix file");
			e.printStackTrace();
		}

		//for each chromosome
		if (verbose) System.out.print("\nScanning chromosomes");
		String[] chromosomes = fetchAllChromosomes();
		for (int i=0; i< chromosomes.length; i++){
			chromosome = chromosomes[i];
			if (chromosome.contains(adapterName)) continue;
			windowScanChromosome();
		}
		if (verbose) System.out.println();
		matrixFileOut.close();
		
		//enough to cluster?
		if (totalNumberReplicas <= 2) cluster = false;

		//convert binomial pvalues to FDRs
		if (verbose) System.out.println("\nCalculating negative binomial p-values and FDRs in R using DESeq (http://www-huber.embl.de/users/anders/DESeq/). This requires patience, 64bit R, and lots of RAM.");
		//call DESeq
		File[] deseqResults = executeDESeq();

		//parse window stats, if no windows then exit
		if (verbose) System.out.println("\nParsing results...");
		if (deseqResults == null || parseDESeqStatResults (deseqResults[0]) == false) return;

		//save window data 
		if (verbose) System.out.println("\nSaving "+totalNumberWindowsPassingThresholds+" serialized window data.");
		swiFile = new File (saveDirectory, "binaryWindowData.swi");
		IO.saveObject(swiFile, smoothingWindowInfo);
		IO.zipAndDelete(swiFile);
		swiFile = new File (saveDirectory, "binaryWindowData.swi.zip");

		//write out bar file graphs, call after saving since Info is modified
		if (verbose) System.out.println("\nWriting bar file window summary graphs.");
		writeBarFileGraphs();

		//cluster
		File clusterFile =null;
		if (cluster){
			//parse and cluster variance corrected count data
			if (verbose) System.out.println("\nHierarchical clustering replicas based on DESeq's variance corrected read counts, experimental.");
			clusterFile = new File (saveDirectory,"clusterPlot.png");
			parseAndClusterDESeqDataResults (deseqResults[1], clusterFile);
		}

		//clean up
		if (deleteTempFiles) {
			matrixFile.delete();
			IO.deleteDirectory(tempRDirectory);
		}
		

	}

	/**Writes stair step window bar graph files*/
	public void writeBarFileGraphs(){
		//make directories
		File windowSummaryTracks = new File (saveDirectory, "WindowSummaryTracks");
		windowSummaryTracks.mkdir();
		File FDR = new File(windowSummaryTracks, "FDR");
		FDR.mkdir();
		File log2Ratio = new File(windowSummaryTracks, "Log2Ratio");
		log2Ratio.mkdir();
		//scores = FDR, Log2Ratio
		//for each chromosome
		for (int i=0; i< smoothingWindowInfo.length; i++){
			Info info = smoothingWindowInfo[i].getInfo();
			SmoothingWindow[] sm = smoothingWindowInfo[i].getSm();
			//convert FDRs sign to match log2Ratio
			for (int x=0; x< sm.length; x++){
				float[] scores = sm[x].getScores();
				if (scores[1]<0){
					scores[0] = -1* scores[0];
				}
			}
			//save graphs
			saveSmoothedHeatMapData (0, sm, info, FDR, "#00FF00", true); //green
			saveSmoothedHeatMapData (1, sm, info, log2Ratio, "#FF0000", true); //red
		}
	}

	/**Saves bar heatmap/ stairstep graph files*/
	public void saveSmoothedHeatMapData (int scoreIndex, SmoothingWindow[] sm, Info info, File dir, String color, boolean posNeg){
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = info.getNotes();		
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		//color red
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, color);
		//what's the source
		String fileNames = Misc.stringArrayToString(IO.fetchFileNames(treatmentPointDirs),",");
		if (controlPointDirs !=null) fileNames = fileNames+" vs "+Misc.stringArrayToString(IO.fetchFileNames(controlPointDirs),",");
		map.put(BarParser.SOURCE_TAG, fileNames);
		//what's window size
		map.put(BarParser.WINDOW_SIZE, windowSize+"");
		//what's the peak shift
		map.put(BarParser.BP_3_PRIME_SHIFT, halfPeakShift+"");
		//what's the unit on the scores
		map.put(BarParser.UNIT_TAG, scoreUnits[scoreIndex]);
		//description
		map.put(BarParser.DESCRIPTION_TAG, scoreDescriptions[scoreIndex]);
		//save in info
		info.setNotes(map);
		//get heatmap positions and values
		PointData pd;		
		if (posNeg){
			HeatMapMakerPosNeg hm = new HeatMapMakerPosNeg(scoreIndex, 0, 0);
			pd = hm.makeHeatMapPositionValues(sm);
		}
		else {
			HeatMapMaker hm = new HeatMapMaker(scoreIndex,0);
			pd = hm.makeHeatMapPositionValues(sm, false);
		}
		pd.setInfo(info);
		pd.writePointData(dir);
		//clean up
		pd.nullPositionScoreArrays();
	}

	/**Window scans a chromosome collecting read count data and calculating binomial p-values.*/
	public void windowScanChromosome(){
		//fetch data
		if (fetchData() == false) return;
		//fetch, shift, and merge all positions from the treatment and control
		int[] positions = fetchShiftStripPointData();
		//make windows using all of the reads
		makeWindows(positions);
		//any windows?
		if (windows.length == 0){
			if (verbose) System.out.println("\n\tSkipping "+chromosome+". No windows found with the minimum number of reads ("+minimumNumberReadsInWindow+") within a window size of "+windowSize);
			return;
		}
		if (verbose) System.out.print(".");
		//scan and write out data
		scan();
	}

	/**Makes a common set of windows using merged positions.*/
	public void makeWindows(int[] positions){
		windows = windowMaker.makeWindows(positions);
		//assign bp positions
		for (int i=0; i< windows.length; i++){
			windows[i][0] = positions[windows[i][0]];
			windows[i][1] = positions[windows[i][1]]+1;	//last base isn't included
		}
		totalNumberWindows += windows.length;
	}

	/**Fetches the names of all the chromosomes in the data.*/
	public String[] fetchAllChromosomes(){
		HashSet<String> c = new HashSet<String>();
		Iterator<String> it = treatmentPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = treatmentMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		if(controlPointDirs != null) {
			it = controlPlusPointData.keySet().iterator();
			while (it.hasNext()) c.add(it.next());
			it = controlMinusPointData.keySet().iterator();
			while (it.hasNext()) c.add(it.next());
		}
		return Misc.hashSetToStringArray(c);
	}


	/**Returns the deseq stats file and the variance stabilized data.*/
	private File[] executeDESeq(){
		File rResultsStats = new File (tempRDirectory, "DESeqResults.txt");
		File rResultsData = new File (tempRDirectory, "DESeqResultsData.txt");
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("library(DESeq)\n");
			sb.append("countsTable = read.delim('"+matrixFile.getCanonicalPath()+"', header=FALSE)\n");
			sb.append("rownames(countsTable) = countsTable$V1\n");
			sb.append("countsTable = countsTable[,-1]\n");
			sb.append("conds = c(");
			for (int i=0; i< numberTreatmentReplicas; i++){
				sb.append("'T',");
			}
			sb.append("'N'");
			for (int i=1; i< numberControlReplicas; i++){
				sb.append(",'N'");
			}
			sb.append(")\n");
			sb.append("cds = newCountDataSet( countsTable, conds)\n");
			sb.append("rm(countsTable)\n");
			sb.append("gc()\n");
			sb.append("cds = estimateSizeFactors( cds )\n");
			if (totalNumberReplicas == 2 || poolReplicas) {
				sb.append("cds = estimateDispersions(cds,method='blind',sharingMode='fit-only')\n");
				if (verbose) System.out.println("\tCalling DESeq's estimateDispersions(cds,method='blind',sharingMode='fit-only') method");
			}
			else if (bypassVarianceOutlierFiltering){
				sb.append("cds = estimateDispersions(cds, sharingMode='fit-only')\n");
				if (verbose) System.out.println("\tCalling DESeq's estimateDispersions(cds, sharingMode='fit-only')");
			}
			else {
				sb.append("cds = estimateDispersions(cds)\n");
				if (verbose) System.out.println("\tCalling DESeq's default estimateDispersions(cds)");
			}
			sb.append("res = nbinomTest( cds, 'N', 'T', pvals_only = FALSE)\n");
			//filter for significant padj
			sb.append("res = res[res$padj <= "+minimumFDR+" ,]\n");
			sb.append("gc()\n");
			//Recalculate log2 ratio with add one, note flip
			sb.append("res[,6] = log2((1+res[,4])/(1+res[,3]))\n");
			//Fred  adjPVal
			sb.append("res[,8] = -10 * log10(res[,8])\n");
			//Parse name, padj, log2ratio
			sb.append("res = res[,c(1,8,6)]\n");
			sb.append("gc()\n");
			//note, the order of the rows is the same as the input
			sb.append("write.table(res, file = '"+rResultsStats.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");
			//pull variance stabilized data for each replica, note these are not flipped
			if (cluster) {
				sb.append("res = getVarianceStabilizedData( cds )\n");
				sb.append("write.table(res, file = '"+rResultsData.getCanonicalPath()+"', quote=FALSE, sep ='\t', row.names = FALSE, col.names = FALSE)\n");
			}
			
			//write script to file
			File scriptFile = new File (tempRDirectory, "RScript.txt");
			File rOut = new File(tempRDirectory, "RScript.txt.Rout");
			IO.writeString(sb.toString(), scriptFile);

			//make command
			String[] command = new String[] {
					fullPathToR.getCanonicalPath(),
					"CMD",
					"BATCH",
					"--no-save",
					"--no-restore",
					scriptFile.getCanonicalPath(),
					rOut.getCanonicalPath()};

			//execute command
			IO.executeCommandLine(command);
			
			if (rResultsStats.exists() == false || (cluster && rResultsData.exists() == false)) {
				//did it fail to fit?
				if (deseqFailedToFitError(rOut)){
					System.err.print("\nWARNING: DESeq's parametric dispersion fit failed.");
					if (poolReplicas == false) {
						System.err.println(" Rerunning with pooled replicas...\n");
						poolReplicas = true;
						return executeDESeq();
					}
					else {
						System.err.println(" Aborting.");
						return null;
					}
				}
				throw new IOException("\nR results file(s) doesn't exist. Check temp files in save directory for error.\n");
			}

			//cleanup
			if (deleteTempFiles) {
				rOut.delete();
				scriptFile.delete();
			}
			

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
			return null;
		}
		return new File[]{rResultsStats,rResultsData};
	}
	
	/**Looks at file's contents for the "dispersion fit" phrase from DESeq when it fails to fit or converge.*/
	public static boolean deseqFailedToFitError(File rOut){
		//look at log and see if it failed to fit
		String[] lines = IO.loadFile(rOut);
		Pattern error = Pattern.compile(".+dispersion fit.+");
		for (String line : lines){
			if (error.matcher(line.toLowerCase()).matches()) return true; 
		}
		return false;
	}

	public boolean parseDESeqStatResults(File results){
		try {
			BufferedReader in = new BufferedReader (new FileReader (results));
			String line;
			String[] tokens;
			Pattern tab = Pattern.compile("\t");
			Pattern coor = Pattern.compile("(.+):(\\d+):(\\d+)");
			float maxFDR = 0;
			ArrayList<SmoothingWindow> alSM = new ArrayList<SmoothingWindow>();
			ArrayList<SmoothingWindowInfo> alSMI = new ArrayList<SmoothingWindowInfo>(); 
			ArrayList<SmoothingWindow> alSMToFixFDR = new ArrayList<SmoothingWindow>();
			String currentChromosome = null;
			int counter = 0;
			while((line = in.readLine())!= null){
				counter++;
				//parse line: name, padj, log2ratio
				tokens = tab.split(line);
				if (tokens.length!=3) Misc.printErrAndExit("One of the DESeq R results rows is malformed -> "+line);
				//parse fdr, watch out for Inf, assign max?
				float fdr = Float.MAX_VALUE;
				if (tokens[1].equals("Inf") == false && tokens[1].equals("NA") == false) {
					fdr = Float.parseFloat(tokens[1]);
					if (fdr > maxFDR) maxFDR = fdr;
				}
				//parse ratio
				float ratio = Float.parseFloat(tokens[2]);
				//parse chromosome:start:stop
				Matcher mat = coor.matcher(tokens[0]);
				if (mat.matches()== false) Misc.printErrAndExit("One of the DESeq R results rows is malformed -> "+line);
				//parse chromosome, start, stop for window
				String chr = mat.group(1);
				int start = Integer.parseInt(mat.group(2));
				int stop = Integer.parseInt(mat.group(3));
				//instantiate new window
				SmoothingWindow sm = new SmoothingWindow(start,stop, new float[]{fdr,ratio});
				//first time through?
				if (currentChromosome == null) currentChromosome = chr;
				//nope new chromosome
				else if (currentChromosome.equals(chr) == false){
					//convert AL to array
					SmoothingWindow[] sma = new SmoothingWindow[alSM.size()];
					alSM.toArray(sma);
					alSM.clear();
					//add new SMI
					alSMI.add(new SmoothingWindowInfo(sma, createInfoObject(currentChromosome)));
					currentChromosome = chr;
					totalNumberWindowsPassingThresholds += sma.length;
				}
				//add window and check fdr
				alSM.add(sm);
				if (fdr == Float.MAX_VALUE) alSMToFixFDR.add(sm);	
			}
			//close last
			SmoothingWindow[] sma = new SmoothingWindow[alSM.size()];
			alSM.toArray(sma);
			alSM = null;
			alSMI.add(new SmoothingWindowInfo(sma, createInfoObject(currentChromosome)));
			totalNumberWindowsPassingThresholds += sma.length;
			//close reader
			in.close();
			//check number
			if (counter == 0){
				if (verbose) System.out.println("\nWARNING: No significant windows found?! FDR threshold is set to "+minimumFDR+" . Try restarting with a less stringent threshold.\n");
				return false;
			}
			if (deleteTempFiles) results.delete();
			//convert Inf to max values * 1%
			maxFDR = maxFDR * 1.01f;
			int numToCovert = alSMToFixFDR.size();
			for (int i=0; i< numToCovert; i++){
				SmoothingWindow swToFix = alSMToFixFDR.get(i);
				float[] scores = swToFix.getScores();
				scores[0] = maxFDR;
			}
			//convert SWI
			smoothingWindowInfo = new SmoothingWindowInfo[alSMI.size()];
			alSMI.toArray(smoothingWindowInfo);
			alSMI = null;	
			return true;
		} catch (Exception e){
			System.err.println("Problem parsing DESeq stats results from R.\n");
			e.printStackTrace();
			System.exit(1);
		}
		return false;
	}

	/**Creates an info object to use for the SmoothingWindowInfo arrays.*/
	private Info createInfoObject(String chrom){
		Info info = new Info();
		HashMap<String,String> notes = new HashMap<String,String>();
		notes.put("totalTreatmentObservations",""+(int)numberTreatmentObservations);
		notes.put("totalControlObservations",""+(int)numberControlObservations);
		notes.put(BarParser.WINDOW_SIZE, windowSize+"");
		notes.put(BarParser.BP_3_PRIME_SHIFT, halfPeakShift+"");
		notes.put(BarParser.DESCRIPTION_TAG, Misc.stringArrayToString(scoreNames, ","));
		notes.put(BarParser.UNIT_TAG, Misc.stringArrayToString(scoreUnits, ","));
		info.setNotes(notes);
		info.setStrand(".");
		info.setChromosome(chrom);
		info.setVersionedGenome(genomeVersion);
		return info;
	}

	
	private void parseAndClusterDESeqDataResults(File results, File clusterFile){
		try {
			//load the counts
			float[][] varCounts = new float[totalNumberWindows][totalNumberReplicas];
			BufferedReader in = new BufferedReader (new FileReader (results));
			String line;
			String[] tokens;
			Pattern tab = Pattern.compile("\t");
			int counter = 0;
			while ((line=in.readLine())!=null){
				//parse line
				tokens = tab.split(line);
				if (tokens.length!=totalNumberReplicas) Misc.printErrAndExit("One of the DESeq data R results rows is malformed -> "+line);
				//parse the scores
				for (int i=0; i< totalNumberReplicas; i++) varCounts[counter][i] = Float.parseFloat(tokens[i]);
				counter++;
			}
			if (counter != totalNumberWindows) Misc.printErrAndExit("\nError, the number of DESeq var corrected count lines ("+counter+") is different than the number of windows ("+totalNumberWindows+") ?!");
			in.close();
			if (deleteTempFiles) results.delete();
			
			//save the float[]s
			File[] files = new File[totalNumberReplicas];
			File tempDir = new File (saveDirectory,".TempClusterDir_"+Passwords.createRandowWord(7));
			tempDir.mkdir();
			tempDir.deleteOnExit();
			Pattern slash = Pattern.compile("/");
			for (int i=0; i< numberTreatmentReplicas; i++){
				String fullPathName = treatmentPointDirs[i].getCanonicalPath();
				fullPathName = slash.matcher(fullPathName).replaceAll("_");
				files[i] = new File(tempDir, fullPathName);
				files[i].deleteOnExit();
			}
			int index = 0;
			for (int i=numberTreatmentReplicas; i< totalNumberReplicas ; i++){
				String fullPathName = controlPointDirs[index++].getCanonicalPath();
				fullPathName = slash.matcher(fullPathName).replaceAll("_");
				files[i] = new File(tempDir, fullPathName);
				files[i].deleteOnExit();
			}
			for (int x=0; x< totalNumberReplicas; x++) IO.saveObject(files[x], varCounts[x]);

			//cluster
			new HierarchicalClustering(files, clusterFile);

			//delete tempDir, not working!
			IO.deleteDirectory(tempDir);
			IO.deleteDirectoryViaCmdLine(tempDir);
			
		} catch (Exception e){
			System.err.println("Problem parsing and clustering DESeq data results from R.\n");
			e.printStackTrace();
			System.exit(1);
		}
	}


	/**Fetchs the data for a particular chromosome.*/
	public boolean fetchData(){		
		//fetch treatments
		treatmentChromPlus = treatmentPlusPointData.get(chromosome);
		treatmentChromMinus = treatmentMinusPointData.get(chromosome);
		//fetch controls
		controlChromPlus = controlPlusPointData.get(chromosome);
		controlChromMinus = controlMinusPointData.get(chromosome);
		//check for nulls
		if (treatmentChromPlus == null || treatmentChromMinus == null || controlChromPlus == null || controlChromMinus == null){
			if (verbose) System.out.println("\n\tWarning: failed to find stranded data from all datasets for "+chromosome+". Skipping!");
			return false;
		}
		//check number
		if (treatmentChromPlus.length != numberTreatmentReplicas || treatmentChromMinus.length != numberTreatmentReplicas || controlChromPlus.length != numberControlReplicas || controlChromMinus.length != numberControlReplicas){
			if (verbose) System.out.println("\n\tWarning: The number of treatment or control replicas differ for chromosome -> "+chromosome+". Skipping!");
			return false;
		}
		//increment number of obs
		for (int i=0; i< numberTreatmentReplicas; i++){
			numberTreatmentReads[i] += treatmentChromPlus[i].getInfo().getNumberObservations();
			numberTreatmentReads[i] += treatmentChromMinus[i].getInfo().getNumberObservations();
		}
		for (int i=0; i< numberControlReplicas; i++){
			numberControlReads[i] += controlChromPlus[i].getInfo().getNumberObservations();
			numberControlReads[i] += controlChromMinus[i].getInfo().getNumberObservations();
		}		//assign genome
		if (genomeVersion == null){
			if (treatmentChromPlus != null) genomeVersion = treatmentChromPlus[0].getInfo().getVersionedGenome();
			else if (treatmentChromMinus != null) genomeVersion = treatmentChromMinus[0].getInfo().getVersionedGenome();
		}
		return true;
	}

	/**For a given loaded chromosome of data collects number of observations under each window. Not stranded. Writes to file for negative binomial*/
	private void scan(){
		//for each window 
		for (int i=0; i< windows.length; i++){

			//fetch stranded alignments for replica
			float[] scores = new float[totalNumberReplicas];

			int index = 0;
			//treatment
			for (int j=0; j< numberTreatmentReplicas; j++){
				if (treatmentChromPlus[j]!= null) scores[index] = treatmentChromPlus[j].sumScoreBP(windows[i][0], windows[i][1]);
				if (treatmentChromMinus[j]!= null) scores[index] += treatmentChromMinus[j].sumScoreBP(windows[i][0], windows[i][1]);
				index++;
			}
			//control
			for (int j=0; j< numberControlReplicas; j++){
				if (controlChromPlus[j]!= null) scores[index] = controlChromPlus[j].sumScoreBP(windows[i][0], windows[i][1]);
				if (controlChromMinus[j]!= null) scores[index] += controlChromMinus[j].sumScoreBP(windows[i][0], windows[i][1]);
				index++;
			}
			//write out scores, chrom:start:stop tabs scores as INTS!
			matrixFileOut.print(chromosome+":"+windows[i][0]+":"+windows[i][1]);
			for (int x =0; x < scores.length; x++){
				matrixFileOut.print("\t");
				matrixFileOut.print((int)scores[x]);
			}
			matrixFileOut.println();
		}
	}

	/**Collects and calculates a bunch of stats re the PointData.*/
	private void calculateReadCountStatistics(){
		//fetch treatment PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (treatmentPointDirs);
		treatmentPlusPointData = PointData.convertArrayList2Array(combo[0]);
		treatmentMinusPointData = PointData.convertArrayList2Array(combo[1]);
		numberTreatmentObservations = PointData.totalObservationsMultiPointData(treatmentPlusPointData);
		numberTreatmentObservations += PointData.totalObservationsMultiPointData(treatmentMinusPointData);
		System.out.println("\t"+(int)numberTreatmentObservations+" Treatment Observations");

		//likewise control data
		combo = PointData.fetchStrandedPointDataNoMerge (controlPointDirs);
		controlPlusPointData = PointData.convertArrayList2Array(combo[0]);
		controlMinusPointData = PointData.convertArrayList2Array(combo[1]);
		numberControlObservations = PointData.totalObservationsMultiPointData(controlPlusPointData);
		numberControlObservations += PointData.totalObservationsMultiPointData(controlMinusPointData);
		System.out.println("\t"+(int)numberControlObservations+" Control Observations");		

		//set number of replicas
		numberTreatmentReplicas = treatmentPointDirs.length;
		numberControlReplicas = controlPointDirs.length;
		totalNumberReplicas = numberTreatmentReplicas + numberControlReplicas;

		numberTreatmentReads = new int[numberTreatmentReplicas];
		numberControlReads = new int[numberControlReplicas];
	}

	/**Shifts the positions halfPeakShift (+ for sense, - for antisense) sets the positions into the data
	 * returns all of the positions after sorting. May replace all scores with 1 if stripScores == true.*/
	private int[] fetchShiftStripPointData(){
		ArrayList<int[]> posAL = new ArrayList<int[]>();
		for (int i=0; i< treatmentChromPlus.length; i++){
			int[] p = treatmentChromPlus[i].getPositions();
			if (halfPeakShift !=0) addShift(p,halfPeakShift);
			posAL.add(p);
			treatmentChromPlus[i].stripScores();
		}
		for (int i=0; i< treatmentChromMinus.length; i++){
			int[] p = treatmentChromMinus[i].getPositions();
			if (halfPeakShift !=0) addShift(p, -1*halfPeakShift);
			posAL.add(p);
			treatmentChromMinus[i].stripScores();
		}
		for (int i=0; i< controlChromPlus.length; i++){
			int[] p = controlChromPlus[i].getPositions();
			if (halfPeakShift !=0) addShift(p,halfPeakShift);
			posAL.add(p);
			controlChromPlus[i].stripScores();
		}
		for (int i=0; i< controlChromMinus.length; i++){
			int[] p = controlChromMinus[i].getPositions();
			if (halfPeakShift !=0) addShift(p, -1*halfPeakShift);
			posAL.add(p);
			controlChromMinus[i].stripScores();
		}
		//merge
		int[][] toMerge = new int[posAL.size()][];
		for (int i=0; i< posAL.size(); i++) toMerge[i] = posAL.get(i);
		int[] merged = Num.collapseIntArray(toMerge);
		Arrays.sort(merged);
		//return
		return merged;
	}

	/**Adds the toAdd to each int.*/
	public static void addShift(int[] positions, int toAdd){
		for (int i=0; i< positions.length; i++){
			positions[i] += toAdd;
			if (positions[i]<0) positions[i] = 0;
		}
	}

	/**Sets score names/ descriptions/ units base on whether control data is present.*/
	public void setScoreStrings(){
		scoreNames = new String[]{
				"FDR",
				"Log2((sumT+1)/(sumC+1))",
		};
		scoreDescriptions = new String[]{
				"Benjamini and Hochberg FDR estimation based on the DESeq's negative binomial p-values",
				"Log2((sumT+1)/(sumC+1)) of linearly scaled sumT and sumC",
		};
		scoreUnits = new String[]{	
				"-10Log10(FDR)",
				"Log2((sumT+1)/(sumC+1))",
		};
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MultipleReplicaScanSeqs(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 't': treatmentPointDirs = IO.extractFiles(args[++i]); break;
					case 'c': controlPointDirs = IO.extractFiles(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'p': peakShift = Integer.parseInt(args[++i]); break;
					case 'w': windowSize = Integer.parseInt(args[++i]); break;
					case 'm': minimumNumberReadsInWindow = Integer.parseInt(args[++i]); break;
					case 'f': minimumFDR = Float.parseFloat(args[++i]); break;
					case 'd': deleteTempFiles = false; break;
					case 'b': bypassVarianceOutlierFiltering = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for point directories
		if (treatmentPointDirs == null || treatmentPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your treatment PointData directories(s)!\n");
		//only one directory look deeper
		if (treatmentPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(treatmentPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) treatmentPointDirs = otherDirs;
		}
		//control data
		if (controlPointDirs != null){
			if (controlPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your control PointData directories(s)!\n");
			//only one directory look deeper
			if (controlPointDirs.length == 1){
				File[] otherDirs = IO.extractOnlyDirectories(controlPointDirs[0]);
				if (otherDirs != null && otherDirs.length > 0) controlPointDirs = otherDirs;
			}
		}
		//set half peak shift and windowSize
		if (peakShift == -1) Misc.printExit("\nPlease enter a peak shift, the results from running the PeakShiftFinder application or use 100-150 if you don't know.\n");
		halfPeakShift = (int)Math.round( ((double)peakShift)/2 );
		if (peakShift == 0 && windowSize == -1) windowSize = peakShift;

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory to save results.\n");
		if (saveDirectory.exists() == false && saveDirectory.mkdir() == false)  Misc.printExit("\nError: could not find nor make the save directory.\n");
		
		//check for R and required libraries
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}
		else {
			String errors = IO.runRCommandLookForError("library(qvalue)", fullPathToR, saveDirectory);
			if (errors == null || errors.length() !=0){
				Misc.printExit("\nError: Cannot find the required R library.  Did you install qvalue " +
						"(http://genomics.princeton.edu/storeylab/qvalue/)?  See the author's websites for installation instructions. Once installed, " +
						"launch an R terminal and type 'library(qvalue)' to see if it is present. R error message:\n\t\t"+errors+"\n\n");
			}
		}
		
		//look for estimateDispersions() function
		DefinedRegionDifferentialSeq.estimateDispersions(fullPathToR, saveDirectory);

		//make matrixFile
		tempRDirectory = new File (saveDirectory, "TempRDir_"+Passwords.createRandowWord(7));
		tempRDirectory.mkdir();
		matrixFile = new File (tempRDirectory, "matrixFile.txt");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Multiple Replica Scan Seqs: May 2012                      **\n" +
				"**************************************************************************************\n" +
				"MRSS uses a sliding window and Ander's DESeq negative binomial pvalue -> Benjamini & \n" +
				"Hochberg FDR statistics to identify enriched and reduced regions in a genome. Both\n" +
				"treatment and control PointData sets are required, one or more biological replicas.\n" +
				"MRSS generates window level differential count tracks for the FDR and normalized\n" +
				"log2Ratio as well as a binary window objec xxx.swi file for downstream use by the\n" +
				"EnrichedRegionMaker. MRSS also makes use of DESeq's variance corrected count data to\n" +
				"cluster your biological replics. Given R's poor memory management, running DESeq\n" +
				"requires lots of RAM, 64bit R, and 1-3 hrs.\n"+

				"\nOptions:\n"+
				"-s Save directory, full path.\n"+
				"-t Treatment replica PointData directories, full path, comma delimited, no spaces,\n" +
				"       one per biological replica. Use the PointDataManipulator app to merge same\n" +
				"       replica and technical replica datasets. Each directory should contain stranded\n" +
				"       chromosome specific xxx_-/+_.bar.zip files. Alternatively, provide one\n" +
				"       directory that contains multiple biological replical PointData directories.\n" +
				"-c Control replica PointData directories, ditto. \n" +
				"-r Full path to 64bit R loaded with DESeq library, defaults to '/usr/bin/R' file, see\n" +
				"       http://www-huber.embl.de/users/anders/DESeq/ . Type 'library(DESeq)' in\n" +
				"       an R terminal to see if it is installed.\n"+
				"-p Peak shift, average distance between + and - strand peaks for chIP-Seq data, see\n" +
				"       PeakShiftFinder or set it to 100bp. For RNA-Seq set it to 0. It will be used\n" +
				"       to shift the PointData by 1/2 the peak shift.\n"+
				"-w Window size, defaults to the peak shift. For chIP-Seq data, a good alternative \n" +
				"       is the peak shift plus the standard deviation, see the PeakShiftFinder app.\n" +
				"       For RNA-Seq data, set this to 100-250.\n"+
				
				"\nAdvanced Options:\n"+
				"-f Minimum FDR threshold for saving windows, defaults to 0.9\n"+
				"-m Minimum number of reads in a window, defaults to 10\n"+
				"-d Don't delete temp files\n" +
				"-b Bypass DESeq variance outlier filtering with estimateDispersions(cds,\n" +
				"       sharingMode='fit-only'), recommended for biologically noisy data.\n\n"+

				"\nExample: java -Xmx4G -jar pathTo/USeq/Apps/MultipleReplicaScanSeqs -t\n" +
				"      /Data/PolIIRep1/,/Data/PolIIRep2/ -c /Data/Input1/,Data/Input2/ -s\n" +
				"      /Data/PolIIResults/ -p 150 -w 250 -b \n\n" +

		"**************************************************************************************\n");

	}

	public int getTotalInterrogatedTreatmentObservations() {
		return totalInterrogatedTreatmentObservations;
	}

	public int getTotalInterrogatedControlObservations() {
		return totalInterrogatedControlObservations;
	}

	public double getNumberTreatmentObservations() {
		return numberTreatmentObservations;
	}

	public double getNumberControlObservations() {
		return numberControlObservations;
	}

	public File getSwiFile() {
		return swiFile;
	}
}
