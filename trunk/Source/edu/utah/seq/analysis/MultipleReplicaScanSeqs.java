package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import trans.tpmap.WindowMaker;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.BarParser;
import edu.utah.seq.useq.apps.Bar2USeq;
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
	private int minimumNumberReadsInWindow = 15;

	//internal fields
	private WindowMaker windowMaker; 
	private int[][] windows;
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
	public MultipleReplicaScanSeqs(File[] treatmentPointDirs, File[] controlPointDirs, File saveDirectory, File fullPathToR, int windowSize, int peakShift, int minimumNumberReadsInWindow, boolean verbose){
		this.treatmentPointDirs = treatmentPointDirs;
		this.controlPointDirs = controlPointDirs;
		this.saveDirectory = saveDirectory;
		this.fullPathToR = fullPathToR;
		this.windowSize = windowSize;
		this.peakShift = peakShift;
		this.minimumNumberReadsInWindow = minimumNumberReadsInWindow;
		this.verbose = verbose;
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
		System.out.println("\t"+minimumNumberReadsInWindow+"\tMinimum number reads in window");

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

		//convert binomial pvalues to AdjPs
		if (verbose) System.out.println("\nCalculating negative binomial p-values and AdjPs in R using DESeq2. This requires patience, 64bit R, and lots of RAM.");
		//call DESeq
		File deseqResults = executeDESeq2();

		//parse window stats, if no windows then exit
		if (verbose) System.out.println("\nParsing results...");
		if (parseDESeq2Results (deseqResults) == false) return;

		//save window data 
		if (verbose) System.out.println("\nSaving "+totalNumberWindowsPassingThresholds+" serialized window data.");
		swiFile = new File (saveDirectory, "binaryWindowData.swi");
		IO.saveObject(swiFile, smoothingWindowInfo);
		IO.zipAndDelete(swiFile);
		swiFile = new File (saveDirectory, "binaryWindowData.swi.zip");

		//write out bar file graphs, call after saving since Info is modified
		if (verbose) System.out.println("\nWriting bar file window summary graphs.");
		writeBarFileGraphs();

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
		File AdjP = new File(windowSummaryTracks, "AdjP");
		AdjP.mkdir();
		File log2Ratio = new File(windowSummaryTracks, "Log2Ratio");
		log2Ratio.mkdir();
		//scores = AdjP, Log2Ratio
		//for each chromosome
		for (int i=0; i< smoothingWindowInfo.length; i++){
			Info info = smoothingWindowInfo[i].getInfo();
			SmoothingWindow[] sm = smoothingWindowInfo[i].getSm();
			//convert AdjPs sign to match log2Ratio
			for (int x=0; x< sm.length; x++){
				float[] scores = sm[x].getScores();
				if (scores[1]<0){
					scores[0] = -1* scores[0];
				}
			}
			//save graphs
			saveSmoothedHeatMapData (0, sm, info, AdjP, "#00FF00", true); //green
			saveSmoothedHeatMapData (1, sm, info, log2Ratio, "#FF0000", true); //red
		}
		//convert to useq format
		new Bar2USeq(AdjP,true);
		new Bar2USeq(log2Ratio,true);
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
	private File executeDESeq2(){
		File rResultsStats = new File (tempRDirectory, "dESeq2Results.txt");
		File clusterFile = new File (saveDirectory,"sampleClusterPlot.pdf");
		try {
			//make R script
			StringBuilder sb = new StringBuilder();
			sb.append("library(DESeq2); library(gplots); library(RColorBrewer)\n");
			sb.append("countTable = read.delim('"+matrixFile.getCanonicalPath()+"', header=FALSE)\n");
			sb.append("rownames(countTable) = countTable[,1]\n");
			sb.append("countTable = countTable[,-1]\n");
			
			//make labels
			StringBuilder groups = new StringBuilder();
			StringBuilder replicas = new StringBuilder();
			for (int i=0; i< numberTreatmentReplicas; i++){
				groups.append("'T',");
				replicas.append("'T"+i+"',");
			}
			groups.append("'C'");
			replicas.append("'C0'");
			for (int i=1; i< numberControlReplicas; i++){
				groups.append(",'C'");
				replicas.append(",'C"+i+"'");
			}
			
			sb.append("sampleInfo = data.frame(condition=as.factor(c("+groups+")))\n");
			sb.append("rownames(sampleInfo) = c("+replicas+")\n");
			sb.append("cds = DESeqDataSetFromMatrix(countData=countTable, colData=sampleInfo, design = ~condition)\n");
			sb.append("cds = DESeq(cds)\n");
			sb.append("rld = rlog(cds)\n");
			sb.append("sampleDists = dist(t(assay(rld)))\n");
			sb.append("sampleDistMatrix <- as.matrix( sampleDists )\n");
			sb.append("colours = colorRampPalette( rev(brewer.pal(9, 'Blues')) )(255)\n");
			sb.append("pdf('"+clusterFile.getCanonicalPath()+"')\n");
			sb.append("heatmap.2( sampleDistMatrix, trace='none', col=colours)\n");
			sb.append("dev.off()\n");
			sb.append("res = results(cds, contrast = c('condition', 'T', 'C'))\n");
			sb.append("res[,6] = -10 * log10(res[,6])\n");
			sb.append("write.table(res, file = '"+rResultsStats.getCanonicalPath()+"', quote=FALSE, sep ='\t')\n");
			
			
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
			
			if (rResultsStats.exists() == false) throw new IOException("\nDESeq2 R results file doesn't exist. Check temp files in save directory for error.\n");

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
		return rResultsStats;
	}

	/**
	 * Parses DESeq2 results:
xxxx	baseMean	log2FoldChange	lfcSE	stat	pvalue	-10Log10(padj)
 0         1              2            3      4        5      6
chr1:137696784:137697035	26.22543451	-4.834447206	0.476843548	-10.13843476	3.73E-24	182.971445
chr1:137696786:137697037	26.35263802	-4.83905225	0.476808555	-10.1488369	3.35E-24	182.971445
chr1:137696787:137697038	26.28706518	-4.834058806	0.477129062	-10.1315539	4.00E-24	182.971445
*/
	public boolean parseDESeq2Results(File results){
		try {
			BufferedReader in = new BufferedReader (new FileReader (results));
			String line;
			String[] tokens;
			Pattern tab = Pattern.compile("\t");
			Pattern coor = Pattern.compile("(.+):(\\d+):(\\d+)");
			float maxAdjP = 0;
			ArrayList<SmoothingWindow> alSM = new ArrayList<SmoothingWindow>();
			ArrayList<SmoothingWindowInfo> alSMI = new ArrayList<SmoothingWindowInfo>(); 
			ArrayList<SmoothingWindow> alSMToFixAdjP = new ArrayList<SmoothingWindow>();
			String currentChromosome = null;
			int counter = 0;
			//skip header
			in.readLine();
			//for each data row
			while((line = in.readLine())!= null){
				counter++;
				//parse line: name, padj, log2ratio
				tokens = tab.split(line);
				if (tokens.length!=7) Misc.printErrAndExit("One of the DESeq2 R results rows is malformed -> "+line);
				//parse -10Log10(adjP), watch out for Inf, assign max?
				float adjP;
				if (tokens[6].equals("Inf")) adjP = Float.MAX_VALUE ;
				else if (tokens[6].equals("NA")) adjP = 0;
				else {
					adjP = Float.parseFloat(tokens[6]);
					if (adjP > maxAdjP) maxAdjP = adjP;
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
				SmoothingWindow sm = new SmoothingWindow(start,stop, new float[]{adjP,ratio});
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
				//add window and check adjP
				alSM.add(sm);
				if (adjP == Float.MAX_VALUE) alSMToFixAdjP.add(sm);	
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
				if (verbose) System.out.println("\nWARNING: No significant windows found?!\n");
				return false;
			}
			if (deleteTempFiles) results.delete();
			//convert Inf to max values * 1%
			maxAdjP = maxAdjP * 1.01f;
			int numToCovert = alSMToFixAdjP.size();
			for (int i=0; i< numToCovert; i++){
				SmoothingWindow swToFix = alSMToFixAdjP.get(i);
				float[] scores = swToFix.getScores();
				scores[0] = maxAdjP;
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
				"AdjP",
				"Log2((sumT+1)/(sumC+1))",
		};
		scoreDescriptions = new String[]{
				"Benjamini and Hochberg adjusted p-value estimation based on the DESeq2's negative binomial",
				"Log2((sumT+1)/(sumC+1)) of linearly scaled sumT and sumC",
		};
		scoreUnits = new String[]{	
				"-10Log10(AdjP)",
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
					case 'd': deleteTempFiles = false; break;
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
		if (windowSize == -1) windowSize = peakShift;

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory to save results.\n");
		if (saveDirectory.exists() == false && saveDirectory.mkdir() == false)  Misc.printExit("\nError: could not find nor make the save directory.\n");
		
		//check for R and required libraries
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}
		else {
			String errors = IO.runRCommandLookForError("library(DESeq2); library(gplots)", fullPathToR, saveDirectory);
			if (errors == null || errors.length() !=0){
				Misc.printExit("\nError: Cannot find the required R libraries?  Did you install DESeq2? gplots? Once installed, " +
						"launch an R terminal and type 'library(DESeq2); library(gplots)' to see if present. R error message:\n\t\t"+errors+"\n\n");
			}
		}

		//make matrixFile
		tempRDirectory = new File (saveDirectory, "TempRDir_"+Passwords.createRandowWord(7));
		tempRDirectory.mkdir();
		matrixFile = new File (tempRDirectory, "matrixFile.txt");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                        Multiple Replica Scan Seqs: May 2014                      **\n" +
				"**************************************************************************************\n" +
				"MRSS uses a sliding window and Ander's DESeq negative binomial pvalue -> Benjamini & \n" +
				"Hochberg AdjP statistics to identify enriched and reduced regions in a genome. Both\n" +
				"treatment and control PointData sets are required, one or more biological replicas.\n" +
				"MRSS generates window level differential count tracks for the AdjP and normalized\n" +
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
				"-m Minimum number of reads in a window, defaults to 15\n"+
				"-d Don't delete temp files\n" +

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
