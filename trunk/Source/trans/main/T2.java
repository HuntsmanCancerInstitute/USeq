package trans.main;
import java.io.*;
import java.util.*;
import java.util.regex.*;

import util.gen.*;

/**Runs entire TiMAT2 suite based on a parameter file.*/
public class T2 {
	//fields
	private File parameterFile;
	private T2Parameter parameters;
	private boolean cluster = false;
	private ChromSet[] replicas;
	private int replicaIndex;
	private File oligoPositions;
	private File[] treatmentReplicas;
	private File[] controlReplicas;
	private File windowDirectory;
	private File enrichedWindows;
	private File reducedWindows = null;
	private File intervalDirectory = null;

	public T2 (String[] args){
		processArgs(args);
		//make a T2Parameters object to hold user options and file locations
		System.out.println("\nParsing your T2 parameter file...");
		parameters = new T2Parameter(parameterFile);
		//convert text cel files to cela files with CelFileConverter?
		if (parameters.getTreatmentCelFiles()[0][0].getName().endsWith("cela") == false){
			System.out.println("\nConverting text xxx.cel files to binary xxx.cela files...");
			convertCelFiles();
			//rename treat and control xxx.cel files to appropriate xxx.cela files in parameters object
			parameters.renameCelFiles();
		}
		else System.out.println("\nCel files already converted to xxx.cela...");
		//normalize files
		System.out.println("\nNormalizing xxx.cela files...");
		if (parameters.isMedianScaleUsingControlSeqs()) System.out.println("\tMedian scaling based on the intensities of the control sequences.");
		else System.out.println("\tMedian scaling based on all of the array intensities.");
		if (parameters.isQuantileNormalizeAll()) System.out.println("\tQuantile normalizing across treatment and control files.");
		else System.out.println("\tQuantile normalizing treatment and control files separately.");
		normalizeCelFiles();
		//make chromosome sets
		System.out.println("\nMaking chromosome sets...");
		makeChromosomeSets();
		//delete split normalized files
		deleteSplitNormalizedFiles();
		//split ChromSet[] into treatment and control replicas
		splitChromSetArray();
		//scan the data
		System.out.println("\nScanning chromosomes...");
		scanChromosomes();
		//windowDirectory = new File (parameters.getResultsDirectory(), "Win");
		System.out.println("\nMerging Window arrays...");
		//merge windows
		mergeWindows();

		//make and process intervals?
		if (parameters.getNumberOfIntervals() != null) {
			//set number interval maker
			System.out.println("\nMaking Interval arrays...");
			setNumberIntervalMaker();
			//load em with intensity info and genomic sequence?
			System.out.println("\nLoading Interval arrays...");
			loadIntervalArrays();
			//find sub binding regions?
			if (parameters.getSubWindowSize() !=0){
				System.out.println("\nFinding Interval peaks...");
				findSubBindingRegions();
			}
			//filter intervals?
			if (parameters.getRepeatRegionFiles() != null){
				System.out.println("\nFlagging Intervals with repetative, low complexity, and simple repeats...");
				filterIntervals();
			}
			//print spreadsheets
			System.out.println("\nPrinting Interval spread sheet reports...");
			printSpreadSheets();
			//print graphs
			System.out.println("\nPrinting Interval graphs...");
			printIntervalGraphs();
		}


		System.out.println("\nT2 is Done!\n");
	}

	/**Returns the number of intervals that remain to be made.*/
	public int numberOfIntervalsToMake(){
		//cal total to make
		String[] numInts = parameters.getNumberOfIntervals().split(",");
		int number = numInts.length;
		if (reducedWindows != null && reducedWindows.exists()) number *= 2;
		//how many made?
		File dir = enrichedWindows.getParentFile();
		String[] fileNames = dir.list();
		int numMade = 0;
		for (int i=0; i< fileNames.length; i++){
			if (fileNames[i].indexOf("Indx") != -1) numMade++;
		}
		return number - numMade;
	}



	/**File management for intervals.*/
	public void makeAndMoveIntervals(){
		intervalDirectory = new File (parameters.getResultsDirectory(), "Intervals");
		if (intervalDirectory.exists() == false) intervalDirectory.mkdir();
		File[] files = enrichedWindows.getParentFile().listFiles();
		//look through all files in window directory and move and rename those that are Intervals into the Interval directory
		Pattern pat = Pattern.compile(".+Indx(\\d+)$");
		for (int i=0; i< files.length; i++){
			String name = files[i].getName();
			//is it an inteval file?
			if (name.indexOf("Indx") != -1){
				String prefix = "en";
				//reduced?
				if (name.startsWith("reduced")) prefix = "rd";
				//parse number on intervals
				Matcher mat = pat.matcher(name);
				mat.matches();
				String numInts = mat.group(1);
				//make new File and move
				File moveTo = new File (intervalDirectory, prefix+parameters.getResultsName()+numInts);
				files[i].renameTo(moveTo);
			}
		}
	}


	/**Creates set number of intervals.*/
	public void setNumberIntervalMaker(){
		if (enrichedWindows.exists() == false || (parameters.isMakeReducedWindows() && reducedWindows.exists() == false)) {
			Misc.printExit("\nError: cannot find merged window array(s)?\n");
		}
		//make general command
		ArrayList command = new ArrayList();
		command.add(parameters.getAppsDirectory()+"/SetNumberIntervalMaker");
		command.add("-a");
		command.add("-s");
		command.add("1");
		command.add("-o");
		command.add(parameters.getMinimumNumberOfOligos()+"");
		command.add("-z");
		command.add(parameters.getSizeOfOligo()+"");
		command.add("-m");
		command.add(parameters.getMaxGap()+"");
		command.add("-n");
		//make number of intervals
		String[] numOfIntsToMake = null;
		//using cluster
		if (cluster){
			numOfIntsToMake = parameters.getNumberOfIntervals().split(",");
		}
		else {
			numOfIntsToMake = new String[]{parameters.getNumberOfIntervals()};
		}

		//for each number of intervals make intervals for enriched windows and reduced windows
		for (int i=0; i< numOfIntsToMake.length; i++){
			//enriched
			ArrayList all = new ArrayList();
			all.addAll(command);
			//num of ints
			all.add(numOfIntsToMake[i]);
			//file
			all.add("-f");
			all.add(enrichedWindows.toString());
			launchJQSub(all, parameters.getMaxMemory());
			//reduced
			if (reducedWindows !=null && reducedWindows.exists()) {
				all.clear();
				all.addAll(command);
				//num of ints
				all.add(numOfIntsToMake[i]);
				//file
				all.add("-f");
				all.add(reducedWindows.toString());
				launchJQSub(all, parameters.getMaxMemory());
			}
		}
		//wait for results
		int timer = 120;	//minutes 
		System.out.print("\t# to make ");
		while (timer != 0){
			int numToMake = numberOfIntervalsToMake();
			System.out.print(numToMake+" ");
			if (numToMake == 0) break;
			timer --;
			Misc.sleep(60);
		}
		System.out.println();
		//did it time out?
		if (timer == 0) Misc.printExit("\nError: timed out on making intervals.\n");
		//make interval folder and move all intervals into it
		else{
			makeAndMoveIntervals();
		}
	}

	/**Creates set number of intervals.*/
	public void loadIntervalArrays(){
		//make general command
		ArrayList command = new ArrayList();
		command.add(parameters.getAppsDirectory()+"/LoadChipSetIntervalOligoInfo");
		if (parameters.getGenomeFastaDirectory() != null){
			command.add("-s");
			command.add(parameters.getGenomeFastaDirectory().toString());
		}
		command.add("-o");
		command.add(oligoPositions.toString());
		command.add("-i");

		File[] intervalFiles = IO.extractOnlyFiles(intervalDirectory);
		int numToLoad = intervalFiles.length;
		String treatmentDirs = IO.concatinateFileFullPathNames(treatmentReplicas,",");
		String controlDirs = IO.concatinateFileFullPathNames(controlReplicas,",");

		//for each number of intervals make intervals for enriched windows and reduced windows
		for (int i=0; i< numToLoad; i++){
			//enriched
			ArrayList all = new ArrayList();
			all.addAll(command);
			//interval file
			all.add(intervalFiles[i].toString());

			String t = "-t";
			String c = "-c";
			//reduced
			if (intervalFiles[i].toString().startsWith("rd")){
				t = "-c";
				c = "-t";
			}
			//treatmentString
			all.add(t);
			all.add(treatmentDirs);
			//controlString
			all.add(c);
			all.add(controlDirs);
			launchJQSub(all, parameters.getMaxMemory());
		}
		//wait for results
		int timer = 240;	//minutes 
		System.out.print("\t# to load ");
		while (timer != 0){
			int num = numToLoad - IO.numberFilesExist(intervalDirectory, "Ld");
			System.out.print(num+" ");
			if (num == 0) break;
			timer --;
			Misc.sleep(60);
		}
		System.out.println();
		//did it time out?
		if (timer == 0) Misc.printExit("\nError: timed out on loading intervals.\n");
		//delete files not ending in Ld, then rename without Ld
		else{
			IO.deleteFilesNotEndingInExtension(intervalDirectory, "Ld");
			IO.removeExtension(intervalDirectory, "Ld");
		}
	}


	/**Creates set number of intervals.*/
	public void findSubBindingRegions(){
		//make general command
		ArrayList command = new ArrayList();
		command.add(parameters.getAppsDirectory()+"/FindSubBindingRegions");
		command.add("-w");
		command.add(parameters.getSubWindowSize()+"");
		command.add("-n");
		command.add(parameters.getMinNumberOligosInSubWin()+"");
		command.add("-s");
		command.add(parameters.getPeakPickerWindowSize()+"");
		command.add("-m");
		command.add(parameters.getMaxNumberPeaks()+"");
		command.add("-i");

		File[] intervalFiles = IO.extractOnlyFiles(intervalDirectory);
		int numToSub = intervalFiles.length;

		//for each number of intervals make intervals for enriched windows and reduced windows
		for (int i=0; i< numToSub; i++){
			//enriched
			ArrayList all = new ArrayList();
			all.addAll(command);
			//interval file
			all.add(intervalFiles[i].toString());
			launchJQSub(all, parameters.getMaxMemory());
		}

		//wait for results
		int timer = 120;	//minutes 
		System.out.print("\t# to pick peaks ");
		while (timer != 0){
			int num = numToSub - IO.numberFilesExist(intervalDirectory, "Sub");
			System.out.print(num+" ");
			if (num == 0) break;
			timer --;
			Misc.sleep(60);
		}
		System.out.println();
		//did it time out?
		if (timer == 0) Misc.printExit("\nError: timed out picking interval peaks.\n");
		//delete files not ending in Sub, then rename without Sub
		else{
			IO.deleteFilesNotEndingInExtension(intervalDirectory, "Sub");
			IO.removeExtension(intervalDirectory, "Sub");
		}
	}


	/**Prints interval reports in a spread sheet via IntervalReportPrinter.*/
	public void printSpreadSheets(){
		//make general command
		ArrayList command = new ArrayList();
		command.add(parameters.getAppsDirectory()+"/IntervalReportPrinter");
		command.add("-i");
		command.add("1");
		command.add("-a");
		command.add("-b"); 
		command.add("-f");

		File[] intervalFiles = IO.extractFiles(intervalDirectory);
		int numToProc = intervalFiles.length;

		//for each number of intervals make intervals for enriched windows and reduced windows
		for (int i=0; i< numToProc; i++){
			//enriched
			ArrayList all = new ArrayList();
			all.addAll(command);
			//interval file
			all.add(intervalFiles[i].toString());
			launchJQSub(all, parameters.getMaxMemory());
		}

		//wait for results
		int timer = 60;	//minutes 
		System.out.print("\t# to print ");
		while (timer != 0){
			int num = numToProc - IO.numberFilesExist(intervalDirectory, "xls");
			System.out.print(num+" ");
			if (num == 0) break;
			timer --;
			Misc.sleep(20);
		}
		System.out.println();
		//did it time out?
		if (timer == 0) Misc.printExit("\nError: timed out printing Intervals spread sheet reports.\n");
		//make directory and move
		else{
			File spreadSheetDir = new File (intervalDirectory, "SpreadSheetReports");
			if (spreadSheetDir.exists() == false) spreadSheetDir.mkdir();
			File[] res = IO.extractFiles(intervalDirectory, ".xls");
			for (int i=0; i< res.length; i++){
				File moved = new File (spreadSheetDir, res[i].getName());
				res[i].renameTo(moved);
			}
		}
	}

	/**Prints interval bed and sgr files via IntervalGraphPrinter.*/
	public void printIntervalGraphs(){
		//make general command
		ArrayList command = new ArrayList();
		command.add(parameters.getAppsDirectory()+"/IntervalGraphPrinter");
		command.add("-f");

		File[] intervalFiles = IO.extractOnlyFiles(intervalDirectory);
		int numToProc = intervalFiles.length;

		//for each interval file
		for (int i=0; i< numToProc; i++){
			//enriched
			ArrayList all = new ArrayList();
			all.addAll(command);
			//interval file
			all.add(intervalFiles[i].toString());
			launchJQSub(all, parameters.getMaxMemory());
		}

		//wait for results
		int timer = 60;	//minutes 
		System.out.print("\t# to graph ");
		while (timer != 0){
			int num = numToProc - IO.numberFilesExist(intervalDirectory, ".bed");
			System.out.print(num+" ");
			if (num == 0) break;
			timer --;
			Misc.sleep(20);
		}
		System.out.println();
		//did it time out?
		if (timer == 0) Misc.printExit("\nError: timed out printing Intervals graphs.\n");
		//make directories and move
		else{
			//sgr files
			File dir = new File (intervalDirectory, "Sgr");
			if (dir.exists() == false) dir.mkdir();
			File[] res = IO.extractFiles(intervalDirectory, ".sgr.zip");
			for (int i=0; i< res.length; i++){
				File moved = new File (dir, res[i].getName());
				res[i].renameTo(moved);
			}
			//bed files
			dir = new File (intervalDirectory, "Bed");
			if (dir.exists() == false) dir.mkdir();
			res = IO.extractFiles(intervalDirectory, ".bed");
			for (int i=0; i< res.length; i++){
				File moved = new File (dir, res[i].getName());
				res[i].renameTo(moved);
			}
		}
	}


	/**Creates set number of intervals.*/
	public void filterIntervals(){
		//make general command
		ArrayList command = new ArrayList();
		command.add(parameters.getAppsDirectory()+"/IntervalFilter");
		command.add("-i");
		command.add("1");
		command.add("-m");
		command.add("1.5"); //don't want to remove any only flag
		command.add("-e");
		command.add(IO.concatinateFileFullPathNames(parameters.getRepeatRegionFiles(),","));
		command.add("-k");

		File[] intervalFiles = IO.extractFiles(intervalDirectory);
		int numToProc = intervalFiles.length;

		//for each number of intervals make intervals for enriched windows and reduced windows
		for (int i=0; i< numToProc; i++){
			//enriched
			ArrayList all = new ArrayList();
			all.addAll(command);
			//interval file
			all.add(intervalFiles[i].toString());
			launchJQSub(all, parameters.getMaxMemory());
		}

		//wait for results
		int timer = 600;	//minutes 
		System.out.print("\t# to filter ");
		while (timer != 0){
			int num = numToProc - IO.numberFilesExist(intervalDirectory, "Good");
			System.out.print(num+" ");
			if (num == 0) break;
			timer --;
			Misc.sleep(60);
		}
		System.out.println();
		//did it time out?
		if (timer == 0) Misc.printExit("\nError: timed out filtering intervals.\n");
		//delete files not ending in Sub, then rename without Sub
		else{
			IO.deleteFilesNotEndingInExtension(intervalDirectory, "Good");
			IO.removeExtension(intervalDirectory, "Good");
		}
	}


	/**Creats one or two window merged window arrays containing all of the chromosomes.*/
	public void mergeWindows(){
		if (windowDirectory.exists() == false) Misc.printExit("\nError: cannot merge windows, the window directory does not exist. "+windowDirectory+"\n");
		//make general command
		ArrayList command = new ArrayList();
		command.add(parameters.getAppsDirectory()+"/MergeWindowArrays");
		command.add("-r");
		command.add("-t");
		command.add("0.1");	//cutoff
		command.add("-i");
		command.add("1");
		command.add("-f");
		command.add(windowDirectory.toString());
		command.add("-n");
		//files
		enrichedWindows = new File (windowDirectory, "enrichedWindows");
		if (enrichedWindows.exists() == false){
			ArrayList all = new ArrayList();
			all.addAll(command);
			all.add(enrichedWindows.toString());
			launchJQSub(all, parameters.getMaxMemory());
		}
		if (parameters.isMakeReducedWindows()) {
			reducedWindows = new File (windowDirectory, "reducedWindows");
			ArrayList all = new ArrayList();
			all.addAll(command);
			all.add(reducedWindows.toString());
			all.add("-m");
			launchJQSub(all, parameters.getMaxMemory());
		}
		//wait for results
		int timer = 60;	//minutes 
		System.out.print("\t# to merge ");
		while (timer != 0){
			int numToMerge = 0;
			if (enrichedWindows.exists() == false) numToMerge++;
			if (parameters.isMakeReducedWindows() && reducedWindows.exists() == false) numToMerge++;
			System.out.print(numToMerge+" ");
			if (numToMerge == 0) break;
			timer --;
			Misc.sleep(60);
		}
		System.out.println();
		//did it time out?
		if (timer == 0) Misc.printExit("\nError: timed out on merging split window files.\n");
	}


	/**Pulls and places replica directories to this.*/
	public void splitChromSetArray(){
		ArrayList t = new ArrayList();
		ArrayList c = new ArrayList();
		for (int i=0; i< replicas.length; i++){
			if (replicas[i].isTreatmentReplica()) t.add(replicas[i].getSaveDirectory());
			else c.add(replicas[i].getSaveDirectory());
		}
		treatmentReplicas = new File[t.size()];
		controlReplicas = new File[c.size()];
		t.toArray(treatmentReplicas);
		c.toArray(controlReplicas);
	}

	public void scanChromosomes(){
		//make check window directory
		windowDirectory = new File (parameters.getResultsDirectory(), "Win");
		if (windowDirectory.exists() == false) windowDirectory.mkdir();

		File[] windowArrays;
		String[] chromsToProcess = null;

		if (cluster){
			//make list of chromosomes to process individually but combine small ones
			String[] chroms = findChromPairs();

			//make File pointers to final window array results, also collect chromosomes for those that don't exist
			windowArrays = new File[chroms.length];
			ArrayList chromsToProcessAL = new ArrayList();
			for (int i=0; i< chroms.length; i++){
				String concat = chroms[i].replaceAll(",","_").replaceAll("chr","");
				windowArrays[i] = new File(windowDirectory,concat+"_Win");
				if (windowArrays[i].exists() == false) chromsToProcessAL.add(chroms[i]);
			}
			chromsToProcess = Misc.stringArrayListToStringArray(chromsToProcessAL);
		}
		else {
			windowArrays = new File[1];
			windowArrays[0] = new File(windowDirectory, "all_Win");
		}

		//do window array already exist?
		int numToScan = windowArrays.length - IO.numberFilesExist(windowArrays);
		if (numToScan ==0){
			System.out.println("\tAll chromosomes scanned. Delete "+windowDirectory+" to rescan.");
			return;
		}

		//make general command
		ArrayList command = new ArrayList();
		//program text
		command.add(parameters.getAppsDirectory()+"/ScanChromosomes");
		//command.add(parameters.getAppsDirectory()+"/ScanChromosomesCNV");	
		//System.out.println("!!!!!!!!CNV analysis!!!!!!");

		//oligo positions
		command.add("-o");
		command.add(oligoPositions.toString());
		//results directory
		command.add("-r");
		command.add(parameters.getResultsDirectory().toString());
		//treatments
		command.add("-t");
		command.add(IO.concatinateFileFullPathNames(treatmentReplicas,","));
		//controls
		command.add("-c");
		command.add(IO.concatinateFileFullPathNames(controlReplicas,","));
		//genome version
		command.add("-v");
		command.add(parameters.getGenomeVersion());
		//strand
		command.add("-k");
		command.add(parameters.getStrand());
		//size oligo
		command.add("-z");
		command.add(parameters.getSizeOfOligo()+"");
		//window size
		command.add("-w");
		command.add(parameters.getWindowSize()+"");
		//minimum number oligos
		command.add("-m");
		command.add(parameters.getMinimumNumberOfOligos()+"");
		//pseudo median?
		if (parameters.isUsePseudoMedian())command.add("-a");
		//random permutation?
		if (parameters.isUseRandomPermutation()){
			//label or position permutation
			if (parameters.getPermutationType().startsWith("l")) command.add("-l");
			else command.add("-u");
			command.add("-n");
			command.add(parameters.getNumberOfRandomPermutations()+"");
		}
		//randomize data!
		if (parameters.isRandomizeData()) command.add("-x");
		//qValues?
		if (parameters.isUseSymmetricNull()){
			command.add("-p");
			command.add("-s");
			command.add(parameters.getSymmetricNullApp().toString());
			command.add("-q");
			command.add(parameters.getRApp().toString());
		}
		//print files
		//-i -e -d -j
		command.add("-j");	//print log2 ratios instead of relative difference ratio in bar files
		command.add("-i");
		command.add("-e");
		command.add("-d");
		if (cluster){
			//particular chromosomes
			command.add("-f");
			//for each chromosome
			for (int i=0; i< chromsToProcess.length; i++){
				ArrayList all = new ArrayList();
				all.addAll(command);
				all.add(chromsToProcess[i]);
				//launch on 32bit, symP and qVal not on 64bit
				launchJQSub(all, parameters.getMaxMemory());
			}
		}
		else launchJQSub(command, parameters.getMaxMemory());

		//wait for results
		int timer = 600;	//minutes
		numToScan = windowArrays.length - IO.numberFilesExist(windowArrays);
		System.out.print("\t# to scan "+numToScan);
		while (timer != 0){
			if (numToScan == 0) break;
			timer --;
			Misc.sleep(60);
			numToScan = windowArrays.length - IO.numberFilesExist(windowArrays);
			System.out.print(" "+numToScan);
		}
		System.out.println();
		//did it time out?
		if (timer == 0) Misc.printExit("\nError: timed out on scanning chromosomes.\n");


	}

	/**Attempts to join smallest chromosomes together so long as the combine file doesn't exceed the size of the biggest.*/
	public String[] findChromPairs(){
		//make FileSize[]
		File[] files = IO.extractFiles(oligoPositions);
		FileSize[] fs = new FileSize[files.length];
		for (int i=0; i< fs.length; i++) fs[i] = new FileSize(files[i]);
		//sort 
		Arrays.sort(fs);
		//join smallest till size just under biggest
		long biggest = fs[fs.length-1].getSize();
		long current = 0;
		int index = 0;
		ArrayList al = new ArrayList();
		for (int i=0; i< fs.length; i++){
			long testSize = fs[i].getSize();
			if ((current+ testSize) <= biggest){
				al.add(fs[i].getName());
				current+= fs[i].getSize();
				index = i;
			}
			else break;
		}
		index++;
		String fusion = Misc.stringArrayListToString(al, ",");
		int sizeArray = fs.length - index +1;
		String[] chromList = new String[sizeArray];
		chromList[0] = fusion;
		int counter = 1;
		for (int i = index; i< fs.length; i++){
			chromList[counter++] = fs[i].getName();
		}
		return chromList;
	}	

	public void deleteSplitNormalizedFiles(){
		for (int i=0; i< replicas.length; i++){
			File[] f = replicas[i].getChipSetDirectories();
			for (int j=0; j< f.length; j++){
				IO.deleteDirectory(f[j]);
			}
		}
	}

	public void makeChromSet (File[][] setRep, boolean treatmentSet){
		//make check PCels directory
		File pCels = new File (parameters.getResultsDirectory(),"PCels");
		if (pCels.exists() == false) pCels.mkdir();
		//for each replica
		for (int r=0; r<setRep[0].length; r++ ){
			File saveDir;
			if (treatmentSet){
				saveDir = new File (pCels, "TreatmentReplica"+(r+1));
				replicas[replicaIndex].setTreatmentReplica(true);
			}
			else{
				saveDir = new File (pCels, "ControlReplica"+(r+1));
				replicas[replicaIndex].setTreatmentReplica(false);
			}
			replicas[replicaIndex].setSaveDirectory(saveDir);

			//process, collect entire chipset for each replica
			File[] set = new File[setRep.length];
			for (int j=0; j< set.length; j++) {
				String name = Misc.removeExtension(setRep[j][r].getName());
				File test = new File (setRep[j][r].getParentFile(), name);
				if (test.exists() == false) Misc.printExit("Missing normalized file! "+test);
				set[j] = test;
			}
			replicas[replicaIndex].setChipSetDirectories(set);

			replicaIndex++;
		}
	}

	public int numberReplicasToMake(){
		int num =0;
		for (int i=0; i< replicas.length; i++) {
			replicas[i].loadNumberOligos();
			if (replicas[i].getNumberOligos() ==0) num++;
		}
		return num;
	}

	public void makeChromosomeSets(){
		int numToProcess = parameters.getTreatmentCelFiles()[0].length + parameters.getControlCelFiles()[0].length;
		replicas = new ChromSet[numToProcess];
		for (int i=0; i< numToProcess; i++) replicas[i] = new ChromSet();
		replicaIndex = 0;
		makeChromSet(parameters.getTreatmentCelFiles(), true);
		makeChromSet(parameters.getControlCelFiles(), false);

		//launch maker
		int numToMake = numberReplicasToMake();	
		if (numToMake !=0){
			for (int i=0; i< replicas.length; i++){
				if (replicas[i].getNumberOligos() !=0) continue;
				//make command
				ArrayList command = new ArrayList();
				//program text
				command.add(parameters.getAppsDirectory()+"/MakeChromosomeSets");
				//files to join
				command.add("-d");
				command.add(IO.concatinateFileFullPathNames(replicas[i].getChipSetDirectories(), ","));
				//save directory
				command.add("-n");
				command.add(replicas[i].getSaveDirectory().toString());
				//skip making oligo positions for all but first one
				if (i !=0) command.add("-s");
				//launch
				launchJQSub(command, parameters.getMaxMemory());
			}
		}
		//wait for results
		int timer = 60;	//minutes
		numToMake = numberReplicasToMake();
		System.out.print("\t# to merge "+numToMake);
		while (timer != 0){
			if (numToMake == 0) break;
			timer --;
			Misc.sleep(60);
			numToMake = numberReplicasToMake();
			System.out.print(" "+numToMake);
		}
		System.out.println();
		//did it time out?
		if (timer == 0) Misc.printExit("\nError: timed out on merging split normalized files.\n");

		//check sizes
		long num =0;
		for (int i=0; i< replicas.length; i++) {
			replicas[i].loadNumberOligos();
			if (num == 0) num = replicas[i].getNumberOligos();
			else if (num != replicas[i].getNumberOligos()) Misc.printExit("\nError: the number of oligos in one of your replicas differ! "+replicas[i].getSaveDirectory()+"\n");
		}
		//move OligoPositions down to base directory
		oligoPositions = new File (parameters.getResultsDirectory(),"OligoPositions");
		File oldOP = new File (replicas[0].getSaveDirectory(),"OligoPositions");
		oldOP.renameTo(oligoPositions);
	}

	public void normalizeCelFiles(){
		//make directory to hold qc
		File qcDir = new File (parameters.getResultsDirectory(), "QC");
		if (qcDir.exists() == false) qcDir.mkdir();

		//for each chip set
		int numChipSets = parameters.getTreatmentCelFiles().length;
		File[][] normalizedFiles = new File[numChipSets][];
		for (int k=0; k< numChipSets; k++){
			//collect all files
			File[] t = parameters.getTreatmentCelFiles()[k];
			File[] c = parameters.getControlCelFiles()[k];
			File[] combine = new File[t.length + c.length];
			int counter =0;
			for (int i=0; i< t.length; i++) combine[counter++] = t[i];
			for (int i=0; i< c.length; i++) combine[counter++] = c[i];

			//check if entire set has been normalized
			boolean allNormalized = true;
			normalizedFiles[k] = new File[combine.length];
			for (int i=0; i< combine.length; i++){
				String name = Misc.removeExtension(combine[i].getName());
				File test = new File (t[0].getParentFile(), name);
				if (test.exists() == false) allNormalized = false;
				normalizedFiles[k][i] = test;
			}
			//normalize
			if (allNormalized == false){
				//make general command
				ArrayList command = new ArrayList();
				//program text
				command.add(parameters.getAppsDirectory()+"/CelProcessor");
				//use control seqs in setting median?
				if (parameters.isMedianScaleUsingControlSeqs()) command.add("-c");
				//tpmap
				command.add("-t");
				command.add(parameters.getTpmapFiles()[k].toString());
				//median target
				command.add("-m");
				command.add(parameters.getTargetMedian()+"");
				//break by chrom
				command.add("-r");
//skip quantile?
//System.out.println("******** WARNING: skipping quantile normalization ********* \n");
//command.add("-q");
				//if combine
				if (parameters.isQuantileNormalizeAll()){
					File clusterName = new File (qcDir, parameters.getResultsName()+"ChipSubSet"+(1+k)+"ClusterPlot.png");
					launchCelNormalization(combine, command, clusterName);
				}
				//seperate treatment and control
				else {
					//treatment
					File clusterName = new File (qcDir, parameters.getResultsName()+"TreatmentChipSubSet"+(1+k)+"ClusterPlot.png");
					launchCelNormalization(t, command, clusterName);
					//control
					clusterName = new File (qcDir, parameters.getResultsName()+"ControlChipSubSet"+(1+k)+"ClusterPlot.png");
					launchCelNormalization(c, command, clusterName);
				}
			}
		}

		//check if they are complete and wait
		int timer = 200;	//minutes
		System.out.print("\t# to normalize ");
		File[] expectedNormed = IO.collapseFileArray(normalizedFiles);
		while (timer != 0){
			timer --;
			int numFound =0;
			for (int i=0; i< expectedNormed.length; i++){
				if (expectedNormed[i].exists()){
					File[] files = IO.extractFiles(expectedNormed[i]);
					if (files.length > 0) numFound++;
				}
			}
			int num = expectedNormed.length - numFound;
			System.out.print(num+" ");
			if (num == 0) break;
			Misc.sleep(60);
		}
		System.out.println();
		//did it time out?
		if (timer == 0) Misc.printExit("\nError: timed out on normalizing xxx.cela files.\n");
	}

	public void launchCelNormalization(File[] celaFiles, ArrayList generalCmd, File clusterName){
		//make ArrayList to hold additional parameters
		ArrayList finalCmd = new ArrayList();
		finalCmd.addAll(generalCmd);
		//files to norm
		finalCmd.add("-d");
		finalCmd.add(IO.concatinateFileFullPathNames(celaFiles, ","));
		//cluster
		//System.out.println("********* Skipping clustering ************");		
		finalCmd.add("-i");
		finalCmd.add(clusterName.toString());
		//launch
		launchJQSub(finalCmd, parameters.getMaxMemory());
	}

	public void convertCelFiles(){
		//make ACels directory
		File aCels = new File (parameters.getResultsDirectory(), "ACels");
		if (aCels.exists() == false) aCels.mkdir();

		//collect txt cel files
		File[] t = IO.collapseFileArray(parameters.getTreatmentCelFiles());
		File[] c = IO.collapseFileArray(parameters.getControlCelFiles());
		File[] combine = new File[t.length+c.length];
		for (int i=0; i< t.length; i++) combine[i] = t[i];
		int counter =0;
		for (int i=t.length; i< combine.length; i++) combine[i] = c[counter++];

		//find those not already converted
		ArrayList toConvert = new ArrayList();
		Pattern p = Pattern.compile("\\.cel$", Pattern.CASE_INSENSITIVE);
		for (int i=0; i< combine.length; i++){
			Matcher mat = p.matcher(combine[i].getName());
			String name = mat.replaceFirst(".cela");
			File test = new File (aCels, name);
			if (test.exists() == false) toConvert.add(combine[i]);
		}

		//convert!
		int numToConvert= toConvert.size();
		for (int i=0; i< numToConvert; i++){
			ArrayList command = new ArrayList();
			//program text
			command.add(parameters.getAppsDirectory()+"/CelFileConverter");
			//params
			command.add("-f");
			//file text
			command.add(((File)toConvert.get(i)).toString());
			//save directory
			command.add("-s");
			command.add(aCels.toString());
			//launch
			launchJQSub(command, parameters.getMaxMemory());
		}

		//check number of files converted
		int timer = 100;
		System.out.print("\t# to convert ");
		File[] converted = null;
		while (timer != 0){
			timer --;
			converted = IO.extractFiles(aCels, "cela");
			int numConverted = 0;
			if (converted != null) numConverted = converted.length;
			int num = combine.length - numConverted;
			System.out.print(num+" ");
			if (num == 0) break;
			Misc.sleep(60);
		}
		System.out.println();
		//did it time out?
		if (timer == 0) Misc.printExit("\nError: timed out on converting xxx.cel files.\n");

		//check size, assumes first file is correct size.
		long size = converted[0].length();
		for (int i=1; i< converted.length; i++) {
			if (converted[i].length() != size) Misc.printExit("\nError: one of your converted xxx.cela files differs in size -> "+converted[i]+" delete and restart.\n");
		}
	}

	public void launchJQSub(ArrayList args, int memory){
		args.add(0, "-jar");
		args.add(0, "-Xmx"+memory+"M");
		args.add(0, parameters.getJava().toString());
		if (cluster) {
			args.add(0, parameters.getAppsDirectory()+"/JQSub");
			args.add(0, "-jar");
			args.add(0, parameters.getJava().toString());

		}
		String[] cmd = Misc.stringArrayListToStringArray(args);
		System.out.println("\tLaunching: "+Misc.stringArrayToString(cmd," "));
		String[] results = IO.executeCommandLineReturnAll(cmd);
		System.out.println("\t"+Misc.stringArrayToString(results,"\n"));
	}





	public static void main(String[] args) {
		if (args.length == 0){
			printDocs();
			System.exit(0);
		}	
		new T2(args);
	}

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': parameterFile = new File(args[i+1]); i++; break;
					case 'c': cluster = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nError: unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		if (parameterFile == null || parameterFile.exists() == false) Misc.printExit("\nCannot find or read your parameter file '"+parameterFile+"'. See the t2ParamFileTemplate.xls for an example.\n");
	}	

	public static void printDocs(){ 
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                T2: April 2007                                    **\n" +
				"**************************************************************************************\n" +
				"T2 launches many of the TiMAT2 applications based on a tab delimited parameter file\n" +
				"(see T2_xxx/Documentation/t2ParamFileTemplate.xls). It converts, normalizes, splits,\n" +
				"and merges txt cel files and then launches ScanChromosomes.  If you wish to make use\n" +
				"of a cluster for parallel processing, configure and test the JQSub TiMAT2 app.\n\n" +

				"-p Full path file text for the tab delimited parameter file.\n"+
				"-c Use the cluster(s) specified by JQSub.\n\n"+

				"Example: java pathTo/T2/Apps/T2 -c -p /my/t2ParamFile.txt\n\n" +

		"**************************************************************************************\n");		
	}


}
