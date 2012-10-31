package edu.utah.seq.data;

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.*;
import util.gen.*;


public class Qseq2Fastq {
	
	//fields
	private File qseqDirectory;
	private File fastqDirectory;
	private boolean printFullHeaders = false;
	private boolean deleteQseqFiles = false;
	private boolean verbose = true;
	private boolean pairsPresent = false;
	private boolean keepAllReads = false;
	private boolean failedParsing = true;
	private LinkedHashMap<String, ArrayList<File>> qseqFiles = new LinkedHashMap<String, ArrayList<File>>();
	private QseqParser[] parsers;

	//constructor
	public Qseq2Fastq(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		if (verbose) System.out.println("Launched...\n");
		
		processArgs(args);
		
		//fetch files
		extractQseqFiles();
		
		//parse files
		parseFiles();
		
		//wait till complete
		monitorThreads();
		
		//check results, be sure to call before deletingQseqFiles!
		checkResults();
		
		//print per lane parsing stats
		if (verbose) printReadStats();
		
		//delete qseq files?
		if (deleteQseqFiles && failedParsing == false) deleteQseqFiles();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/(1000*60);
		if (verbose) System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}

	//methods
	public void deleteQseqFiles(){
		Iterator<String> it = qseqFiles.keySet().iterator();
		while (it.hasNext()){
			ArrayList<File> al = qseqFiles.get(it.next());
			int num = al.size();
			for (int i=0; i< num; i++) al.get(i).deleteOnExit();
		}
	}
	
	public void printReadStats(){
		System.out.println("Lane parsing statistics:\n");
		String p = "";
		if (pairsPresent) p= " paired";
		for (int i=0; i< parsers.length; i++) {
			System.out.println("\tLane "+parsers[i].getSaveDirectory().getName());
			System.out.println("\t\t"+parsers[i].getTotalNumberReads()+"\tTotal number"+p+" reads");
			if (keepAllReads == false) {
				double total = parsers[i].getTotalNumberReads();
				double fractPass = (double)parsers[i].getNumberPassingReads()/ total;
				double fractFail = (double)parsers[i].getNumberFailedReads()/ total;
				System.out.println("\t\t"+parsers[i].getNumberPassingReads()+"\tNumber passing ("+Num.formatPercentOneFraction(fractPass)+")");
				System.out.println("\t\t"+parsers[i].getNumberFailedReads()+"\tNumber failed ("+Num.formatPercentOneFraction(fractFail)+")");
			}
			System.out.println("\t\t"+parsers[i].getNumberMalformedQseqLines()+"\tNumber malformed qseq lines");
		}
	}
	
	public void checkResults(){
		failedParsing = false;
		for (int i=0; i< parsers.length; i++) {
			if (parsers[i].isPassed() == false) {
				failedParsing = true;
				break;
			}
		}
		if (failedParsing) Misc.printErrAndExit("Error: failed to convert!");
	}
	
	
	public void monitorThreads(){
		while (true){
			//wait 10 secs
			try {
				Thread.sleep(10000L);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			//check threads if complete
			int numComplete = 0;
			for (int i=0; i< parsers.length; i++) {
				if (parsers[i].isComplete()) numComplete++;
			}
			if (numComplete == parsers.length) break;
		}
	}
	
	public void parseFiles(){
		try{
	        //for each lane, parse in a new thread
	        Iterator<String> it = qseqFiles.keySet().iterator();
	        int numberLanes = qseqFiles.size();
	        parsers = new QseqParser[numberLanes];
	        for (int i=0; i< numberLanes; i++){
	        	String lane = it.next();
	        	File laneDir = new File (fastqDirectory, lane);
	        	laneDir.mkdir();
	        	parsers[i] = new QseqParser(laneDir, qseqFiles.get(lane), printFullHeaders, pairsPresent, keepAllReads);
	        	parsers[i].start();
	        }
	        
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing single qseq files.\n");
		}
	}
	
	/**Fetches qseq files and performs a variety of checks.*/
	public void extractQseqFilesDirectlyNOTWorking(){
		//maps to hold files by read number (1st, 2nd, 3rd)
		HashMap<String, File> readOne = new HashMap<String, File>();
		HashMap<String, File> readTwo = new HashMap<String, File>();
		HashMap<String, File> readThree = new HashMap<String, File>();
		
		HashMap<String, QseqLane> lanes = new HashMap<String, QseqLane>();
		HashSet<String> tiles = new HashSet<String>();
		
		//fetch files
		File[][] tot = new File[6][];
		tot[0] = IO.extractFiles(qseqDirectory,".qseq");
		tot[1] = IO.extractFiles(qseqDirectory,".qseq.gz");
		tot[2] = IO.extractFiles(qseqDirectory,".qseq.zip");
		tot[3] = IO.extractFiles(qseqDirectory,".qseq.txt");
		tot[4] = IO.extractFiles(qseqDirectory,".qseq.txt.gz");
		tot[5] = IO.extractFiles(qseqDirectory,".qseq.txt.zip");
		File[] dataFiles = IO.collapseFileArray(tot);
		
		//check names
		Pattern fileName = Pattern.compile("s_(\\d)_(\\d)_(\\d+)_qseq.*");
		LinkedHashMap<String,File> namesHM = new LinkedHashMap<String,File>();
		for (int i=0; i< dataFiles.length; i++){
			Matcher mat = fileName.matcher(dataFiles[i].getName());
			if (mat.matches() == false) {
				Misc.printErrAndExit("\nError: the following qseq file did not fit the Illumina naming convention (e.g. s_5_1_0062_qseq...)\n\n\t-> "+dataFiles[i].getName()+"\n");
			}
			//laneNumber readNumber tileNumber
			String laneNumber = mat.group(1);
			String readNumber = mat.group(2);
			String tileNumber = mat.group(3);
			
			//namesHM.put(mat.group(1)+"_"+mat.group(2)+"_"+mat.group(3), dataFiles[i]);
			
			//laneNumber, tileNumber
			if (readNumber.equals("1")) readOne.put(laneNumber+"_"+tileNumber, dataFiles[i]);
			else if (readNumber.equals("2")) readTwo.put(laneNumber+"_"+tileNumber, dataFiles[i]);
			else if (readNumber.equals("3")) readThree.put(laneNumber+"_"+tileNumber, dataFiles[i]);
			else Misc.printErrAndExit("\nError: found a qseq file with an unexpected read number (not 1, 2, or 3)\n\t-> "+dataFiles[i].getName()+"\n");
			
			if (lanes.containsKey(laneNumber) == false) lanes.put(laneNumber,new QseqLane(laneNumber));
			tiles.add(tileNumber);
		}
		
		//check numbers
		if (readThree.size() !=0){
			if (readThree.size() != readTwo.size() || readTwo.size() != readOne.size()) Misc.printErrAndExit("\nError: found an unequal number of qseq files for each of the three reads.\n");
		}
		else if (readTwo.size() !=0){
			if (readTwo.size() != readOne.size()) Misc.printErrAndExit("\nError: found an unequal number of qseq files for each of the two reads.\n");
		}
		else if (readOne.size() == 0) Misc.printErrAndExit("\nError: no qseq files were found?\n");
		
		//build QseqFileSet and put in QseqLane
		
		
		
		//look for qualified files
		fetchQualifiedFiles(namesHM);
	}	
	
	/**Fetches paired or unpaired checked */
	public void extractQseqFiles(){
		
		//fetch files
		File[][] tot = new File[6][];
		tot[0] = IO.extractFiles(qseqDirectory,".qseq");
		tot[1] = IO.extractFiles(qseqDirectory,".qseq.gz");
		tot[2] = IO.extractFiles(qseqDirectory,".qseq.zip");
		tot[3] = IO.extractFiles(qseqDirectory,".qseq.txt");
		tot[4] = IO.extractFiles(qseqDirectory,".qseq.txt.gz");
		tot[5] = IO.extractFiles(qseqDirectory,".qseq.txt.zip");
		File[] dataFiles = IO.collapseFileArray(tot);
		
		//check names
		Pattern fileName = Pattern.compile("s_(\\d)_(\\d)_(\\d+)_qseq.*");
		LinkedHashMap<String,File> namesHM = new LinkedHashMap<String,File>();
		for (int i=0; i< dataFiles.length; i++){
			Matcher mat = fileName.matcher(dataFiles[i].getName());
			if (mat.matches() == false) {
				Misc.printErrAndExit("\nError: the following qseq file did not fit the Illumina naming convention (e.g. s_5_1_0062_qseq...)\n\n\t-> "+dataFiles[i].getName()+"\n");
			}
			//laneNumber readNumber tileNumber
			namesHM.put(mat.group(1)+"_"+mat.group(2)+"_"+mat.group(3), dataFiles[i]);
		}
		
		//look for qualified files
		fetchQualifiedFiles(namesHM);
	}
	
	/**Returns a File[] containing paired 1st then 2nd qseq files or single files.
	 * Sets the pairsPresent boolean
	 * Checks for unpaired pairs and other nonsense.*/
	public void fetchQualifiedFiles(LinkedHashMap<String,File> nameFile){
		Pattern fileName = Pattern.compile("(\\d)_(\\d)_(\\d+)");
		Iterator<String> it = nameFile.keySet().iterator();
		LinkedHashMap<String, ArrayList<File>> singles = new LinkedHashMap<String, ArrayList<File>>();
		LinkedHashMap<String, ArrayList<File>> pairs = new LinkedHashMap<String, ArrayList<File>>();
		ArrayList<File> al;
		int numberPairs =0;
		int numberSingles =0;
		while (it.hasNext()){
			String test = it.next();
			Matcher mat = fileName.matcher(test);
			mat.matches();
			//2nd pair found?
			if (mat.group(2).equals("2")){
				//look for first file
				String firstFileName = mat.group(1)+"_1_"+mat.group(3);
				if (nameFile.containsKey(firstFileName)== false){
					Misc.printErrAndExit("\nError: found an unpaired 2nd pair qseq file.  Cannot locate -> s_"+firstFileName+"_qseq.....\n");
				}
				else {
					//fetch ArrayList based on lane number
					al = pairs.get(mat.group(1));
					if (al == null) {
						al = new ArrayList<File>();
						pairs.put(mat.group(1), al);
					}
					//add first and second
					al.add(nameFile.get(firstFileName));
					al.add(nameFile.get(test));
					numberPairs++;
				}
			}
			else {
				al = singles.get(mat.group(1));
				if (al == null) {
					al = new ArrayList<File>();
					singles.put(mat.group(1), al);
				}
				al.add(nameFile.get(test));
				numberSingles++;
			}
		}

		//any pairs?
		if (numberPairs != 0) {
			//check that all are paired up
			if (numberPairs != numberSingles) {
				Misc.printErrAndExit("\nError: found some paired and some unpaired qseq files?! Aborting.\n");
			}
			pairsPresent = true;
			qseqFiles = pairs;			
		}
		//no pairs
		else {
			pairsPresent = false;
			qseqFiles = singles;
		}
	}
	
	public static void main(String[] args) {
		if (args.length<1){
			printDocs();
			System.exit(0);	
		}
		new Qseq2Fastq(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'q': qseqDirectory = new File(args[++i]); break;
					case 'f': fastqDirectory = new File(args[++i]); break;
					case 'p': printFullHeaders = true; break;
					case 'd': deleteQseqFiles = true; break;
					case 'a': keepAllReads = true; break;
					case 's': verbose = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nError, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nError, something doesn't look right with this parameter request: -"+test);
				}
			}
		}
		
		if (qseqDirectory == null || qseqDirectory.isDirectory() == false || qseqDirectory.canRead() == false){
			Misc.printErrAndExit("\nError: cannot find, read, or identify a directory containing qseq files?! -> "+ qseqDirectory);
		}
		//make the save directory
		if (fastqDirectory == null) fastqDirectory = new File (qseqDirectory, "Qseq2Fastq");
		if (fastqDirectory.exists() == false){
			if (fastqDirectory.mkdirs() == false) Misc.printErrAndExit("\nError: failed to make the fastq directory for saving your converted files?! -> "+fastqDirectory);
		}
		if (fastqDirectory.isDirectory() == false || fastqDirectory.canWrite() == false) Misc.printErrAndExit("\nError: cannot write to or identify a the fastq save directory?! -> "+fastqDirectory);

	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                  Qseq2Fastq: Aug 2010                            **\n" +
				"**************************************************************************************\n" +
				"Parses, filters out reads failing QC, compresses and converts single and paired read\n" +
				"qseq files to Illumina fastq format. Does not concatinate tiles.\n" +

				"\nRequired Parameters:\n"+
				"-q Qseq directory. This should contain all of the qseq files for a sequencing run,\n" +
				"       multiple lanes, paired and single reads. (e.g. s_5_1_0025_qseq.txt(.gz/zip OK))\n"+
				
				"\nOptional Parameters:\n"+
				"-f Fastq save directory. Defaults to the qseq directory.\n"+
				"-a Keep all reads. Defaults to removing those failing the QC flag. Paired reads are\n" +
				"       only removed if both reads fail QC.\n"+
				"-p Print full fastq headers. Defaults to using read count.\n"+
				"-d Delete qseq files upon successfull parsing of all files. Be carefull!\n"+
				"-s Silence non error output.\n"+

				"\n"+
				"Example: java -Xmx2G -jar pathTo/USeq/Apps/Qseq2Fastq -f /Runs/7/Fastq -q\n" +
				"      /Runs/100726_SN141_0265_A207D4ABXX/Data/Intensities/BaseCalls \n\n" +

		"**************************************************************************************\n");		
	}
}
