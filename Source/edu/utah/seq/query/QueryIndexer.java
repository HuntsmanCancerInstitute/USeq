package edu.utah.seq.query;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.useq.data.RegionScoreText;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndex;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class QueryIndexer {

	//user params
	private File dataDir;
	private File chrNameLength;
	private File indexDir;
	private boolean verbose = true;
	private File bgzip = null;
	private File tabix = null;
	private int numberThreads = 0;
	private String[] skipDirs = null;

	//internal fields
	private String email = "bioinformaticscore@utah.edu";
	private File[] dataFilesToParse;
	private String[] fileExtToIndex = {"vcf.gz", "bed.gz", "bedgraph.gz", "maf.txt.gz"};
	private HashSet<File> skippedDataSources = new HashSet<File>();
	private static HashMap<String, int[]> extensionsStartStopSub = new HashMap<String, int[]>();
	private TreeMap<String, Integer> chrLengths;
	private HashMap<String, ArrayList<File>> chrFiles = new HashMap<String, ArrayList<File>>();
	private HashMap<File, Integer> fileId = new HashMap<File, Integer>();
	private Integer[] ids = null;
	private long totalRecordsProcessed = 0;
	
	//prior index
	private int toTruncatePoint = -1;
	private HashMap<String, Integer> priorTruncFileNameId = null;
	private HashMap<String, Long> priorTruncFileNameSize = null;

	//per chr fields
	private ArrayList<IndexRegion>[] workingIndex = null; 
	private int bpBlock = 250000000; //max size hg19 chr1 is 249,250,621
	private ArrayList<File> workingFilesToParse = new ArrayList<File>();
	private ArrayList<File> workingBedFilesToMerge = new ArrayList<File>();
	private int workingFilesToParseIndex = 0;
	private String workingChr = null;
	private int workingStartBp;
	private int workingStopBp;
	private long workingParsed = 0;
	private PrintWriter out = null;

	//constructor
	public QueryIndexer(String[] args) {
		long startTime = System.currentTimeMillis();

		processArgs(args);
		
		loadPrior();
		parseSkipDirs(skipDirs);
		setKnownFileTypeSSS();
		parseChromLengthFile();
		
		parseDataSources();
		
		if (priorTruncFileNameId!= null) contrastPriorWithCurrent();
		else createFileIdHash();
		
		createChrFiles();
		
		createFileIdArray();

		//for each chromosome
		IO.pl("\nIndexing records by chr...\n\tChrBlock\t#Parsed");
		for (String chr: chrLengths.keySet()) {
			workingChr = chr;
			parseChr();
		}
 
		//bgzip and tabix
		IO.pl("\nCompressing master index...");
		bgzipAndTabixIndex();

		IO.pl("\nSaving file objects...");
		saveFileIds();

		String diffTime = Num.formatNumberOneFraction(((double)(System.currentTimeMillis() -startTime))/1000/60);
		IO.pl("\n"+ diffTime+" Min to parse ~"+ totalRecordsProcessed+" records and build the query index");
	}

	private void createFileIdArray() {
		int maxValue = 0;
		for (int i: fileId.values()) if (i>maxValue) maxValue = i;
		ids = new Integer[maxValue+1];
		for (Integer i: fileId.values()) ids[i] = i;
	}

	private void createFileIdHash() {
		IO.pl("\nCreating file id hash...");
		for (int i=0; i< dataFilesToParse.length; i++) {
			Integer id = new Integer(i);
			fileId.put(dataFilesToParse[i], id);
		}
	}
	
	private void contrastPriorWithCurrent() {
		IO.pl("\nComparing prior index with current data sources...");
		
		int numOldDataSources = 0;
		int numNewDataSources = 0;
		int numNewDataSourcesDiffSize = 0;
		int numOldDataSourcesMissingInNew = 0;

		//walk current files to index
		HashSet<String> currentTrimmedDataSourceNames = new HashSet<String>();
		for (int i=0; i< dataFilesToParse.length; i++) {
			
			//add to fileId hash
			Integer id = new Integer(i);
			fileId.put(dataFilesToParse[i], id);
			
			//was it already parsed?
			String trimmedName = dataFilesToParse[i].toString().substring(toTruncatePoint);
			currentTrimmedDataSourceNames.add(trimmedName);
			Integer idTest = priorTruncFileNameId.get(trimmedName);

			if (idTest == null) numNewDataSources++;
			else {
				//same size?
				Long size = priorTruncFileNameSize.get(trimmedName);
				if (dataFilesToParse[i].length() == size) numOldDataSources++;
				else numNewDataSourcesDiffSize++;
			}
		}
		
		//walk old Index and see which are missing in new and should be deleted
		for (String oldTN: priorTruncFileNameId.keySet()) {
			if (currentTrimmedDataSourceNames.contains(oldTN) == false) numOldDataSourcesMissingInNew++;
		}
		
		
		if (verbose) {
			//old datasets already parsed and in index - skip
			IO.pl("\t"+numOldDataSources +" Datasets present in index");
			//entirely new datasets - to parse
			IO.pl("\t"+numNewDataSources +" New datasets to add to index");
			//new datasets with preexisting same named file but diff size, - delete old from index then parse
			IO.pl("\t"+numNewDataSourcesDiffSize +" New datasets to replace an existing data source");
			//old datasets to delete, either getting replace or were deleted from current datasource list
			IO.pl("\t"+numOldDataSourcesMissingInNew+" Datasets missing from index ");
		}
		
		//check if any work to do
		if (numNewDataSources==0 && numNewDataSourcesDiffSize==0 && numOldDataSourcesMissingInNew==0)  Misc.printExit("\nIndex is up to date. Exiting.");	
		else IO.pl("\t\tRebuilding Index");
	}

	private void loadPrior() {
		try {
			
			//look for the fileIds.obj, fileSizes.obj
			File ids = new File(indexDir, "fileIds.obj");
			File sizes = new File(indexDir, "fileSizes.obj");
			
			if (ids.exists() == false || sizes.exists() == false) return;
			
			IO.pl("Loading prior file objects...");
			
			priorTruncFileNameId = (HashMap<String, Integer>)IO.fetchObject(ids);
			priorTruncFileNameSize = (HashMap<String, Long>)IO.fetchObject(sizes);
			File[] priorChromIndexFiles = IO.extractFiles(indexDir, ".bed.gz");
			
			if (priorChromIndexFiles == null || priorChromIndexFiles.length == 0) throw new IOException("\nFailed to find your chrXXX.bed.gz index files in your index directory?");
			if (verbose) {
				IO.pl("\t"+priorTruncFileNameId.size()+"\tIndexed data sources ");
				IO.pl("\t"+priorChromIndexFiles.length+"\tChromosome indexes");
			}

		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: opening prior index objects, aborting\n");
		}
		
	}

	public void parseChr(){
		try {
			workingFilesToParse.clear();
			workingFilesToParse.addAll(chrFiles.get(workingChr));
			workingBedFilesToMerge.clear();

			//remove any that belong to a skip dir
			removeSkipDirs(workingFilesToParse);
			
			//any work to do?
			if (workingFilesToParse.size() == 0 ) {
				IO.pl("\t"+workingChr+"\tNothing to do");
				return;
			}
			
			//start io
			File queryIndexFile = new File(indexDir, workingChr+".qi.bed");
			out = new PrintWriter(new FileWriter((queryIndexFile)));
			
			//for each block
			int chromLength = chrLengths.get(workingChr)+2;
			for (int i=0; i<chromLength; i+=bpBlock) {
				workingStartBp = i;
				workingStopBp = workingStartBp+bpBlock;
				if (workingStopBp > chromLength) workingStopBp = chromLength;
				IO.p("\t"+workingChr+":"+workingStartBp+"-"+workingStopBp);
				
				//prep for new chr seg
				workingIndex = new ArrayList[workingStopBp-workingStartBp+1];
				workingParsed = 0;
				parseChrBlock();
			}

			out.close();

		} catch (IOException e){
			e.printStackTrace();
			Misc.printErrAndExit("\nFATAL error with loading "+workingChr+", aborting.");
		}
	}


	private void parseChrBlock() throws IOException {
		
		//load new datasets?
		if (workingFilesToParse.size() !=0) {
			workingFilesToParseIndex = 0;
			//try to make a loader for each file
			int numToMake= workingFilesToParse.size();
			if (numToMake > numberThreads) numToMake = numberThreads;			
			QueryIndexFileLoader[] loader = new QueryIndexFileLoader[numToMake];			
			ExecutorService executor = Executors.newFixedThreadPool(numToMake);
			for (int i=0; i< loader.length; i++){
				loader[i] = new QueryIndexFileLoader(this, workingChr, workingStartBp, workingStopBp);
				executor.execute(loader[i]);
			}
			executor.shutdown();

			//spins here until the executer is terminated, e.g. all threads complete
			while (!executor.isTerminated()) {}

			//check loaders 
			for (QueryIndexFileLoader c: loader) {
				if (c.isFailed()) throw new IOException("ERROR: File Loader issue! \n"+c);
			}
		}
		

		//save the index for this chrom block
		if (workingParsed !=0) saveWorkingChrBlock();

		//stats
		IO.pl("\t"+workingParsed);
		totalRecordsProcessed+= workingParsed;
	}

	/*Removes any files that are found in a skip dir*/
	private void removeSkipDirs(ArrayList<File> al) {
		try {
			if (skipDirs != null){
				ArrayList<File> bad = new ArrayList<File>();
				//for each file
				for (File f: al){
					String conPath = f.getCanonicalPath();
					//for each skip dir
					for (String sd: skipDirs){
						if (conPath.startsWith(sd)) {
							bad.add(f);
							break;
						}
					}
				}
				al.removeAll(bad);
				skippedDataSources.addAll(bad);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} 
	}
	
	synchronized void addRegions (ArrayList<IndexRegion> regions) {
		for (IndexRegion region: regions){
			//add start
			int sIndex = region.start;
			if (workingIndex[sIndex] == null) workingIndex[sIndex] = new ArrayList<IndexRegion>();
			workingIndex[sIndex].add(region);

			//add stop 
			sIndex = region.stop;
			if (workingIndex[sIndex] == null) workingIndex[sIndex] = new ArrayList<IndexRegion>();
			workingIndex[sIndex].add(region);
		}
		workingParsed+= regions.size();
		regions.clear();
	}

	/**Gets a file to parse containing records that match a particular chr. Thread safe.*/
	synchronized File getFileToParse(){
		if (workingFilesToParseIndex < workingFilesToParse.size()) return workingFilesToParse.get(workingFilesToParseIndex++);
		return null;
	}

	private void createChrFiles() {
		IO.pl("\nChecking tbi files for chr content...");
		
		//load hash to hold files with a particular chr
		for (String chr: chrLengths.keySet()) chrFiles.put(chr, new ArrayList<File>());

		try {
			//for each file
			for (File f: dataFilesToParse){
				//make an index
				File i = new File (f+".tbi");
				TabixIndex ti = new TabixIndex(i);
				//for each chromosome
				for (String chr: chrLengths.keySet()){
					if (ti.containsChromosome(chr) || ti.containsChromosome("chr"+chr)) chrFiles.get(chr).add(f);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem with testing indexes for particular chroms\n");
		}
	}
	
	public void saveWorkingChrBlock(){

		HashSet<IndexRegion> openRegions = new HashSet<IndexRegion>();
		int startPos = -1;
		ArrayList<RegionToPrint> toSave = new ArrayList<RegionToPrint>();

		for (int i=0; i< workingIndex.length; i++){
			if (workingIndex[i]== null) continue;

			//first in block?
			if (openRegions.size() == 0){
				openRegions.addAll(workingIndex[i]);
				startPos = i;
				if (toSave.size() != 0) saveRegions(toSave, out);
			}

			//must be adding a new so need to save old
			else {
				toSave.add(new RegionToPrint(startPos, i, fetchIds(openRegions)));
				startPos = i;

				//go through each of the current index regions
				for (IndexRegion ir: workingIndex[i]){
					//already in open then this must be an end so remove it
					if (openRegions.contains(ir)) openRegions.remove(ir);
					//not present so this is a beginning, add it
					else openRegions.add(ir);
				}
			}
		}
		//save last?
		if (toSave.size() != 0) saveRegions(toSave, out);
	}
	
	private void saveRegions(ArrayList<RegionToPrint> al, PrintWriter out){
		//set first
		RegionToPrint rtp = al.get(0);
		
		//walk looking for same ids in following regions and merge em
		int num = al.size();
		for (int i=1; i< num; i++){
			RegionToPrint next = al.get(i);
			//same file id's
			if (next.ids.equals(rtp.ids)) rtp.stop = next.stop;
			else {
				StringBuilder sb = new StringBuilder(workingChr);
				rtp.appendInfo(sb, workingStartBp);
				out.println(sb);
				rtp = next;
			}
		}
		//print last
		StringBuilder sb = new StringBuilder(workingChr);
		rtp.appendInfo(sb, workingStartBp);
		out.println(sb);
		al.clear();
		
	}

	private class RegionToPrint{
		int start;
		int stop;
		String ids;
		
		public RegionToPrint(int start, int stop, String ids){
			this.start = start;
			this.stop = stop;
			this.ids = ids;
		}

		public void appendInfo(StringBuilder sb, int blockStartPosition) {
			sb.append("\t");
			sb.append(start+blockStartPosition);
			sb.append("\t");
			sb.append(stop+blockStartPosition);
			sb.append("\t");
			sb.append(ids);
		}
	}
	
	
	private static String fetchIds(HashSet<IndexRegion> al) {
		//just one?
		if (al.size() == 1) return new Integer (al.iterator().next().fileId).toString();
		
		//hash ids and sort, might have dups
		TreeSet<Integer> ids = new TreeSet<Integer>();
		Iterator<IndexRegion> it = al.iterator();
		while (it.hasNext())ids.add(it.next().fileId);
	
		//create string
		Iterator<Integer>iit = ids.iterator();
		StringBuilder sb = new StringBuilder(iit.next().toString());
		while (iit.hasNext()){
			sb.append(",");
			sb.append(iit.next().toString());
		}
		
		return sb.toString();
	}

	public void bgzipAndTabixIndex() {
		File[] beds = IO.extractFiles(indexDir, ".qi.bed");
		for (File bed: beds){
			//compress with bgzip, this will replace any existing compressed file and delete the uncompressed
			String[] cmd = { bgzip.toString(), "-f", "--threads", numberThreads+"", bed.toString()};
			String[] output = IO.executeViaProcessBuilder(cmd, false);
			File compBed = new File (bed+".gz");
			if (output.length != 0 || bed.exists() == true || compBed.exists() == false){
				Misc.printErrAndExit("\nERROR: Failed to bgzip compress "+bed+"\nMessage: "+Misc.stringArrayToString(output, "\n"));
			}

			//tabix
			//must use -0 --sequence 1 --begin 2 --end 3; -p bed doesn't work with java tabix!!!!
			cmd = new String[]{ tabix.toString(), "-0", "--sequence", "1", "--begin", "2", "--end", "3", compBed.toString()};

			output = IO.executeViaProcessBuilder(cmd, false);
			File indexBed = new File (compBed+".tbi");
			if (output.length != 0 || indexBed.exists() == false){
				Misc.printErrAndExit("\nERROR: Failed to tabix index "+compBed+"\nMessage: "+Misc.stringArrayToString(output, "\n"));
			}
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new QueryIndexer (args);
	}	

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		File tabixBinDirectory = null;
		skipDirs = null;
		IO.pl("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': chrNameLength = new File(args[++i]); break;
					case 'd': dataDir = new File(args[++i]).getCanonicalFile(); break;
					case 'i': indexDir = new File(args[++i]); break;
					case 's': skipDirs = Misc.COMMA.split(args[++i]); break;
					case 'q': verbose = false; break;
					case 't': tabixBinDirectory = new File(args[++i]); break;
					case 'n': numberThreads = Integer.parseInt(args[++i]); break;
					case 'b': bpBlock = Integer.parseInt(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//pull tabix and bgzip executables
		if (tabixBinDirectory == null) Misc.printExit("\nError: please point to the dir containing the tabix and bgzip HTSlib executibles (e.g. /Users/Clinton/BioApps/HTSlib/1.3/bin/ )\n");
		bgzip = new File (tabixBinDirectory, "bgzip");
		tabix = new File (tabixBinDirectory, "tabix");
		if (bgzip.canExecute() == false || tabix.canExecute() == false) Misc.printExit("\nCannot find or execute bgzip or tabix executables from "+bgzip+" "+tabix);

		if (chrNameLength == null) Misc.printErrAndExit("\nError: please provide a bed file of chromosome and their max lengths to index. e.g. X 0 155270560\n" );
		if (dataDir == null || dataDir.isDirectory() == false) Misc.printErrAndExit("\nERROR: please provide a directory containing gzipped and tabix indexed bed, vcf, maf.txt, and bedGraph files to index." );
		if (indexDir == null ) Misc.printErrAndExit("\nERROR: please provide a directory in which to write the master query index." );
		
		if (indexDir.exists() == false) indexDir.mkdirs();

		//threads to use
		int numAvail = Runtime.getRuntime().availableProcessors();
		if (numberThreads < 1 || numberThreads > numAvail) numberThreads =  numAvail - 1;
		IO.pl(numAvail +" available processors, using "+numberThreads);
		
		toTruncatePoint = dataDir.getParentFile().toString().length()+1;
	}	

	private void parseSkipDirs(String[] skipDirs) {
		try {
			if (skipDirs != null){
				String indexDirCP = dataDir.getCanonicalPath();
				for (int i=0; i< skipDirs.length; i++){
					File f = new File(skipDirs[i]);
					if (f.exists() == false) Misc.printErrAndExit("\nError: failed to find one of your skip directories? "+skipDirs[i]);
					if (f.isDirectory() == false) Misc.printErrAndExit("\nError: is this skip directory not a directory? "+skipDirs[i]);
					String conPath = f.getCanonicalPath();
					if (conPath.startsWith(indexDirCP) == false) Misc.printErrAndExit("\nError: this skip directory doesn't appear to be within the data directory? "+skipDirs[i]);
					skipDirs[i] = conPath;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void parseChromLengthFile() {
		HashMap<String, RegionScoreText[]> chrLen = Bed.parseBedFile(chrNameLength, true, false);
		chrLengths = new TreeMap<String, Integer>();
		for (String chr: chrLen.keySet()){
			//find max
			RegionScoreText[] regions = chrLen.get(chr);
			int max = -1;
			for (RegionScoreText r: regions){
				if (r.getStop()> max) max = r.getStop();
			}
			//drop chr
			if (chr.startsWith("chr")) chr = chr.substring(3);

			//any already present?
			Integer priorLen = chrLengths.get(chr);
			if (priorLen == null || priorLen.intValue() < max) chrLengths.put(chr, max);
		}
	}

	private void setKnownFileTypeSSS() {
		//0 index start stop column info for "bed.gz", "bedGraph.gz", "maf.txt.gz"
		//the last is the number to subtract from the start to convert to interbase
		//note vcf.gz is it's own beast and handled separately
		extensionsStartStopSub.put("bed.gz", new int[]{1,2,0});
		extensionsStartStopSub.put("bedgraph.gz", new int[]{1,2,0});
		extensionsStartStopSub.put("maf.txt.gz", new int[]{5,6,1});
	}
	
	synchronized int[] getSSS(File f) {
		int[] startStopSubtract = null;
		try {
			//pull known
			String name = f.getName().toLowerCase();
			if (name.endsWith(".bed.gz")) startStopSubtract = extensionsStartStopSub.get("bed.gz");
			else if (name.endsWith(".maf.txt.gz")) startStopSubtract =  extensionsStartStopSub.get("maf.txt.gz");
			else if (name.endsWith(".bedgraph.gz")) startStopSubtract =  extensionsStartStopSub.get("bedgraph.gz");
			
			//pull from tbi
			//I'm suspicious of the flags field so the subtract may be wrong.  Best to manually set for each file type.
			else {
				File index = new File (f.getCanonicalPath()+".tbi");
				TabixIndex ti = new TabixIndex(index);
				TabixFormat tf = ti.getFormatSpec();
				int startIndex = tf.startPositionColumn- 1;
				int stopIndex = tf.endPositionColumn- 1;
				int sub = 1;
				if (tf.flags == 65536) sub = 0;
				startStopSubtract = new int[]{startIndex, stopIndex, sub};
				//IO.pl("New data type! ");
				//Misc.printArray(startStopSubtract);
			}
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("Error loading tbi index, aborting.");
		}
		return startStopSubtract;
	}

	
	private void saveFileIds(){
		try {

			//find all the parsed files
			HashSet<File> parsedFiles = new HashSet<File>();
			for (ArrayList<File> al: chrFiles.values()) parsedFiles.addAll(al);

			//save ids with truncated paths
			int toSkip = dataDir.getParentFile().toString().length()+1;

			HashMap<String, Integer> fileStringId = new HashMap<String, Integer>();
			HashMap<String, Long> fileStringSize = new HashMap<String, Long>();
			TreeSet<String> skippedFileNames = new TreeSet<String>();

			for (File f: parsedFiles){
				//Hmm, this only works with actual files, not soft links!
				String trimmedName = f.toString().substring(toSkip);

				//save trmmed name : id  and size
				Integer id = fileId.get(f);
				if (id == null) throw new Exception("\nFailed to find an id for "+f);
				fileStringId.put(trimmedName, id);
				fileStringSize.put(trimmedName, f.length());
				
				//a skipped file?
				if (skippedDataSources.contains(f)) skippedFileNames.add(trimmedName);
			}
			
			//save the trunk file name : id
			File ids = new File(indexDir, "fileIds.obj");
			IO.saveObject(ids, fileStringId);
			
			//save the trunk file name : size
			File sizes = new File(indexDir, "fileSizes.obj");
			IO.saveObject(sizes, fileStringSize);
			
			//save the skipped files, might be empty
			File skippedFiles = new File(indexDir, "skippedSources.obj");
			IO.saveObject(skippedFiles, skippedFileNames);
			

		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: saving file objects, aborting\n");
		}
	}
	
	
	/*  Don't save these, it'll explode, just pull from live file
		private String[] stringHeaderPatternStarts = {"^[#/<@].+", "^browser.+", "^track.+", "^color.+", "^url.+", "^Hugo_Symbol.+"};
	Pattern[] headerStarts = new Pattern[stringHeaderPatternStarts.length];
	for (int i =0; i< stringHeaderPatternStarts.length; i++) headerStarts[i] = Pattern.compile(stringHeaderPatternStarts[i]); 
	private ArrayList<String> parseHeader(File f, Pattern[] headerLineStarts) throws IOException {
		ArrayList<String> header = new ArrayList<String>();
		BufferedReader in = IO.fetchBufferedReader(f);
		String line;
		while ((line = in.readLine()) != null){
			//skip blanks
			line = line.trim();
			if (line.length()==0) continue;
			//scan patterns
			boolean found = false;
			for (Pattern p: headerLineStarts) {
				Matcher m = p.matcher(line);
				if (m.matches()) {
					header.add(line);
					found = true;
					break;
				}
			}
			if (found == false) break;
		}
		in.close();
		return header;
	}*/
	
	private void parseDataSources(){
		IO.pl("\nSearching for tabix indexes and bgzipped data sources...");
		try {
			ArrayList<File> goodDataSources = new ArrayList<File>();
			ArrayList<File> tbiMissingDataSources = new ArrayList<File>();
			ArrayList<File> unrecognizedDataSource = new ArrayList<File>();

			//find .tbi files and check
			File[] tbis = IO.fetchFilesRecursively(dataDir, ".gz.tbi");
			for (File tbi: tbis){
				String path = tbi.getCanonicalPath();

				//look for data file
				File df = new File(path.substring(0, path.length()-4));
				if (df.exists() == false){
					tbiMissingDataSources.add(tbi);
					continue;
				}
				
				//is it a known type
				boolean recognized = false;
				for (String knownExt: fileExtToIndex){
					if (df.getName().toLowerCase().endsWith(knownExt)){
						recognized = true;
						break;
					}
				}
				if (recognized) goodDataSources.add(df); 
				else unrecognizedDataSource.add(df);
			}

			int numFiles = goodDataSources.size() + unrecognizedDataSource.size();

			//print messages
			if (numFiles == 0) Misc.printErrAndExit("\nERROR: failed to find any data sources with xxx.tbi indexes in your dataDir?!\n");
			IO.pl("\t"+goodDataSources.size()+" Data sources with known formats ("+ Misc.stringArrayToString(fileExtToIndex, ", ")+")");
			if (tbiMissingDataSources.size() !=0){
				IO.pl("\t"+tbiMissingDataSources.size()+" WARNING: The data source file(s) for the following tbi index(s) could not be found, skipping:");
				for (File f: tbiMissingDataSources) IO.pl("\t\t"+f.getCanonicalPath());
			}
			if (unrecognizedDataSource.size() !=0){
				IO.pl("\t"+unrecognizedDataSource.size()+" WARNING: Data sources with unknown format(s). The format of the "
						+ "following files will be set using info from the tabix index and may be incorrect. Contact "+email+" to add.");
				for (File f: unrecognizedDataSource) IO.pl("\t\t"+f.getCanonicalPath());
			}
			
			//make final set
			goodDataSources.addAll(unrecognizedDataSource);
			dataFilesToParse = new File[goodDataSources.size()];
			goodDataSources.toArray(dataFilesToParse);
			Arrays.sort(dataFilesToParse);
			
		} catch (IOException e){
			e.printStackTrace();
			Misc.printErrAndExit("\nError examining data source files. \n");
		}
	}

	public static File[] returnFilesWithTabix(File[] tabixFiles) {
		ArrayList<File> goodFiles = new ArrayList<File>();
		for (File tb: tabixFiles){
			File index = new File (tb.toString()+".tbi");
			if (index.exists()) goodFiles.add(tb);
		}
		//any files?
		if (goodFiles.size()==0) return null;
		File[] toReturn = new File[goodFiles.size()];
		goodFiles.toArray(toReturn);
		return toReturn;
	}


	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                               Query Indexer: Jan 2019                            **\n" +
				"**************************************************************************************\n" +
				"Builds index files for Query by recursing through a data directory looking for bgzip\n"+
				"compressed and tabix indexed genomic data files (e.g. vcf, bed, maf, and custom).\n"+
				"Interval trees are built containing regions that overlap with one or more data sources.\n"+
				"These are used by the Query REST service to rapidly identify which data files contain\n"+
				"records that overlap user's ROI. This app is threaded for simultanious file loading\n"+
				"and requires >30G RAM to run on large data collections so use a big analysis server.\n"+
				"Note, relative file paths are saved. So long as the structure of the Data Directory is\n"+
				"preserved, the QueryIndexer and Query REST service don't need to run on the same file\n"+
				"system.\n"+


				"\nRequired Params:\n"+
				"-c A bed file of chromosomes and their lengths (e.g. chr21 0 48129895) to use to \n"+
				"     building the intersection index. Exclude those you don't want to index. For\n"+
				"     multiple builds and species, add all, duplicates will be collapsed taking the\n"+
				"     maximum length. Any 'chr' prefixes are ignored when indexing and searching.\n"+
				"-d A data directory containing bgzipped and tabix indexed data files. Known file types\n"+
				"     include xxx.vcf.gz, xxx.bed.gz, xxx.bedGraph.gz, and xxx.maf.txt.gz. Others will\n"+
				"     be parsed using info from the xxx.gz.tbi index. Be sure to normalize and\n"+
				"     decompose_blocksub all VCF records, see http://genome.sph.umich.edu/wiki/Vt.\n"+
				"     Files may be hard linked but not soft.\n"+
				"-t Full path directory containing the compiled bgzip and tabix executables. See\n" +
				"     https://github.com/samtools/htslib\n"+
				"-i A directory in which to save the index files\n"+

				"\nOptional Params:\n"+
				"-s One or more directory paths, comma delimited no spaces, to skip when building\n"+
				"     interval trees but make available for data source record retrieval. Useful for\n"+
				"     whole genome gVCFs and read coverage files that cover large genomic regions.\n"+
				"-q Quiet output, no per record warnings.\n"+
				"-u Update query index, defaults to building it anew. Only works when adding or\n"+
				"     deleting files from the data dir. Any other changes requires a full rebuild.\n"+
				"-b BP block to process, defaults to 250000000. Reduce if out of memory issues occur.\n"+
				"-n Number cores to use, defaults to all\n"+
				
				"\nWARNING: for bed file tabix indexing don't use -p/--preset option, it doesn't\n"+
				"     work with the htsjdk java classes. Use '-0 --sequence 1 --begin 2 --end 3'\n"+

				"\nExample for generating the test index using the GitHub Query/TestResources files see\n"+
				"https://github.com/HuntsmanCancerInstitute/Query\n\n"+
				
				"d=/pathToYourLocalGitHubInstalled/Query/TestResources\n"+
				"java -Xmx120G -jar pathToUSeq/Apps/QueryIndexer -c $d/b37Chr20-21ChromLen.bed -d $d/Data\n"+
				"-i $d/Index -t ~/BioApps/HTSlib/1.3/bin/ -s $d/Data/Public/B37/GVCFs \n\n" +

				"**************************************************************************************\n");
	}

	public int getWorkingIndexLength() {
		return workingIndex.length;
	}

	public boolean isVerbose() {
		return verbose;
	}

	public HashMap<File, Integer> getFileId() {
		return fileId;
	}

	public Integer[] getIds() {
		return ids;
	}

}
