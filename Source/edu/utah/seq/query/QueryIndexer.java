package edu.utah.seq.query;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.json.JSONObject;
import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
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
	private File[] dataFilesToParse;
	private String[] fileExtToIndex = {"vcf.gz", "bed.gz", "bedgraph.gz", "maf.txt.gz"};
	private HashSet<File> skippedDataSources = new HashSet<File>();
	private static HashMap<String, int[]> extensionsStartStopSub = new HashMap<String, int[]>();
	private TreeMap<String, Integer> chrLengths;
	private String[] stringHeaderPatternStarts = {"^[#/<@].+", "^browser.+", "^track.+", "^color.+", "^url.+", "^Hugo_Symbol.+"};

	private HashMap<String, ArrayList<File>> chrFiles = new HashMap<String, ArrayList<File>>();
	private HashMap<File, Integer> fileId = new HashMap<File, Integer>();
	private long totalRecordsProcessed = 0;


	//per chr fields
	private ArrayList<IndexRegion>[] workingIndex = null;
	private ArrayList<File> workingFilesToParse = null;
	private long workingPassed = 0;
	private long workingFailed = 0;
	private String workingChr = null;

	//constructor
	@SuppressWarnings("unchecked")
	public QueryIndexer(String[] args) {
		long startTime = System.currentTimeMillis();

		processArgs(args);

		System.out.println("Creating file id hash...");
		createFileIdHash();

		System.out.println("Checking tbi files for chr content...");
		createChrFiles();

		//for each chromosome
		if (verbose) System.out.println("Indexing records by chr...\n\tchr\t#pass\t#fail\tmemUsed");
		else System.out.print("Indexing records by chr...\n\t"); 
		for (String chr: chrLengths.keySet()) {
			workingChr = chr;
			parseChr();
		}

		if (verbose == false) System.out.println(": "+IO.memory());

		//bgzip and tabix
		System.out.println("Compressing master index...");
		bgzipAndTabixIndex();

		System.out.println("Saving file headers...");
		saveFileHeadersAndIds();

		String diffTime = Num.formatNumberOneFraction(((double)(System.currentTimeMillis() -startTime))/1000/60);
		System.out.println("\n"+ diffTime+" Min to parse "+ totalRecordsProcessed+" records and build the query index");
	}

	public void parseChr(){
		try {
			//any files?
			workingFilesToParse = new ArrayList<File>();
			workingFilesToParse.addAll(chrFiles.get(workingChr));
			
			//remove any that belong to a skip dir
			removeSkipDirs(workingFilesToParse);
			
			if (workingFilesToParse.size() == 0) return; 

			//prep for new chr
			workingIndex = new ArrayList[chrLengths.get(workingChr)+1];
			workingPassed = 0;
			workingFailed = 0;

			//try to make a loader for each file
			int numToMake= workingFilesToParse.size();
			if (numToMake > numberThreads) numToMake = numberThreads;			
			QueryIndexFileLoader[] loader = new QueryIndexFileLoader[numToMake];			
			ExecutorService executor = Executors.newFixedThreadPool(numToMake);
			for (int i=0; i< loader.length; i++){
				loader[i] = new QueryIndexFileLoader(this, workingChr);
				executor.execute(loader[i]);
			}
			executor.shutdown();

			//spins here until the executer is terminated, e.g. all threads complete
			while (!executor.isTerminated()) {}

			//check loaders 
			for (QueryIndexFileLoader c: loader) {
				if (c.isFailed()) throw new IOException("ERROR: File Loader issue! \n"+c);
			}

			//save the index for this chrom
			saveIndex();

			//cleanup
			if (verbose) System.out.println("\t"+workingChr+"\t"+workingPassed+"\t"+workingFailed+"\t"+IO.memory());
			else System.out.print(workingChr+" ");
			totalRecordsProcessed+= workingPassed;
			totalRecordsProcessed+= workingFailed;
		} catch (IOException e){
			e.printStackTrace();
			Misc.printErrAndExit("\nFATAL error with loading "+workingChr+", aborting.");
		}
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
	}

	/**Gets a file to parse containing records that match a particular chr. Thread safe.*/
	synchronized File getFileToParse(){
		if (workingFilesToParse.size() == 0) return null;
		return workingFilesToParse.remove(0);
	}

	public synchronized void incrementPassFail(long pass, long fail){
		workingPassed+= pass;
		workingFailed+= fail;
	}

	private void createChrFiles() {
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


	private void createFileIdHash() {
		for (int i=0; i< dataFilesToParse.length; i++) {
			fileId.put(dataFilesToParse[i], i);
		}
		//contains 0 thru final no interruptions
	}


	public void saveIndex(){
		try {
			File queryIndexFile = new File(indexDir, workingChr+".qi.bed");
			PrintWriter out = new PrintWriter(new FileWriter((queryIndexFile)));
			HashSet<IndexRegion> openRegions = new HashSet<IndexRegion>();
			int startPos = -1;
			ArrayList<RegionToPrint> toSave = new ArrayList<RegionToPrint>();
			
			for (int i=0; i< workingIndex.length; i++){
				if (workingIndex[i]== null) continue;
				//System.out.println(" xxxxxxxxxxxxxxxxxxxxxxxxxxxxx Index "+i);
				
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
			
			//close it
			out.close();
		} catch (IOException e){
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem saving the index for "+workingChr+", aborting.\n");
		}
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
				rtp.appendInfo(sb);
				out.println(sb);
				rtp = next;
			}
		}
		//print last
		StringBuilder sb = new StringBuilder(workingChr);
		rtp.appendInfo(sb);
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

		public void appendInfo(StringBuilder sb) {
			sb.append("\t");
			sb.append(start);
			sb.append("\t");
			sb.append(stop);
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
			String[] cmd = { bgzip.toString(), "-f", bed.toString()};
			String[] output = IO.executeViaProcessBuilder(cmd, false);
			File compBed = new File (bed+".gz");
			if (output.length != 0 || bed.exists() == true || compBed.exists() == false){
				Misc.printErrAndExit("\nERROR: Failed to bgzip compress "+bed+"\nMessage: "+Misc.stringArrayToString(output, "\n"));
			}

			//tabix
			cmd = new String[]{ tabix.toString(), "-f", "-p", "bed", compBed.toString() };
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
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
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
		
		//remove contents of indexDir
		IO.deleteDirectory(indexDir);
		indexDir.mkdirs();
		
		parseSkipDirs(skipDirs);
		setKnownFileTypeSSS();
		parseDataSources();
		parseChromLengthFile();

		//threads to use
		int numAvail = Runtime.getRuntime().availableProcessors();
		if (numberThreads < 1) numberThreads =  numAvail - 1;
		System.out.println(numAvail +" available processors, using "+numberThreads+" threaded loaders");
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

	/*Methods to modify with hard coded file extensions*/

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
				//System.out.println("New data type! ");
				//Misc.printArray(startStopSubtract);
			}
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("Error loading tbi index, aborting.");
		}
		return startStopSubtract;
	}

	private void saveFileHeadersAndIds(){
		try {
			//make Patterns to scan
			Pattern[] headerStarts = new Pattern[stringHeaderPatternStarts.length];
			for (int i =0; i< stringHeaderPatternStarts.length; i++) headerStarts[i] = Pattern.compile(stringHeaderPatternStarts[i]); 

			//find all the parsed files
			HashSet<File> parsedFiles = new HashSet<File>();
			for (ArrayList<File> al: chrFiles.values()) parsedFiles.addAll(al);

			//save ids with truncated paths
			int toSkip = dataDir.getParentFile().toString().length()+1;

			HashMap<String, Integer> fileStringId = new HashMap<String, Integer>();
			HashMap<String, String> fileHeaders = new HashMap<String, String>();
			TreeSet<String> skippedFileNames = new TreeSet<String>();

			for (File f: parsedFiles){
				//Hmm, this only works with actual files, not soft links!
				String trimmedName = f.toString().substring(toSkip);

				//save trmmed name : id ?
				Integer id = fileId.get(f);
				if (fileStringId.containsKey(trimmedName)){
					Misc.printErrAndExit("Duplicate trimmedName! "+trimmedName);
				}
				else fileStringId.put(trimmedName, id);

				//save trimmed name: json header, have to use string rep since JSONObject is not serializable
				ArrayList<String> header = parseHeader(f, headerStarts);
				JSONObject jo = new JSONObject();
				jo.put("source", trimmedName);
				jo.put("header", header);
				fileHeaders.put(trimmedName, jo.toString(1));
				
				//a skipped file?
				if (skippedDataSources.contains(f)) skippedFileNames.add(trimmedName);
			}

			//save the headers for all
			File h = new File(indexDir, "fileHeaders.obj");
			IO.saveObject(h, fileHeaders);
			
			//save the trunk file name : id
			File ids = new File(indexDir, "fileIds.obj");
			IO.saveObject(ids, fileStringId);
			
			//save the skipped files, might be empty
			File skippedFiles = new File(indexDir, "skippedSources.obj");
			IO.saveObject(skippedFiles, skippedFileNames);
			

		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: parsing headers, aborting\n");
		}
	}

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
	}
	
	private void parseDataSources(){
		try {
			System.out.println("Searching for tabix indexes and bgzipped data sources...");
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
			System.out.println("\t"+goodDataSources.size()+" Data sources with known formats ("+ Misc.stringArrayToString(fileExtToIndex, ", ")+")");
			if (tbiMissingDataSources.size() !=0){
				System.out.println("\t"+tbiMissingDataSources.size()+" WARNING: The data source file(s) for the following tbi index(s) could not be found, skipping:");
				for (File f: tbiMissingDataSources) System.out.println("\t\t"+f.getCanonicalPath());
			}
			if (unrecognizedDataSource.size() !=0){
				System.out.println("\t"+unrecognizedDataSource.size()+" WARNING: Data sources with unknown format(s). The format of the "
						+ "following files will be set using info from the tabix index and may be incorrect. Contact bioinformaticscore@utah.edu to add.");
				for (File f: unrecognizedDataSource) System.out.println("\t\t"+f.getCanonicalPath());
			}
			
			//make final set
			goodDataSources.addAll(unrecognizedDataSource);
			dataFilesToParse = new File[goodDataSources.size()];
			goodDataSources.toArray(dataFilesToParse);
			
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
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Query Indexer: June 2019                           **\n" +
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


				"\nExample for generating the test index using the GitHub Query/TestResources files see\n"+
				"https://github.com/HuntsmanCancerInstitute/Query\n\n"+
				
				"d=/pathToYourLocalGitHubInstalled/Query/TestResources\n"+
				"java -Xmx10G -jar pathToUSeq/Apps/QueryIndexer -c $d/b37Chr20-21ChromLen.bed -d $d/Data\n"+
				"-i $d/Index -t ~/BioApps/HTSlib/1.3/bin/ -s $d/Data/Public/B37/GVCFs \n\n" +

				"**************************************************************************************\n");
	}

	public ArrayList<IndexRegion>[] getWorkingIndex() {
		return workingIndex;
	}

	public boolean isVerbose() {
		return verbose;
	}

	public HashMap<File, Integer> getFileId() {
		return fileId;
	}

}
