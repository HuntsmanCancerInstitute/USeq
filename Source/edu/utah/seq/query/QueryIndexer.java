package edu.utah.seq.query;

import java.io.BufferedReader;
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
import org.json.JSONObject;
import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
import edu.utah.seq.useq.data.RegionScoreText;
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
	private String[] fileExtToIndex = null;
	private String[] supportedFileExtensions = {"vcf.gz", "bed.gz", "bedGraph.gz", "maf.txt.gz"};
	private static HashMap<String, int[]> extensionsStartStopSub = new HashMap<String, int[]>();
	private TreeMap<String, Integer> chrLengths;
	private String[] stringHeaderPatternStarts = {"^#.+", "^browser.+", "^track.+", "^Hugo_Symbol.+"};

	private HashMap<String, ArrayList<File>> chrFiles = new HashMap<String, ArrayList<File>>();
	private HashMap<File, Integer> fileId = new HashMap<File, Integer>();
	private long totalRecordsProcessed = 0;


	//per chr fields
	private TreeSet<File>[] workingIndex = null;
	private HashMap<String, TreeSet<File>> fileNames2TheirPointers = new HashMap<String, TreeSet<File>>();
	private ArrayList<File> workingFilesToParse = null;
	private long workingPassed = 0;
	private long workingFailed = 0;

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
		for (String chr: chrLengths.keySet()) parseChr(chr);

		if (verbose == false) System.out.println(": "+IO.memory());

		//bgzip and tabix
		System.out.println("Compressing master index...");
		bgzipAndTabixIndex();

		System.out.println("Saving file headers...");
		saveFileHeadersAndIds();

		System.out.print("Building and saving interval trees...\n\t");
		saveIntervalTrees();

		String diffTime = Num.formatNumberOneFraction(((double)(System.currentTimeMillis() -startTime))/1000/60);
		System.out.println("\n"+ diffTime+" Min to parse "+ totalRecordsProcessed+" records and build the query index trees");
	}

	public void parseChr(String workingChr){
		try {
			//any files?
			workingFilesToParse = new ArrayList<File>();
			workingFilesToParse.addAll(chrFiles.get(workingChr));
			
			//remove any that belong to a skip dir
			removeSkipDirs(workingFilesToParse);
			
			if (workingFilesToParse.size() == 0) return; 

			//prep for new chr
			workingIndex = new TreeSet[chrLengths.get(workingChr)+1];
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
			saveIndex(workingChr);

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
			}
		} catch (IOException e) {
			e.printStackTrace();
		} 
	}

	/**Gets a file to parse contining records that match a particular chr. Thread safe.*/
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
					if (ti.containsChromosome(chr)) chrFiles.get(chr).add(f);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem with testing indexes for particular chroms\n");
		}
	}


	private void createFileIdHash() {
		int counter = 0;
		for (File f: dataFilesToParse) fileId.put(f, counter++);
	}

	public void saveIntervalTrees(){
		try {
			File[] chrBed = IO.extractFiles(indexDir, ".qi.bed.gz");
			HashMap<String, IntervalST<int[]>> chrTrees = new HashMap<String, IntervalST<int[]>>();
			HashMap<String, int[]> idsFileIndex = new HashMap<String, int[]>();
			//for each file
			for (File f: chrBed){
				BufferedReader in = IO.fetchBufferedReader(f);
				IntervalST<int[]> st = new IntervalST<int[]>();
				String line;
				while ((line = in.readLine())!= null){
					String[] t = Misc.TAB.split(line);
					int start = Integer.parseInt(t[1]);
					int stop = Integer.parseInt(t[2]);
					int[] ids = idsFileIndex.get(t[3]);
					if (ids == null){
						ids = stringToInts(t[3], Misc.COMMA);
						idsFileIndex.put(t[3], ids);
					}
					//the end is included in IntervalST so sub 1 from end
					st.put(new Interval1D(start, stop-1), ids);
				}
				in.close();

				//save tree
				String chr = f.getName().replace(".qi.bed.gz", "");

				chrTrees.put(chr, st);
				System.out.print(chr+" ");
			}
			System.out.println(": "+IO.memory());

			File t = new File (indexDir, "its.obj");
			IO.saveObject(t, chrTrees);

		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem saving interval trees, aborting.\n");
		}
	}

	/**Given a String of ints delimited by something, will parse or return null.*/
	public static int[] stringToInts(String s, Pattern pat){
		String[] tokens = pat.split(s);
		int[] num = new int[tokens.length];
		try {
			for (int i=0; i< tokens.length; i++){
				num[i] = Integer.parseInt(tokens[i]);
			}
			return num;
		} catch (Exception e){
			return null;
		}
	}


	public void saveIndex(String workingChr){
		try {
			File queryIndexFile = new File(indexDir, workingChr+".qi.bed");
			PrintWriter out = new PrintWriter(new FileWriter((queryIndexFile)));
			//collect bps with the same set of files
			//start first entry, should always contain one or more files
			int start = Integer.MIN_VALUE;
			String priorIds = "";
			int stop = Integer.MIN_VALUE;

			for (int i=0; i< workingIndex.length; i++){
				if (workingIndex[i] == null) {
					//in a block? write out old if not the first
					if (start != Integer.MIN_VALUE){
						out.print(workingChr); out.print("\t");
						out.print(start); out.print("\t");
						out.print(i); out.print("\t");
						out.println(priorIds);

						//reset start 
						priorIds = "";
						start = Integer.MIN_VALUE;
					}
				}

				//nope hash present
				else {
					//save last good position
					stop = i;
					//different ids?
					String currentIds = fetchSortedIds(workingIndex[i]);
					if (priorIds.equals(currentIds) == false){
						//write out old if not the first
						if (start != Integer.MIN_VALUE){
							out.print(workingChr); out.print("\t");
							out.print(start); out.print("\t");
							out.print(stop); out.print("\t");
							out.println(priorIds);
						}
						//start new
						priorIds = currentIds;
						start = i;
					}
				}
			}
			//save last?
			if (start != Integer.MIN_VALUE){
				out.print(workingChr); out.print("\t");
				out.print(start); out.print("\t");
				out.print(stop+1); out.print("\t");
				out.println(priorIds);
			}
			//close it
			out.close();
		} catch (IOException e){
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem saving the index for "+workingChr+", aborting.\n");
		}
	}



	/**Returns the int id concat rep for all the files in the sorted TreeSet*/
	private String fetchSortedIds(TreeSet<File> files) {
		int[] ids = new int[files.size()];
		int counter = 0;
		for (File f: files) ids[counter++] = fileId.get(f);
		return Misc.intArrayToString(ids, ",");
	}



	/**Need to run through some gymnastics here since hits to the current TreeSet[] will be out of order and coming in from multiple threads.
	 * Using this method as the bottleneck to keep the additions current. */
	synchronized void addRef(int bp, File f) {

		Integer fId = fileId.get(f);
		//nothing entered just yet?
		if (workingIndex[bp] == null) {
			//any existing tree
			TreeSet<File> ts = fileNames2TheirPointers.get(fId.toString());
			if (ts == null){
				ts = new TreeSet<File>();
				ts.add(f);
				fileNames2TheirPointers.put(fId.toString(), ts);
			}
			workingIndex[bp] = ts;
		}

		//does the current hash already contain the file due to overlapping regions?
		else if (workingIndex[bp].contains(f) == false){
			//ok something there but it doesn't contain this file
			//attempt to find an existing hash with these files
			String key = fetchCompositeIdKey(workingIndex[bp], fId);
			TreeSet<File> ts = fileNames2TheirPointers.get(key);
			if (ts == null){
				ts = new TreeSet<File>();
				ts.addAll(workingIndex[bp]);
				ts.add(f);
				fileNames2TheirPointers.put(key, ts);
			}
			workingIndex[bp] = ts;
		}
		//yes, so don't do anything it's already there
	}


	/**This is expensive in computing time. Returns a String of sorted ints representing the file ids.*/
	private String fetchCompositeIdKey(TreeSet<File> set, Integer newId) {
		int numOld = set.size();

		//make int[]
		int[] ids = new int[numOld+1];
		Iterator<File> it = set.iterator();
		for (int i=0; i< numOld; i++) ids[i] = fileId.get(it.next());
		ids[numOld] = newId;

		//sort
		Arrays.sort(ids);

		//concat
		StringBuilder sb = new StringBuilder();
		sb.append(ids[0]);
		for (int i=1; i< ids.length; i++){
			sb.append(",");
			sb.append(ids[i]);
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
					case 'd': dataDir = new File(args[++i]); break;
					case 'i': indexDir = new File(args[++i]); break;
					case 'e': fileExtToIndex = Misc.COMMA.split(args[++i]); break;
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
		indexDir.mkdirs();
		
		parseSkipDirs(skipDirs);
		parseFileExtensions();
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

	private void parseFileExtensions() {
		if (fileExtToIndex == null) fileExtToIndex = supportedFileExtensions;
		else {
			HashSet<String> ok = new HashSet<String>();
			for (String e : supportedFileExtensions) ok.add(e);
			for (String t: fileExtToIndex){
				if (ok.contains(t) == false) Misc.printErrAndExit("\nError: one or more of your requested data source extensions isn't supported, choose from these "+ok+"\n");
			}
		}
		//0 index start stop column info for "bed.gz", "bedGraph.gz", "maf.txt.gz"
		//the last is the number to subtract from the start to convert to interbase
		//note vcf.gz is it's own beast and handled separately
		extensionsStartStopSub.put("bed.gz", new int[]{1,2,0});
		extensionsStartStopSub.put("bedGraph.gz", new int[]{1,2,0});
		extensionsStartStopSub.put("maf.txt.gz", new int[]{5,6,1});
	}

	public int[] getSSS(File f) {
		String name = f.getName();
		if (name.endsWith(".bed.gz")) return extensionsStartStopSub.get("bed.gz");
		if (name.endsWith(".maf.txt.gz")) return extensionsStartStopSub.get("maf.txt.gz");
		if (name.endsWith(".bedGraph.gz")) return extensionsStartStopSub.get("bedGraph.gz");
		return null;
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

			for (File f: parsedFiles){
				String trimmedName = f.toString().substring(toSkip);

				//save trmmed name : id
				Integer id = fileId.get(f);
				fileStringId.put(trimmedName, id);

				//save trimmed name: json header, have to use string rep since JSONObject is not serializable
				ArrayList<String> header = parseHeader(f, headerStarts);
				JSONObject jo = new JSONObject();
				jo.put("source", trimmedName);
				jo.put("header", header);
				fileHeaders.put(trimmedName, jo.toString(1));				
			}

			//save them
			File ids = new File(indexDir, "fileIds.obj");
			IO.saveObject(ids, fileStringId);
			File h = new File(indexDir, "fileHeaders.obj");
			IO.saveObject(h, fileHeaders);

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
		System.out.println("Searching for tabix data sources...");

		File[][] toCombine = new File[fileExtToIndex.length][];
		for (int i=0; i< fileExtToIndex.length; i++) toCombine[i] = IO.fetchFilesRecursively(dataDir, fileExtToIndex[i]);

		dataFilesToParse = IO.collapseFileArray(toCombine);
		int numFiles = dataFilesToParse.length;

		dataFilesToParse = returnFilesWithTabix(dataFilesToParse);
		int numFilesWithIndex = dataFilesToParse.length;

		if (numFilesWithIndex == 0) Misc.printErrAndExit("\nERROR: failed to find any data sources with xxx.tbi indexes in your dataDir?!\n");
		System.out.println("\t"+numFiles+" Data sources ("+ Misc.stringArrayToString(fileExtToIndex, ", ")+")");
		System.out.println("\t"+numFilesWithIndex+" Data sources with xxx.tbi indexes");
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
				"**                                Query Indexer: Nov 2016                           **\n" +
				"**************************************************************************************\n" +
				"Builds index files for Query by recursing through a data directory looking for vcf,\n"+
				"maf, bed, and bedGraph files that have been bgzip compressed and tabix indexed.\n"+
				"Interval trees are build containing regions that overlap with one or more data sources.\n"+
				"These are used by the Query REST service to rapidly identify which data files contain\n"+
				"records that overlap each users query. This app typically needs >20G RAM to run.\n"+


				"\nRequired Params:\n"+
				"-c A bed file of chromosomes and their lengths (e.g. 21 0 48129895) to use to \n"+
				"     building the intersection index. Exclude those you don't want to index. For\n"+
				"     multiple builds and species, add all, duplicates will be collapsed taking the\n"+
				"     maximum length. Any 'chr' prefixes are ignored when indexing and searching.\n"+
				"-d A data directory containing gzipped tabixed xxx.vcf.gz, xxx.bed.gz, xxx.bedGraph.gz,\n"+
				"     and xxx.maf.txt.gz files. Recurses through all sub directories. Except for gVCF,\n"+
				"     Be sure to normalize and decompose_blocksub all VCF records, see \n"+
				"     http://genome.sph.umich.edu/wiki/Vt\n"+
				"-t Full path directory containing the compiled bgzip and tabix executables. See\n" +
				"     https://github.com/samtools/htslib\n"+
				"-i A directory in which to save the index files\n"+

				"\nOptional Params:\n"+
				"-e Comma delimited list, no spaces, of file type extensions to index. Defaults to\n"+
				"     vcf.gz,bed.gz,bedGraph.gz,maf.txt.gz\n"+
				"-s One or more directory paths, comma delimited no spaces, to skip when building\n"+
				"     interval trees but make available for data source record retrieval. Useful for\n"+
				"     whole genome gVCFs and read coverage files that cover large genomic regions.\n"+
				"-q Quiet output, no per record warnings.\n"+


				"\nExample for generating the test index using the GitHub Query/TestResources files see\n"+
				"https://github.com/HuntsmanCancerInstitute/Query\n\n"+
				"d=/pathToYourLocalGitHubInstalled/Query/TestResources\n"+
				"java -Xmx50G -jar pathToUSeq/Apps/QueryIndexer -c $d/b37Chr20-21ChromLen.bed -d $d/Data\n"+
				"-i $d/Index -t ~/BioApps/HTSlib/1.3/bin/ -s $d/Data/B37/GVCFs \n\n" +

				"**************************************************************************************\n");
	}

	public TreeSet<File>[] getWorkingIndex() {
		return workingIndex;
	}

	public boolean isVerbose() {
		return verbose;
	}

}
