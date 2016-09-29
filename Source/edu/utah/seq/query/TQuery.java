package edu.utah.seq.query;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.useq.data.RegionScoreText;
import htsjdk.tribble.readers.TabixReader;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;


/** Tool for intersecting regions with a collection of vcf and bed data. */
public class TQuery {

	//fields
	private File chrLengthFile;
	private File[] vcfDataFiles;
	private File[] bedDataFiles;
	private int numberThreads = 0;
	private int numberQueriesInChunk = 1000;
	private boolean printWarnings = true;
	private boolean printStats = true;
	private boolean fetchData = true;
	private boolean printJson = true;
	private int bpPadding = 0;

	//internal
	private long recordsLoaded = 0;
	private long recordsSkipped = 0;
	private boolean scoresSet = false;
	private TQueryFilter filter; 
	private ArrayList<TabixChunk> tabixChunks = new ArrayList<TabixChunk>();
	private Iterator<TabixChunk> chunkIterator = null;
	public static final Pattern END_POSITION = Pattern.compile(".*END=(\\d+).*", Pattern.CASE_INSENSITIVE);

	//interbase coordinates, 0 is first base, last base in any region is not included.
	private HashMap<String, HashSet<File>[]> chromFileIndex = new HashMap<String, HashSet<File>[]>();

	/*Contains a File pointer to a data source and the associated TabixQueries that intersect it.
	 * These are to be used in a tabix search.*/
	private HashMap<File, ArrayList<TabixQuery>> fileTabixQueries = new HashMap<File, ArrayList<TabixQuery>>();

	/*Contains all the TabixQueries split by chromosome. Some may contain results.*/
	private HashMap<String, TabixQuery[]> chrTabixQueries = new HashMap<String, TabixQuery[]>();

	/*Container for active TabixReaders*/
	private HashMap<File, TabixReader> tabixReaders = new HashMap<File, TabixReader>();

	//constructor
	public TQuery (String[] args) {
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			buildEngine();

			//print some stats on building the engine
			String diffTime = Num.formatNumberOneFraction(((double)(System.currentTimeMillis() -startTime))/1000);
			int numFiles = vcfDataFiles.length + bedDataFiles.length;
			System.err.println("\n"+ diffTime+" Sec to build using "+IO.memory()+ " of RAM");
			System.err.println("\t"+numFiles+"\tData sources loaded");
			System.err.println("\t"+recordsLoaded+"\tRecords indexed");
			System.err.println("\t"+recordsSkipped+"\tRecords skipped\n");

			System.err.println(filter.fetchSummary());

			queryFilesFromCmdLine();

			//close the readers
			closeTabixReaders();

		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("\nProblem with executing the TQuery!");
		}
	}

	private void queryFilesFromCmdLine() throws IOException {
		while (true){
			System.err.println("\nProvide a bed or vcf file to intersect (or blank to exit):");//plus any key=val1,val2,... (reset=, fetchData=, dataType=, fileType=, dataSource= , no spaces) options:\n");

			//parse input
			String in = (new BufferedReader(new InputStreamReader(System.in))).readLine().trim();
			if (in == null || in.length()==0) {
				System.err.println();
				break;
			}
			File inputFile = parseCmdLineInput(in);
			if (inputFile == null) continue;

			//find dataSource files that intersect regions in the input file bed or vcf, will return false if bed or vcf not found
			if (intersectFileRegionsWithIndex(inputFile)){

				//filter which files to fetch data from base on the parent dir name, the actual file, the file extension (bed or vcf)
				//these might not be set so just skipped
				filter.filter(fileTabixQueries);

				//fetch the underlying data (slow) or just indicate which file hit which region (fast)
				if (fetchData) loadTabixQueriesWithData();
				else loadTabixQueriesWithFileSources();

				//print results to stdout in json format
				if (printJson) printAllTabixQueries();
			}
		}
	}

	private File parseCmdLineInput(String cmd){

		//split on white space
		String[] args = Misc.WHITESPACE.split(cmd);

		//load hash
		File bedOrVcfFile = null;
		HashMap<String, String> options = new HashMap<String, String>();
		for (int i=0; i< args.length; i++) {

			//attempt to split on = sign
			String[] keyVal = Misc.PATTERN_EQUALS.split(args[i]);

			//nope must be a cmd or file
			if (keyVal.length !=2) {
				String lcArg = args[i].toLowerCase();
				if (lcArg.contains("reset")){
					filter.reset();
					fetchData = true;
					printJson = true;
					return null;
				}
				else if (lcArg.contains("filter")){
					System.err.println("\nCurrent Filters:\n"+filter.getFilters("\n"));
					return null;
				}
				//is it the file?
				else {
					bedOrVcfFile = new File(args[i]);
					if (bedOrVcfFile.exists() == false) {
						System.err.println("\nHmm can't split this key=value or understand the single cmd argument or find the file?! -> "+args[i]);
						return null;
					}
				}
			}
			else options.put(keyVal[0].toLowerCase(), keyVal[1]);
		}

		//walk through each option
		for (String key: options.keySet()){

			//fetchData?
			if (key.equals("fetchdata")){
				String bool = options.get("fetchdata").toLowerCase();
				if (bool.startsWith("t")) fetchData = true;
				else fetchData = false;
			}

			//print results
			else if (key.equals("printjson")){
				String bool = options.get("printjson").toLowerCase();
				if (bool.startsWith("t")) printJson = true;
				else printJson = false;
			}

			//filter on dataType?
			else if (key.equals("datatypes")){
				String[] dataTypes = Misc.COMMA.split(options.get("datatypes"));
				TreeSet<String> dtToReturn = new TreeSet<String>();
				//for each
				for (String dt: dataTypes){
					if (filter.getAvailableDataTypes().contains(dt)) dtToReturn.add(dt);
					else {
						System.err.println("\nHmm can't find this data type?! -> "+dt+ " in "+filter.fetchDataTypes(", "));
						return null;
					}
				}
				filter.setDataTypesToReturn(dtToReturn);
				filter.setFilterOnDataTypes(true);
			}

			//filter on fileType?
			else if (key.equals("filetypes")){
				String[] fileTypes = Misc.COMMA.split(options.get("filetypes"));
				TreeSet<String> toReturn = new TreeSet<String>();
				//for each
				for (String dt: fileTypes){
					if (filter.getAvailableFileTypes().contains(dt)) toReturn.add(dt);
					else {
						System.err.println("\nHmm can't find this file type?! -> "+dt+ " in "+filter.fetchFileTypes(", "));
						return null;
					}
				}

				filter.setFileTypesToReturn(toReturn);
				filter.setFilterOnFileTypes(true);
			}

			//filter on data source/ actual files?
			else if (key.equals("datasources")){
				String[] dataFiles = Misc.COMMA.split(options.get("datasource"));
				TreeSet<File> toReturn = new TreeSet<File>();
				//for each
				for (String fs: dataFiles){
					File f = new File(filter.getPathToTrimmedFile()+fs);
					if (filter.getAvailableDataFiles().contains(f)) toReturn.add(f);
					else {
						System.err.println("\nHmm can't find this data soure file?! -> "+fs+ " in "+filter.fetchDataFilesRelative(", "));
						return null;
					}
				}
				filter.setDataFilesToReturn(toReturn);
				filter.setFilterOnDataFiles(true);
			}


			else {
				System.err.println("\nSorry, this key=value isn't recognized -> "+ key+"="+options.get(key));
				return null;
			}
		}

		return bedOrVcfFile;

	}



	/**Synchronized method for getting or setting a TabixReader.*/
	public synchronized TabixReader getSetTabixReader(File tabixFile, TabixReader reader) {
		//do they want an active reader?
		if (reader == null) return tabixReaders.remove(tabixFile);
		//nope, want to put one back, do so if it doesn't exist
		else {
			//save it?
			if (tabixReaders.containsKey(tabixFile) == false) {
				tabixReaders.put(tabixFile, reader);
			}
			//nope, close it
			else reader.close();
			return null;
		}
	}

	/**Method for freeing up the file handles.*/
	private void closeTabixReaders(){
		for (TabixReader tr: tabixReaders.values()) tr.close();
	}

	private void loadTabixQueriesWithFileSources() {
		//for each data source file
		for (File tabixFile: fileTabixQueries.keySet()){
			//pull queries that intersect it
			ArrayList<TabixQuery> toFetch = fileTabixQueries.get(tabixFile);
			//for each query, add the data source
			for (TabixQuery tq: toFetch) tq.getSourceResults().put(tabixFile, null);

		}

	}

	private void printAllTabixQueries() {
		for (String chr: chrTabixQueries.keySet()){
			TabixQuery[] tqs = chrTabixQueries.get(chr);
			if (tqs.length > 1) {
				System.out.println("[");
				System.out.print(tqs[0].toJson(scoresSet, filter.getPathToTrimmedFile()));
				for (int i=1; i< tqs.length; i++){
					System.out.println(",");
					System.out.print(tqs[i].toJson(scoresSet, filter.getPathToTrimmedFile()));
				}
				System.out.println("\n]");
			}
			else System.out.println(tqs[0].toJson(scoresSet, filter.getPathToTrimmedFile()));
		}
		//flush the stream to while waiting for another query set
		System.out.flush();
	}

	private int countChunks(int size){
		int chunks = 0;
		for (File f: fileTabixQueries.keySet()){
			ArrayList<TabixQuery> al = fileTabixQueries.get(f);
			int numTQs = al.size();
			if (numTQs <= size) chunks++;
			else {
				//add whole chunks
				chunks+=  numTQs/size;
				//any remainder?
				int r = numTQs % size;
				if (r != 0) chunks++;
			}
		}
		return chunks;
	}

	/**Returns the next chunk or null.*/
	public synchronized TabixChunk getChunk(){
		if (chunkIterator.hasNext()) return chunkIterator.next();
		return null;
	}

	private void makeTabixChunks() {
		//clear prior
		tabixChunks.clear();

		//walk through files
		for (File f: fileTabixQueries.keySet()){
			ArrayList<TabixQuery> al = fileTabixQueries.get(f);
			int numTQs = al.size();
			if (numTQs <= numberQueriesInChunk) tabixChunks.add( new TabixChunk(f, al));
			else {
				//walk it
				ArrayList<TabixQuery> tqSub = new ArrayList<TabixQuery>();
				int counter = 0;
				for (TabixQuery tq : al){
					tqSub.add(tq);
					//hit max?
					if (tqSub.size() == numberQueriesInChunk) {
						tabixChunks.add(new TabixChunk(f, tqSub));
						counter = 0;
						tqSub = new ArrayList<TabixQuery>();
					}
				}
				//add last?
				if (counter !=0) tabixChunks.add(new TabixChunk(f, tqSub));
			}
		}

		//create iterator for Loaders to pull from
		chunkIterator = tabixChunks.iterator();

		if (printStats) System.err.println( tabixChunks.size()+"\tQuery chunks created");
	}

	/**This takes File sources that have been previously identified to contain overlapping TabixQuery objects and performs a tabix lookup to load the associated data into each TabixQuery.
	 * Its threaded for speed since the tabix look up is relatively slow. */
	private void loadTabixQueriesWithData() throws IOException {
		long startTime = System.currentTimeMillis();

		//make TabixChunks
		makeTabixChunks();

		//make one loader per thread
		TabixLoader[] loader = new TabixLoader[numberThreads];
		ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
		for (int i=0; i< loader.length; i++){
			loader[i] = new TabixLoader(this);
			executor.execute(loader[i]);
		}
		executor.shutdown();

		//spins here until the executer is terminated, e.g. all threads complete
		while (!executor.isTerminated()) {}

		//check loaders and tabulate results
		long numFetchedDataLines = 0;
		int numQueriesWithDataLines = 0;
		for (TabixLoader c: loader) {
			if (c.isFailed()) throw new IOException("\nERROR: TabixLoader issue! \n"+c);
			numFetchedDataLines+= c.getNumberRetrievedResults();
			numQueriesWithDataLines+= c.getNumberQueriesWithResults();
		}

		if (printStats) {
			System.err.println( (System.currentTimeMillis() -startTime)+"\tMillisec to load queries with data");
			System.err.println(numQueriesWithDataLines+"\tQueries found with tabix records");
			System.err.println(numFetchedDataLines+"\tRecords loaded");
		}
	}



	public void buildEngine() throws IOException {
		createChromIndex();
		filter = new TQueryFilter(this);
		System.err.println("\nParsing data sources...");
		System.err.println("\nFileName\tRecordsParsed\tRecordsSkipped");
		loadChromIndexWithVcf();
		loadChromIndexWithBed();
	}

	/**This takes a bed file, parses each region into a TabixQuery object and identifies which file data sources contain data that intersects it.*/
	public boolean intersectFileRegionsWithIndex(File inputFile) throws IOException {

		String name = inputFile.getName();
		if (name.endsWith(".bed.gz") || name.endsWith(".bed")){

			//load bed file of regions they want to intersect
			HashMap<String, RegionScoreText[]> chrRegions = Bed.parseBedFile(inputFile, true, false);

			//are scores set in the bed?
			scoresSet = Bed.scoresSet(inputFile);

			//convert them to TabixQuery and set in this
			chrTabixQueries = convert2TabixQuery(chrRegions);
		}

		else if (name.endsWith(".vcf.gz") || name.endsWith(".vcf")){
			chrTabixQueries = convertVcfToTabixQuery(inputFile);
		}

		else {
			if (printWarnings) System.err.println("\nERROR: user input files must be either bed or vcf (.gz OK)");
			return false;
		}
		//perform the file search
		queryFileIndex();	
		return true;
	}

	
	private void queryFileIndex() {
		long startTime = System.currentTimeMillis();
		long numQueries = 0;
		long numSkippedQueries = 0;
		long numQueriesWithHits = 0;
		long numberHits = 0;
		fileTabixQueries.clear();

		//for each chromosome of regions to query
		for (String chr: chrTabixQueries.keySet()){

			//check that chr exists in index
			if (chromFileIndex.containsKey(chr) == false){
				int numSkipped = chrTabixQueries.get(chr).length;
				if (printWarnings) System.err.println("\nWARNING: chromosome '"+chr+"' not found in index? Skipping "+numSkipped+" query regions.");
				numSkippedQueries += numSkipped; 
			}

			else {
				HashSet<File>[] index = chromFileIndex.get(chr);
				//for each region
				TabixQuery[] regions = chrTabixQueries.get(chr);
				numQueries += regions.length;
				for (TabixQuery tq: regions){
					HashSet<File> fileHits = intersect(index, tq);
					if (fileHits.size() !=0) {
						numberHits += fileHits.size();
						numQueriesWithHits++;
					}
					addHits(fileHits, fileTabixQueries, tq);
				}
			}
		}
		if (printStats){
			long diffTime = System.currentTimeMillis() -startTime;
			System.err.println("\nQuery Stats Pre Filtering:");
			System.err.println(numQueries+ "\tNum index queries");
			System.err.println(numSkippedQueries+ "\tNum skipped index queries");
			System.err.println(numQueriesWithHits+ "\tNum queries with hits");
			System.err.println(numberHits+ "\tNum total hits for tabix retrieval");
			System.err.println(diffTime+"\tMillisec to complete index search");
		}
	}

	/**Adds the TabixQuery to an ArrayList associated with a file resource to fetch the data from.*/
	private void addHits(HashSet<File> hits, HashMap<File, ArrayList<TabixQuery>> toQuery, TabixQuery tq) {
		for (File fHit: hits){
			ArrayList<TabixQuery> al = toQuery.get(fHit);
			if (al == null){
				al = new ArrayList<TabixQuery>();
				toQuery.put(fHit, al);
			}
			al.add(tq);
		}
	}

	/**Checks the begin and end for out of bounds then uses a hash to collapse all the Files found to have region that overlaps the query.*/
	private HashSet<File> intersect(HashSet<File>[] index, TabixQuery tq) {
		HashSet<File> hits = new HashSet<File>();
		int begin = tq.getStart();
		if (begin < 0) begin = 0;
		int end = tq.getStop();
		if (end >= index.length) end = index.length;
		for (int i=begin; i<end; i++) {
			if (index[i] != null) hits.addAll(index[i]);
		}
		return hits;
	}

	public HashMap<String, TabixQuery[]> convert2TabixQuery(HashMap<String, RegionScoreText[]> chrRegions) throws IOException {
		HashMap<String, TabixQuery[]> chrTQ = new HashMap<String, TabixQuery[]>();

		for (String chr: chrRegions.keySet()){
			RegionScoreText[] regions = chrRegions.get(chr);
			TabixQuery[] tq = new TabixQuery[regions.length];
			for (int i=0; i< regions.length; i++) {
				if (bpPadding !=0) {
					int start = regions[i].getStart()-bpPadding;
					if (start < 0) start = 0; 
					regions[i].setStart(start);
					regions[i].setStop(regions[i].getStart()+bpPadding);
				}
				tq[i] = new TabixQuery(chr, regions[i]);
			}
			chrTQ.put(chr, tq);
		}
		return chrTQ;
	}

	private void loadChromIndexWithBed() throws IOException {

		for (int i=0; i< bedDataFiles.length; i++){

			//load filter with file, extension, parent dir name
			addFileToFilter(bedDataFiles[i]);

			BufferedReader in = IO.fetchBufferedReader(bedDataFiles[i]);

			String[] t;
			String line;
			String currChrom = "";
			HashSet<File>[] currIndex = null;
			long numLoadedRecords = 0;
			while ((line = in.readLine()) != null){
				if (line.startsWith("#") == false){
					t = Misc.TAB.split(line);
					//diff chrom?
					if (currChrom.equals(t[0]) == false){
						currIndex = chromFileIndex.get(t[0]);
						if (currIndex == null) throw new IOException("\nError: Failed to find a chromosome for '"+t[0]+ "' from "+bedDataFiles[i]+" aborting on line "+line);
						currChrom = t[0];
					}

					//parse start and stop, note interbase coordinates
					int start = Integer.parseInt(t[1]);
					int stop = Integer.parseInt(t[2]);

					//add in references to source file over the covered bases, stop isn't covered.
					for (int j= start; j< stop; j++){
						if (currIndex[j] == null) currIndex[j] = new HashSet<File>(1);
						currIndex[j].add(bedDataFiles[i]);
					}
					numLoadedRecords++;
				}
			}
			//clean up and stat incrementing
			in.close();
			System.err.println(bedDataFiles[i].getParentFile().getName()+ File.separator+ bedDataFiles[i].getName()+ "\t"+ numLoadedRecords+"\t0");
			recordsLoaded += numLoadedRecords;
		}
	}

	public HashMap<String, TabixQuery[]> convertVcfToTabixQuery(File vcf) throws IOException{
		BufferedReader in = IO.fetchBufferedReader(vcf);
		String[] t;
		String line;
		ArrayList<Bed> bedAl = new ArrayList<Bed>();
		while ((line = in.readLine()) != null){
			if (line.startsWith("#") == false){
				t = Misc.TAB.split(line);

				//check chrom
				if (chromFileIndex.containsKey(t[0]) == false){
					if (printWarnings) System.err.println("\tWARNING: Failed to find a chromosome for "+vcf+" skipping "+line);
				}

				else {
					//fetch start and stop for effected bps.
					int[] startStop = fetchEffectedBps(t);
					if (startStop != null) {
						//pad it?
						if (bpPadding !=0) {
							int start = startStop[0]-bpPadding;
							if (start < 0) start = 0; 
							startStop[0] = start;
							startStop[1]+= bpPadding;
						}
						bedAl.add(new Bed(t[0], startStop[0], startStop[1], line.replace("\t", "\\t"), 0, '.'));
					}
				}
			}
		}
		//clean up and stat incrementing
		in.close();

		//covert to TQ
		HashMap<String, TabixQuery[]> tqs = new HashMap<String, TabixQuery[]>();
		Bed[] bed = new Bed[bedAl.size()];
		bedAl.toArray(bed);
		Arrays.sort(bed);
		HashMap<String, Bed[]> bedHash = Bed.splitBedByChrom(bed);
		for (String chr: bedHash.keySet()){
			Bed[] b = bedHash.get(chr);
			TabixQuery[] tq = new TabixQuery[b.length];
			for (int i=0; i< b.length; i++) tq[i] = new TabixQuery(b[i]);
			tqs.put(chr, tq);
		}
		return tqs;
	}


	private void loadChromIndexWithVcf() throws IOException {

		for (int i=0; i< vcfDataFiles.length; i++){

			//load filter with file, extension, parent dir name
			addFileToFilter(vcfDataFiles[i]);

			BufferedReader in = IO.fetchBufferedReader(vcfDataFiles[i]);

			String[] t;
			String line;
			String currChrom = "";
			HashSet<File>[] currIndex = null;
			long numLoadedRecords = 0;
			long numRecordsSkipped = 0;
			while ((line = in.readLine()) != null){
				if (line.startsWith("#") == false){
					t = Misc.TAB.split(line);

					//diff chrom? pull index
					if (currChrom.equals(t[0]) == false){
						currIndex = chromFileIndex.get(t[0]);
						if (currIndex == null) throw new IOException("\nError: Failed to find a chromosome for '"+t[0]+ "' from "+vcfDataFiles[i]+" aborting on line "+line);
						currChrom = t[0];
					}

					//fetch start and stop for effected bps.
					int[] startStop = fetchEffectedBps(t);
					if (startStop == null) numRecordsSkipped++;
					else {
						//add in references to source file over the covered bases, stop isn't covered.
						for (int j= startStop[0]; j< startStop[1]; j++){
							if (currIndex[j] == null) currIndex[j] = new HashSet<File>(1);
							currIndex[j].add(vcfDataFiles[i]);
						}
						numLoadedRecords++;
					}
				}
			}

			//clean up and stat incrementing
			in.close();
			System.err.println(vcfDataFiles[i].getParentFile().getName()+ File.separator+ vcfDataFiles[i].getName()+"\t"+numLoadedRecords+"\t"+numRecordsSkipped);
			recordsLoaded += numLoadedRecords;
			recordsSkipped += numRecordsSkipped;
		}
	}

	/**Loads the filter with the file, extension, and parent dir name*/
	private void addFileToFilter(File file) {
		//add to data sources
		filter.getAvailableDataFiles().add(file);

		//add parent dir name
		filter.getAvailableDataTypes().add(file.getParentFile().getName());

		//add extension minus .gz
		Matcher mat = TQueryFilter.FILE_TYPE.matcher(file.getName());
		mat.matches();
		filter.getAvailableFileTypes().add(mat.group(1));
	}

	/**Returns the interbase start stop region of effected bps for simple SNV and INDELs. 
	 * SNV=iPos,iPos+LenRef; INS=iPos,iPos+lenRef+1; DEL=iPos+1,iPos+lenRef; iPos=pos-1.
	 * For multi alts, returns min begin and max end of all combinations.
	 * For alts with < indicating a CNV or trans, attempts to parse the END=number from the INFO column. */
	public int[] fetchEffectedBps(String[] vcf) throws IOException{
		//CHROM	POS	ID	REF	ALT	
		//  0    1   2   3   4  
		//put into interbase coordinates
		int iPos = Integer.parseInt(vcf[1]) - 1;
		String ref= vcf[3];
		String alt= vcf[4];

		//any commas/ multi alts?
		if (alt.contains(",") == false) return fetchEffectedBpsNoMultiAlt(iPos, ref, alt, vcf);

		//OK commas present, these need to be tested for max effect
		//There is complexity with multi alts, best to deconvolute and left justify! 
		String[] alts = Misc.COMMA.split(alt);
		int begin = Integer.MAX_VALUE;
		int end = -1;
		for (int i=0; i< alts.length; i++){
			int[] ss = fetchEffectedBpsNoMultiAlt(iPos, ref, alts[i], vcf);
			if (ss == null) return null;
			if(ss[0] < begin) begin = ss[0];
			if(ss[1]> end) end = ss[1];
		}

		return new int[]{begin, end};
	}

	private int[] fetchEffectedBpsNoMultiAlt(int iPos, String ref, String alt, String[] vcf) throws IOException{
		int begin = -1;
		int end = -1;
		int lenRef = ref.length();
		int lenAlt = alt.length();

		//watch out for < in the alt indicative of a CNV or structural var
		if (alt.contains("<")){
			//CHROM	POS	ID	REF	ALT QUAL FILTER INFO	
			//  0    1   2   3   4   5      6     7
			Matcher mat = END_POSITION.matcher(vcf[7]);
			if (mat.matches()) end = Integer.parseInt(mat.group(1));
			else {
				if (printWarnings) System.err.println("\n\tWARNING: found a < or > containing alt, failed to parse END=number position, skipping -> "+Misc.stringArrayToString(vcf, "\t"));
				return null;
			}
			begin = iPos;
		}

		//single or multi adjacent snp? return just the changed bps,  GC->AT or G->A
		else if (lenRef == lenAlt) {
			begin = iPos;
			end = iPos+ lenRef;
		}
		//ins? return the bases on either side of the insertion GC->GATTA or G->ATTA
		else if (lenAlt > lenRef) {
			begin = iPos;
			end = iPos+ lenRef +1;
			if (ref.charAt(0) != alt.charAt(0)) {
				if (printWarnings) System.err.println("\n\tWARNING: Odd INS vcf record, the first base in the ref and alt must be the same, use vt to normalize your variants, skipping -> "+Misc.stringArrayToString(vcf, "\t"));
				return null;
			}

		}
		//del? return the bps that are deleted, AT->A, ATTCG->ACC
		else if (lenRef > lenAlt) {
			begin = iPos+1;
			end = iPos + lenRef;
			if (ref.charAt(0) != alt.charAt(0)) {
				if (printWarnings) System.err.println("\n\tWARNING: Odd DEL vcf record, the first base in the ref and alt must be the same, use vt to normalize your variants, skipping -> "+Misc.stringArrayToString(vcf, "\t"));
				return null;
			}
		}
		//odd, shouldn't hit this
		else throw new IOException("\nError: Contact admin! Odd vcf record, can't parse effected bps for -> "+Misc.stringArrayToString(vcf, "\t"));

		return new int[]{begin, end};
	}


	private void createChromIndex() throws IOException {
		HashMap<String, RegionScoreText[]> chrLen = Bed.parseBedFile(chrLengthFile, true, false);

		System.err.println("Building empty chrom index...");
		for (String chr: chrLen.keySet()){
			RegionScoreText[] regions = chrLen.get(chr);
			if (regions.length !=1) throw new IOException("\nError: there can be only one bed region for each chromosome, see "+chr);
			HashSet<File>[] indexes = new HashSet[regions[0].getStop()+1];
			chromFileIndex.put(chr, indexes);
			System.err.println("\t"+chr+"\t"+regions[0].getStop());
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new TQuery (args);
	}	

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		System.err.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		File tabixDataDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': chrLengthFile = new File(args[++i]); break;
					case 'd': tabixDataDir = new File(args[++i]); break;
					case 'n': numberThreads = Integer.parseInt(args[++i]); break;
					case 'p': bpPadding = Integer.parseInt(args[++i]); break;
					case 'i': printWarnings = false; break;
					case 's': printStats = false; break;
					case 'f': fetchData = false; break;
					case 'q': this.numberQueriesInChunk = Integer.parseInt(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (chrLengthFile == null) Misc.printErrAndExit("\nError: please provide a bed file of chromosome and their lengths, e.g. X 0 155270560\n" );
		if (tabixDataDir == null || tabixDataDir.isDirectory() == false) Misc.printErrAndExit("\nError: please provide a directory containing tabix indexed xxx.vcf.gz and xxx.bed.gz files with their associated xxx.gz.tbi indexes" );

		vcfDataFiles = IO.fetchFilesRecursively(tabixDataDir, "vcf.gz");
		bedDataFiles = IO.fetchFilesRecursively(tabixDataDir, "bed.gz");
		if (vcfDataFiles.length == 0 && bedDataFiles.length == 0) Misc.printErrAndExit("\nError: failed to find any xxx.bed.gz or xxx.vcf.gz tabix files in your tabixDataDir -> "+tabixDataDir);
		//check for index
		lookForTabixIndex(vcfDataFiles);
		lookForTabixIndex(bedDataFiles);


		//threads to use
		int numAvail = Runtime.getRuntime().availableProcessors();
		if (numberThreads < 1) numberThreads =  numAvail - 1;
		System.err.println(numAvail +" Available processors, using "+numberThreads+" threaded loaders\n");


	}	

	private void lookForTabixIndex(File[] tabixFiles) {
		ArrayList<String> badFiles = new ArrayList<String>();
		for (File tb: tabixFiles){
			File index = new File (tb.toString()+".tbi");
			if (index.exists() == false) badFiles.add(tb.toString());
		}
		if (badFiles.size() !=0) Misc.printErrAndExit("\nError: the following files are missing their xxx.gz.tbi Tabix indexes?\n"+ Misc.stringArrayListToString(badFiles, "\n"));
	}

	public static void printDocs(){
		System.err.println("\n" +
				"**************************************************************************************\n" +
				"**                                    TQuery: Sept 2016                             **\n" +
				"**************************************************************************************\n" +
				"TQ returns bed and vcf records that overlap user supplied regions in bed or vcf format."+
				"Strict adherence to bed format is assumed (1st base is 0, last base is not included,\n"+
				"last base > first; vcf is in 1 base so subtract 1 from pos to convert to bed).\n"+
				"This app needs > 27G of RAM for human/ mouse/ plant. A two step search is performed\n"+
				"using an in memory file : bp index to find intersecting regions followed by a tabix\n"+
				"query to pull the data.  Multiple data source filters can be applied to limit output.\n"+

				"\nRequired Params:\n"+
				"-c A bed file of chromosomes and their lengths (e.g. chr21 0 48129895) to use to \n"+
				"     building the intersection index. Note, the chr name must match across all \n"+
				"     datasets, no mixing of chr21 and 21.\n"+
				"-d A data directory containing gzipped tabixed (https://github.com/samtools/htslib)\n"+
				"     vcf and bed files. Recurses through all sub directories. Be sure to normalize\n"+
				"     the vcf records with a package like Vt, see http://genome.sph.umich.edu/wiki/Vt\n"+

				"\nDefault Params:\n"+
				"-p Pad input bed/vcf regions +/- x bps, defaults to 0.\n"+
				"-i Silence warnings\n"+
				"-s Don't print query statistics\n"+
				"-f Don't fetch records, just print the files that intersect each region, very fast.\n"+
				"-n Number of processors to use, defaults to all-1.\n"+
				"-q Number of queries in each lookup chunk, defaults to 1000\n"+

				"\nExample: java -Xmx64G -jar pathToUSeq/Apps/TQuery -c b37ChrLen.bed \n"+
				"     -t TabixDataFiles/ | gzip > results.json.gz  \n\n" +

				"**************************************************************************************\n");
	}

	public boolean isPrintWarnings() {
		return printWarnings;
	}

	public boolean isFetchData() {
		return fetchData;
	}

	public boolean isPrintJson() {
		return printJson;
	}


}






