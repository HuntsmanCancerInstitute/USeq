package edu.utah.seq.query;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;
import edu.utah.seq.useq.data.RegionScoreText;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class QueryRequest {

	//fields
	private TQueryFilter filter; 
	private TQuery tQuery;
	private DataSources dataSources;
	private QueryIndex queryIndex;
	private ArrayList<TabixChunk> tabixChunks = new ArrayList<TabixChunk>();
	private Iterator<TabixChunk> chunkIterator = null;
	private boolean printWarnings;
	private long numberRetrievedResults = 0;
	private int numberQueriesWithResults = 0;
	private long numberSavedResults = 0;
	private long time2Load;

	/*Contains all the TabixQueries split by chromosome. Some may contain results.*/
	private HashMap<String, TabixQuery[]> chrTabixQueries = new HashMap<String, TabixQuery[]>();
	
	/*Contains a File pointer to a data source and the associated TabixQueries that intersect it.
	 * These are to be used in a tabix search.*/
	private HashMap<File, ArrayList<TabixQuery>> fileTabixQueries = new HashMap<File, ArrayList<TabixQuery>>();

	
	public QueryRequest(String in, TQuery tQuery) throws IOException {
		this.tQuery = tQuery;
		printWarnings = tQuery.isPrintWarnings();
		this.dataSources = tQuery.getDataSources();
		this.queryIndex = tQuery.getQueryIndex();
		
		//create a new file filter and modify it with the input
		filter = new TQueryFilter();
		
		//parse the input, if no file do nothing
		File inputFile = parseCmdLineInput(in);
		if (inputFile == null) return;

		//find dataSource files that intersect regions in the input file bed or vcf, will return false if bed or vcf not found
		if (intersectFileRegionsWithIndex(inputFile)){

			//filter which files to fetch data from base on the parent dir name, the actual file, the file extension 
			filter.filter(fileTabixQueries);

			//any left?
			if (fileTabixQueries.size()!=0){
				//fetch the underlying data (slow) or just indicate which file hit which region (fast)
				if (filter.isFetchData()) {
					makeTabixChunks();
					tQuery.getQueryLoader().loadTabixQueriesWithData(this);
					//print the results
					if (tQuery.isPrintStats()) {
						System.err.println(time2Load+"\tMillisec to load queries with data");
						System.err.println(numberQueriesWithResults+"\tQueries found with tabix records");
						System.err.println(numberRetrievedResults+"\tRecords retrieved");
						System.err.println(numberSavedResults+"\tRecords returned");
					}
				}
				else loadTabixQueriesWithFileSources();

				//print results to stdout in json format
				if (filter.isPrintJson()) printTabixQueriesWithResults();
			}
			else {
				System.err.println("No intersecting files after filtering!");
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
				if (lcArg.contains("filter")){
					System.err.println("\nCurrent Filters:\n"+filter.getCurrentFilters("\n"));
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
				if (bool.startsWith("t")) filter.setFetchData(true);
				else filter.setFetchData(false);
			}
			//print results
			else if (key.equals("printjson")){
				String bool = options.get("printjson").toLowerCase();
				if (bool.startsWith("t")) filter.setPrintJson(true);
				else filter.setPrintJson(false);
			}
			//print results
			else if (key.equals("matchvcf")){
				String bool = options.get("matchvcf").toLowerCase();
				if (bool.startsWith("t")) filter.setMatchVcf(true);
				else filter.setMatchVcf(false);
			}
			//filter on primaryDataType? parent
			else if (key.equals("primarydatatypes")){
				String[] dataTypes = Misc.COMMA.split(options.get("primarydatatypes"));
				TreeSet<String> dtToReturn = new TreeSet<String>();
				//for each
				for (String dt: dataTypes){
					if (dataSources.getAvailablePrimaryDataTypes().contains(dt)) dtToReturn.add(dt);
					else {
						System.err.println("\nHmm can't find this data type?! -> "+dt+ " in "+dataSources.fetchPrimaryDataTypes(", "));
						return null;
					}
				}
				filter.setPrimaryDataTypesToReturn(dtToReturn);
				filter.setFilterOnPrimaryDataTypes(true);
			}
			//filter on secondaryDataType? grand parent
			else if (key.equals("secondarydatatypes")){
				String[] dataTypes = Misc.COMMA.split(options.get("secondarydatatypes"));
				TreeSet<String> dtToReturn = new TreeSet<String>();
				//for each
				for (String dt: dataTypes){
					if (dataSources.getAvailableSecondaryDataTypes().contains(dt)) dtToReturn.add(dt);
					else {
						System.err.println("\nHmm can't find this data type?! -> "+dt+ " in "+dataSources.fetchSecondaryDataTypes(", "));
						return null;
					}
				}
				filter.setSecondaryDataTypesToReturn(dtToReturn);
				filter.setFilterOnSecondaryDataTypes(true);
			}
			//filter on fileType?
			else if (key.equals("filetypes")){
				String[] fileTypes = Misc.COMMA.split(options.get("filetypes"));
				TreeSet<String> toReturn = new TreeSet<String>();
				//for each
				for (String dt: fileTypes){
					if (dataSources.getAvailableFileTypes().contains(dt)) toReturn.add(dt);
					else {
						System.err.println("\nHmm can't find this file type?! -> "+dt+ " in "+dataSources.fetchFileTypes(", "));
						return null;
					}
				}
				filter.setFileTypesToReturn(toReturn);
				filter.setFilterOnFileTypes(true);
			}
			//filter on data source/ actual files?
			else if (key.equals("datasources")){
				String[] dataFiles = Misc.COMMA.split(options.get("datasources"));
				TreeSet<File> toReturn = new TreeSet<File>();
				//for each
				for (String fs: dataFiles){
					File f = new File(dataSources.getPathToTrimmedFile()+fs);
					if (dataSources.getAvailableDataFiles().contains(f)) toReturn.add(f);
					else {
						System.err.println("\nHmm can't find this data soure file?! -> "+fs+ " in "+dataSources.fetchDataFilesRelative(", "));
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
	
	/**This takes a bed file, parses each region into a TabixQuery object and identifies which file data sources contain data that intersects it.*/
	public boolean intersectFileRegionsWithIndex(File inputFile) throws IOException {

		String name = inputFile.getName();
		if (name.endsWith(".bed.gz") || name.endsWith(".bed")){

			//load bed file of regions they want to intersect
			HashMap<String, RegionScoreText[]> chrRegions = Bed.parseBedFile(inputFile, true, false);

			//are scores set in the bed?
			filter.setScoresSet(Bed.scoresSet(inputFile));

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
		queryIndex.queryFileIndex(this);	
		return true;
	}
	
	public HashMap<String, TabixQuery[]> convert2TabixQuery(HashMap<String, RegionScoreText[]> chrRegions) throws IOException {
		HashMap<String, TabixQuery[]> chrTQ = new HashMap<String, TabixQuery[]>();
		int bpPadding = filter.getBpPadding();
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
	
	private void loadTabixQueriesWithFileSources() {
		//for each data source file
		for (File tabixFile: fileTabixQueries.keySet()){
			//pull queries that intersect it
			ArrayList<TabixQuery> toFetch = fileTabixQueries.get(tabixFile);
			//for each query, add the data source
			for (TabixQuery tq: toFetch) tq.getSourceResults().put(tabixFile, null);
		}
	}
	
	public HashMap<String, TabixQuery[]> convertVcfToTabixQuery(File vcf) throws IOException{
		HashMap<String, HashSet<File>[]> chromFileIndex = queryIndex.getChromFileIndex();
		int bpPadding = filter.getBpPadding();
		BufferedReader in = IO.fetchBufferedReader(vcf);
		String[] t;
		String line;
		ArrayList<Bed> bedAl = new ArrayList<Bed>();
		while ((line = in.readLine()) != null){
			if (line.startsWith("#") == false){
				t = Misc.TAB.split(line);

				//check chrom
				if (chromFileIndex.containsKey(t[0]) == false){
					if (printWarnings) System.err.println("\tWARNING: Failed to find a chromosome for '"+t[0]+ "' skipping -> "+line);
				}

				else {
					//fetch start and stop for effected bps.
					int[] startStop = QueryIndex.fetchEffectedBps(t, printWarnings);
					if (startStop != null) {
						//pad it?
						if (bpPadding !=0) {
							int start = startStop[0]-bpPadding;
							if (start < 0) start = 0; 
							startStop[0] = start;
							startStop[1]+= bpPadding;
						}
						bedAl.add(new Bed(t[0], startStop[0], startStop[1], line, 0, '.'));
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
			for (int i=0; i< b.length; i++) {
				tq[i] = new TabixQuery(b[i]);
				tq[i].parseVcf();
			}
			tqs.put(chr, tq);
		}
		return tqs;
	}


	private void printAllTabixQueries() {
		String pathToTrimmed = dataSources.getPathToTrimmedFile();
		boolean scoresSet = filter.isScoresSet();
		for (String chr: chrTabixQueries.keySet()){
			TabixQuery[] tqs = chrTabixQueries.get(chr);
			if (tqs.length > 1) {
				System.out.println("[");
				System.out.print(tqs[0].toJson(scoresSet, pathToTrimmed));
				for (int i=1; i< tqs.length; i++){
					System.out.println(",");
					System.out.print(tqs[i].toJson(scoresSet, pathToTrimmed));
				}
				System.out.println("\n]");
			}
			else System.out.println(tqs[0].toJson(scoresSet, pathToTrimmed));
		}
		//flush the stream to while waiting for another query set
		System.out.flush();
	}
	
	private void printTabixQueriesWithResults() {
		String pathToTrimmed = dataSources.getPathToTrimmedFile();
		boolean scoresSet = filter.isScoresSet();
		
		for (String chr: chrTabixQueries.keySet()){
			TabixQuery[] tqs = chrTabixQueries.get(chr);
			ArrayList<String> res = new ArrayList<String>();
			for (int i=0; i< tqs.length; i++){
				if (tqs[i].getSourceResults().size() !=0) res.add(tqs[i].toJson(scoresSet, pathToTrimmed));
			}
			int num = res.size();
			if (num!=0) {
				System.out.println("[");
				System.out.print(res.get(0));
				for (int i=1; i< num; i++){
					System.out.println(",");
					System.out.print(res.get(i));
				}
				System.out.println("\n]");
			}
		}
		//flush the stream to while waiting for another query set
		System.out.flush();
	}

	/*
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
	}*/

	/**Returns the next chunk or null.*/
	public TabixChunk getChunk(){
		if (chunkIterator.hasNext()) return chunkIterator.next();
		return null;
	}

	private void makeTabixChunks() {
		int numberQueriesInChunk = tQuery.getNumberQueriesInChunk();

		//walk through files
		for (File f: fileTabixQueries.keySet()){
			ArrayList<TabixQuery> al = fileTabixQueries.get(f);
			int numTQs = al.size();			

			if (numTQs <= numberQueriesInChunk) tabixChunks.add( new TabixChunk(f, al, this));
			else {
				//walk it
				ArrayList<TabixQuery> tqSub = new ArrayList<TabixQuery>();
				for (TabixQuery tq : al){
					tqSub.add(tq);
					//hit max?
					if (tqSub.size() == numberQueriesInChunk) {
						tabixChunks.add(new TabixChunk(f, tqSub, this));
						tqSub = new ArrayList<TabixQuery>();
					}
				}
				//add last?
				if (tqSub.size() !=0) tabixChunks.add(new TabixChunk(f, tqSub, this));
			}
		}

		//create iterator for Loaders to pull from
		chunkIterator = tabixChunks.iterator();

		if (tQuery.isPrintStats()) System.err.println( tabixChunks.size()+"\tQuery chunks created");
	}
	
	public synchronized void incrementNumRetrievedResults(long n){
		numberRetrievedResults += n;
	}
	public synchronized void incrementNumQueriesWithResults(int n){
		numberQueriesWithResults += n;
	}
	public synchronized void incrementNumSavedResults(long n){
		numberSavedResults += n;
	}
	public TQueryFilter getFilter() {
		return filter;
	}
	public TQuery gettQuery() {
		return tQuery;
	}
	public ArrayList<TabixChunk> getTabixChunks() {
		return tabixChunks;
	}
	public Iterator<TabixChunk> getChunkIterator() {
		return chunkIterator;
	}
	public HashMap<String, TabixQuery[]> getChrTabixQueries() {
		return chrTabixQueries;
	}
	public HashMap<File, ArrayList<TabixQuery>> getFileTabixQueries() {
		return fileTabixQueries;
	}
	public long getNumberRetrievedResults() {
		return numberRetrievedResults;
	}
	public int getNumberQueriesWithResults() {
		return numberQueriesWithResults;
	}
	public long getNumberSavedResults() {
		return numberSavedResults;
	}
	public void setTime2Load(long l) {
		time2Load = l;
	}

}
