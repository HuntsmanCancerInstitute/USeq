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
import edu.utah.seq.vcf.GatkRunnerChunk;
import htsjdk.tribble.readers.TabixReader;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;


/** Test of using a big genome index and tabix lookup for a variant store.
 * 23669MB to load empty gIndex. 
 * 25526MB to load cosmic, 3 10K vcfs, dbSNP into gIndex.
 * 790ms to int entire exome with gIndex.
 * 116,745ms to fetch overlap vars for entire exome (257,726 regions intersected  2,221,581 records)
 * 15,855 ms for 14111 short vars pulling 13841 Records
 * So 15.8 sec to fetch all intersecting records with 14K T-N exome variants using Tabix readers.  Hmm the tabix retrieval is slow, need to parallelize! Doubt could get this to < 2 Sec.
 * Better to use a NoSQL db, key: data.
 * */
public class TQuery {

	//fields
	private File chrLenBedFile;
	private File[] vcfFiles;
	private File[] bedFiles;
	private Gzipper results;
	private int numberThreads = 0;
	private boolean printWarnings = true;
	private boolean printStats = true;
	private boolean printFindings = false;

	//internal
	private long recordsLoaded = 0;
	private long recordsSkipped = 0;

	//interbase coordinates, 0 is first base, last base in any region is not included.
	private HashMap<String, HashSet<File>[]> chromFileIndex = new HashMap<String, HashSet<File>[]>();

	//per group query search

	/*Contains a File pointer to a datasource and the associatted TabixQueries that intersect it.
	 * These are to be used in a tabix search.*/
	private HashMap<File, ArrayList<TabixQuery>> fileTabixQueries = new HashMap<File, ArrayList<TabixQuery>>();

	/*Contains all the TabixQueries split by chromosome. Some may contain results*/
	private HashMap<String, TabixQuery[]> chrTabixQueries = new HashMap<String, TabixQuery[]>();

	public TQuery (String[] args) {
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			buildEngine();

			String diffTime = Num.formatNumberOneFraction(((double)(System.currentTimeMillis() -startTime))/1000);
			System.out.println("\n"+ diffTime+" Sec to build using "+IO.memory()+ " of RAM");
			System.out.println("\t"+recordsLoaded+"\tRecords indexed");
			System.out.println("\t"+recordsSkipped+"\tRecords skipped");

			queryBedFilesFromCmdLine();


		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("\nProblem with executing the TQuery!");
		}
	}

	private void queryBedFilesFromCmdLine() throws IOException {
		while (true){
			System.out.print("\nProvide a bed file of regions to query against the engine (or blank to exit):\n");
			String bedFileString = (new BufferedReader(new InputStreamReader(System.in))).readLine().trim();
			if (bedFileString == null || bedFileString.length()==0) {
				System.out.println();
				break;
			}
			File bed = new File (bedFileString);
			if (bed.exists() == false) {
				System.err.println("\nHmm can't seem to find that file?!");
				continue;
			}

			intersectBedFileRegionsWithIndex(bed);

			//filter which files to fetch data from, then retrieve it
			//filterFiles(hits);

			loadTabixQueriesWithData();

			if (printFindings){
				printTabixQueriesWithHitsBySource();
				printAllTabixQueries();
			}
		}
	}

	private synchronized void appendResults(String res) {
		try {
			results.println(res);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void printAllTabixQueries() {
		System.out.println("\nPrinting All Queries:");
		for (String chr: chrTabixQueries.keySet()){
			TabixQuery[] tqs = chrTabixQueries.get(chr);
			for (TabixQuery tq : tqs){
				System.out.println(tq.getInterbaseCoordinates());
				HashMap<File, ArrayList<String>> sourceRes = tq.getSourceResults();
				for (File s: sourceRes.keySet()){
					System.out.println("\t"+s);
					for (String r: sourceRes.get(s)){
						System.out.println("\t\t"+r);
					}
				}
			}
		}
	}

	private void printTabixQueriesWithHitsBySource() {
		System.out.println("\nPrinting Queries With Hits:");
		for (File source: fileTabixQueries.keySet()){
			System.out.println(source);
			ArrayList<TabixQuery> al = fileTabixQueries.get(source);
			for (TabixQuery tq: al){
				System.out.println("\t"+tq.getInterbaseCoordinates());
				ArrayList<String> res = tq.getSourceResults().get(source);
				for (String s: res) System.out.println("\t\t"+s);
			}
			System.out.println();
		}
	}

	/**This takes File sources that have been previously identified to contain overlapping TabixQuery objects and performs a tabix lookup to load the associated data into each TabixQuery.
	 * Its threaded for speed since the tabix look up is relatively slow. */
	private void loadTabixQueriesWithData() throws IOException {

		//TODO NEED to optimize better use of threads when one file source has many TQs to look up. Might also help to create a pool of TabixReaders instead of instantiating them for every query.

		long startTime = System.currentTimeMillis();
		//dumb division, make one loader per file
		TabixLoader[] loader = new TabixLoader[fileTabixQueries.size()];
		ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
		int index = 0;
		for (File tabixFile: fileTabixQueries.keySet()){
			ArrayList<TabixQuery> toFetch = fileTabixQueries.get(tabixFile);
			loader[index] = new TabixLoader(tabixFile, toFetch, printWarnings);
			executor.execute(loader[index++]);
		}
		executor.shutdown();
		
		//spins here until the executer is terminated, e.g. all threads complete
		while (!executor.isTerminated()) {
		}
		//check loaders
		for (TabixLoader c: loader) {
			if (c.isFailed()) throw new IOException("\nERROR: Failed TabixLoader, aborting! \n"+c);
		}

		if (printStats) System.out.println( (System.currentTimeMillis() -startTime)+"\tMillisec to load queries with data");
	}

	public void buildEngine() throws IOException {
		createChromIndex();
		loadChromIndexWithVcf();
		loadChromIndexWithBed();
	}

	/**This takes a bed file, parses each region into a TabixQuery object and identifies which file data sources contain data that intersects it.*/
	public void intersectBedFileRegionsWithIndex(File bed) throws IOException {
		//load bed file of regions they want to intersect
		HashMap<String, RegionScoreText[]> chrRegions = Bed.parseBedFile(bed, true, false);

		//convert them to TabixQuery and set in this
		chrTabixQueries = convert2TabixQuery(chrRegions);

		//perform the file search
		queryFileIndex();	
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
			System.out.println("\nQuery Stats:");
			System.out.println(numQueries+ "\tNum index queries");
			System.out.println(numSkippedQueries+ "\tNum skipped index queries");
			System.out.println(numQueriesWithHits+ "\tNum queries with hits");
			System.out.println(numberHits+ "\tNum total hits for tabix retrieval");
			System.out.println(diffTime+"\tMillisec to complete index search");
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

	public static HashMap<String, TabixQuery[]> convert2TabixQuery(HashMap<String, RegionScoreText[]> chrRegions) throws IOException {
		HashMap<String, TabixQuery[]> chrTQ = new HashMap<String, TabixQuery[]>();
		for (String chr: chrRegions.keySet()){
			RegionScoreText[] regions = chrRegions.get(chr);
			TabixQuery[] tq = new TabixQuery[regions.length];
			for (int i=0; i< regions.length; i++) tq[i] = new TabixQuery(chr, regions[i]);
			chrTQ.put(chr, tq);
		}
		return chrTQ;
	}

	private void loadChromIndexWithBed() throws IOException {

		for (int i=0; i< bedFiles.length; i++){
			BufferedReader in = IO.fetchBufferedReader(bedFiles[i]);
			System.out.println("Loading bed "+bedFiles[i].getName());
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
						if (currIndex == null) throw new IOException("\nError: Failed to find a chromosome for '"+t[0]+ "' from "+bedFiles[i]+" aborting on line "+line);
						currChrom = t[0];
					}

					//parse start and stop, note interbase coordinates
					int start = Integer.parseInt(t[1]);
					int stop = Integer.parseInt(t[2]);

					//add in references to source file over the covered bases, stop isn't covered.
					for (int j= start; j< stop; j++){
						if (currIndex[j] == null) currIndex[j] = new HashSet<File>(1);
						currIndex[j].add(bedFiles[i]);
					}
					numLoadedRecords++;
				}
			}
			//clean up and stat incrementing
			in.close();
			System.out.println("\t"+numLoadedRecords+" regions indexed");
			recordsLoaded += numLoadedRecords;
		}
	}


	private void loadChromIndexWithVcf() throws IOException {

		for (int i=0; i< vcfFiles.length; i++){
			BufferedReader in = IO.fetchBufferedReader(vcfFiles[i]);
			System.out.println("\nParsing vcf "+vcfFiles[i].getName());

			String[] t;
			String line;
			String currChrom = "";
			HashSet<File>[] currIndex = null;
			long numLoadedRecords = 0;
			long numRecordsSkipped = 0;
			while ((line = in.readLine()) != null){
				if (line.startsWith("#") == false){
					t = Misc.TAB.split(line);

					//diff chrom?
					if (currChrom.equals(t[0]) == false){
						currIndex = chromFileIndex.get(t[0]);
						if (currIndex == null) throw new IOException("\nError: Failed to find a chromosome for '"+t[0]+ "' from "+vcfFiles[i]+" aborting on line "+line);
						currChrom = t[0];
					}
					//does ALT have any <> characters indicating it is not a simple snv or indel? if so skip
					// TODO:  should parse out the end positions and not skip!
					if (t[4].contains("<")){
						if (printWarnings) System.err.println("\tWARNING: Only SNV/INDEL variants supported, skipping -> "+line);
						numRecordsSkipped++;
					}
					else {
						//fetch start and stop for effected bps.
						int[] startStop = fetchEffectedBps(t);
						if (startStop == null) numRecordsSkipped++;
						else {
							//add in references to source file over the covered bases, stop isn't covered.
							for (int j= startStop[0]; j< startStop[1]; j++){
								if (currIndex[j] == null) currIndex[j] = new HashSet<File>(1);
								currIndex[j].add(vcfFiles[i]);
							}
							numLoadedRecords++;
						}
					}
				}
			}
			//clean up and stat incrementing
			in.close();
			System.out.println("\t"+numLoadedRecords+" variants indexed, "+numRecordsSkipped+" skipped");
			recordsLoaded += numLoadedRecords;
			recordsSkipped += numRecordsSkipped;
		}
	}

	/**Returns the interbase start stop region of effected bps for simple SNV and INDELs. 
	 * SNV=iPos,iPos+LenRef; INS=iPos,iPos+lenRef+1; DEL=iPos+1,iPos+lenRef; iPos=pos-1.
	 * For multi alts, returns min begin and max end of all combinations.
	 * @throws an exception if alts with < or cases where the first base differ between ref and alt for INDELS. */
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
		//TODO: implement!
		if (alt.contains("<")){
			throw new IOException("\nError: < or > containing alts not parsed, please remove them your "
					+ "vcf files, offending record -> "+Misc.stringArrayToString(vcf, "\t"));
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
				if (printWarnings) System.err.println("\tWARNING: Odd INS vcf record, the first base in the ref and alt should be the same, skipping, see -> "+Misc.stringArrayToString(vcf, "\t"));
				return null;
			}

		}
		//del? return the bps that are deleted, AT->A, ATTCG->ACC
		else if (lenRef > lenAlt) {
			begin = iPos+1;
			end = iPos + lenRef;
			if (ref.charAt(0) != alt.charAt(0)) {
				if (printWarnings) System.err.println("\tWARNING: Odd DEL vcf record, the first base in the ref and alt should be the same, skipping, see -> "+Misc.stringArrayToString(vcf, "\t"));
				return null;
			}
		}
		//odd, shouldn't hit this
		else throw new IOException("\nError: Contact admin! Odd vcf record, can't parse effected bps for -> "+Misc.stringArrayToString(vcf, "\t"));

		return new int[]{begin, end};
	}


	private void createChromIndex() throws IOException {
		HashMap<String, RegionScoreText[]> chrLen = Bed.parseBedFile(chrLenBedFile, true, false);

		System.out.println("Building empty chrom index...");
		for (String chr: chrLen.keySet()){
			RegionScoreText[] regions = chrLen.get(chr);
			if (regions.length !=1) throw new IOException("\nError: there can be only one bed region for each chromosome, see "+chr);
			HashSet<File>[] indexes = new HashSet[regions[0].getStop()+1];
			chromFileIndex.put(chr, indexes);
			System.out.println("\t"+chr+"\t"+regions[0].getStop());
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
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		File tabixDataDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': chrLenBedFile = new File(args[++i]); break;
					case 'd': tabixDataDir = new File(args[++i]); break;
					case 'n': numberThreads = Integer.parseInt(args[++i]); break;
					case 'i': printWarnings = false; break;
					case 's': printStats = false; break;
					case 'p': printFindings = true; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (chrLenBedFile == null) Misc.printErrAndExit("\nError: please provide a bed file of chromosome and their lengths, e.g. X 0 155270560\n" );
		if (tabixDataDir == null || tabixDataDir.isDirectory() == false) Misc.printErrAndExit("\nError: please provide a directory containing tabix indexed xxx.vcf.gz and xxx.bed.gz files with their associated xxx.gz.tbi indexes" );

		vcfFiles = IO.fetchFilesRecursively(tabixDataDir, "vcf.gz");
		bedFiles = IO.extractFiles(tabixDataDir, "bed.gz");
		if (vcfFiles.length == 0 && bedFiles.length == 0) Misc.printErrAndExit("\nError: failed to find any xxx.bed.gz or xxx.vcf.gz tabix files in your tabixDataDir -> "+tabixDataDir);
		
		//threads to use
		int numAvail = Runtime.getRuntime().availableProcessors();
		if (numberThreads < 1) numberThreads =  numAvail - 1;
		System.out.println(numAvail +" Available processors, using "+numberThreads+"\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                    TQuery: Sept 2016                             **\n" +
				"**************************************************************************************\n" +
				"TQ returns bed and vcf records that overlap lists of regions.\n"+
				"Strict adherence to bed format is assumed (1st base is 0, last base is not included,\n"+
				"last base > first; vcf is in 1 base so subtract 1 from pos to convert to bed).\n"+
				"This app needs > 27G of RAM for human/ mouse/ plant. A two step search is performed\n"+
				"using an in memory file : bp index to find intersecting regions followed by a tabix\n"+
				"query to pull the data.\n"+

				"\nRequired Params:\n"+
				"-c A bed file of chromosomes and their lengths (e.g. chr21 0 48129895) to use to \n"+
				"     building the intersection index. Note, the chr name must match across all \n"+
				"     datasets, no mixing of chr21 and 21.\n"+
				"-d A data directory containing gzipped tabixed (https://github.com/samtools/htslib)\n"+
				"     vcf and bed files. Recurses through all sub directories. Be sure to normalize\n"+
				"     the vcf records with a package like Vt, see http://genome.sph.umich.edu/wiki/Vt\n"+

				"\nDefault Params:\n"+
				"-i Ignore warnings, defaults to printing to stdout\n"+
				"-s Don't print query statistics\n"+
				"-p Print query results to stdout\n"+
				"-n Number of processors to use, defaults to all-1.\n"+

				"\nExample: java -Xmx64G -jar pathToUSeq/Apps/TQuery -c b37ChrLen.bed \n"+
				"     -t TabixDataFiles/ \n\n" +

				"**************************************************************************************\n");
	}
}






