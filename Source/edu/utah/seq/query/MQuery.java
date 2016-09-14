package edu.utah.seq.query;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.UnknownHostException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import com.mongodb.BasicDBObject;
import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.DBObject;
import com.mongodb.MongoClient;
import edu.utah.seq.useq.data.RegionScoreText;
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
 * So 15.8 sec to fetch all intersecting records with 14K T-N exome variants using Mongo readers.  Hmm the tabix retrieval is slow, need to parallelize! Doubt could get this to < 2 Sec.
 * Better to use a NoSQL db, key: data.
 * */
public class MQuery {

	//user defined fields
	private File chrLenBedFile;
	private File[] vcfFiles;
	//private File[] bedFiles;
	private Gzipper results;
	private int numberThreads = 0;
	private String genomeBuild = "b37";
	private boolean loadDb = true;
	private boolean printWarnings = true;
	private boolean printStats = true;
	private boolean printFindings = false;

	//internal
	private long recordsLoaded = 0;
	private long recordsSkipped = 0;
	/* To restrict the fetching to particular sources then instatiate this Hash, otherwise every hit will be returned */
	private HashSet<String> sourcesToLoad = null;

	//mongo
	MongoClient mongoClient = null;
	DBCollection vcfCollection = null;
	//DBCollection bedCollection = null;

	/* Chrom: BpIndex[] of ArrayList<String>, interbase coordinates, 0 is first base, last base in any region is not included. */
	private HashMap<String, ArrayList<String>[]> chromFileIndex = new HashMap<String, ArrayList<String>[]>();

	/*Contains a data source and the MongoQueries that intersect it.*/
	private HashMap<String, ArrayList<MongoQuery>> sourceMongoQueries = new HashMap<String, ArrayList<MongoQuery>>();
	private ArrayList<MongoQuery> queriesWithHits = new ArrayList<MongoQuery>();

	/*Contains all the MongoQueries split by chromosome. Some may contain results*/
	private HashMap<String, MongoQuery[]> chrMongoQueries = new HashMap<String, MongoQuery[]>();


	public MQuery (String[] args) {
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			buildEngine();

			String diffTime = Num.formatNumberOneFraction(((double)(System.currentTimeMillis() -startTime))/(1000*60));
			System.out.println("\n"+ diffTime+"Min to build using "+IO.memory()+ " of RAM");
			System.out.println("\t"+recordsLoaded+"\tRecords loaded");
			System.out.println("\t"+recordsSkipped+"\tRecords skipped");

			queryBedFilesFromCmdLine();

			//shut down db connection
			mongoClient.close();

		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("\nProblem with executing the GQuery!");
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
			
			//make resultsWriter
			File output = new File(bed.getParentFile(), Misc.removeExtension(bed.getName())+"_gq.txt.gz");
			results = new Gzipper(output);

			intersectBedFileRegionsWithIndex(bed);

			loadMongoQueriesWithData();

			if (printFindings){
				printMongoQueriesWithHitsBySource();
				printAllMongoQueries();
			}
			
			//close writer
			results.close();
		}
	}
	
	private synchronized void appendResults(String res) {
		try {
			results.println(res);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void printAllMongoQueries() {
		System.out.println("\nPrinting All Queries...");
		for (String chr: chrMongoQueries.keySet()){
			MongoQuery[] tqs = chrMongoQueries.get(chr);
			for (MongoQuery tq : tqs){
				System.out.println(tq.getInterbaseCoordinates());
				HashMap<String, ArrayList<DBObject>> sourceRes = tq.getSourceResults();
				if (sourceRes !=null){
					for (String s: sourceRes.keySet()){
						System.out.println("\t"+s);
						for (DBObject r: sourceRes.get(s)){
							System.out.println("\t\t"+r);
						}
					}
				}
			}
		}
	}

	private void printMongoQueriesWithHitsBySource() {
		System.out.println("\nPrinting Queries With Hits by Data Source...");
		for (String source: sourceMongoQueries.keySet()){
			System.out.println(source);
			ArrayList<MongoQuery> al = sourceMongoQueries.get(source);
			for (MongoQuery mq: al){
				System.out.println("\t"+mq.getInterbaseCoordinates());
				ArrayList<DBObject> res = mq.getSourceResults().get(source);
				for (DBObject s: res) System.out.println("\t\t"+s.get("record"));
			}
			System.out.println();
		}
	}

	/**This takes the MongoQueries that intersect a resource, chunks into batches, then fires a thread on each batch to execute the mongo find(). */
	private void loadMongoQueriesWithData() throws IOException {
		long startTime = System.currentTimeMillis();

		//anything to query?
		if (queriesWithHits.size() == 0) return; 

		//chunk
		ArrayList<MongoLoader> loaders = makeLoaders(vcfCollection);
		int numChunks = loaders.size();		

		//just one? execute it
		if (numChunks == 1) loaders.get(0).run();

		//nope, create executor and run it.
		else {

			ExecutorService executor = Executors.newFixedThreadPool(numChunks);
			for (MongoLoader ml: loaders) executor.execute(ml);
			executor.shutdown();

			//spins here until the executer is terminated, e.g. all threads complete
			while (!executor.isTerminated()) {}
		}

		if (printStats) System.out.println( (System.currentTimeMillis() -startTime)+"\tMillisec to load queries with data");
	}

	public ArrayList<MongoLoader> makeLoaders(DBCollection col){
		
		int numOfQuer = queriesWithHits.size();

		//figure out how many queries to put in each loader
		int numToLoad = (int)Math.round((double)numOfQuer/(double)numberThreads);

		//just a small number? then use one thread.
		if (numToLoad < 20) numToLoad = numOfQuer;

		ArrayList<MongoLoader> loaders = new ArrayList<MongoLoader>();
		ArrayList<MongoQuery> chunk = new ArrayList<MongoQuery>(numToLoad);

		for (int i=0; i< numOfQuer; i++){
			//add it
			chunk.add(queriesWithHits.get(i));
			//made count?
			if (chunk.size() == numToLoad){
				loaders.add(new MongoLoader(col, chunk));
				chunk = new ArrayList<MongoQuery>(numToLoad);
			}
		}
		//add last?
		if (chunk.size() > 0) loaders.add(new MongoLoader(vcfCollection, chunk));
		return loaders;
	}

	public void buildEngine() throws IOException {
		connectToMongoDb();

		createChromIndex();

		loadChromIndexWithVcf();
		//loadChromIndexWithBed();
	}

	/**This takes a bed file, parses each region into a MongoQuery object and identifies which file data sources contain data that intersects it.*/
	public void intersectBedFileRegionsWithIndex(File bed) throws IOException {
		//load bed file of regions they want to intersect
		HashMap<String, RegionScoreText[]> chrRegions = Bed.parseBedFile(bed, true, false);

		//convert them to MongoQuery and set in this
		chrMongoQueries = convert2MongoQuery(chrRegions);

		//perform the file search
		queryFileIndex();	
	}

	private void queryFileIndex() {
		long startTime = System.currentTimeMillis();
		int numQueries = 0;
		int numSkippedQueries = 0;
		long numberHits = 0;
		sourceMongoQueries.clear();
		queriesWithHits.clear();

		//for each chromosome of regions to query
		for (String chr: chrMongoQueries.keySet()){

			//check that chr exists in index
			if (chromFileIndex.containsKey(chr) == false){
				int numSkipped = chrMongoQueries.get(chr).length;
				if (printWarnings) System.err.println("\nWARNING: chromosome '"+chr+"' not found in index? Skipping "+numSkipped+" query regions.");
				numSkippedQueries += numSkipped; 
			}

			else {
				ArrayList<String>[] index = chromFileIndex.get(chr);

				//for each query region
				MongoQuery[] regions = chrMongoQueries.get(chr);
				numQueries += regions.length;

				for (MongoQuery q: regions){
					//fetch sources that intersect this region
					HashSet<String> sourceHits = intersect(index, q);

					if (sourceHits.size() !=0) {
						numberHits += sourceHits.size();
						addHits(sourceHits, q);
						queriesWithHits.add(q);
					}
				}
			}
		}
		if (printStats){
			long diffTime = System.currentTimeMillis() -startTime;
			System.out.println("\nQuery Stats:");
			System.out.println(numQueries+ "\tNum index queries");
			System.out.println(numSkippedQueries+ "\tNum skipped index queries");
			System.out.println(queriesWithHits.size()+ "\tNum queries with hits");
			System.out.println(numberHits+ "\tNum total MongoDb searches.");
			System.out.println(diffTime+"\tMillisec to complete index search");
		}
	}

	private void filterSources(HashSet<String> sourceHits) {
		if (sourcesToLoad != null){
			ArrayList<String> toRemove = new ArrayList<String>();
			for (String s: sourceHits){
				if (sourcesToLoad.contains(s) == false) toRemove.add(s);
			}
			sourceHits.removeAll(toRemove);
		}
	}

	/**Adds the MongoQuery to an ArrayList associated with a resource to fetch the data from.*/
	private void addHits(HashSet<String> hits, MongoQuery mq) {

		//create results container for MQ
		HashMap<String, ArrayList<DBObject>> sourceResults = new HashMap<String, ArrayList<DBObject>>();

		//for each source
		for (String source: hits){
			//add to global hash
			ArrayList<MongoQuery> al = sourceMongoQueries.get(source);
			if (al == null){
				al = new ArrayList<MongoQuery>();
				sourceMongoQueries.put(source, al);
			}
			al.add(mq);

			//add to MongoQuery so the loader knows what to search for
			sourceResults.put(source, null);
		}

		//save it in MQ
		mq.setSourceResults(sourceResults);
	}

	/**Checks the begin and end for out of bounds then uses a hash to collapse all the Sources found to have region that overlaps the query.
	 * Also filters against the sourcesToLoad if set.*/
	private HashSet<String> intersect(ArrayList<String>[] index, MongoQuery tq) {
		HashSet<String> hits = new HashSet<String>();
		int begin = tq.getStart();
		if (begin < 0) begin = 0;
		int end = tq.getStop();
		if (end >= index.length) end = index.length;
		for (int i=begin; i<end; i++) {
			if (index[i] != null) hits.addAll(index[i]);
		}

		//filter against those to search?
		if (hits.size() !=0) filterSources(hits);

		return hits;
	}

	public static HashMap<String, MongoQuery[]> convert2MongoQuery(HashMap<String, RegionScoreText[]> chrRegions) throws IOException {
		HashMap<String, MongoQuery[]> chrTQ = new HashMap<String, MongoQuery[]>();
		for (String chr: chrRegions.keySet()){
			RegionScoreText[] regions = chrRegions.get(chr);
			MongoQuery[] tq = new MongoQuery[regions.length];
			for (int i=0; i< regions.length; i++) tq[i] = new MongoQuery(chr, regions[i]);
			chrTQ.put(chr, tq);
		}
		return chrTQ;
	}

	/*private void loadChromIndexWithBed() throws IOException {

		for (int i=0; i< bedFiles.length; i++){

			//for mongo
			String source = vcfFiles[i].toString();
			ArrayList<DBObject> dbo = new ArrayList<DBObject>();

			BufferedReader in = IO.fetchBufferedReader(bedFiles[i]);
			System.out.println("Loading bed "+bedFiles[i].getName());
			String[] t;
			String line;
			String currChrom = "";
			ArrayList<String>[] currIndex = null;
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
						if (currIndex[j] == null) currIndex[j] = new ArrayList<String>(1);
						currIndex[j].add(bedFiles[i]);
					}
					numLoadedRecords++;

					// TODO ASSUMING one bed! Rewrite to collect all source chr pos
					MongoBed mb = new MongoBed(source, new String[]{line});
					dbo.add(mb.createDBObject());
				}
			}
			//add records to db in batch
			bedCollection.insert(dbo);

			System.out.println("\t"+numLoadedRecords+" regions loaded");
			in.close();
		}
	}*/

	/**This walks through each vcf file, loads each record, and for each base covered by the variant, it makes a reference to it in the chr:String[] hash.
	 * It also if so indicated, load the data into mongo.  At present, the "source" reference is just File.toString(). */
	private void loadChromIndexWithVcf() throws IOException {

		//for each vcf file
		for (int i=0; i< vcfFiles.length; i++){

			//for mongo
			String source = vcfFiles[i].toString();
			ArrayList<DBObject> dbo = new ArrayList<DBObject>();

			BufferedReader in = IO.fetchBufferedReader(vcfFiles[i]);
			System.out.print("Parsing vcf "+vcfFiles[i].getName());

			String[] t;
			String line;
			String currChrom = "";
			ArrayList<String>[] currIndex = null;
			long numLoadedRecords = 0;
			long numRecordsSkipped = 0;
			int counter = 0;

			//for each vcf record
			while ((line = in.readLine()) != null){

				//not a header line?
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
							//currently just using File.toString() to pull the unique name of the source
							for (int j= startStop[0]; j< startStop[1]; j++){
								if (currIndex[j] == null) currIndex[j] = new ArrayList<String>(1);
								currIndex[j].add(vcfFiles[i].toString());
							}
							//make db obj
							if (loadDb) {
								//ready to insert chunk?
								if (loadDb && ++counter == 50000) {
									vcfCollection.insert(dbo);
									counter = 0;
									dbo.clear();
									System.out.print(".");
								}
								MongoVcf mv = new MongoVcf(source, startStop, line);
								dbo.add(mv.createDBObject());
							}
							numLoadedRecords++;
						}
					}
				}
			}
			//add records to db in batch
			if (loadDb && dbo.size() !=0) vcfCollection.insert(dbo);

			System.out.println("\n\t"+numLoadedRecords+" variants loaded, "+numRecordsSkipped+" skipped\n");
			recordsLoaded += numLoadedRecords;
			recordsSkipped += numRecordsSkipped;
			in.close();
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


	/**This creates a HashMap of chr: ArrayList<String>[] where each position in the [] represents a base and contains all the sources that overlap. */
	private void createChromIndex() throws IOException {
		HashMap<String, RegionScoreText[]> chrLen = Bed.parseBedFile(chrLenBedFile, true, false);

		System.out.println("Building empty chrom index...");
		for (String chr: chrLen.keySet()){
			RegionScoreText[] regions = chrLen.get(chr);
			if (regions.length !=1) throw new IOException("\nError: there can be only one bed region for each chromosome, see "+chr);
			ArrayList<String>[] indexes = new ArrayList[regions[0].getStop()+1];
			chromFileIndex.put(chr, indexes);
			System.out.println("\t"+chr+"\t"+regions[0].getStop());
		}
	}

	private void connectToMongoDb() throws UnknownHostException{
		// Connect to mongodb server
		mongoClient = new MongoClient( "localhost" , 27017 );

		// Delete old db if it exists
		if (loadDb) mongoClient.dropDatabase(genomeBuild);

		// Connect to/ create the db
		DB db = mongoClient.getDB( genomeBuild );

		//create collection
		if (loadDb){
			vcfCollection = db.getCollection("vcf");
			vcfCollection.createIndex(new BasicDBObject("source", 1));
			vcfCollection.createIndex(new BasicDBObject("chr", 1));
			vcfCollection.createIndex(new BasicDBObject("start", 1));
			vcfCollection.createIndex(new BasicDBObject("stop", 1));
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MQuery (args);
	}	

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		File dataDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': chrLenBedFile = new File(args[++i]); break;
					case 'd': dataDir = new File(args[++i]); break;
					case 'n': numberThreads = Integer.parseInt(args[++i]); break;
					case 'g': genomeBuild = args[++i]; break;
					case 'i': printWarnings = false; break;
					case 'e': loadDb = false; break;
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
		if (dataDir == null || dataDir.isDirectory() == false) Misc.printErrAndExit("\nError: please provide a directory containing tabix indexed xxx.vcf.gz and xxx.bed.gz files with their associated xxx.gz.tbi indexes" );

		vcfFiles = IO.fetchFilesRecursively(dataDir, "vcf.gz");
		//bedFiles = IO.extractFiles(tabixDataDir, "bed.gz");
		//if (vcfFiles.length == 0 && bedFiles.length == 0) Misc.printErrAndExit("\nError: failed to find any xxx.bed.gz or xxx.vcf.gz tabix files in your tabixDataDir -> "+tabixDataDir);
		if (vcfFiles.length == 0) Misc.printErrAndExit("\nError: failed to find any xxx.bed.gz or xxx.vcf.gz files in your DataDir -> "+dataDir);

		//threads to use
		int numAvail = Runtime.getRuntime().availableProcessors();
		if (numberThreads < 1) numberThreads =  numAvail - 1;
		System.out.println(numAvail +" Available threads, using "+numberThreads+"\n");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                    MQuery: Sept 2016                             **\n" +
				"**************************************************************************************\n" +
				"MQ returns vcf records that overlap a provided list of regions in bed format.\n"+
				"Strict adherence to bed format is assumed (1st base is 0, last base is not included,\n"+
				"last base > first; vcf is in 1 base so subtract 1 from pos to convert to bed).\n"+
				"This app needs > 27G of RAM for human/mouse. Requires a running instance of MongoDB. \n"+

				"\nRequired Params:\n"+
				"-c A bed file of chromosomes and their lengths to use to building the intersection\n"+
				"     index. Note, the chr name must match across all datasets, no mixing of chrX and X\n"+
				"-d A directory containing gzipped vcf data files to index. Recurses through\n"+
				"     all sub directories looking for xxx.vcf.gz. Be sure to normalize\n"+
				"     the vcf records with a package like Vt but don't decompose, see \n"+
				"     http://genome.sph.umich.edu/wiki/Vt .\n"+

				"\nDefault Params:\n"+
				"-g Genome build - data base name, defaults to 'b37'\n"+
				"-e Use existing db, defaults to dropping and reloading\n"+
				"-i Ignore warnings, defaults to printing to stdout\n"+
				"-s Don't print query statistics\n"+
				"-p Print query results to stdout\n"+
				"-n Number of processors to use, defaults to all-1. By default, mongo uses 10, so if \n"+
				"     running mongo on the same server set this to all-11 .\n"+

				"\nExample: java -Xmx64G -jar pathToUSeq/Apps/MQuery -c b37ChrLen.bed \n"+
				"     -d ~/VCFFiles/ \n\n" +

				"**************************************************************************************\n");
	}

	public String getGenomeBuild() {
		return genomeBuild;
	}
}






