package edu.utah.seq.query;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.useq.data.RegionScoreText;
import htsjdk.tribble.readers.TabixReader;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;


/** Test of using a big genome index and tabix lookup for a variant store.
 * 23669MB to load empty gIndex. 
 * 25526MB to load cosmic, 3 10K vcfs, dbSNP into gIndex.
 * 790ms to int entire exome with gIndex.
 * 116,745ms to fetch overlap vars for entire exome (257,726 regions intersected  2,221,581 records)
 * 15,855 ms for 14111 short vars pulling 13841 Records
 * So 15.8 sec to fetch all intersecting records with 14K T-N exome variants using Tabix readers.  Hmm the tabix retrieval is slow, need to parallelize! Doubt could get this to < 2 Sec.
 * Better to use a NoSQL db, key: data.
 * */
public class IntersectionEngine {

	//fields
	private File chrLenBedFile;
	private File[] vcfFiles;
	private File[] bedFiles;
	private File queryRegionsBedFile;
	private int numberThreads = 0;
	private boolean verbose = true;

	//internal
	private HashMap<String, HashSet<File>[]> chromFileIndex = new HashMap<String, HashSet<File>[]>();


	public IntersectionEngine (String[] args) {
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			buildEngine();

			long diffTime = System.currentTimeMillis() -startTime;
			System.out.println("\n"+ diffTime+" MS to instatiate using "+IO.memory()+ " of RAM");

			while (true){
				System.out.print("\nProvide a bed file of regions to query against the engine (or blank to exit):\n");
				String bedFileString = (new BufferedReader(new InputStreamReader(System.in))).readLine();
				if (bedFileString == null || bedFileString.length()==0) {
					System.out.println();
					break;
				}
				queryRegions(new File (bedFileString));
			}


		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("\nProblem with IntersectionEngine!");
		}
	}

	public void buildEngine() throws IOException {
		createChromIndex();
		loadChromIndexWithVcf();
		loadChromIndexWithBed();
	}

	public void queryRegions(File bed) {
		//load as bed
		HashMap<String, RegionScoreText[]> chrRegions = Bed.parseBedFile(bed, true, false);
		//convert to Tabix
		HashMap<String, TabixQuery[]> chrTQ = convert2TabixQuery(chrRegions);
		//fetch files to retrieve actual data.
		HashMap<File, ArrayList<TabixQuery>> files = queryFileIndex(chrTQ);	
	}

	private HashMap<File, ArrayList<TabixQuery>> queryFileIndex(HashMap<String, TabixQuery[]> chrTQ) {
		long startTime = System.currentTimeMillis();
		long numQueries = 0;
		long numSkippedQueries = 0;
		long numQueriesWithHits = 0;
		long numberHits = 0;

		//want to create a hashmap of File:ArrayList<TabixQuery> to query
		HashMap<File, ArrayList<TabixQuery>> toQuery = new HashMap<File, ArrayList<TabixQuery>>();

		//for each chromosome of regions to query
		for (String chr: chrTQ.keySet()){

			//check that chr exists in index
			if (chromFileIndex.containsKey(chr) == false){
				int numSkipped = chrTQ.get(chr).length;
				System.err.println("\nWARNING: chromosome '"+chr+"' not found in index? Skipping "+numSkipped+" query regions.");
				numSkippedQueries += numSkipped; 
			}

			else {
				HashSet<File>[] index = chromFileIndex.get(chr);
				//for each region
				TabixQuery[] regions = chrTQ.get(chr);
				numQueries += regions.length;
				for (TabixQuery tq: regions){
					HashSet<File> hits = intersect(index, tq);
					if (hits.size() !=0) {
						numberHits += hits.size();
						numQueriesWithHits++;
					}
					addHits(hits, toQuery, tq);
				}
			}
		}
		if (verbose){
			long diffTime = System.currentTimeMillis() -startTime;
			System.out.println("\nFile Query Stats:");
			System.out.println(numQueries+ "\tNum executed queries");
			System.out.println(numSkippedQueries+ "\tNum skipped queries");
			System.out.println(numQueriesWithHits+ "\tNum queries with hits");
			System.out.println(numberHits+ "\tNum total hits");
			System.out.println(diffTime+"\tMillisec to complete");
			System.out.println("\nDataSource\tHits");
			for (File f: toQuery.keySet()) {
				System.out.println(f.getName()+"\t"+TabixQuery.getCoordinates(toQuery.get(f)));
			}

		}
		return toQuery;
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

	private static HashMap<String, TabixQuery[]> convert2TabixQuery(HashMap<String, RegionScoreText[]> chrRegions) {
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
			System.out.println("\t"+numLoadedRecords+" regions loaded");
			in.close();
		}
	}


	private void loadChromIndexWithVcf() throws IOException {

		for (int i=0; i< vcfFiles.length; i++){
			BufferedReader in = IO.fetchBufferedReader(vcfFiles[i]);
			System.out.println("Loading vcf "+vcfFiles[i].getName());
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
					if (t[4].contains("<")){
						if (numRecordsSkipped ==0) System.err.println("\tWARNING: Only SNV/INDEL variants supported, skipping lines like "+line);
						numRecordsSkipped++;
					}

					else {
						//check HashSet at the genomic position, note interbase coordinates so subtract 1 from vcf pos
						int pos = Integer.parseInt(t[1]) - 1;

						//define the stop, for snvs and indels its pos + len of ref plus one since the stop isn't included in interbase
						int stop = pos + t[3].length() +1;

						//add in references to source file over the covered bases, stop isn't covered.
						for (int j= pos; j< stop; j++){
							if (currIndex[j] == null) currIndex[j] = new HashSet<File>(1);
							currIndex[j].add(vcfFiles[i]);
						}
						numLoadedRecords++;
					}

				}
			}
			System.out.println("\t"+numLoadedRecords+" variants loaded, "+numRecordsSkipped+" skipped");
			in.close();
		}
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
		System.out.println("\t"+IO.memory()+"\n");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new IntersectionEngine (args);
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
					case 't': tabixDataDir = new File(args[++i]); break;
					case 'q': queryRegionsBedFile = new File(args[++i]); break;
					case 'n': numberThreads = Integer.parseInt(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (numberThreads < 1) numberThreads = Runtime.getRuntime().availableProcessors() - 1;
		if (chrLenBedFile == null) Misc.printErrAndExit("\nError: please provide a bed file of chromosome and their lengths, e.g. X 0 155270560\n" );
		if (tabixDataDir == null || tabixDataDir.isDirectory() == false) Misc.printErrAndExit("\nError: please provide a directory containing tabix indexed xxx.vcf.gz and xxx.bed.gz files with their associated xxx.gz.tbi indexes" );

		vcfFiles = IO.extractFiles(tabixDataDir, "vcf.gz");
		bedFiles = IO.extractFiles(tabixDataDir, "bed.gz");
		if (vcfFiles.length == 0 && bedFiles.length == 0) Misc.printErrAndExit("\nError: failed to find any xxx.bed.gz or xxx.vcf.gz tabix files in your tabixDataDir -> "+tabixDataDir);

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Intersection Engine: Aug 2016                        **\n" +
				"**************************************************************************************\n" +
				"Beta. Needs > 25G of RAM.\n"+

				"\nRequired Options:\n"+
				"-c A bed file of chromosomes and their lengths to use to build the index, note the chr\n"+
				"     name must match all datasets, no mixing of chrX and X, e.g. chrX 0 155270560\n"+
				"-t A directory containing gzipped and tabix indexed vcf and bed files.\n"+
				"-n Number of processors to use, defaults to all available.\n"+

				"\nExample: java -Xmx64G -jar pathToUSeq/Apps/IntersectionEngine -c b37ChrLen.bed \n"+
				"     -t TabixDataFiles/ -n 40 \n\n" +

				"**************************************************************************************\n");

	}
}






