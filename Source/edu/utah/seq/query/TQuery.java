package edu.utah.seq.query;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;


/** Tool for intersecting regions with a collection of vcf and bed data. */
public class TQuery {

	//fields
	private File chrLengthFile;
	private File[] vcfDataFiles = null;
	private File[] mafDataFiles = null;
	private File[] bedDataFiles = null;
	private int numberThreads = 0;
	private int numberQueriesInChunk = 1000;
	private boolean printWarnings = true;
	private boolean printStats = true;
	
	//internal
	private DataSources dataSources;
	private QueryIndex queryIndex;
	private QueryLoader queryLoader;

	//constructor
	public TQuery (String[] args) {
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);
			
			queryIndex = new QueryIndex(this);
			queryLoader = new QueryLoader(this);

			//print some stats on building the engine
			String diffTime = Num.formatNumberOneFraction(((double)(System.currentTimeMillis() -startTime))/1000);
			int numFiles = vcfDataFiles.length + bedDataFiles.length + mafDataFiles.length;
			System.err.println("\n"+ diffTime+" Sec to build using "+IO.memory()+ " of RAM");
			System.err.println("\t"+numFiles+"\tData sources loaded");
			System.err.println("\t"+dataSources.getRecordsLoaded()+"\tRecords indexed");
			System.err.println("\t"+dataSources.getRecordsSkipped()+"\tRecords skipped\n");
			
			//print summary of available filters 
			System.err.println(dataSources.fetchSummary());

			queryFilesFromCmdLine();

			//release file handles
			queryLoader.closeTabixReaders();

		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("\nProblem with executing the TQuery!");
		}
	}

	
	private void queryFilesFromCmdLine() throws IOException {
		while (true){
			System.err.println("\nProvide a bed or vcf file to intersect (or blank to exit):");

			//parse input
			String in = (new BufferedReader(new InputStreamReader(System.in))).readLine().trim();
			if (in == null || in.length()==0) {
				System.err.println();
				break;
			}
			
			//create new QueryRequest
			try {
				new QueryRequest(in, this);
			} catch (IOException e) {
				System.err.println("\nERROR executing this query\n");
				e.printStackTrace();
			}
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
					case 'q': numberQueriesInChunk = Integer.parseInt(args[++i]); break;
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

		//pull data sources
		vcfDataFiles = IO.fetchFilesRecursively(tabixDataDir, "vcf.gz");
		bedDataFiles = IO.fetchFilesRecursively(tabixDataDir, "bed.gz");
		mafDataFiles = IO.fetchFilesRecursively(tabixDataDir, "maf.txt.gz");
		if (vcfDataFiles.length == 0 && bedDataFiles.length == 0 && mafDataFiles.length == 0) Misc.printErrAndExit("\nError: failed to find any xxx.bed.gz, xxx.vcf.gz, or xxx.maf.txt.gz tabix files in your tabixDataDir -> "+tabixDataDir);

		//check for index
		lookForTabixIndex(vcfDataFiles);
		lookForTabixIndex(bedDataFiles);
		lookForTabixIndex(mafDataFiles);

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
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                    TQuery: Oct 2016                             **\n" +
				"**************************************************************************************\n" +
				"TQ returns vcf, maf, and bed records that overlap user's regions in bed or vcf format.\n"+
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
				"     vcf, bed, and maf files. Recurses through all sub directories. Be sure to normalize\n"+
				"     the vcf records with a package like Vt, see http://genome.sph.umich.edu/wiki/Vt .\n"+
				"     Use the MafParser to sort and tabix your TCGA xxx.maf.txt files.\n"+

				"\nDefault Params:\n"+
				"-n Number of processors to use, defaults to all-1.\n"+
				"-q Number of queries in each lookup chunk, defaults to 1000\n"+

				"\nExample: java -Xmx64G -jar pathToUSeq/Apps/TQuery -c b37ChrLen.bed \n"+
				"     -t TabixDataFiles/ | gzip > results.json.gz  \n\n" +

				"**************************************************************************************\n");
	}

	public boolean isPrintWarnings() {
		return printWarnings;
	}

	public File getChrLengthFile() {
		return chrLengthFile;
	}

	public void setChrLengthFile(File chrLengthFile) {
		this.chrLengthFile = chrLengthFile;
	}

	public File[] getVcfDataFiles() {
		return vcfDataFiles;
	}

	public void setVcfDataFiles(File[] vcfDataFiles) {
		this.vcfDataFiles = vcfDataFiles;
	}

	public File[] getMafDataFiles() {
		return mafDataFiles;
	}

	public void setMafDataFiles(File[] mafDataFiles) {
		this.mafDataFiles = mafDataFiles;
	}

	public File[] getBedDataFiles() {
		return bedDataFiles;
	}

	public void setBedDataFiles(File[] bedDataFiles) {
		this.bedDataFiles = bedDataFiles;
	}

	public int getNumberThreads() {
		return numberThreads;
	}

	public void setNumberThreads(int numberThreads) {
		this.numberThreads = numberThreads;
	}

	public int getNumberQueriesInChunk() {
		return numberQueriesInChunk;
	}

	public void setNumberQueriesInChunk(int numberQueriesInChunk) {
		this.numberQueriesInChunk = numberQueriesInChunk;
	}

	public boolean isPrintStats() {
		return printStats;
	}

	public void setPrintStats(boolean printStats) {
		this.printStats = printStats;
	}

	public DataSources getDataSources() {
		return dataSources;
	}

	public void setDataSources(DataSources dataSources) {
		this.dataSources = dataSources;
	}

	public void setPrintWarnings(boolean printWarnings) {
		this.printWarnings = printWarnings;
	}

	public QueryIndex getQueryIndex() {
		return queryIndex;
	}

	public QueryLoader getQueryLoader() {
		return queryLoader;
	}


}






