package edu.utah.seq.barcodes;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.parsers.PairedAlignmentChrParser;
import htsjdk.samtools.*;
import util.bio.annotation.Bed;
import util.bio.annotation.Coordinate;
import util.gen.*;

/**
 * @author Nix
 * */
public class Consensus {

	//user defined fields
	private File bamFile;
	private File saveDirectory;
	private int numberThreads = 0;
	
	//for chunking
	private int numRecordsInChunk = 200000;
	private int maxRecordsInChunk = 300000;
	private int minBasePairGap = 500;
	private int halfGap = minBasePairGap/2;
	
	//for the BarcodeClusterEngine for clustering reads based on their barcode
	private int minBarcodeBaseQuality = 20;
	private double minNumPassingBarcodeBases = 6;
	private double minFractionBarcodeIdentity = 0.75; 
	private int maxAlignmentsToCluster = 5000;
	
	//for the ConsensusEngine for collapsing clustered reads
	private int minReadBaseQuality = 13;
	private double minReadBaseFractionIdentity = 0.66;
	
	//internal 
	protected SamReaderFactory readerfactory;
	protected SAMFileWriterFactory writerFactory;
	private ArrayList<Coordinate> chunks = new ArrayList<Coordinate>();
	private int numHardBreaks;
	private HashMap<String, HashMap<String, SAMRecord>[]> chromDiscordAlign = new HashMap<String, HashMap<String, SAMRecord>[]>();
	protected static final String firstOfPair = "/1";
	protected static final String secondOfPair = "/2";
	
	//stats collected from workers
	private long numRecords = 0;
	private long numFastqPairsOut = 0;
	private long numFastqUnpairedOut = 0;
	private long numPassingBamOut = 0;
	private long numFailingBamOut = 0;
	private long numBasesSubsampled = 0;
	private long numPropPairMatesNotFound;
	private Histogram famSizeHist = new Histogram(0, 100, 100);

	/**For stand alone app.*/
	public Consensus(String[] args){
		try {
			
			//start clock
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			printParams();

			makeFactories();
			
			System.out.print("Loading discordant alignments and chunking genome...\n\t");
			fetchChromsAndLoadDiscordantRecords();
			
			//just for testing
			Coordinate[] c = new Coordinate[chunks.size()];
			chunks.toArray(c);
			Coordinate.writeToFile(c, new File (saveDirectory,"chunks.bed"));
			
			writeOutUnmapped();

			System.out.print("\nClustering and calling consensus on "+chunks.size()+" chunks\n\t");
			parseAlignmentFile();
			
			printStats();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");
			
		} catch (Exception e){
			System.err.println("Problem parsing "+bamFile+"\n");
			e.printStackTrace();
		}
	}
	
	private void printStats() {
		//stats
		System.out.println("\n\nProcessing statistics:");
		System.out.println("Total records:       "+numRecords);
		System.out.println("Paired fastq:        "+numFastqPairsOut);
		System.out.println("Unpaired fastq:      "+numFastqUnpairedOut);
		System.out.println("Passing bam records: "+numPassingBamOut);
		System.out.println("Failing bam records: "+numFailingBamOut);
		System.out.println("Num hard chunk breaks: "+numHardBreaks);
		System.out.println("Num prop pair mates not found: "+numPropPairMatesNotFound);
		System.out.println("Num bases subsampled w/ too many alignments: "+ numBasesSubsampled);
		//histo
		System.out.println("\nFamily size histogram:");
		famSizeHist.printScaledHistogram();
		
	}

	public void fetchChromsAndLoadDiscordantRecords(){
		try {
			//get reader 
			SamReader samReader = readerfactory.open(bamFile);
			//fetch header
			//samHeader = samReader.getFileHeader().getTextHeader().trim();
			//samHeader = samHeader+"\n@PG\tID:MergePairedAlignments\tCL: "+programArguments;
			//get chromosomes
			List<SAMSequenceRecord> chrs = samReader.getFileHeader().getSequenceDictionary().getSequences();
			//chromosomes.add("21");
			samReader.close();
			
			//load discordant alignments
			loadDiscordantAlignments(chrs);
			
		} catch (Exception e) {
			System.err.println("\nError loading chromsomes from "+bamFile.getName());
			e.printStackTrace();
			System.exit(0);
		} 
	}
	/*This is going to be a memory hog.  But it's needed since the picard queryMate() method takes seconds to complete.
	 * I've been forced to use it for discordant pairs so will try to see if they can be squeezed into memory for a faster lookup.
	 * This creates a HashMap with the keys a String rep of each chromosome with data.  If there is no chrom then there is no discordant alignments.
	 * The Value is a HashMap[] where each index represents a potential alignment start.
	 * The HashMap for from a give index contains the readName/readOrder as the key and the value the assoicated SAMRecord.
	 * Again, if no entry then no discordant pair is available.
	 * Don't think any of this needs to be synchronized since no structural changes are made to the hash once it is created. 
	 * It's OK to null the associated name and SAMRecord, but not remove the entry.*/
	private void loadDiscordantAlignments(List<SAMSequenceRecord> chrs) throws IOException {
		SamReader samReader = readerfactory.open(bamFile);
		
		//for each chromosome
		for (SAMSequenceRecord r : chrs) {
			long numDis = 0;
			//create entry for the hash hash
			String chromName = r.getSequenceName();
			int length = r.getSequenceLength();
			HashMap<String, SAMRecord>[] posHash = new HashMap[length+150];
			boolean foundDiscordantRecords = false;
			
			//for bed creation
			int start = -1;
			int lastUnclippedStart = 0;
			int numRec = 0;
			
			//iterate through records for this chrom
			SAMRecordIterator it = samReader.queryOverlapping(chromName, 0, 0);
			
			while (it.hasNext()){
				SAMRecord sam = it.next();
				//is this unmapped or secondary?
				if (sam.getReadUnmappedFlag() || sam.isSecondaryOrSupplementary()) continue;
				
				int unclippedStart = sam.getUnclippedStart();
				//first position?
				if (start == -1) {
					start = unclippedStart-halfGap;
					if (start < 0) start = 0;
				}
				
				//check chunk
				if (numRec > numRecordsInChunk){
					//exceeded max? this is bad, don't want to do it! how else?
					if (numRec > maxRecordsInChunk){
						//make hard break
						numHardBreaks++;
						int end = unclippedStart;
						chunks.add(new Coordinate (chromName, start, end));
						start = unclippedStart;
						numRec = 0;
					}
					else {
						//OK, have collected enough records, now look for big enough gap
						int diff = unclippedStart-lastUnclippedStart;
						if (diff > minBasePairGap){
							//make bed
							int end = lastUnclippedStart + halfGap;
							chunks.add(new Coordinate (new String (chromName), start, end));
							start = unclippedStart-halfGap;
							numRec = 0;
						}
					}
				}
				numRec++;
				lastUnclippedStart = unclippedStart;
				
				//is it discordant? or mate unmapped?
				if (sam.getProperPairFlag() || sam.getMateUnmappedFlag()) continue;
				
				//OK, this sam is discordant, mapped, not secondary or supp, and has a mapped mate, save it!
				//fetch the hashmap at the aligned position (using this one since it's what is recorded in the mate position)
				int alignmentStart = sam.getAlignmentStart();
				if (posHash[alignmentStart] == null){
					posHash[alignmentStart] = new HashMap<String, SAMRecord>();
				}
				//create a name / 1 or 2 to represent firstOfPair or secondOfPair 
				String nameOrder;
				if (sam.getFirstOfPairFlag()) nameOrder = sam.getReadName()+ firstOfPair;
				else nameOrder = sam.getReadName()+ secondOfPair;
				//add it
				posHash[alignmentStart].put(nameOrder, sam);
				foundDiscordantRecords = true;
				numDis++;
			}
			it.close();
			
			//make last chunk?
			if (numRec !=0) chunks.add(new Coordinate (chromName, start, lastUnclippedStart + halfGap));
			
			//save it?
			if (foundDiscordantRecords) {
				chromDiscordAlign.put(chromName, posHash);
				System.out.print(chromName+":"+numDis+" ");
			}
		}
		System.out.println();
		samReader.close();
	}

	private void writeOutUnmapped() throws IOException {
		SamReader samReader = readerfactory.open(bamFile);
		
		String name = Misc.removeExtension( bamFile.getName() );
		File bamFile = new File (saveDirectory, "failingNotMerged_Unmapped_"+name+".bam");
		SAMFileWriter bamWriter = writerFactory.makeBAMWriter(samReader.getFileHeader(), false, bamFile);
		
		SAMRecordIterator it = samReader.queryUnmapped();
		while (it.hasNext()) {
			bamWriter.addAlignment(it.next());
			numFailingBamOut++;
		}
		numRecords = numFailingBamOut;
		it.close();
		bamWriter.close();
		samReader.close();
	}

	private void makeFactories() throws IOException {
		readerfactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		writerFactory = new SAMFileWriterFactory();
		writerFactory.setTempDirectory(bamFile.getParentFile());
		//check bam
		SamReader reader = readerfactory.open(bamFile);
		boolean indexPresent = reader.hasIndex();
		reader.close();
		if ( indexPresent == false) Misc.printErrAndExit("\nError: cannot find the index for your bam file? "+bamFile);
	}





	public void parseAlignmentFile(){
		try {
			//chunks.clear();
			//chunks.add(new Coordinate("21", 39047560, 40036261));
			//chunks.add(new Coordinate ("21", 42514616, 42885579));
			
			ThreadPoolExecutor executor = (ThreadPoolExecutor) Executors.newFixedThreadPool(4);
			//for each chunk lauch with the executor
			for (int i=0; i< chunks.size(); i++) {
				Coordinate region = chunks.get(i);
				executor.execute( new BarcodeChromMerger(new String(region.getChromosome()), region.getStart(), region.getStop(), this));
			}
			
			executor.shutdown();
			while (!executor.isTerminated()) {
	        }
			
		} catch (Exception e){
			System.err.println("\nError in parsing alignment file!");
			e.printStackTrace();
			System.exit(1);
		}

	}

	/**Runs each parser in its own thread to max threads defined by user.*/
	/*
	private void execute(BarcodeChromMerger[] mergers) {
		while (true){ 
			try {
				//how many running and get AL of those that could be started
				int numRunning = 0;
				int numNull = 0;
				ArrayList<BarcodeChromMerger> toStart = new ArrayList<BarcodeChromMerger>();
				for (BarcodeChromMerger p : mergers){
					if (p == null) {
						numNull++;
					}
					else {
						if (p.isWorkStarted() && p.isWorkComplete()==false) {
							System.out.println(p+" "+p.getNumRecords() +" "+p.getNumWorkingRecords());
							numRunning++;
						}
						if (p.isWorkStarted() == false) toStart.add(p);
					}
				}
				//any to start?
				int numToStart = numberThreads - numRunning;
				System.out.println("NumRunning: "+numRunning+" NumRemaining: "+toStart.size() +" "+ IO.memoryUsed());
				if (numToStart > 0){
					int stop = numToStart;
					if (stop > toStart.size()) stop = toStart.size();
					for (int i=0; i< stop; i++){
						BarcodeChromMerger p = toStart.get(i);
						p.start();
						//System.out.println("Starting "+p);
					}
				}
				boolean allComplete = true;
				for (BarcodeChromMerger p : mergers){
					if (p != null) {
						if (p.isWorkComplete() == false) {
							allComplete = false;
							break;
						}
						//null it and call garbage collection
						else {
							//System.out.println("Nulling "+p);
							p = null;
							System.gc();
						}
					}
				}
				if (allComplete) break;
				//sleep
				Thread.sleep(1000);
			} catch (InterruptedException ie) {
				ie.printStackTrace();
				Misc.printErrAndExit("\nProblem running threads!");
			}
		}
	} */

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new Consensus(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bamFile = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 't': numberThreads = Integer.parseInt(args[++i]); break;
					case 'q': minBarcodeBaseQuality = Integer.parseInt(args[++i]); break;
					case 'n': minNumPassingBarcodeBases = Integer.parseInt(args[++i]); break;
					case 'x': maxAlignmentsToCluster = Integer.parseInt(args[++i]); break;
					case 'f': minFractionBarcodeIdentity = Double.parseDouble(args[++i]); break;
					case 'u': minReadBaseQuality = Integer.parseInt(args[++i]); break;
					case 'r': minReadBaseFractionIdentity = Double.parseDouble(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		


		//check bam
		if (bamFile == null || bamFile.canRead() == false) Misc.printExit("\nError: cannot find your xxx.bam file!\n");
		
		//any save directory?
		if (saveDirectory == null){
			String name = bamFile.getName();
			name = "Consensus_"+Misc.removeExtension(name);
			saveDirectory = new File (bamFile.getParentFile(), name);
		}
		saveDirectory.mkdirs();
		if (saveDirectory.canWrite() == false) Misc.printExit("\nError: cannot create or write to the save directory -> "+saveDirectory);
		
		//number of threads to use?
		if (numberThreads == 0) numberThreads = Runtime.getRuntime().availableProcessors()/2;

	}	
	
	public void printParams() throws IOException{
		StringBuilder sb = new StringBuilder();
		sb.append("Settings:");
		sb.append("\n\tBam file:\t"); sb.append(bamFile.getCanonicalPath());
		sb.append("\n\tSave dir:\t"); sb.append(saveDirectory.getCanonicalPath());
		sb.append("\n\tNum threads:\t"); sb.append(numberThreads);

		sb.append("\nBarcode similarity clustering:");
		sb.append("\n\tMin barcode base qual:\t"); sb.append(minBarcodeBaseQuality);
		sb.append("\n\tMin num pass barcode bases:\t"); sb.append((int)minNumPassingBarcodeBases);
		sb.append("\n\tMin fract barcode identity:\t"); sb.append(minFractionBarcodeIdentity);
		sb.append("\n\tMax num align to cluster:\t"); sb.append(maxAlignmentsToCluster);

		sb.append("\nConsensus calling on clustered alignments:");
		sb.append("\n\tMin read base qual for con calling:\t"); sb.append(minReadBaseQuality);
		sb.append("\n\tMin read base fract identity:\t"); sb.append(minReadBaseFractionIdentity);
		sb.append("\n");
		
		System.out.println(sb);
	}

	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                  Consensus : Aug 2015                            **\n" +
				"**************************************************************************************\n" +
				"Consensus clusters alignments sharing the same unclipped start position and molecular\n"+
				"barcode. It then calls consensus on the clustered alignments outputing fastq for\n"+
				"realignment and unmodified bam records. After running, align the fastq files and merge\n"+
				"the new bams with those in the save directory. ~100G of memory is needed for an exome\n"+
				"sized dataset. \n\n" +

				"Required arguments:\n"+
				"-b Path to a coordinate sorted, indexed, paired and molecular barcoded bam file.\n" +
				"     Use the FastqBarcodeTagger to attach barcodes to fastq headers before aligning.\n"+
				
				"\nOptional Arguments:\n"+
				"-s Path to a directory to save the results, defaults to a derivative of the\n"+
				"     bam file.\n"+
				"-t Number concurrent threads to run, defaults to the max available to the jvm / 2.\n"+
				"-x Maximum number of alignments to cluster before subsampling, defaults to 5000.\n"+
				"-q Minimum barcode base quality, defaults to 20, anything less is assigned an N.\n"+
				"-n Minimum number of non N barcode bases, defaults to 6, anything less is tossed.\n"+
				"-f Minimum fraction barcode identity for inclusion in a cluster, defaults to 0.75 .\n"+
				"-u Minimum read base quality for inclusion in consensus calling, defaults to 13.\n"+
				"-r Minimum read base fraction identity to call a consensus base, defaults to 0.66 .\n"+
				"     Anything less is assigned an N.\n"+

				"\nExample: java -Xmx100G -jar pathTo/USeq/Apps/Consensus -b /Data/res.sorted.bam \n\n"+

				"**************************************************************************************\n");

	}

	public SamReaderFactory getReaderfactory() {
		return readerfactory;
	}

	public SAMFileWriterFactory getWriterFactory() {
		return writerFactory;
	}

	public File getSaveDirectory() {
		return saveDirectory;
	}

	public File getBamFile() {
		return bamFile;
	}

	public int getMinBarcodeBaseQuality() {
		return minBarcodeBaseQuality;
	}

	public double getMinNumPassingBarcodeBases() {
		return minNumPassingBarcodeBases;
	}

	public double getMinFractionBarcodeIdentity() {
		return minFractionBarcodeIdentity;
	}

	public int getMinReadBaseQuality() {
		return minReadBaseQuality;
	}

	public double getMinReadBaseFractionIdentity() {
		return minReadBaseFractionIdentity;
	}

	public int getMaxAlignmentsToCluster() {
		return maxAlignmentsToCluster;
	}

	public HashMap<String, HashMap<String, SAMRecord>[]> getChromDiscordAlign() {
		return chromDiscordAlign;
	}
	
	public void incrementNumRecords(long n){
		numRecords +=n;
	}
	public void incrementNumFastqPairsOut(long n){
		numFastqPairsOut +=n;
	}
	public void incrementNumFastqUnpairedOut(long n){
		numFastqUnpairedOut +=n;
	}
	public void incrementNumPassingBamOut(long n){
		numPassingBamOut +=n;
	}
	public void incrementNumFailingBamOut(long n){
		numFailingBamOut +=n;
	}
	public void incrementNumBasesSubsampled(long n){
		numBasesSubsampled +=n;
	}
	public void incrementNumPropPairMatesNotFound(long n){
		numPropPairMatesNotFound +=n;
	}
	public void incrementFamSizeHist(Histogram h) throws Exception{
		famSizeHist.addCounts(h);
	}


}
