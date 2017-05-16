package edu.utah.seq.barcodes;

import java.io.*;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.*;

/** This app takes the passing bam file from MatchMates after sorting and clusters barcodes associated with alignments with the same unclipped start.
 * Consensus is then called on each stack of clustered reads and written to fastq for subsequent realignment.
 * @author Nix
 * */
public class Consensus {

	//user defined fields
	private File bamFile;
	private File saveDirectory;
	private int numberThreads = 0;
	
	//for chunking
	private int numRecordsInChunk = 1000000;
	
	//for the BarcodeClusterEngine for clustering reads based on their barcode
	private int minBarcodeBaseQuality = 13;
	private double minNumPassingBarcodeBases = 7;
	private double minFractionBarcodeIdentity = 0.875; 
	private int maxAlignmentsToCluster = 20000;
	
	//for the ConsensusEngine for collapsing clustered reads
	private int minReadBaseQuality = 13;
	private double minReadBaseFractionIdentity = 0.66;
	
	//internal 
	protected static final String firstOfPair = "/1";
	protected static final String secondOfPair = "/2";
	private long startTime;
	private BarcodeLoader barcodeLoader;

	/**For stand alone app.*/
	public Consensus(String[] args){
		try {
			startTime = System.currentTimeMillis();

			processArgs(args);
			
			printParams();
			
			//make a BarcodeLoader, wait for initialization completion, this is in it's own thread
			barcodeLoader = new BarcodeLoader(this);
			while (barcodeLoader.isCompletedInitialization() == false) Thread.sleep(1000);
			
			//make mergers, these self feed from the BarcodeLoader
			System.out.println("Clustering and calling consensus on "+barcodeLoader.getChunks().size()+" chunks (Thread Region RamUsed Time):");
			ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
			BarcodeChromMerger[] bcm = new BarcodeChromMerger[numberThreads];
			for (int i=0; i< numberThreads; i++) {
				bcm[i] = new BarcodeChromMerger(this, i+1);
				executor.execute(bcm[i]);
			}
			executor.shutdown();
			//spins here until the executer is terminated, e.g. all threads complete
	        while (!executor.isTerminated()) {
	        }
	        
			//close the loader
			barcodeLoader.completeTasks();
			
			//look for errors
			for (BarcodeChromMerger b :bcm) {
				if (b.isWorkFailed()) throw new Exception ("One of the threads encountered an issue while executing. ");
			}
			
			System.out.println("\nConcatinating temp gz files...");
			concatinateFiles(barcodeLoader, bcm);
			
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");
			
		} catch (Exception e){
			System.err.println("Problem processing "+bamFile+"\n");
			e.printStackTrace();
			System.exit(1);
		}
	}

	/**Concatinates the gzipped temp files.
	 * @throws IOException */
	private void concatinateFiles(BarcodeLoader bc, BarcodeChromMerger[] bcm) throws IOException {
		ArrayList<File> pass  = new ArrayList<File>();
		ArrayList<File> fail  = new ArrayList<File>();
		ArrayList<File> first  = new ArrayList<File>();
		ArrayList<File> second  = new ArrayList<File>();
		ArrayList<File> single  = new ArrayList<File>();
		
		pass.add(bc.getPassingSamFile());
		fail.add(bc.getFailingSamFile());
		
		for (BarcodeChromMerger b: bcm){
			pass.add(b.getPassingSamFile());
			fail.add(b.getFailingSamFile());
			first.add(b.getFastqFileReadOne());
			second.add(b.getFastqFileReadTwo());
			single.add(b.getFastqFileUnpaired());
		}
		
		IO.concatinateFiles(pass, new File (saveDirectory, "passing.sam.gz"));
		IO.concatinateFiles(fail, new File (saveDirectory, "failing.sam.gz"));
		IO.concatinateFiles(first, new File (saveDirectory, "paired_1.fastq.gz"));
		IO.concatinateFiles(second, new File (saveDirectory, "paired_2.fastq.gz"));
		IO.concatinateFiles(single, new File (saveDirectory, "unpaired.fastq.gz"));
		
	}

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
					case 'c': numRecordsInChunk = Integer.parseInt(args[++i]); break;
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
		if (bamFile == null || bamFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find your xxx.bam file!\n");
		
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
		sb.append("\n\tNum align per chunk:\t"); sb.append(numRecordsInChunk);

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
				"**                                  Consensus : March 2017                          **\n" +
				"**************************************************************************************\n" +
				"Consensus clusters alignments sharing the same unclipped start position and molecular\n"+
				"barcode. It then calls consensus on the clustered alignments outputing fastq for\n"+
				"realignment and unmodified bam records. After running, align the fastq files and merge\n"+
				"the new bams with those in the save directory. \n\n "+

				"Required arguments:\n"+
				"-b Path to the mate matched bam file created by FastqBarcodeTagger | cutadapt | bwa |\n"+
				"     MatchMates.  See FBT and MM for details. \n"+
				
				"\nOptional Arguments:\n"+
				"-s Path to a directory to save the results, defaults to a derivative of the\n"+
				"     bam file.\n"+
				"-t Number concurrent threads to run, defaults to the max available to the jvm / 2.\n"+
				"-c Number of alignments to process in one chunk, defaults to 1,000,000. Adjust for the\n"+
				"     availible RAM.\n"+
				"-x Maximum number of alignments to cluster before subsampling, defaults to 20000.\n"+
				"-q Minimum barcode base quality, defaults to 13, anything less is assigned an N.\n"+
				"-n Minimum number of non N barcode bases, defaults to 7, anything less is tossed.\n"+
				"-f Minimum fraction barcode identity for inclusion in a cluster, defaults to 0.875 .\n"+
				"-u Minimum read base quality for inclusion in consensus calling, defaults to 13.\n"+
				"-r Minimum read base fraction identity to call a consensus base, defaults to 0.66 .\n"+
				"     Anything less is assigned an N.\n"+

				"\nExample: java -Xmx100G -jar pathTo/USeq/Apps/Consensus -b MM/passingMM.sorted.bam \n\n"+

				"**************************************************************************************\n");

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
	public int getNumberThreads() {
		return numberThreads;
	}
	public int getNumRecordsInChunk() {
		return numRecordsInChunk;
	}
	public BarcodeLoader getBarcodeLoader() {
		return barcodeLoader;
	}
}
