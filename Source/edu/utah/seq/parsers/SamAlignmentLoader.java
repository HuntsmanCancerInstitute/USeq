package edu.utah.seq.parsers;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import util.bio.annotation.Bed;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class SamAlignmentLoader {
	
	//fields
	private File bam = null;
	private File bed = null;
	private File resultsDir = null;
	private QueryInterval[] regions = null;
	private int regionIndex = 0;
	private int chunkSize = 200;
	private String samHeader = null;
	private File headerFile = null;
	private File[] parsedBams = null;
	
	
	public SamAlignmentLoader(File bam, File bed, File workingDir, int numberThreads){
		try {
			this.bam = bam;
			this.bed = bed;

			//make the dir and delete
			resultsDir = new File(workingDir, Misc.getRandomString(10)+"_"+Misc.removeExtension(bam.getName()));
			resultsDir.mkdir();
			parseRegions();
			System.out.println("\tParsed "+regions.length+ " regions...");
			
			//create loaders
			System.out.println("\tLaunching "+numberThreads +" loaders...");
			BamLoader[] loaders = new BamLoader[numberThreads];
			for (int i=0; i< numberThreads; i++) loaders[i] = new BamLoader(this, i+1);
			chunkSize = regions.length/ numberThreads;
			if (chunkSize < 10) chunkSize = 10;

			//launch em!
			ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
			for (BamLoader l: loaders) executor.execute(l);
			executor.shutdown();
			while (!executor.isTerminated()) {}  //wait here until complete

			//check loaders 
			for (BamLoader l: loaders) {
				if (l.isFailed()) throw new IOException("ERROR: Failed to extract alignments from "+bam.getName());
			}
			
			//write out header
			headerFile = new File(resultsDir, "0_header.sam.gz");
			Gzipper header = new Gzipper( headerFile );
			header.println(samHeader);
			header.close();
			
			parsedBams = combineAndSortResults();
			
			System.out.println("\tComplete");
			
		} catch (Exception e){
			IO.deleteDirectory(resultsDir);
			IO.deleteDirectoryViaCmdLine(resultsDir);
			e.printStackTrace();
		}
	}
	
	private File[] combineAndSortResults() throws Exception {
		//cat results R1
		System.out.println("\tCombinging temp files...");
		ArrayList<File> al = new ArrayList<File>();
		al.add(headerFile);
		for (File f: IO.extractFiles(resultsDir, "_R1.sam.gz")) al.add(f);
		File toSortR1= new File(resultsDir,"toSortR1.sam.gz");
		IO.concatinateFiles(al, toSortR1);
		al.clear();
		
		//cat results R2
		al.add(headerFile);
		for (File f: IO.extractFiles(resultsDir, "_R2.sam.gz")) al.add(f);
		File toSortR2= new File(resultsDir,"toSortR2.sam.gz");
		IO.concatinateFiles(al, toSortR2);
		
		//sort em
		System.out.println("\tSorting results...");
		File r1Bam = new File (resultsDir, "sortedR1.bam");
		File r2Bam = new File (resultsDir, "sortedR2.bam");
		File[] bams = new File[]{r1Bam, r2Bam};
		sortSams(new File[]{toSortR1, toSortR2}, bams, 2);
		
		//delete sams
		IO.deleteFiles(resultsDir, "sam.gz");
		
		return bams;
		
	}
	
	/**Launches a threaded version of Picard's sort sam.*/
	public static void sortSams(File[] toSort, File[] bamOutput, int maxNumThreads) throws IOException{
		SamSorter[] ss = new SamSorter[toSort.length];
		for (int i=0; i< toSort.length; i++) ss[i] = new SamSorter(toSort[i], bamOutput[i]);
		int numThreads = ss.length;
		if (numThreads > maxNumThreads) numThreads = maxNumThreads;
		
		//launch em!
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		for (SamSorter l: ss) executor.execute(l);
		executor.shutdown();
		while (!executor.isTerminated()) {}  //wait here until complete
		
		//check loaders 
		for (SamSorter l: ss) if (l.isFailed()) throw new IOException("ERROR: Failed to sort alignments!");
		
	}

	public synchronized boolean loadRegions(ArrayList<QueryInterval> al){
		int count = 0;
		for (; regionIndex< regions.length; regionIndex++){
			al.add(regions[regionIndex]);
			if (++count >= chunkSize) break;
		}
		if (al.size() != 0) return true;
		return false;
	}
	
	private void parseRegions() throws Exception {
		//parse regions
		Bed[] bedRegions = Bed.parseFile(bed, 0, 0);
		regions = new QueryInterval[bedRegions.length];
				
		//pull indexes
		SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		SamReader samReader = factory.open(bam);
		SAMSequenceDictionary sd = samReader.getFileHeader().getSequenceDictionary();
		
		
		//create QueryIntervals
		for (int i=0; i< bedRegions.length; i++) {
			int index = sd.getSequenceIndex(bedRegions[i].getChromosome());
			if (index == -1) throw new Exception("Failed to find a chromosome index for this bed line "+bedRegions[i].toString() +" in "+bam.getName());
			regions[i] = new QueryInterval(index, bedRegions[i].getStart(), bedRegions[i].getStop());
		}	
		Arrays.sort(regions);
		//samReader.getFileHeader().setProgramRecords(new ArrayList<SAMProgramRecord>());
		samHeader = samReader.getFileHeader().getTextHeader().trim();
		
		samReader.close();
	}
	
	/**For testing.*/
	public static void main (String[] args){
		if (args.length == 0) System.out.println("\nUSAGE: bamFile, bedFile, workingDir, numberThreads, bpPadding\n");
		else new SamAlignmentLoader(new File(args[0]), new File(args[1]), new File(args[2]), 5);
	}


	public File getBam(){
		return bam;
	}
	public File getResultsDir() {
		return resultsDir;
	}

	public File[] getParsedBams() {
		return parsedBams;
	}

}


