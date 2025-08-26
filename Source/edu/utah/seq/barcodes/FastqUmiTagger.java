
package edu.utah.seq.barcodes;

import java.io.*;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.*;
import util.gen.*;

/**@author david.nix@hci.utah.edu*/
public class FastqUmiTagger{
	
	//user defined fields
	private File firstReadFastq;
	private File secondReadFastq;
	private File umiFastq;
	private int umiLength = -1;
	private int bpToTrim = -1;
	private int numberThreads = 1;
	private int chunkSize = 1000;
	private int maxNsInUMI = -1;
	
	//internal fields
	private BufferedReader firstFastqIn;
	private BufferedReader secondFastqIn;
	private BufferedReader umiFastqIn;
	private long numberFastqLoaded = 0;
	private long numberFastqWithFailedUmi = 0;
	private boolean allRead = false;

	//constructors
	public FastqUmiTagger(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);

		try {
			makeReaders();

			parse();
			
			closeIO();
			
		} catch (IOException e) {
			System.err.println("\n\nProblems found parsing fastq data with the FastqUmiTagger!");
			Misc.printErrAndExit(e.getMessage());
		} 
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		String res = "FastqUmiTagger2 Done "+Math.round(diffTime)+" min\n\t"+
		   numberFastqLoaded+" fastq records loaded using "+numberThreads+" threads\n\t"+
		   numberFastqWithFailedUmi+" UMIs with too many N's\n\t"+	
		   (numberFastqLoaded - numberFastqWithFailedUmi)+" fastq records printed to stream";
		System.err.println(res);
	}
	

	
	synchronized boolean load3Reads(ArrayList<String[]> reads) throws Exception {
		//all done, tell the worker to shut down
		if (allRead) return false;

		for (int x=0; x< chunkSize; x++) {
			//attempt to load 4 lines
			int i=0;
			String[] first = new String[4];
			String[] second = new String[4];
			String[] umi = new String[4];
			boolean failed = false;

			for (; i<4; i++){
				first[i] = firstFastqIn.readLine();
				second[i] = secondFastqIn.readLine();
				umi[i] = umiFastqIn.readLine();
				//check all are good
				if (first[i] == null || second[i] == null || umi[i] == null) {
					failed = true;
					if (first[i] == null && second[i] == null && umi[i] == null) allRead = true;
					break;
				}
			}
			
			//good load?
			if (failed == false) {
				reads.add(first);
				reads.add(second);
				reads.add(umi);
				numberFastqLoaded++;
			}
			//nope if failed
			else {
				closeIO();
				//all failed so this is the last of the fastq
				if (allRead) return true;
				//nope truncation!
				else Misc.printErrAndExit("\nError: looks like one of your fastq files is truncated! Only loaded part of a four line fastq record.");
			}
		}
		return true;
	}
	
	synchronized boolean load2Reads(ArrayList<String[]> reads) throws Exception {
		//all done, tell the worker to shut down
		if (allRead) return false;

		for (int x=0; x< chunkSize; x++) {
			//attempt to load 4 lines
			int i=0;
			String[] first = new String[4];
			String[] second = new String[4];
			boolean failed = false;

			for (; i<4; i++){
				first[i] = firstFastqIn.readLine();
				second[i] = secondFastqIn.readLine();
				//check all are good
				if (first[i] == null || second[i] == null) {
					failed = true;
					if (first[i] == null && second[i] == null) allRead = true;
					break;
				}
			}
			
			//good load?
			if (failed == false) {
				reads.add(first);
				reads.add(second);
				numberFastqLoaded++;
			}
			//nope if failed
			else {
				closeIO();
				//all failed so this is the last of the fastq
				if (allRead) return true;
				//nope truncation!
				else Misc.printErrAndExit("\nError: looks like one of your fastq files is truncated!");
			}
		}
		return true;
	}

	public void parse() throws IOException {
		//make a UmiWorkers
		UmiWorker[] workers = new UmiWorker[numberThreads];
		
		//make workers, these self feed 
		IO.el("Launching "+numberThreads+" workers to load and process fastq...");
		ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
		for (int i=0; i< numberThreads; i++) {
			workers[i] = new UmiWorker(this);
			executor.execute(workers[i]);
		}
		executor.shutdown();
		//spins here until the executer is terminated, e.g. all threads complete
        while (!executor.isTerminated()) {}
		
		//look for errors and sum bad umis
		for (UmiWorker w :workers) {
			if (w.isFailed()) Misc.printErrAndExit("\nError: one of the FastqUmiTagger workers encountered an issue while running.");
			numberFastqWithFailedUmi += w.getNumUmiNFailures();
		}

	}
	
	synchronized void printReads(ArrayList<String> interlacedReads) throws IOException {
		for (String s: interlacedReads) System.out.println(s);
	}


	public void makeReaders() throws IOException{
		//readers
		firstFastqIn = IO.fetchBufferedReader(firstReadFastq);
		secondFastqIn = IO.fetchBufferedReader(secondReadFastq);
		if (umiFastq != null) umiFastqIn = IO.fetchBufferedReader(umiFastq);
	}

	public void closeIO() throws IOException{
		//close readers
		firstFastqIn.close();
		secondFastqIn.close();
		if (umiFastqIn != null) umiFastqIn.close();
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new FastqUmiTagger(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': firstReadFastq = new File(args[++i]); break;
					case 's': secondReadFastq = new File(args[++i]); break;
					case 'b': umiFastq = new File(args[++i]); break;
					case 'l': umiLength = Integer.parseInt(args[++i]); break;
					case 't': bpToTrim = Integer.parseInt(args[++i]); break;
					case 'c': numberThreads = Integer.parseInt(args[++i]); break;
					case 'n': maxNsInUMI = Integer.parseInt(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e) {
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		String params = "\n"+IO.fetchUSeqVersion()+" FastqUmiTagger2 Arguments: "+Misc.stringArrayToString(args, " ")+"\n";
		System.err.println(params);
		
		//check fields
		if (firstReadFastq == null || firstReadFastq.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your first read fastq file -> "+firstReadFastq);
		if (secondReadFastq == null || secondReadFastq.canRead() == false) Misc.printErrAndExit("\nError: cannot find or read your second read fastq file -> "+secondReadFastq);
		
		//umi file?
		if (umiFastq != null) {
			if (umiLength != -1 || bpToTrim != -1) Misc.printErrAndExit("\nError: cannot provide a UMI file AND set inline UMI length info");
			if (umiFastq.canRead()==false) Misc.printErrAndExit("\nError: cannot find or read your UMI fastq file -> "+umiFastq);
		}
		//must be inline
		else {
			if (umiLength == -1) Misc.printErrAndExit("\nError: if a UMI file is not provided, then please set a UMI length value for inline UMI extraction.");
			if (bpToTrim == -1) bpToTrim = umiLength;
		}
	}

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                            Fastq Umi Tagger: July 2025                           **\n" +
				"**************************************************************************************\n" +
				"Slightly faster version of the FastqBarcodTagger. Only outputs interlaced fastq.\n"+
				"Takes 2 or 3 fastq files (paired end reads and possibly a third containing unique \n"+
				"molecular indexes), appends the UMI and its quality scores to the fastq header and \n"+
				"writes out interlaced fastq records to standard out. For in line UMI datasets, provide\n"+
				"a length to extract the UMI from the beginning of each fastq. \n"+

				"\nOptions:\n"+
				"-f First fastq file, .gz/.zip OK.\n" +
				"-s Second fastq file, .gz/.zip OK.\n" +
				"-b UMI fastq file, .gz/.zip OK.    Skip for in line UMI extraction.\n" + 
				"-l For in line, length of UMI to parse from the beginning of each fastq read.\n"+
				"-t For in line, length to trim the fastq reads, defaults to -l\n"+
				"-c Number of cpu to use in parsing. Defaults to 1. More isn't necessarily better!\n"+
				"-n Maximum number of Ns in UMI to save.\n"+
				

				"\nExample: java -Xmx1G -jar pathToUSeq/Apps/FastqUmiTagger -f lob_1.fastq.gz\n" +
				"     -s lob_2.fastq.gz -l 3 -t 5 -c 1 -n 0 | bwa mem -p /ref/hg19.fa -c 5 \n\n" +

		"**************************************************************************************\n");

	}

	public int getChunkSize() {
		return chunkSize;
	}

	public File getUmiFastq() {
		return umiFastq;
	}

	public int getUmiLength() {
		return umiLength;
	}

	public int getBpToTrim() {
		return bpToTrim;
	}

	public int getMaxNsInUMI() {
		return maxNsInUMI;
	}
}
