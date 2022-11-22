package edu.utah.seq.cnv.cfdna;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.*;
import htsjdk.samtools.*;
import htsjdk.samtools.cram.ref.ReferenceSource;
import trans.cel.QuantileNormalization;
import trans.main.WilcoxonSignedRankTest;
import util.bio.annotation.Bed;
import util.gen.*;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream; //needed by cram

/** 
 * @author Nix
 * */
public class LiquidBiopsyCA {

	//fields
	private File bedFile = null;
	private LiquidBiopsySample[] samples = null;
	private File resultsDirectory;
	private File fastaFile;
	private int minMappingQuality = 13;
	private SamReaderFactory samFactory;
	private LinkedHashMap<String, ArrayList<Bed>> groupNameBed = null;
	private int numCpu = 0;
	private ArrayList<LookupJob[]> jobs = new ArrayList<LookupJob[]>();
	private int[][] targetSampleCounts = null;
	private int totalNumberTargets = 0;

	//constructor for stand alone use
	public LiquidBiopsyCA(String[] args){
		try {
			long startTime = System.currentTimeMillis();

			processArgs(args);
			
			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Num.formatNumber(diffTime, 2)+" minutes");
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR running the copy analysis.");
		}
	}
	
	public void doWork() throws Exception {
		//create reader factory, try without the reference, no need to add in the sequence
		if (fastaFile!=null) samFactory = SamReaderFactory.makeDefault().referenceSource(new ReferenceSource(fastaFile)).validationStringency(ValidationStringency.SILENT);
		else samFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		
		//load and split the bed into gene name blocks
		parseBed();
		
		//make jobs
		makeLookupJobs();
		
		//make and execute loaders
		IO.pl("\nLaunching "+numCpu+" loaders...\n\tJobName\t#IntFragments\t#NonIntFragments");
		ExecutorService executor = Executors.newFixedThreadPool(numCpu);
		LookupJobRunner[] loaders = new LookupJobRunner[numCpu];
		for (int i=0; i< numCpu; i++) {
			loaders[i] = new LookupJobRunner(this);
			executor.execute(loaders[i]);
		}
		executor.shutdown();
		while (!executor.isTerminated()) {}
		
		//close sam file readers
		closeSamReaders();
		
		//check loaders fetch files
		for (LookupJobRunner l: loaders) {
			if (l.isFailed()) throw new IOException("\nERROR: Lookup Loader issue, aborting!");	
		}
		
//saveCountTable();
		saveCountsForR();
		
		//normalizeCounts();
		
		//wilcoxonTestPairedData();

	}
	
	
	private void wilcoxonTestPairedData() {
		IO.pl("\nTesting each gene for copy differences...");
		LookupJob[][] germlines = new LookupJob[samples.length][];
		int normIndex = 0;
		for (LiquidBiopsySample s: samples) {
			IO.pl(s.getSampleName());
			LookupJob[] cf = s.getCfJobs();
			LookupJob[] norm = s.getGermlineJobs();
			germlines[normIndex++] = norm;
			for (int i=0; i<cf.length; i++) {
				String geneName = cf[i].getRegions().get(0).getName();
				IO.pl("\t"+geneName);
				IO.pl("\t\tCF  : "+Num.floatArrayToString(cf[i].getNormalizedCounts()," "));
				IO.pl("\t\tNorm: "+Num.floatArrayToString(norm[i].getNormalizedCounts()," "));
				WilcoxonSignedRankTest w = new WilcoxonSignedRankTest(cf[i].getNormalizedCounts(), norm[i].getNormalizedCounts());
				IO.pl("\t\tPVal: "+w.getTransformedPValue());
				
			}
		}
		IO.pl("\nTesting each gene for mock copy differences between germline samples...");
		for (int k=0; k< germlines.length; k++) {
			LookupJob[] normA = germlines[k];
			for (int j=k+1; j< germlines.length; j++) {
				LookupJob[] normB = germlines[j];
				for (int i=0; i<normA.length; i++) {
					String geneName = normA[i].getRegions().get(0).getName();
					IO.p("\t"+geneName);
					//IO.pl("\t\tCF  : "+Num.floatArrayToString(cf[i].getGermlineizedCounts()," "));
					//IO.pl("\t\tNorm: "+Num.floatArrayToString(norm[i].getGermlineizedCounts()," "));
					WilcoxonSignedRankTest w = new WilcoxonSignedRankTest(normA[i].getNormalizedCounts(), normB[i].getNormalizedCounts());
					IO.pl("\tPVal: "+w.getTransformedPValue());
					
				}
			}
		}
		
	}

	private void normalizeCounts() {
		IO.pl("\nQuantile germlineizing target counts...");
		// TODO Auto-generated method stub
		float[][] intensities = new float[samples.length*2][];
		int index = 0;
		for (LiquidBiopsySample s: samples) {
			intensities[index++] = s.getCfCounts(totalNumberTargets);
			intensities[index++] = s.getGermlineCounts(totalNumberTargets);
		}
		QuantileNormalization qn = new QuantileNormalization (intensities);
		intensities = qn.getIntensities();
		addGermlineizedCountsToLookupJobs(intensities);
		//printIntensityStats(intensities);
		
	}
	
	private void addGermlineizedCountsToLookupJobs(float[][] intensities) {
		int index = 0;
		//add back to the LookupJob
		for (LiquidBiopsySample s: samples) {
			LookupJob[] cf = s.getCfJobs();
			LookupJob[] norm = s.getGermlineJobs();
			float[] normCf = intensities[index++];
			float[] normNorm = intensities[index++];
			int counter = 0;
			for (int i=0; i< cf.length; i++) {
				int numTargetsInJob = cf[i].getRegions().size();
				int stop = counter+numTargetsInJob;
				float[] nCf = Arrays.copyOfRange(normCf, counter, stop);
				float[] nNorm = Arrays.copyOfRange(normNorm, counter, stop);
				cf[i].setNormalizedCounts(nCf);
				norm[i].setNormalizedCounts(nNorm);
				counter = stop;
			}
		}
	}

	/**Prints statistics about the current intensity float[][] arrays*/
	public void printIntensityStats(float[][] intensities){
		for (int i=0; i<intensities.length; i++){
			Num.statFloatArray(intensities[i], false);
		}
	}

	private void saveCountTable() {
		IO.pl("\nSaving count table...");
		//create rowNames
		String[] rowNames = new String[totalNumberTargets];
		int counter = 0;
		for (String blockName: groupNameBed.keySet()) {
			for (Bed bed: groupNameBed.get(blockName)) {
				rowNames[counter++] = blockName+":"+bed.toStringCoordinates();
			}
		}
		
		//create columnNames
		String[] columnNames = new String[(samples.length*2) + 1];
		counter = 0;
		columnNames[counter++] = "Target:Coordinates";
		for (LiquidBiopsySample l: samples) {
			columnNames[counter++] = l.getSampleName()+":Somatic";
			columnNames[counter++] = l.getSampleName()+":Germline";
		}
		
		//create int[][] table
		targetSampleCounts = new int[totalNumberTargets][samples.length*2];
		
		//for each targetBlock
		int numGroup = groupNameBed.keySet().size();	
		int rowCounter = 0;
		for (int t = 0; t< numGroup; t++) {
			//for each sample set
			int numRowsInBlock = samples[0].getCfJobs()[t].getRegions().size();
			int columnCounter = 0;
			for (int s = 0; s< samples.length; s++) {
				int blockRowCounter = rowCounter;
				LookupJob tumor = samples[s].getCfJobs()[t];
				LookupJob germline = samples[s].getGermlineJobs()[t];
				for (int r=0; r< numRowsInBlock; r++) {
					targetSampleCounts[blockRowCounter][columnCounter] = tumor.getCounts()[r];
					targetSampleCounts[blockRowCounter++][columnCounter+1] = germline.getCounts()[r];
				}
				columnCounter+=2;
			}
			rowCounter+= numRowsInBlock;
		}
		
		//print it
		StringBuilder sb = new StringBuilder();
		
		//column names
		for (String cn: columnNames) {
			sb.append(cn);
			sb.append("\t");
		}
		sb.append("\n");
		
		//rows
		for (int i=0; i< totalNumberTargets; i++) {
			sb.append(rowNames[i]);
			sb.append("\t");
			
			int[] rowCounts = targetSampleCounts[i];
			for (int rc: rowCounts) {
				sb.append(rc);
				sb.append("\t");
			}
			sb.append("\n");
		}
		
		//IO.pl(sb);
		
	}
	
	private void saveCountsForR() throws IOException {
		IO.pl("\nSaving Counts For R...");
		File forR = new File(resultsDirectory, "countsForR.txt");
		PrintWriter out = new PrintWriter( new FileWriter(forR));
		for (LiquidBiopsySample s: samples) {
			s.addCfCounts(out);
			s.addGermlineCounts(out);
		}
		out.close();
	}

	private void closeSamReaders() throws IOException {
		for (LiquidBiopsySample s: samples) {
			s.getCfReader().close();
			s.getGermlineReader().close();
			s.getCfBedOut().close();
			s.getGermlineBedOut().close();
		}
	}

	public synchronized LookupJob[] fetchNextJob() {
		if (jobs.size() == 0) return null;
		return jobs.remove(0);
	}

	private void makeLookupJobs() throws IOException {
		for (LiquidBiopsySample s: samples) {
			jobs.add(s.buildCFJobs(groupNameBed, samFactory, minMappingQuality));
			jobs.add(s.buildGermlineJobs(groupNameBed, samFactory, minMappingQuality));
		}
		
	}

	private void parseBed() {
		IO.pl("Parsing the bed file...");
		Bed[] targetRegions = Bed.parseFile(bedFile, 0, 0);
		totalNumberTargets = targetRegions.length;
		Arrays.sort(targetRegions);
		groupNameBed = new LinkedHashMap<String, ArrayList<Bed>>();
		for (Bed b: targetRegions) {
			ArrayList<Bed> beds = groupNameBed.get(b.getName());
			if (beds == null) {
				beds = new ArrayList<Bed>();
				groupNameBed.put(b.getName(), beds);
			}
			beds.add(b);
		}
		IO.pl("\tCopyGroup\tNumberRegionsInGroup");
		for (String name: groupNameBed.keySet()) {
			int numBeds = groupNameBed.get(name).size(); 
			IO.pl("\t"+ name + "\t" + numBeds);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new LiquidBiopsyCA(args);
	}		
	
	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args) throws Exception{
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 's': forExtraction = new File(args[++i]).getCanonicalFile(); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'r': fastaFile = new File(args[++i]); break;
					case 'o': resultsDirectory = new File(args[++i]); break;
					case 'm': minMappingQuality = Integer.parseInt(args[++i]); break;
					case 'p': numCpu = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//create sample objects, will throw exceptions 
		if (forExtraction == null || forExtraction.exists() == false || forExtraction.isDirectory()== false) Misc.printErrAndExit("\nError: please enter a path to directory containing one or more sample specific sub directories.\n");
		File[] sampleDirectories = IO.extractOnlyDirectories(forExtraction);
		samples = new LiquidBiopsySample[sampleDirectories.length];
		for (int i=0; i< samples.length; i++) samples[i] = new LiquidBiopsySample(sampleDirectories[i], resultsDirectory);
		
		//if (fastaFile == null || fastaFile.canRead() == false)  Misc.printErrAndExit("\nError: please provide an reference genome fasta file and it's index.");		
		if (bedFile == null ||  bedFile.canRead() == false) Misc.printErrAndExit("\nError: please provide a file of regions in bed format.");
		if (resultsDirectory == null ) Misc.printErrAndExit("\nError: please provide a directory for saving the output results");
		resultsDirectory.mkdirs();
			
		//number of workers
		int numProc = Runtime.getRuntime().availableProcessors();
		if (numCpu == 0 || numCpu > numProc) numCpu = numProc;
		int numJobs = samples.length*2;
		if (numCpu> numJobs) numCpu = numJobs;

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             LiquidBiopsyCopyAnalysis:  Nov 2022                  **\n" +
				"**************************************************************************************\n" +
				"\n"+
				
				"\nRequired Options:\n"+
				"-s Path to a directory containing sample specific sub directories with two aligment files representing the cfDNA and matched germline cram or bam files with their indexes. The germline file name should contain the 'norm' string.\n"+
				"-b Path to a bed file of non overlapping regions to estimate copy ratio differences, xxx.bed.gz/.zip OK. The name column will be used to group each region into a larger area for merged copy ratio differences.\n"+
				"-r Path to the reference fasta with xxx.fai index used for alignment.\n"+
				"-o Path to the output results directory.\n"+

				"\nDefault Options:\n"+
				"-m Minimum alignment mapping quality, defaults to 13\n"+
				"-p Number processors to use, defaults to all, reduce if out of memory errors occur.\n"+


				"\nExample: java -Xmx100G -jar pathTo/USeq/Apps/BamPileup \n\n" +

				"**************************************************************************************\n");
	}

}
