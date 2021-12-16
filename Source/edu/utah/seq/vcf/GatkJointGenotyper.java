package edu.utah.seq.vcf;

import java.io.*;
import java.util.regex.*;
import util.bio.annotation.Bed;
import util.gen.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * @author david.nix@hci.utah.edu 
 * Takes a bed file of desired regions, splits it, writes out each, then calls GATK GenomicsDBImport and GenotypeGVCFs, finally merges the results.
 * GATK really should thread it's own tools!!!!
 **/
public class GatkJointGenotyper {

	//user defined fields
	private File bedFile = null;
	private File gvcfDirectory = null;
	private File tmpDirectory = null;
	private File gatkExecutable = null;
	private File fasta = null;
	private File vcfResults = null;
	private double numberGBPerChunk = 10.0;
	private boolean deleteTmpDir = true;
	
	//internal fields
	private int numberThreads = 0;
	private GatkJointGenotyperChunk[] runners;
	private File sampleMap = null;
	private File[] splitBed = null;
	private int splitBedIndex = -1;
	private int numChunks = 250;
	private int maxBpOfRegion = 1000;
	private boolean keepRunning = true;
	
	//constructors
	public GatkJointGenotyper(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running GatkJointGenotyper!");
		}
	}

	public void doWork() throws Exception{
		
		//parse regions 
		System.out.println("Parsing, spliting large regions, chunking bed files...");
		splitBed();
		
		//write out a sample map
		System.out.println("Writing sample map...");
		writeSampleMap();
		
		//create runners and start
		System.out.println("Launching...\n");
		runners = new GatkJointGenotyperChunk[numberThreads];
		ExecutorService executor = Executors.newFixedThreadPool(runners.length);
		for (int i=0; i< runners.length; i++) {
			runners[i] = new GatkJointGenotyperChunk(this, ""+i);
			executor.execute(runners[i]);
		}
		executor.shutdown();
		//spins here until the executer is terminated, e.g. all threads complete or one errors
        while (!executor.isTerminated() && keepRunning) {}
		
		//check runners
        for (GatkJointGenotyperChunk c: runners) {
			if (c.isFailed()) {
				File failedTmp = new File(gvcfDirectory.getParentFile(), "FailedTmpDir");
				IO.copyDirectoryRecursive(tmpDirectory, failedTmp, null);
				Misc.printErrAndExit("\nERROR: Failed runner, copying over TmpDir to "+failedTmp+", aborting! \n"+c.getCmd());
			}
		}
		
		//merge vcf
		System.out.println("\nMerging split vcfs...");
		mergeSplitVcfs();
		
		//delete tmp
		if (deleteTmpDir) IO.deleteDirectory(tmpDirectory);
	}
	
	public synchronized File getNextBed() {
		splitBedIndex++;
		if (splitBedIndex < splitBed.length) return splitBed[splitBedIndex];
		return null;
	}

	private void splitBed() {
		Bed[] regions = Bed.parseFile(bedFile, 0, 0);
		
		//sort them just in case the user provides an unsorted list
		Arrays.sort(regions);
	
		int numInitialRegions = regions.length;
		
		//split big regions, this will come back sorted
		regions = Bed.splitBigRegions(regions, maxBpOfRegion);
		
		int numRegionsPerChunk = (int)Math.round((double)regions.length/(double)numChunks)+1;
		
		ArrayList<Bed[]> splitBedAL = Bed.splitByNumber(regions, numRegionsPerChunk);
		splitBed = new File[splitBedAL.size()];
		
		File bedDir = new File(tmpDirectory, "SplitBeds");
		bedDir.mkdirs();
		for (int i=0; i< splitBed.length; i++) {
			//write out bed
			splitBed[i] = new File (bedDir, i+".bed");
			Bed.writeToFile(splitBedAL.get(i), splitBed[i]);
		}
		IO.pl("\t"+numInitialRegions+" initial regions split into "+regions.length+" regions and "+splitBed.length+" job files with ~"+numRegionsPerChunk+" each.");
		regions = null;
	}

	private void writeSampleMap() throws IOException {
		//pull full path names
		File[] gvcfFiles = IO.extractFiles(gvcfDirectory, ".vcf.gz");
		if (gvcfFiles == null || gvcfFiles.length == 0) throw new IOException("\nERROR: failed to find one or more xxx.vcf.gz files in "+gvcfDirectory);
		String[] trimmedNames = new String[gvcfFiles.length];
		for (int x=0; x<gvcfFiles.length; x++) trimmedNames[x] = gvcfFiles[x].getName();
		trimmedNames = Misc.trimCommon(trimmedNames);
		sampleMap = new File (tmpDirectory, "sampleMap.txt");
		PrintWriter out = new PrintWriter (new FileWriter(sampleMap));
		int num = 0;
		for (int i=0; i< gvcfFiles.length; i++) {
			String line = trimmedNames[i]+"\t"+gvcfFiles[i].getCanonicalPath();
			out.println(line);
			num++;
		}
		out.close();
		IO.pl("\t"+num+" gVCFs");
	}

	private void mergeVcfs() throws Exception{
		ArrayList<VCFRecord> recordsAL = new ArrayList<VCFRecord>();
		VCFParser vp = null;
		HashSet<String> unique = new HashSet<String>();
		for (GatkJointGenotyperChunk c: runners) {
			for (File v:  c.getVcfs()) {
				vp = new VCFParser(v, true, false, false);
				VCFRecord[] records = vp.getVcfRecords();
				if (records == null) continue;
				for (VCFRecord r: records) {
					String key = r.getChromosome()+"_"+r.getPosition()+"_"+r.getReference()+"_"+r.getAlternate();
					if (unique.contains(key) == false) {
						recordsAL.add(r);
						unique.add(key);
					}
				}
			}
		}
		VCFRecord[] all = new VCFRecord[recordsAL.size()];
		recordsAL.toArray(all);
		Arrays.sort(all);

		Gzipper out = new Gzipper(vcfResults);

		//write out header
		ArrayList<String> header = vp.getVcfComments().getOriginalComments();
		for (String h : header) out.println(h);

		//write out records
		for (VCFRecord r : all) out.println(r.getOriginalRecord());
		out.close();
		System.out.println("\t"+unique.size()+"\tvariants found");
	}
	
	private void mergeSplitVcfs() throws Exception {
		Gzipper out = new Gzipper(vcfResults);
		HashSet<String> prior = new HashSet<String>();
		long uniqueRecords = 0;
		long totalRecords = 0;
		for (int i=0; i< splitBed.length; i++) {
			HashSet<String> current = new HashSet<String>();
			String name = splitBed[i].getName().replace(".bed", "_temp.vcf");
			File vcfFile = new File (tmpDirectory, name);
			if (vcfFile.exists() == false) throw new IOException("ERROR: failed to find the tmp vcf for chunk "+vcfFile);
			BufferedReader in = new BufferedReader( new FileReader(vcfFile));
			String line = null;
			String[] t = null;
			while ((line = in.readLine())!= null) {
				if (line.startsWith("#")) {
					if (i==0) out.println(line);
				}
				else {
					totalRecords++;
					t = Misc.TAB.split(line);
					//chr = t[0]; pos = t[1]; ref = t[3]; alt = t[4];
					String key = t[0]+"_"+t[1]+"_"+t[3]+"_"+t[4];
					if (prior.contains(key) == false) {
						out.println(line);
						current.add(key);
						uniqueRecords++;
					}
				}
			}
			in.close();
			prior = current;
		}
		
		out.close();
		System.out.println("\t"+uniqueRecords+"\tunique variants found in "+totalRecords+ " processed records.");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new GatkJointGenotyper(args);
	}		
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		File userTmpDir = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bedFile = new File(args[++i]); break;
					case 'g': gvcfDirectory = new File(args[++i]); break;
					case 'v': vcfResults = new File(args[++i]); break;
					case 'r': numberGBPerChunk = Double.parseDouble(args[++i]); break;
					case 't': userTmpDir = new File(args[++i]); break;
					case 'e': gatkExecutable = new File(args[++i]); break;
					case 'c': numChunks = Integer.parseInt(args[++i]); break;
					case 'f': fasta = new File(args[++i]); break;
					case 'd': deleteTmpDir = false; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//threads to use
		double gbAvail = ((double)Runtime.getRuntime().maxMemory())/ 1073741824.0 * 0.9;
		int numChunks = (int)(gbAvail/numberGBPerChunk);
		int procAvail = Runtime.getRuntime().availableProcessors();
		if (numChunks > procAvail) numberThreads = procAvail;
		else numberThreads = numChunks;
		
		System.out.println("Core usage:\n\tGB available to Java * 0.9 :\t"+ Num.formatNumber(gbAvail, 1));
		System.out.println("\tAvailable cores:\t"+procAvail);
		System.out.println("\tNumber jobs to launch with a minimum "+numberGBPerChunk+"GB/Job -> "+numberThreads+"\n");
		
		//check bed
		if (bedFile == null || bedFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find your bed file of regions to interrogate?\n"+bedFile);
		//check gvcfDirectory
		if (gvcfDirectory == null || gvcfDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: cannot find your directory of GVCF files to interrogate?\n"+gvcfDirectory);	
		//check vcfResults
		if (vcfResults == null ) Misc.printErrAndExit("\nError: cannot find your final output vcf file\n"+vcfResults);	
		//check userTmpDir
		if (userTmpDir == null ) Misc.printErrAndExit("\nError: cannot find your tmp directory?\n"+userTmpDir);
		userTmpDir.mkdirs();
		tmpDirectory = new File(userTmpDir, "GatkJointGenotyperTmpDir");
		tmpDirectory.mkdirs();	
		//check gatkExecutable
		if (gatkExecutable == null || gatkExecutable.canExecute() == false) Misc.printErrAndExit("\nError: cannot find the gatk executable?\n"+gatkExecutable);
		//check fasta
		if (fasta == null || fasta.canRead() == false) Misc.printErrAndExit("\nError: cannot find your indexed fasta reference?\n"+fasta);
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             Gatk Joint Genotyper: Dec 2021                       **\n" +
				"**************************************************************************************\n" +
				"The GJG takes a bed file of target regions, splits it into 250 jobs, executes the GATK\n"+
				"GenomicsDBImporter and GenotypeGVCFs on each, and merges the results. Set the -Xmx for\n"+
				"this app to the maximum available on the machine to enable correct thread usage.\n"+

				"\nOptions:\n"+
				"-b Bed file of sorted regions (chr, start, stop,... gz/zip OK), to call variants. No\n"+
				"      overlaps. Run the USeq MergeRegions if unsure.\n"+
				"-g Directory containing GATK xxx.g.vcf.gz files to call, bgzippped, tabix indexed.\n"+
				"-v Path to write the final results joint genotyped vcf file, must end with vcf.gz\n"+
				"-t Tmp directory.\n"+
				"-e Path to the gatk executable.\n"+
				"-f Path to the indexed fasta reference.\n"+
				"-r GB RAM per worker thread, defaults to 10, increase if memory errors occur.\n"+
				"-c Number of chunks to split bed file, defaults to 100.\n"+
				"-d Don't delete the tmp directory.\n"+

				"\nExample: java -Xmx128G -jar pathToUSeq/Apps/GatkJointGenotyper\n"+
				"-b hg38NimIdtMergedPad150bp.bed.gz -g ToGenotype/ -v jg.vcf.gz -r 20 \n"+
				"-t /scratch/local/u0028003/4239066/ -e ~/BioApps/GATK/gatk-4.2.3.0/gatk\n"+
				"-f ~/TNRunner/GATKResourceBundleAug2021/Homo_sapiens_assembly38.fasta\n\n"+

				"**************************************************************************************\n");

	}

	public File getTmpDirectory() {
		return tmpDirectory;
	}

	public File getGatkExecutable() {
		return gatkExecutable;
	}

	public File getFasta() {
		return fasta;
	}

	public File getSampleMap() {
		return sampleMap;
	}

	public double getNumberGBPerChunk() {
		return numberGBPerChunk;
	}

	public synchronized void setKeepRunning(boolean keepRunning) {
		this.keepRunning = keepRunning;
	}
}
