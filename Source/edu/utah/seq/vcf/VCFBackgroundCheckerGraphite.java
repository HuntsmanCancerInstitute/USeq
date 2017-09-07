package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import edu.utah.seq.parsers.mpileup.MpileupSample;
import edu.utah.seq.parsers.mpileup.MpileupTabixLoaderGraphite;
import edu.utah.seq.vcf.GraphiteZScore.GraphiteResult;
import edu.utah.seq.vcf.xml.SimpleVcf;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Given a vcf file compares each record against a set of normal samples from a tabix indexed mpileup file and marks those with evidence of the variant in the normals.
 * tabix-0.2.6/bgzip ~/repro3Ct8Norm.mpileup
 * tabix-0.2.6/tabix -s 1 -b 2 -e 2 ~/repro3Ct8Norm.mpileup.gz
 * @author Nix*/
public class VCFBackgroundCheckerGraphite {
	
	//user defined fields
	private File[] vcfFiles;
	private File mpileup;
	private File saveDir;
	private File indexedFasta;
	private String optionalBashCmds = null;
	private int minBaseQuality = 20;
	private int minReadCoverage = 20;
	private int minNumSamples = 3;
	private boolean verbose = false;
	private boolean removeNonZScoredRecords = false;
	private double minimumZScore = 0;
	private double maxSampleAF = 0.3;
	private boolean replaceQualScore = false;
	private int numberThreads = 0;
	private File graphite = null;
	private File[] graphiteBams = null;
	private int bpProximalSnv = 125;
	
	//internal
	private static final int numVcfToProc = 100;
	private MpileupTabixLoaderGraphite[] loaders;
	private boolean deleteTempFiles = false;

	//working
	private File vcfFile;
	private File vcfForGraphite;
	private File vcfForMPileUp;
	private BufferedReader vcfIn;
	private PrintWriter vcfOut;
	private int numRecords = 0;
	private int numNotScored = 0;
	private int numFailingZscore = 0;
	private int numSaved = 0;
	private ArrayList<String> tooFewSamples = new ArrayList<String>();
	private HashMap<String, String> shortLongVcfGraphite = new HashMap<String, String>();
	
	//constructor
	public VCFBackgroundCheckerGraphite(String[] args) {
		long startTime = System.currentTimeMillis();
		try {	
			processArgs(args);
			
			//make Loaders
			loaders = new MpileupTabixLoaderGraphite[numberThreads];
			for (int i=0; i< loaders.length; i++) loaders[i] = new MpileupTabixLoaderGraphite(mpileup, this);
			

			//for each vcf file
			System.out.println("\nParsing vcf files...");
			for (int i=0; i< vcfFiles.length; i++){
				vcfFile = vcfFiles[i];
				resetWorkingParams();
				String name = Misc.removeExtension(vcfFile.getName());
				System.out.println(name);
				
				//split into isolated snvs for mpileup and others for graphite
				vcfForGraphite = new File(saveDir, name+"ForGraphiteTemp.vcf");
				vcfForMPileUp = new File(saveDir, name+"ForMpileUpTemp.vcf");
				int[] goodBad = parseOutIsolatedSnvVars(vcfForMPileUp, vcfForGraphite);
				
				File tempVcf = new File(saveDir, name+"_BKZed.vcf");
				if (deleteTempFiles) tempVcf.deleteOnExit();
				vcfOut = new PrintWriter(new FileWriter(tempVcf));
				createReaderSaveHeader(vcfForMPileUp);
				
				System.out.println("\tProcessing "+goodBad[0]+" SNVs with mpileup...");
				ExecutorService executor = Executors.newFixedThreadPool(numberThreads);
				for (MpileupTabixLoaderGraphite l: loaders) executor.execute(l);
				executor.shutdown();

				//spins here until the executer is terminated, e.g. all threads complete
				while (!executor.isTerminated()) {}

				//check loaders 
				for (MpileupTabixLoaderGraphite l: loaders) {
					if (l.isFailed()) throw new IOException("ERROR: File Loader issue! \n");
				}

				vcfIn.close();
				
				//run graphite 
				System.out.println("\tProcessing "+goodBad[1]+" INDELS with graphite (slow)...");
				GraphiteZScore g = new GraphiteZScore(vcfForGraphite, graphite, graphiteBams, indexedFasta, optionalBashCmds, minReadCoverage, minNumSamples, maxSampleAF);
				addGraphiteResults(g);
				
				//close out
				vcfOut.close();
				
				//sort, in memory so hopefully it doesn't blow up 
				File sorted = new File(tempVcf.getParentFile(), tempVcf.getName()+".gz");
				SimpleVcf.sortVcf(tempVcf, sorted);
				
				//print summary stats
				if (tooFewSamples.size() !=0){
					System.err.println("WARNING, the following did not have any mpileup lines (e.g. outside the compiled regions?) or enough background "
							+ "samples passing thresholds to calculate a z-score. ");
					System.err.println(Misc.stringArrayListToString(tooFewSamples, "\n"));
				}
				System.out.println("\n#Rec="+numRecords+" #Saved="+numSaved+" #NotScored="+numNotScored+" #FailingBKZ="+numFailingZscore+"\n");
			}
			
			//delete temp files
			if (deleteTempFiles)  {
				if (vcfForGraphite != null) vcfForGraphite.delete();
				if (vcfForMPileUp != null) vcfForMPileUp.delete();
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			//shut down loaders
			for (MpileupTabixLoaderGraphite m : loaders) m.getTabixReader().close();
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
			
		}
	}
	
	private void addGraphiteResults(GraphiteZScore g) throws IOException {
		HashMap<String, GraphiteResult> vcfGraphiteResults = g.getVcfGraphiteResults();
		
		if (vcfGraphiteResults.size() != shortLongVcfGraphite.size()) throw new IOException ("The number of graphite results isn't equal to the number of input vcf records?!");
		for (String shortVcf: vcfGraphiteResults.keySet()){
			GraphiteResult gr = vcfGraphiteResults.get(shortVcf);
			String oriVcf = shortLongVcfGraphite.get(shortVcf);
			String[] oriFields = Misc.TAB.split(oriVcf);
			if (oriVcf == null) throw new IOException ("Failed to find an original vcf record for this short graphite result?! "+shortVcf);
			//too few samples?
			if (gr == null) {
				tooFewSamples.add(oriVcf);
				numNotScored++;
				if (removeNonZScoredRecords == false) printScoredVcf(oriFields, oriVcf);
			}
			
			//append min zscore and save
			else {
				//fetch AFs
				String bkafs = fetchFormattedAFs(gr.getBkAFs());
				String bkzString = Num.formatNumberNoComma(gr.getzScore(), 2);
				oriFields[7] = "BKZ="+bkzString+";BKAF="+bkafs+";"+oriFields[7];
				if (replaceQualScore) oriFields[5] = bkzString;
				String modRecord = Misc.stringArrayToString(oriFields, "\t");
				//threshold it?
				if (minimumZScore == 0) {
					vcfOut.println(modRecord);
					numSaved++;
				}
				else {
					if (gr.getzScore() < minimumZScore) numFailingZscore++;
					else {
						vcfOut.println(modRecord);
						numSaved++;
					}
				}
			}
		}
	}
	
	public static String fetchFormattedAFs(double[] afs) {
		StringBuilder sb = new StringBuilder(Num.formatNumber(afs[0], 4));
		for (int i=1; i< afs.length; i++){
			sb.append(",");
			sb.append(Num.formatNumber(afs[i], 4));
		}
		return sb.toString();
	}
	
	private void printScoredVcf(String[] fields, String record) throws IOException{
		numSaved++;
		if (replaceQualScore) {
			fields[5] = "0";
			vcfOut.println(Misc.stringArrayToString(fields, "\t"));
		}
		else vcfOut.println(record);
	}

	private int[] parseOutIsolatedSnvVars(File forMpileup, File forGraphite) throws IOException{
		PrintWriter mpileupOut = new PrintWriter (new FileWriter(forMpileup));
		PrintWriter graphiteOut = new PrintWriter (new FileWriter(forGraphite));
		VCFParser p = new VCFParser(vcfFile, true, false, true);
		graphiteOut.println("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		mpileupOut.println(Misc.stringArrayListToString(p.getVcfComments().getOriginalComments(), "\n"));
		
		p.splitVCFRecordsByChromosome();
		HashMap<String, VCFLookUp> chrVcf = p.getChromosomeVCFRecords();
		int numGoodSnvs = 0;
		int numOthers = 0;
		//for each chromosome
		for (String chr: chrVcf.keySet()){
			VCFLookUp vlu = chrVcf.get(chr);
			//for each record
			VCFRecord[] records = vlu.getVcfRecord();
			for (VCFRecord r: records){
				if (r.isSNP()){
					int start = r.getPosition()-bpProximalSnv;
					if (start < 0) start = 0;
					VCFRecord[] ints = vlu.fetchVCFRecords(start, r.getPosition()+bpProximalSnv);
					//just itself?
					if (ints == null || ints.length <=1) {
						mpileupOut.println(r.getOriginalRecord());
						numGoodSnvs++;
					}
					else {
						String shortVcf = r.getTruncatedRecord();
						graphiteOut.println(shortVcf);
						shortLongVcfGraphite.put(shortVcf, r.getOriginalRecord());
						numOthers++;
					}
				}
				else {
					String shortVcf = r.getTruncatedRecord();
					graphiteOut.println(shortVcf);
					shortLongVcfGraphite.put(shortVcf, r.getOriginalRecord());
					numOthers++;
				}
			}
		}
		mpileupOut.close();
		graphiteOut.close();
		numRecords = numGoodSnvs + numOthers;
		return new int[]{numGoodSnvs, numOthers};
	}

	private void resetWorkingParams() {
		numRecords = 0;
		numNotScored = 0;
		numFailingZscore = 0;
		numSaved = 0;
		tooFewSamples.clear();
		shortLongVcfGraphite.clear();
	}

	/**Just prints out records.*/
	public synchronized void saveModifiedVcf(ArrayList<String> vcf) throws IOException{
		for (String s: vcf) vcfOut.println(s);
		numSaved+= vcf.size();
	}

	/**Loads the AL with vcf lines up to the numVcfToProc, if none could be read, returns false, otherwise true.*/
	public synchronized boolean loadVcfRecords(ArrayList<String> chunk) throws IOException{
		String record;
		int count = 0;
		while ((record = vcfIn.readLine()) != null){
			record = record.trim();
			if (record.length() == 0) continue;
			chunk.add(record);
			if (count++ > numVcfToProc) break;
		}
		
		if (chunk.size() == 0) return false;
		return true;
	}
	
	public synchronized void update(int numNotScored, int numFailingZscore, ArrayList<String> tooFewSamples) {
		this.numNotScored += numNotScored;
		this.numFailingZscore+= numFailingZscore;
		this.tooFewSamples.addAll(tooFewSamples);
	}
	
	private void createReaderSaveHeader(File forMpileupVcf) throws Exception {
		vcfIn = IO.fetchBufferedReader(forMpileupVcf);
		String record;
		boolean addInfo = true;
		boolean endFound = false;
		while ((record = vcfIn.readLine()) != null){
			record = record.trim();
			if (record.length() == 0) continue;
			if (record.startsWith("#")){
				if (addInfo && record.startsWith("##INFO=")) {
					vcfOut.println("##INFO=<ID=BKZ,Number=1,Type=Float,Description=\"Smallest AF z-score calculated from background AFs over effected bases. "
							+ "Values < ~4 are suspicous, non reference observations are likely present in the background samples.\">");
					vcfOut.println("##INFO=<ID=BKAF,Number=1,Type=String,Description=\"Non-reference AFs from the background samples used to calculate the BKZ.\">");
					addInfo = false;
				}
				vcfOut.println(record);
				if (record.startsWith("#CHROM")){
					endFound = true;
					break;
				}
			}
			else {
				endFound = true;
				break;
			}
		}
		if (endFound == false) Misc.printErrAndExit("\nERROR: failed to find the #CHROM line, aborting\n");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFBackgroundCheckerGraphite(args);
	}		
	
	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
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
					case 'v': forExtraction = new File(args[++i]); break;
					case 's': saveDir = new File(args[++i]); break;
					case 'm': mpileup = new File(args[++i]); break;
					case 'g': graphite = new File(args[++i]); break;
					case 'i': indexedFasta = new File(args[++i]); break;
					case 'o': optionalBashCmds = args[++i]; break;
					case 'n': graphiteBams = IO.extractFiles(new File(args[++i]), ".bam"); break;
					case 'f': maxSampleAF = Double.parseDouble(args[++i]); break;
					case 'z': minimumZScore = Double.parseDouble(args[++i]); break;
					case 'q': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'b': bpProximalSnv = Integer.parseInt(args[++i]); break;
					case 'c': minReadCoverage = Integer.parseInt(args[++i]); break;
					case 'a': minNumSamples = Integer.parseInt(args[++i]); break;
					case 'e': removeNonZScoredRecords = true; break;
					case 'd': verbose = true; break;
					case 'r': replaceQualScore = true; break;
					case 't': numberThreads = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//pull vcf files
		if (forExtraction == null || forExtraction.exists() == false) Misc.printErrAndExit("\nError: please enter a path to a vcf file or directory containing such.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.fetchFilesRecursively(forExtraction, ".vcf");
		tot[1] = IO.fetchFilesRecursively(forExtraction,".vcf.gz");
		tot[2] = IO.fetchFilesRecursively(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz OK) file(s)!\n");

		//save dir?
		if (saveDir == null) Misc.printExit("\nError: please provide a directory in which to save the marked vcf files\n");
		saveDir.mkdirs();
		
		//tabix indexed mpileup?
		if (mpileup == null || mpileup.canRead() == false) Misc.printExit("\nError: please provide a path to a bgzipped tabix indexed mpileup file.\n");
		File index = new File (mpileup.toString()+".tbi");
		if (index.exists() == false) Misc.printExit("\nError: cannot find the '"+index.getName()+"' index file corresponding to this indexed mpileup file "+mpileup.getName());	
		
		//threads to use
		int numAvail = Runtime.getRuntime().availableProcessors();
		if (numberThreads < 1 || numberThreads > numAvail) numberThreads =  numAvail - 1;
		
		//graphite
		if (graphite == null || graphite.canExecute() == false)  Misc.printExit("\nError: cannot find or execute your provided graphite app?! "+graphite);
		if (graphiteBams == null || graphiteBams.length == 0)  Misc.printExit("\nError: failed to find any background normal bam files for use with graphite?");
		OverdispersedRegionScanSeqs.lookForBaiIndexes(graphiteBams, false);
		if (indexedFasta == null || indexedFasta.exists() == false) Misc.printExit("\nError: cannot find your indexed fasta file used to align the normal bams?! "+indexedFasta);
		
		printSettings();
	}	
	
	public void printSettings(){
		System.out.println("Settings:\nNormal mpileup\t"+mpileup);
		System.out.println("Save dir\t"+saveDir);
		System.out.println("Normal bam dir\t"+ graphiteBams[0].getParent());
		System.out.println(minReadCoverage+"\tMin mpileup sample read coverage");
		System.out.println(minBaseQuality+"\tMin mpileup sample base quality");
		System.out.println(maxSampleAF+"\tMax mpileup sample AF");
		System.out.println(minNumSamples+"\tMin # samples for z-score calc");
		System.out.println(minimumZScore+"\tMin vcf AF z-score to save");
		System.out.println(bpProximalSnv + "\tBP intersection radius for pushing SNVs to graphite AF estimation ");
		System.out.println(numberThreads+"\tCPUs");
		System.out.println(removeNonZScoredRecords+ "\tExclude vcf records that could not be z-scored");
		System.out.println(verbose+"\tVerbose");
		System.out.println(replaceQualScore+"\tReplace QUAL score with z-score and set non scored records to 0");
	}

	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                         VCF Background Checker : April 2016                      **\n" +
				"**************************************************************************************\n" +
				"VBC calculates non-ref allele frequencies (AF) from a normal background multi-sample \n"+
				"mpileup file over each snv or using normal bams and graphite for indel and multi alts.\n"+
				"It then calculates a z-score for the vcf AF and appends it to the INFO field. Z-scores\n"+
				"< ~4 are indicative of non reference bps in the background samples. Note, VBC requires\n"+
				"an AF tag in the INFO field of each record. \n"+

				"\nRequired:\n"+
				"-v Path to a xxx.vcf(.gz/.zip OK) file or directory containing such.\n" +
				"-m Path to a bgzip compressed and tabix indexed multi-sample mpileup file. e.g.:\n"+
				"      1) Mpileup: 'echo \"#SampleOrder: \"$(ls *bam) > bkg.mpileup; samtools mpileup\n"+
				"             -B -q 13 -d 1000000 -f $fastaIndex -l $bedFile *.bam >> bkg.mpileup'\n"+
				"      2) Bgzip: 'tabix-0.2.6/bgzip bkg.mpileup'\n"+
                "         Tabix: 'tabix-0.2.6/tabix -s 1 -b 2 -e 2 bkg.mpileup.gz'\n"+
				"-g Path to graphite executable, see https://github.com/dillonl/graphite.\n"+
				"-n Directory containing normal background indexed bams for graphite AF estimation.\n"+
				"-i Indexed fasta file used to align the normal bams.\n"+
				"-s Path to a directory to save the modified vcf file(s)\n"+
						
				"\nOptional:\n" +
				"-z Minimum vcf z-score, defaults to 0, no filtering\n"+
				"-q Minimum mpileup sample bp quality, defaults to 20\n"+
				"-c Minimum sample read coverge, defaults to 20\n"+
				"-f Maximum sample AF, defaults to 0.3\n"+
				"-a Minimum # samples for z-score calculation, defaults to 3\n" +
				"-e Exclude vcf records that could not be z-scored\n"+
				"-r Replace QUAL value with z-score. Un scored vars will be assigned 0\n"+
				"-d Print verbose debugging output\n" +
				"-t Number of threads to use, defaults to all\n"+
				"-o Additional bash shell commands before executing graphite,e.g. 'module load gcc/4.9'\n"+
				"-b BP intersection radius for pushing SNVs to graphite for AF estimation, default 125\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/VCFBackgroundChecker -v SomaticVcfs/ -z 3\n"+
				"-m bkg.mpileup.gz -s BkgFiltVcfs/ -b 50 -q 13 -e -r -g ~/Graphite/bin/graphite -n\n"+
				"~/BkgrndNormalBams/ -o 'module load gcc/4.9.2' -i ~/B37/b37.fa \n\n"+

		        "**************************************************************************************\n");
	}

	//getters and setters
	public int getMinBaseQuality() {
		return minBaseQuality;
	}

	public int getMinReadCoverage() {
		return minReadCoverage;
	}

	public int getMinNumSamples() {
		return minNumSamples;
	}

	public boolean isRemoveNonZScoredRecords() {
		return removeNonZScoredRecords;
	}

	public double getMinimumZScore() {
		return minimumZScore;
	}

	public double getMaxSampleAF() {
		return maxSampleAF;
	}

	public boolean isReplaceQualScore() {
		return replaceQualScore;
	}

	public boolean isVerbose() {
		return verbose;
	}

}
