package edu.utah.seq.parsers.mpileup.concordance;

import java.io.*;
import java.util.regex.*;
import util.bio.annotation.Bed;
import util.gen.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class SampleConcordance {

	//user defined fields
	private File bedFile;
	private File bamPileupFile;
	private File commonSnvBed;
	private File gender;
	private int minSnvDP = 25;
	private double minAFForHom = 0.95;
	private double minAFForMatch = 0.90;
	private double minAFForHis = 0.05;
	private double minFracSim = 0.85;
	private File jsonOutputFile = null;
	private String toIgnoreForCall = "RNA";
	private double minMale = 2.5;
	private double maxFemale = 1.5;
	

	//internal fields
	private File tempDirectory;
	private int numberThreads = 0;
	private double maxIndel = 0.01;
	private Histogram[] afHist;
	private Histogram[] chrXAfHist;
	private SimilarityBamPileup[] similarities;
	private String[] sampleNames = null;
	private String[] similarityForJson = null;
	private String[] genderForJson = null;
	private String[] bamFileNamesForJson = null;
	private String bamNames = null;
	private String concordanceCheck = "";
	private String genderCheck = "";
	private Boolean clinicalGenderMatchesFemale = null;
	private String clinicalGenderBCComparison = "NA";
	
	//internal fields
	ConcordanceChunkBamPileup[] runners;
	

	
	//constructors
	public SampleConcordance(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);
			
			loadSampleNames();

			printThresholds();
			
			doWork();

			printStats();
			
			printGenderRatios();

			printHist();
			
			saveJson();
			
			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem running Bam Concordance Calculator!");
		}
	}

	private void loadSampleNames() throws IOException {
		//sample names
		/*
		# MinMapQual    13
		# MinBaseQual   10
		# IncludeOverlappingBpCounts false
		# PrintAll false
		# Bed /uufs/chpc.utah.edu/common/PE/hci-bioinformatics1/TNRunner/Bed/AllExonHg38Bed8April2020/hg38AllGeneExonsPad175bp.bed.gz
		# BamCram       0       /scratch/general/pe-nfs1/u0028003/Avatar/AJobs/1121013_13-0021350a/Alignment/1121013_13-0021350a_TumorRNA/Bam/1121013_13-0021350a_TumorRNA_Hg38.bam
		# BamCram       1       /scratch/general/pe-nfs1/u0028003/Avatar/AJobs/1121013_13-0021350a/Alignment/1121013_TumorDNA/Bam/1121013_TumorDNA_Hg38_final.bam
		# Chr   1BasePos        Ref     A,C,G,T,N,Del,Ins,FailBQ
		*/
		BufferedReader in = IO.fetchBufferedReader(bamPileupFile);
		String line = null;
		String[] fields = null;
		ArrayList<String> namesAL = new ArrayList<String>();
		while ((line = in.readLine())!= null) {
			line = line.trim();
			if (line.startsWith("# Chr")) break;
			if (line.startsWith("# BamCram")){
				fields = Misc.WHITESPACE.split(line);
				if (fields.length != 4) throw new IOException("ERROR: failed to parse the sample names from your bam pileup file "+bamPileupFile);
				namesAL.add(fields[3]);
			}
		}
		sampleNames = new String[namesAL.size()];
		bamFileNamesForJson = new String[sampleNames.length];
		StringBuilder sb = new StringBuilder();
		for (int i=0; i< sampleNames.length; i++) {
			String fullPathName = namesAL.get(i);
			fields = Misc.FORWARD_SLASH.split(fullPathName);
			String justName = fields[fields.length-1];
			if (justName.endsWith(".bam")) sampleNames[i] = justName.replace(".bam", "");
			else sampleNames[i] = justName.replace(".cram", "");
			bamFileNamesForJson[i] = justName;
			sb.append(fullPathName); sb.append(" ");
		}
		bamNames = sb.toString();
	}

	public void doWork() throws Exception{
		
		//attempt to pull gender
		if (gender != null) {
			clinicalGenderMatchesFemale = parseGenderFile(gender);
			IO.pl("\nClinically reported gender is female -> "+clinicalGenderMatchesFemale);
		}
		
		System.out.println("\nParsing and spliting bed files...");
		
		//parse regions and sort
		Bed[] regions = Bed.parseFile(bedFile, 0, 0);
		Arrays.sort(regions);
		int numRegionsPerChunk = 1+ (int)Math.round( (double)regions.length/ (double)numberThreads );
		
		//split regions
		ArrayList<Bed[]> splitBed = Bed.splitByNumber(regions, numRegionsPerChunk);
		
		//create runners and start
		System.out.println("\nLaunching "+splitBed.size() +" jobs...");
		runners = new ConcordanceChunkBamPileup[splitBed.size()];
		ExecutorService executor = Executors.newFixedThreadPool(runners.length);
		
		for (int i=0; i< runners.length; i++) {
			runners[i] = new ConcordanceChunkBamPileup(splitBed.get(i), this, ""+i);
			executor.execute(runners[i]);
		}
		executor.shutdown();
		//spins here until the executer is terminated, e.g. all threads complete
        while (!executor.isTerminated()) {}
		
		//check runners and pull bed gzippers
        ArrayList<File> toMerge = new ArrayList<File>();
        for (ConcordanceChunkBamPileup c: runners) {
			if (c.isFailed()) Misc.printErrAndExit("\nERROR: Failed runner, aborting! \n"+c.getCmd());
			toMerge.add(c.getMisMatchBed().getGzipFile());
		}
        
        //write out mismatch bed file
        File wd = new File (System.getProperty("user.dir"));
        
        File misMatchBed = new File(wd, "misMatch.bed.gz");
        IO.concatinateFiles(toMerge, misMatchBed);
        
        aggregateStats(runners);
        
		//delete temp dir and any remaining files
		IO.deleteDirectory(tempDirectory);
	}

	private void aggregateStats(ConcordanceChunkBamPileup[] runners) throws Exception {
        //collect stats
    	afHist = runners[0].getAfHist(); 	
    	chrXAfHist = runners[0].getChrXAfHist();
    	similarities = runners[0].getSimilarities();
    	    	
    	//for each runner
    	for (int i=1; i< runners.length; i++){
    		
    		//histograms
    		Histogram[] a = runners[i].getAfHist();
    		Histogram[] b = runners[i].getChrXAfHist();
    		if (a == null || b == null) continue;
    		for (int x=0; x< afHist.length; x++) {
    			if (a[x] != null) afHist[x].addCounts(a[x]);
    			if (b[x] != null) chrXAfHist[x].addCounts(b[x]);
    		}
    	
    		//similarities
    		SimilarityBamPileup[] s = runners[i].getSimilarities();
    		for (int x=0; x< s.length; x++) similarities[x].add(s[x]);
    		
    		runners[i].getMisMatchBed().getGzipFile().delete();
    	}
	}

	public void printThresholds(){
		IO.pl("Settings:");
		IO.pl(minSnvDP+"\tMinSnvDP");
		IO.pl(minAFForHom+"\tMinAFForHom");
		IO.pl(minAFForMatch+"\tminAFForMatch");
		IO.pl(minAFForHis+"\tMinAFForHis");
		IO.pl(maxIndel+"\tMaxIndel");
		IO.pl(minFracSim+"\tMinFractionSimilarity");
		IO.pl(maxFemale+"\tMaxFemale");
		IO.pl(minMale+"\tMinMale");
		IO.pl(toIgnoreForCall+"\tSample name to ignore in passing");
		if (commonSnvBed!= null) IO.pl(commonSnvBed.getName()+ "\tCommon SNV exclusion file");
		else IO.pl("All SNVs counted, common and uncommon.");
	}

	public void printStats(){
		similarityForJson = new String[similarities.length];
		//calculate max match and sort
		for (int i=0; i< similarities.length; i++) similarities[i].calculateMaxMatch();
		Arrays.sort(similarities);
		System.out.println("\nStats:");
		boolean failConcordance = false;
		boolean warnConcordance = false;
		for (int i=0; i< similarities.length; i++) {
			System.out.println(similarities[i].toString(sampleNames));
			similarityForJson[i] = similarities[i].toStringShort(sampleNames);
			String comp = similarities[i].passSimilarity(sampleNames);
			if (comp.equals("FAIL")) failConcordance = true;
			else if (comp.equals("WARNING")) warnConcordance = true; 
		}
		//set concordance check
		if (failConcordance) concordanceCheck = "FAIL";
		else if (warnConcordance) concordanceCheck = "WARNING";
		else concordanceCheck = "PASS";
		
		IO.pl(concordanceCheck+"\tConcordance Check (warning samples containing '"+toIgnoreForCall+"' that fail)\n");
	}

	public void printGenderRatios() throws IOException{
		genderForJson = new String[afHist.length];
		IO.pl("Het/Hom histogram AF count ratios for AllChrs, ChrX, log2(All/X)");
		boolean maleFound = false;
		boolean femaleFound = false;
		boolean indeterminant = false;
		for (int i=0; i< afHist.length; i++){
			double allChr = ratioCenterVsLast(afHist[i]);
			double chrX = ratioCenterVsLast(chrXAfHist[i]);
			double lgrto = Num.log2(allChr/chrX);
			//check gender?
			String genderCall = "\tSKIPPED";
			if (sampleNames[i].contains(toIgnoreForCall)==false && Double.isNaN(lgrto)==false && Double.isInfinite(lgrto)==false){
				if (lgrto > minMale) {
					maleFound = true;
					genderCall = "\tMALE";
				}
				else if (lgrto < maxFemale) {
					femaleFound = true;
					genderCall = "\tFEMALE";
				}
				else {
					indeterminant = true;
					genderCall = "\tINDETERMINANT";
				}
			}
			String g = sampleNames[i] +"\t"+Num.formatNumber(allChr, 3)+"\t"+Num.formatNumber(chrX, 3)+"\t"+Num.formatNumber(lgrto, 3)+genderCall;
			IO.pl(g);
			genderForJson[i] = Misc.TAB.matcher(g).replaceAll(" ");
		}
		//set for whole sample set
		if (indeterminant) genderCheck = "INDETERMINANT";
		else if (maleFound==true && femaleFound==true) genderCheck = "CONFLICTING";
		else if (maleFound) genderCheck = "MALE";
		else if (femaleFound) genderCheck = "FEMALE";
		else genderCheck = "NA";
		
		//compare to clinically reported? by default clinicalGenderBCComparison is NA
		if (clinicalGenderMatchesFemale != null) {
			if (genderCheck.equals("FEMALE")){
				if (clinicalGenderMatchesFemale) {
					clinicalGenderBCComparison = "PASS";
				}
				else clinicalGenderBCComparison = "FAIL";
			}
			else if (genderCheck.equals("MALE")) {
				if (clinicalGenderMatchesFemale == false) {
					clinicalGenderBCComparison = "PASS";
				}
				else clinicalGenderBCComparison = "FAIL";
			}
		}
		
		IO.pl("\n"+genderCheck+"\tGender check matches reported gender\t"+ clinicalGenderBCComparison);
	}
	/**
	 * Looks at each line in the file .gz/.zip OK for the word gender, case-insensitive ci
	 * Once found, looks to match female, male, f, m each subsequently, all ci*/
	public static Boolean parseGenderFile(File gender) throws IOException {
		//any gender file provided?
		if (gender == null) return null;
		
		Boolean clinicalMatchesFemale = null;
		
		//walk each line looking for female and male
		Pattern genderPat = Pattern.compile(".*gender.*", Pattern.CASE_INSENSITIVE);
		Pattern sexPat = Pattern.compile(".*sex.*", Pattern.CASE_INSENSITIVE);
		Pattern femalePat = Pattern.compile(".*female.*", Pattern.CASE_INSENSITIVE);
		Pattern malePat = Pattern.compile(".*male.*", Pattern.CASE_INSENSITIVE);
		Pattern fPat = Pattern.compile(".*f.*", Pattern.CASE_INSENSITIVE);
		Pattern mPat = Pattern.compile(".*m.*", Pattern.CASE_INSENSITIVE);
		
		BufferedReader in = IO.fetchBufferedReader(gender);
		String line;
		Matcher mat;
		Matcher sex;
		while ((line = in.readLine())!=null) {
			line = line.trim();
			if (line.startsWith("#")) continue;
			
			//gender
			mat = genderPat.matcher(line);
			sex = sexPat.matcher(line);
			if (mat.matches() || sex.matches()) {
				mat = femalePat.matcher(line);
				if (mat.matches()) clinicalMatchesFemale = true;
				else {
					mat = malePat.matcher(line);
					if (mat.matches()) clinicalMatchesFemale = false;
					else {
						mat = fPat.matcher(line);
						if (mat.matches()) clinicalMatchesFemale = true;
						else {
							mat = mPat.matcher(line);
							if (mat.matches()) clinicalMatchesFemale = false;
							else throw new IOException("Found a gender line but could not parse the actual gender from -> "+line);
						}
					}
				}
				break;
			}
		}
		in.close();
		return clinicalMatchesFemale;
	}

	/**Sums the num of middle and end bins, returns the ratio. Hard coded!*/
	public double ratioCenterVsLast(Histogram his){
		long[] counts = his.getBinCounts();
		
		//0.895 - 1.0
		double end = 0;
		end+= counts[22];
		end+= counts[23];
		end+= counts[24];
		if (end ==0) end++;
		
		//0.396 - 0.588
		double middle = 0;
		middle+= counts[9];
		middle+= counts[10];
		middle+= counts[11];
		middle+= counts[12];
		
		return middle/end;
	}
	
	public void printHist(){
		IO.pl("\nHistograms of AFs observed:");
		for (int i=0; i< afHist.length; i++){
			System.out.println("Sample\t"+sampleNames[i]);
			System.out.println("All NonRef AFs >= "+minAFForHis+":");
			afHist[i].printScaledHistogram();
			System.out.println("ChrX NonRef AFs >= "+minAFForHis+" - gender check:");
			chrXAfHist[i].printScaledHistogram();
			System.out.println();
		}
		IO.pl("\nHistograms of the AFs for mis matched homozygous variants - contamination check:");
		for (int i=0; i< similarities.length; i++) {
			similarities[i].toStringMismatch(sampleNames);
			System.out.println();
		}
	}
	
	@SuppressWarnings("unchecked")
	private void saveJson() {
		if (jsonOutputFile == null) return;
		
		try {
			//output simple json, DO NOT change the key names without updated downstream apps that read this file!
			Gzipper gz = new Gzipper(jsonOutputFile);
			gz.println("{");
			gz.printJson("sampleNames", sampleNames, true);
			gz.printJson("bamFileNames", bamFileNamesForJson, true);
			gz.printJson("similarities", similarityForJson, true);
			gz.printJson("concordanceCheck", concordanceCheck, true);
			gz.printJson("genderChecks", genderForJson, true);
			gz.printJson("genderCall", genderCheck, true);
			gz.printJson("genderComparison", clinicalGenderBCComparison, false);
			gz.println("}");
			gz.close();
			
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem writing json file! "+jsonOutputFile);
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SampleConcordance(args);
	}		
	
	/**This method will process each argument and assign new varibles
	 * @throws IOException */
	public void processArgs(String[] args) throws IOException{
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'r': bedFile = new File(args[++i]); break;
					case 'c': commonSnvBed = new File(args[++i]); break;
					case 'g': gender = new File(args[++i]); break;
					case 'b': bamPileupFile = new File(args[++i]); break;
					case 'd': minSnvDP = Integer.parseInt(args[++i]); break; 
					case 'a': minAFForHom = Double.parseDouble(args[++i]); break;
					case 'x': maxFemale = Double.parseDouble(args[++i]); break;
					case 'y': minMale = Double.parseDouble(args[++i]); break;
					case 'm': minAFForMatch = Double.parseDouble(args[++i]); break;
					case 'j': jsonOutputFile = new File(args[++i]); break;
					case 'n': minFracSim = Double.parseDouble(args[++i]); break;
					case 'e': toIgnoreForCall = args[++i]; break;
					case 't': numberThreads = Integer.parseInt(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//check bed
		if (bedFile == null || bedFile.canRead() == false) Misc.printErrAndExit("\nError: cannot find your bed file of regions to interrogate? "+bedFile);
		//check bam pileup
		if (bamPileupFile == null || bamPileupFile.getName().endsWith(".gz") == false) Misc.printErrAndExit("\nError: cannot find your tabix indexed and bgzip compressed multi sample BamPileup -> "+bamPileupFile);
		
		//threads to use
		int numPossThreads = Runtime.getRuntime().availableProcessors();
		if (numberThreads == 0) numberThreads = 20;
		if (numberThreads > numPossThreads) numberThreads = numPossThreads -1;
		
		System.out.println("Number threads to use:\t"+numberThreads+"\n");
		
				
		//make temp dir
		tempDirectory = new File ("TempDirBCC_"+Misc.getRandomString(6));
		tempDirectory.mkdir();
	
	}	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Sample Concordance: August 2021                        **\n" +
				"**************************************************************************************\n" +
				"SC calculates sample level concordance and gender based on uncommon homozygous SNVs in\n"+
				"alignment files. Those from the same person will show high similarity (>0.9). To run\n"+
				"SC, first create a multi-sample pileup file with the USeq BamPileup and BamPileupMerger\n"+
				"apps. Try SC with related bams (e.g tumor & normal exomes, tumor RNASeq) plus an\n"+
				"unrelated bam for omparison. Mismatches passing filters are written to file. BC also\n"+
				"generates a variety of AF histograms for checking gender and sample contamination.\n\n"+
				
				"Note re FFPE derived RNASeq data, a fair bit of systematic error is found in these\n"+
				"datasets.  As such, the RNA-> DNA contrasts are low. Yet the DNA->RNA are > 0.9\n\n"+

				"Options:\n"+
				"-r Path to a bed file of regions to interrogate.\n"+
				"-b Path to a multi-sample bam/ cram pileup file, see the USeq BamPileup and\n"+
				"       BamPileupMerger apps.\n"+
				"-c Path to a tabix indexed bed file of common dbSNPs. Download 00-common_all.vcf.gz \n"+
				"       from ftp://ftp.ncbi.nih.gov/snp/organisms/, grep for 'G5;' containing lines, \n"+
				"       run VCF2Bed, bgzip and tabix it with https://github.com/samtools/htslib,\n"+
				"       defaults to no exclusion from calcs. \n"+
				"-g Path to a txt file containing a line with the word 'gender' followed by 'male',\n"+
				"       'female', 'm', or 'f', surrounded by whitespace, case insensitive, gz/zip OK.\n"+
				"       This is compared to the BC gender estimate.\n"+
				"-d Minimum SNV read depth, defaults to 25.\n"+
				"-a Minimum allele frequency to count as a homozygous variant, defaults to 0.95\n"+
				"-m Minimum allele frequency to count a homozygous match, defaults to 0.9\n"+
				"-n Minimum fraction similarity to pass sample set, defaults to 0.85\n"+
				"-x Maximum log2Rto score for calling a sample female, defaults to 1.5\n"+
				"-y Minimum log2Rto score for calling a sample male, defaults to 2.5\n"+
				"-e Sample name to ignore in scoring similarity and gender, defaults to 'RNA'\n"+
				"-j Write gzipped summary stats in json format to this file.\n"+
				"-t Number of threads to use.  If not set, determines this based on the number of\n"+
				"      threads and memory available to the JVM so set the -Xmx value to the max.\n\n"+

				"Example: java -Xmx100G -jar pathTo/USeq/Apps/SampleConcordance -r ~/exomeTargets.bed\n"+
				"      -b ~/PatientZ.bp.txt.gz -d 30 -a 0.9 -m 0.8 -c ~/B38/b38ComSnps.bed.gz -j\n"+
				"      bc.json.gz -g gender.json.gz\n\n" +

				"**************************************************************************************\n");

	}

	public int getMinSnvDP() {
		return minSnvDP;
	}

	public double getMinAFForHom() {
		return minAFForHom;
	}

	public double getMinAFForMatch() {
		return minAFForMatch;
	}

	public double getMinAFForHis() {
		return minAFForHis;
	}

	public File getTempDirectory() {
		return tempDirectory;
	}

	public double getMaxIndel() {
		return maxIndel;
	}

	public File getCommonSnvBed() {
		return commonSnvBed;
	}

	public String getBamNames() {
		return bamNames;
	}

	public String getToIgnoreForCall() {
		return toIgnoreForCall;
	}

	public double getMinFracSim() {
		return minFracSim;
	}

	public File getBamPileupFile() {
		return bamPileupFile;
	}

	public String[] getSampleNames() {
		return sampleNames;
	}

	
}
