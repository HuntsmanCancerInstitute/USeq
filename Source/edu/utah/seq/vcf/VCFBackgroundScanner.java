package edu.utah.seq.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.parsers.mpileup.MpileupLine;
import edu.utah.seq.parsers.mpileup.MpileupSample;
import edu.utah.seq.query.QueryIndexFileLoader;
import htsjdk.tribble.readers.TabixReader;
import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Given a vcf file compares each record against a set of normal samples from a tabix indexed mpileup file and marks those with evidence of the variant in the normals.
 * tabix-0.2.6/bgzip ~/repro3Ct8Norm.mpileup
 * tabix-0.2.6/tabix -s 1 -b 2 -e 2 ~/repro3Ct8Norm.mpileup.gz
 * @author Nix*/
public class VCFBackgroundScanner {
	
	//user defined fields
	private File[] vcfFiles;
	private File mpileup;
	private File saveDir;
	private HashSet<Integer> sampleIndexesToExamine = null;
	//private File fullPathToR = new File ("/usr/bin/R");
	private int bpBuffer = 0;
	private int minBaseQuality = 20;
	private int minReadCoverage = 20;
	private int minNumSamples = 3;
	private boolean verbose = false;
	private boolean removeNonZScoredRecords = false;
	private double minimumZScore = 0;
	private double maxSampleAF = 0.3;
	
	//internal
	private TabixReader reader = null;
	private static final Pattern DP = Pattern.compile(".*DP=(\\d+).*");
	private static final Pattern AF = Pattern.compile(".*AF=([\\d\\.]+).*");
	private double zscoreForInfinity = 100;
	private NumberFormat fourDecimalMax = NumberFormat.getNumberInstance();
	
	
	//working
	private File vcfFile;
	private Gzipper vcfOut;
	private int numRecords = 0;
	private int numNotScored = 0;
	private int numFailingZscore = 0;
	private int numSaved = 0;
	private ArrayList<String> tooFewSamples = new ArrayList<String>();
	
	//constructor
	public VCFBackgroundScanner(String[] args){
		try {	
			processArgs(args);
			
			//tabix reader
			reader = new TabixReader(mpileup.toString());

			//for each vcf file
			System.out.println("\nParsing vcf files...");
			for (int i=0; i< vcfFiles.length; i++){
				vcfFile = vcfFiles[i];
				numRecords = 0;
				numNotScored = 0;
				numFailingZscore = 0;
				numSaved = 0;
				tooFewSamples.clear();
				System.out.println(vcfFile.getName());
				String name = Misc.removeExtension(vcfFile.getName());
				vcfOut = new Gzipper(new File(saveDir, name+"_BKZed.vcf.gz"));
				scoreWorkingVcfFile();
				vcfOut.close();
				
				//print summary stats
				if (tooFewSamples.size() !=0){
					System.err.println("WARNING, the following did not have any mpileup lines (e.g. outside the compiled regions?) or enough background "
							+ "samples passing thresholds to calculate a z-score. ");
					System.err.println(Misc.stringArrayListToString(tooFewSamples, "\n"));
				}
				System.out.println("\t#Rec="+numRecords+" #Saved="+numSaved+" #NotScored="+numNotScored+" #FailingBKZ="+numFailingZscore+"\n");
			}

		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			reader.close();
		}
	}

	private void scoreWorkingVcfFile() throws Exception {
		BufferedReader in = IO.fetchBufferedReader(vcfFile);
		String record;
		boolean addInfo = true;
		while ((record = in.readLine()) != null){
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
			}
			else {
				numRecords++;
				score(record);
			}
		}
		in.close();
	}

	private void score(String record) throws Exception {
		//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
		String[] fields = Misc.TAB.split(record);
		
		//fetch DP and AF
		Matcher mat = AF.matcher(fields[7]);
		if (mat.matches() == false) throw new IOException ("Failed to parse AF= number from the INFO field in this variant:\n"+record);
		double freq = Double.parseDouble(mat.group(1));
		//mat = DP.matcher(fields[7]);
		//if (mat.matches() == false) throw new IOException ("Failed to parse DP= number from the INFO field in this variant:\n"+record);
		//double depth = Double.parseDouble(mat.group(1));
		
		//interbase coor of effected bps
		int[] startStop = QueryIndexFileLoader.fetchEffectedBps(fields, true);
		if (startStop == null) {
			System.err.println("Failed to parse the effected bps.");
			if (removeNonZScoredRecords == false) {
				vcfOut.println(record);
				numSaved++;
			}
			tooFewSamples.add(record);
			numNotScored++;
			return;
		}
		
		//pull mpileup records, if none then warn and save vcf record
		int start = startStop[0]+1-bpBuffer;
		if (start < 1) start = 1;
		String tabixCoor = fields[0]+":"+start+"-"+(startStop[1]+bpBuffer);
		TabixReader.Iterator it = fetchInteratorOnCoordinates(tabixCoor);
		if (it == null) {
			if (removeNonZScoredRecords == false) {
				vcfOut.println(record);
				numSaved++;
			}
			tooFewSamples.add(record);
			numNotScored++;
			
			return;
		}
		
		if (verbose) System.out.println("VCF\t"+record);
		
		//for each mpileup record, find lowest z-score
		String mpileupLine = null;
		double minZScore = Double.MAX_VALUE;
		MpileupSample[] minZScoreSamples = null;
		while ((mpileupLine = it.next()) != null){
			
			//parse mpilup line and pull filtered set of samples
			MpileupLine ml = new MpileupLine(mpileupLine, minBaseQuality);
			if (ml.getChr() == null) throw new IOException ("Failed to parse the mpileup line:\n"+mpileupLine);
			MpileupSample[] toExamine = fetchMpileupSamples(ml);
			if (toExamine.length < minNumSamples) continue;
			
			//calculate zscore
			double[] meanStd = calcMeanStdev(toExamine);
			double zscore = (freq-meanStd[0])/meanStd[1];
			if (Double.isInfinite(zscore)) zscore = zscoreForInfinity;
			if (zscore < minZScore) {
				minZScore = zscore;
				minZScoreSamples = toExamine;
			}
			if (verbose) printPileupInfo(ml, toExamine, zscore);			
		}
		
		//was a z-score calculated?
		if (minZScore == Double.MAX_VALUE) {
			if (removeNonZScoredRecords == false) {
				vcfOut.println(record);
				numSaved++;
			}
			tooFewSamples.add(record);
			numNotScored++;
		}
		
		//append min zscore and save
		else {
			//fetch AFs
			String bkafs = fetchFormattedAFs(minZScoreSamples);
			fields[7] = "BKZ="+Num.formatNumber(minZScore, 3)+";BKAF="+bkafs+";"+fields[7];
			String modRecord = Misc.stringArrayToString(fields, "\t");
			//threshold it?
			if (minimumZScore == 0) {
				vcfOut.println(modRecord);
				numSaved++;
			}
			else {
				if (minZScore < minimumZScore) numFailingZscore++;
				else {
					vcfOut.println(modRecord);
					numSaved++;	
				}
			}
		}
	}
	
	private String fetchFormattedAFs(MpileupSample[] toExamine) {
		StringBuilder sb = new StringBuilder(Num.formatNumber(toExamine[0].getAlleleFreqNonRefPlusIndels(), 4));
		for (int i=1; i< toExamine.length; i++){
			sb.append(",");
			sb.append(Num.formatNumber(toExamine[i].getAlleleFreqNonRefPlusIndels(), 4));
		}
		return sb.toString();
	}

	private void printPileupInfo(MpileupLine ml, MpileupSample[] toExamine, double zscore) {
		StringBuilder sb = new StringBuilder();
		sb.append("\tMpileup\t"+ml.getChr()+"\t"+(ml.getZeroPos()+1) +"\t"+ml.getRef()+"\t"+zscore+"\n");
		//generate sample info
		for (MpileupSample ms: toExamine){
				sb.append("\t\tSample\t");
				ms.appendCounts(sb, true);
				sb.append("\t");
				sb.append(ms.getAlleleFreqNonRefPlusIndels());
				sb.append("\n");
		}
		System.out.println(sb.toString());
	}

	private double[] calcMeanStdev(MpileupSample[] samples) {
		double[] d = new double[samples.length];
		double mean = 0;
		for (int i=0; i< samples.length; i++) {
			d[i] = samples[i].getAlleleFreqNonRefPlusIndels();
			mean += d[i];
		}
		mean = mean/(double)samples.length;
		
		double stdev = Num.standardDeviation(d, mean);
		return new double[]{mean, stdev};
	}

	private MpileupSample[] fetchMpileupSamples (MpileupLine ml) throws Exception{

		//fetch samples to Examine
		ArrayList<MpileupSample> al = new ArrayList<MpileupSample>();
		MpileupSample[] allSamples = ml.getSamples();
		for (int i=0; i< allSamples.length; i++){
			//check depth
			if (allSamples[i].getReadCoverageAll() < minReadCoverage) continue;
			//check af for each base and ins del, if any exceed then skip, likely germline
			if (allSamples[i].findMaxSnvAF() > maxSampleAF || allSamples[i].getAlleleFreqINDEL() > maxSampleAF) continue;
			//select samples?
			if (sampleIndexesToExamine != null){
				if (sampleIndexesToExamine.contains(i)) al.add(allSamples[i]);
			}
			else al.add(allSamples[i]);
		}
		
		MpileupSample[] toExamine = new MpileupSample[al.size()];
		al.toArray(toExamine);
		return toExamine;
	}
	
	private TabixReader.Iterator fetchInteratorOnCoordinates(String coordinates) {
		TabixReader.Iterator it = null;
		//watch out for no retrieved data error from tabix
		try {
			it = reader.query(coordinates);
		} catch (ArrayIndexOutOfBoundsException e){}
		return it;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFBackgroundScanner(args);
	}		
	


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		int[] sampleIndexes = null;
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
					case 'i': sampleIndexes = Num.parseInts(args[++i], Misc.COMMA); break;
					case 'z': minimumZScore = Double.parseDouble(args[++i]); break;
					case 'b': bpBuffer = Integer.parseInt(args[++i]); break;
					case 'q': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'c': minReadCoverage = Integer.parseInt(args[++i]); break;
					case 'a': minNumSamples = Integer.parseInt(args[++i]); break;
					case 'e': removeNonZScoredRecords = true; break;
					case 'd': verbose = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (sampleIndexes!= null){
			sampleIndexesToExamine = new HashSet<Integer>(sampleIndexes.length);
			for (Integer i: sampleIndexes) sampleIndexesToExamine.add(i);
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
		
		//check for R and required libraries
		//if (fullPathToR == null || fullPathToR.canExecute()== false) {
			//Misc.printErrAndExit("\nError: Cannot find or execute your R application -> "+fullPathToR+"\n");
		//}
		/*
		String errors = IO.runRCommandLookForError("library(DESeq2); library(gplots)", fullPathToR, saveDir);
		if (errors == null || errors.length() !=0){
			Misc.printErrAndExit("\nError: Cannot find the required R libraries (DESeq2 and gplots).  Did you install DESeq2 " +
					"(http://www-huber.embl.de/users/anders/DESeq2/)?  See the author's websites for installation instructions. Once installed, " +
					"launch an R terminal and type 'library(DESeq2)' to see if it is present. R error message:\n\t\t"+errors+"\n\n");
		}*/
		
		printSettings();
		fourDecimalMax.setMaximumFractionDigits(4);
	}	
	
	public void printSettings(){
		System.out.println("Settings:\nBackground\t"+mpileup);
		System.out.println("Save dir\t"+saveDir);
		System.out.println(bpBuffer+"\tBP buffer");
		System.out.println(minReadCoverage+"\tMin mpileup sample read coverage");
		System.out.println(minBaseQuality+"\tMin mpileup sample base quality");
		System.out.println(maxSampleAF+"\tMax mpileup sample AF");
		System.out.println(minNumSamples+"\tMin # samples for z-score calc");
		System.out.println(removeNonZScoredRecords+ "\tExclude vcf records that could not be z-scored");
		System.out.println(verbose+"\tVerbose");
		System.out.println(minimumZScore+"\tMin vcf AF z-score to save");
		if (sampleIndexesToExamine!=null) System.out.println(Misc.hashSetToString(sampleIndexesToExamine, ",")+"\tMpileup samples to examine");
		else System.out.println("All\tMpileup samples to examine");
	}

	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           VCF Background Scanner : Jan 2016                      **\n" +
				"**************************************************************************************\n" +
				"VBS calculates non-reference allele frequencies (AF) from a background multi-sample \n"+
				"mpileup file over each vcf record. It then calculates a z-score for the vcf AF and \n"+
				"appends it to the INFO field. If multiple bps are affected (e.g. INDELs) or bp padding\n"+
				"provided, the lowest bp z-score is appended. Z-scores < ~4 are indicative of non\n"+
				"reference bps in the background samples. Note, VBS requires an AF tag in the INFO\n"+
				"field of each record. \n"+

				"\nRequired:\n"+
				"-v Path to a xxx.vcf(.gz/.zip OK) file or directory containing such.\n" +
				"-m Path to a bgzip compressed and tabix indexed multi-sample mpileup file. e.g.:\n"+
				"      1) Mpileup: 'echo \"#SampleOrder: \"$(ls *bam) > bkg.mpileup; samtools mpileup\n"+
				"             -B -q 13 -d 1000000 -f $fastaIndex -l $bedFile *.bam >> bkg.mpileup'\n"+
				"      2) Bgzip: 'tabix-0.2.6/bgzip bkg.mpileup'\n"+
                "         Tabix: 'tabix-0.2.6/tabix -s 1 -b 2 -e 2 bkg.mpileup.gz'\n"+
				"-s Path to directory in which to save the modified vcf file(s)\n"+
						
				"\nOptional:\n" +
				"-i Comma delimited list (zero is 1st sample, no spaces) of mpileup sample indexes to\n"+
				"     examine, defaults to all.\n"+
				"-b BP padding, defaults to 0\n"+
				"-z Minimum vcf z-score, defaults to 0, no filtering\n"+
				"-q Minimum mpileup sample bp quality, defaults to 20\n"+
				"-c Minimum mpileup sample read coverge, defaults to 20\n"+
				"-f Maximum mpileup sample AF, defaults to 0.3\n"+
				"-a Minimum # mpileup samples for z-score calculation, defaults to 3\n" +
				"-e Exclude vcf records that could not be z-scored\n"+
				"-d Print verbose debugging output\n" +

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/VCFBackgroundScanner -v SomaticVcfs/ -z 3\n"+
				"-m bkg.mpileup.gz -s BkgFiltVcfs/ -b 1 -q 13 -e \n\n"+

		        "**************************************************************************************\n");
	}

}
