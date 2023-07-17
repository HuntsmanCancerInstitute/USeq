
package edu.utah.seq.parsers.mpileup;
import java.io.*;
import java.util.ArrayList;
import java.util.regex.*;

import util.gen.*;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class EstimateErrorRates{
	//fields
	private File mpileupFile;
	private File countFile;
	private int minReadCoverage = 100;
	private double maxIndelAlleleFreq = 0.1;
	private double maxNonRefAlleleFreq = 0.1;
	private double maxPoorQualAlleleFreq = 0.5;
	private int flank = 3;
	private int minBaseQuality = 20;
	private int[] samplesToProcess = null;

	//internal
	private int windowSize =0;
	private long parsedLines = 0;
	private long scoredSamples = 0;
	private Histogram histGATC = null;
	private PrintWriter outCounts = null;

	//for calculating error rates on passing bases
	private long gBases = 0;
	private long aBases = 0;
	private long tBases = 0;
	private long cBases = 0;
	private long totalBases = 0;
	private long totalBasesPlusIndels = 0;
	private long insertions = 0;
	private long deletions = 0;

	//constructors
	public EstimateErrorRates(String[] args) {
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);
			printSettings();

			//any count table?
			if (outCounts != null) outCounts.println("#Chr\tPos(1based)\tRef\tSample base counts for passing GATC forward, GATC reverse, ins, del, failed minQual, nonRef, nonRefFreq");

			System.out.println("\nParsing...");
			parse();

			System.out.println("\n"+parsedLines+"\tParsed mpileup records");
			System.out.println(scoredSamples+"\tScored mpileup samples");

			printErrorRates();

			System.out.println("\nHistogram of per base nonRef GATC AlleleFreq");
			histGATC.printScaledHistogram();

			if (outCounts != null) {
				//add a line that can be tabixed
				outCounts.print("observations\t1\tError Observations: G, A, T, C, totalSnv, INS, DEL, totalWithIndel, nonRef, nonRefAF\t");
				outCounts.print(gBases); outCounts.print(",");
				outCounts.print(aBases); outCounts.print(",");
				outCounts.print(tBases); outCounts.print(",");
				outCounts.print(cBases); outCounts.print(",");
				outCounts.print(totalBases); outCounts.print(",");
				outCounts.print(insertions); outCounts.print(",");
				outCounts.print(deletions); outCounts.print(",");
				outCounts.println(totalBasesPlusIndels);
				outCounts.close();
			}

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");
		} catch (Exception e){
			e.printStackTrace();
		}
	}



	private void printErrorRates() {
		System.out.println("\nError rates:\nReferenceBase\t#NonRef\t#Total\t%Error");
		System.out.println(fetchErrorRateLine("G", gBases, totalBases));
		System.out.println(fetchErrorRateLine("A", aBases, totalBases));
		System.out.println(fetchErrorRateLine("T", tBases, totalBases));
		System.out.println(fetchErrorRateLine("C", cBases, totalBases));
		System.out.println(fetchErrorRateLine("GATC", (gBases+aBases+tBases+cBases), totalBases));
		
		System.out.println("\nIndelType\t#Observed\t#Total\t%Error");
		System.out.println(fetchErrorRateLine("INS", insertions, totalBasesPlusIndels));
		System.out.println(fetchErrorRateLine("DEL", deletions, totalBasesPlusIndels));
	}

	public static String fetchErrorRateLine(String name, long baseCount, long totalBaseCount){
		double af = 100.0 * ((double)baseCount / (double) totalBaseCount);
		StringBuilder sb = new StringBuilder(name);
		sb.append("\t");
		sb.append(baseCount); 
		sb.append("\t");
		sb.append(totalBaseCount); 
		sb.append("\t");
		sb.append(af);
		return sb.toString();
	}



	/**The goal here is to identify bases in non problematic regions (good read depth, low indels, low poor qual bps, ref not N) two good bp flanks, so working xxTxx */
	private void parse() {
		try{
			BufferedReader in = IO.fetchBufferedReader(mpileupFile);
			String line;

			ArrayList<MpileupSample> workingSamples = new ArrayList<MpileupSample>();

			while ((line=in.readLine())!= null){
				if (line.startsWith("#") || line.trim().length() == 0) continue;
				MpileupLine ml = new MpileupLine(line, minBaseQuality);
				
				//watch out for missing chrom names and if the ref base is N
				if (ml.getChr() == null || ml.getRef().equals("N")) continue;
				if (outCounts != null) saveCounts(ml);

				//merge sample counts?
				MpileupSample[] toMerge = ml.getSamples();
				if (toMerge == null || toMerge.length == 0) {
					workingSamples.clear();
					continue;
				}

				if (samplesToProcess != null){
					MpileupSample[] sub = new MpileupSample[samplesToProcess.length];
					int counter = 0;
					for (int index: samplesToProcess) {
						sub[counter++] = toMerge[index];
					}
					toMerge = sub;
				}
				MpileupSample sample = toMerge[0];
				if (toMerge.length > 1) sample = MpileupSample.mergeSampleCounts(toMerge);

				//check quality
				checkQuality(sample);
				
				//score indel error?
				if (sample.isPass()){
					scoreIndel(sample);
					//now see if this is good for snv calling
					if (sample.getAlleleFreqINDEL() > maxIndelAlleleFreq) sample.setPass(false);
				}
				
				//poor so clear array, need windowSize good ones to score center base
				if (sample.isPass() == false) workingSamples.clear();
				//good so add it
				else {
					//add it to the end
					workingSamples.add(sample);
					//have five good ones?
					if (workingSamples.size() == windowSize){
						//look to see if all adjacent
						MpileupSample toScore = processWindow(workingSamples);
						if (toScore != null) score(toScore); 
						//remove from beginning
						workingSamples.remove(0);
					}
				}
				parsedLines++;
			}
			in.close();

		}catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem paring mpileup file! Aborting."+mpileupFile);
		}

	}

	private void saveCounts(MpileupLine ml) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append(ml.getChr()); sb.append("\t");
		sb.append(ml.getZeroPos()+1); sb.append("\t");
		sb.append(ml.getRef()); 

		for (MpileupSample ms: ml.getSamples()){
			sb.append("\t");
			ms.appendCounts(sb, false);
		}
		outCounts.println(sb.toString());
	}

	private void scoreIndel(MpileupSample s){
		//does it look like a snv or indel is present?
		double nrAF = s.getAlleleFreqNonRefPlusIndels();
		if (nrAF > maxNonRefAlleleFreq) return;
		//increment
		totalBasesPlusIndels += s.getReadCoverageAll();
		insertions += s.getInsertions();
		deletions += s.getDeletions();
	}

	private void score(MpileupSample s) {
		double nrAF = s.getAlleleFreqNonRef();
		//does it look like a snv?
		if (nrAF > maxNonRefAlleleFreq) return;
		//process
		scoredSamples++;
		incrementNonRefBaseCounts(s);
	}

	public void incrementNonRefBaseCounts(MpileupSample s){
		//totals
		double total = s.getReadCoverageForwardBases() + s.getReadCoverageReverseBases();
		totalBases+= total;

		//individual bases
		int[] gatcFor = s.getForwardGATC();
		int[] gatcRev = s.getReverseGATC();
		int g = 0;
		int a = 0;
		int t = 0;
		int c = 0;

		char ref = s.getRecord().getRef().charAt(0);
		if (ref == 'G') {
			//add A,T,C forward
			a+= gatcFor[MpileupSample.A_INDEX];
			a+= gatcRev[MpileupSample.A_INDEX];
			t+= gatcFor[MpileupSample.T_INDEX];
			t+= gatcRev[MpileupSample.T_INDEX];
			c+= gatcFor[MpileupSample.C_INDEX];
			c+= gatcRev[MpileupSample.C_INDEX];
		}
		else if (ref == 'A'){
			//add G,T,C
			g+= gatcFor[MpileupSample.G_INDEX];
			t+= gatcFor[MpileupSample.T_INDEX];
			c+= gatcFor[MpileupSample.C_INDEX];
			g+= gatcRev[MpileupSample.G_INDEX];
			t+= gatcRev[MpileupSample.T_INDEX];
			c+= gatcRev[MpileupSample.C_INDEX];
		}
		else if (ref == 'T'){
			//add G,A,C
			g+= gatcFor[MpileupSample.G_INDEX];
			a+= gatcFor[MpileupSample.A_INDEX];
			c+= gatcFor[MpileupSample.C_INDEX];
			g+= gatcRev[MpileupSample.G_INDEX];
			a+= gatcRev[MpileupSample.A_INDEX];
			c+= gatcRev[MpileupSample.C_INDEX];
		}
		else if (ref == 'C'){
			//add G,T,A
			g+= gatcFor[MpileupSample.G_INDEX];
			t+= gatcFor[MpileupSample.T_INDEX];
			a+= gatcFor[MpileupSample.A_INDEX];
			g+= gatcRev[MpileupSample.G_INDEX];
			t+= gatcRev[MpileupSample.T_INDEX];
			a+= gatcRev[MpileupSample.A_INDEX];
		}
		else {
			s.debug();
			Misc.printErrAndExit("ERROR: unrecognized reference base, see above.\n");
		}
		//increment total
		gBases+=g;
		aBases+=a;
		tBases+=t;
		cBases+=c;

		//af
		double af = (double)(g+a+t+c) / total;
		histGATC.count(af);
	}




	/**Looks to see if same chr and adjacent positions*/
	private MpileupSample processWindow(ArrayList<MpileupSample> workingSamples) {
		//pull chr and pos
		MpileupSample first = workingSamples.get(0);
		String chr = first.getRecord().getChr();
		int pos = first.getRecord().getZeroPos();
		//scan others
		for (int i=1; i< windowSize; i++){
			MpileupSample test = workingSamples.get(i);
			if (test.getRecord().getChr().equals(chr) == false) return null;
			pos++;
			if (test.getRecord().getZeroPos() != pos) return null;
		}
		return workingSamples.get(flank);
	}

	private void checkQuality(MpileupSample sample) {
		//low read coverage?
		if (sample.getReadCoverageAll() < minReadCoverage) sample.setPass(false);
		//poor qual bps
		else if (sample.getAlleleFreqPoorQual() > maxPoorQualAlleleFreq) sample.setPass(false);
		//ref is N?
		else if (sample.getRecord().getRef().equals("N"))  sample.setPass(false);
		//good!
		else sample.setPass(true);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'm': mpileupFile = new File(args[++i]); break;
					case 'b': minBaseQuality = Integer.parseInt(args[++i]); break;
					case 'r': minReadCoverage = Integer.parseInt(args[++i]); break;
					case 'i': maxIndelAlleleFreq = Double.parseDouble(args[++i]); break;
					case 'n': maxNonRefAlleleFreq = Double.parseDouble(args[++i]); break;
					case 'p': maxPoorQualAlleleFreq = Double.parseDouble(args[++i]); break;
					case 'f': flank = Integer.parseInt(args[++i]); break;
					case 's': samplesToProcess = Num.parseInts(args[++i], Misc.COMMA);break;
					case 'c': countFile = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (mpileupFile == null || mpileupFile.exists() == false) Misc.printErrAndExit("\nERROR: Can't find your mpileup file? "+mpileupFile);

		//set window size
		windowSize = (1+flank*2);

		//make histograms
		histGATC = new Histogram(0, maxNonRefAlleleFreq, 25);

		if (countFile != null)
			try {
				outCounts = new PrintWriter( new FileWriter(countFile));
			} catch (Exception e) {
				e.printStackTrace();
			} 
	}

	public void printSettings(){
		System.out.println("Parsing "+mpileupFile);
		System.out.println("Count table "+countFile);
		System.out.println(minBaseQuality+"\tMin Base Quality");
		System.out.println(minReadCoverage+"\tMin Read Coverage");
		System.out.println(flank+"\tBP buffer (window size = "+windowSize+")");
		System.out.println(maxIndelAlleleFreq+"\tMax Indel Allele Freq");
		System.out.println(maxNonRefAlleleFreq+"\tMax Non Ref Allele Freq");
		System.out.println(maxPoorQualAlleleFreq+"\tMax Poor Qual Allele Freq");
		if (samplesToProcess!=null) System.out.println(Num.intArrayToString(samplesToProcess, ",")+"\tSamples to process");
		else System.out.println("All\tSamples to process");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new EstimateErrorRates(args);
	}


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Estimate Error Rates: May 2023                        **\n" +
				"**************************************************************************************\n" +
				"EER scans an mpileup file looking for short windows of adjacent bps (default 7) where\n"+
				"1) each base exceeds a minimum read depth of high quality bases (>100)\n"+
				"2) shows little evidence of indels (<0.1), and \n"+
				"3) the fraction of poor quality bps isn't excessive (<0.5). \n"+
				"The non reference snv observations are then tabulated for the center base in each\n"+
				"window, if low (<0.1), they are assumed error and saved. For indel error calculations,\n"+
				"each bp is scored as above sans the indel filter and window requirement. Insertions\n"+
				"are counted once regardless of the size, where as deletions are counted for every base\n"+
				"affected. Run this app on samples where real snvs and indels are expected to have an\n"+
				"allele frequency of  > ~0.5 , e.g. normal or pure single clone somatic.\n"+

				"\nRequired Options:\n"+
				"-m Path to a normal sample mpileup file (gz/zip OK), 'samtools mpileup -B -q 20 -d \n"+
				"     1000000 -f $fastaIndex -l $bedFile *.bam | gzip > mpileup.gz' Multiple samples\n" +
				"     in the file are merged.\n"+

				"\nDefault Options:\n"+
				"-b Minimum base quality, default 20\n"+
				"-r Minimum good base coverage, default 100\n"+
				"-i Maximum INDEL allele freq for snv counting, default 0.1\n"+
				"-n Maximum non reference allelic freq, default 0.1\n"+
				"-p Maximum failing base allele freq, default 0.5\n"+
				"-f Number flanking bp to define scorable region, default 3\n"+
				"-s Comma delimited list (zero is 1st sample, no spaces) of sample indexes to merge,\n"+
				"     defaults to all.\n"+
				"-c File path to save a count table of parsed observations, defaults to none.\n"+
				//"-r Full path to R (version 3+) loaded with DESeq2, samr, and gplots defaults to\n"+
				//"       '/usr/bin/R' file, see http://www.bioconductor.org . Type 'library(DESeq2);\n"+
				//"       library(samr); library(gplots)' in R to see if they are installed. \n"+

				"\nExample: java -Xmx4G -jar pathToUSeq/Apps/EstimateErrorRates -m normExo.mpileup.gz\n" +
				"     -r 200 -i 0.15 -f 2 -s 0,3,4 -c countTable.txt\n" +

				"**************************************************************************************\n");

	}
}
