
package edu.utah.seq.parsers.mpileup;
import java.awt.Window;
import java.io.*;
import java.util.ArrayList;
import java.util.regex.*;

import javax.naming.Reference;

import util.gen.*;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class EstimateErrorRates{
	//fields
	private File mpileupFile;
	private int minReadCoverage = 100;
	private double maxIndelAlleleFreq = 0.1;
	private double maxNonRefAlleleFreq = 0.1;
	private double maxPoorQualAlleleFreq = 0.5;
	private int flank = 3;
	private int minBaseQuality = 20;

	//internal
	private int windowSize =0;
	private long parsedLines = 0;
	private long scoredSamples = 0;
	private Histogram histGATC = null;
	
	//for calculating error rates on passing bases
	private long gBases = 0;
	private long aBases = 0;
	private long tBases = 0;
	private long cBases = 0;
	private long totalBases = 0;

	//constructors
	public EstimateErrorRates(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		printSettings();
		
		System.out.println("\nParsing...");
		parse();
		
		System.out.println("\n"+parsedLines+"\tParsed mpileup records");
		System.out.println(scoredSamples+"\tScored mpileup samples");
		
		printErrorRates();
		
		System.out.println("\nHistogram of per base nonRef GATC AlleleFreq");
		histGATC.printScaledHistogram();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");


	}



	private void printErrorRates() {
		System.out.println("\nError rates (Base, #NonRef, #Total, %Error)");
		System.out.println(fetchErrorRateLine("G", gBases, totalBases));
		System.out.println(fetchErrorRateLine("A", aBases, totalBases));
		System.out.println(fetchErrorRateLine("T", tBases, totalBases));
		System.out.println(fetchErrorRateLine("C", cBases, totalBases));
		System.out.println(fetchErrorRateLine("GATC", (gBases+aBases+tBases+cBases), totalBases));
		
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
				MpileupLine ml = new MpileupLine(line, minBaseQuality);
				if (ml.getChr() == null) continue;
				
				//merge sample counts, not working!
				MpileupSample sample = MpileupSample.mergeSampleCounts(ml.getSamples());

				//check quality
				checkQuality(sample);
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
		}

	}

	private void score(MpileupSample s) {

		double nrAF = s.getAlleleFreqNonRef();
		//if (nrAF > 0.05) {
			//s.printMinInfo();
			//s.debug();
		//}
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
		//indels?
		else if (sample.getAlleleFreqINDEL() > maxIndelAlleleFreq) sample.setPass(false);
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
	}
	
	public void printSettings(){
		System.out.println("Parsing "+mpileupFile);
		System.out.println(minBaseQuality+"\tMin Base Quality");
		System.out.println(minReadCoverage+"\tMin Read Coverage");
		System.out.println(flank+"\tBP buffer (window size = "+windowSize+")");
		System.out.println(maxIndelAlleleFreq+"\tMax Indel Allele Freq");
		System.out.println(maxNonRefAlleleFreq+"\tMax Non Ref Allele Freq");
		System.out.println(maxPoorQualAlleleFreq+"\tMax Poor Qual Allele Freq");
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
				"**                            Estimate Error Rates: Dec 2016                        **\n" +
				"**************************************************************************************\n" +
				"EER scans an mpileup file looking for short windows of adjacent bps (default 7) where\n"+
				"1) each base exceeds a minimum read depth of high quality bases (>100)\n"+
				"2) shows little evidence of INDELs (<0.1), and \n"+
				"3) the fraction of poor quality bps isn't excessive (< 0.5). \n"+
				"The non reference base observations are then tabulated for the center base in each\n"+
				"window, if low (<0.1), they are assumed error and saved. Run this app on samples where \n"+
				"real snvs are expected to be > 0.5 allele freq, e.g. normal or pure single clone\n"+
				"somatic.\n"+

				"\nRequired Options:\n"+
				"-m Path to a normal sample mpileup file (gz/zip OK), 'samtools mpileup -B -q 20 -d \n"+
				"     1000000 -f $fastaIndex -l $bedFile *.bam | gzip > mpileup.gz' Multiple samples\n" +
				"     in the file are merged.\n"+

				"\nDefault Options:\n"+
				"-b Minimum base quality, default 20\n"+
				"-r Minimum good base coverage, default 100\n"+
				"-i Maximum INDEL allele freq, default 0.1\n"+
				"-n Maximum non reference allelic freq, default 0.1\n"+
				"-p Maximum failing base allele freq, default 0.5\n"+
				"-f Number flanking bp to define scorable region, default 3\n"+

				"\nExample: java -Xmx4G -jar pathToUSeq/Apps/EstimateErrorRates -m normExo.mpileup.gz\n" +
				"     -r 200 -i 0.15 -f 2  \n" +

				"**************************************************************************************\n");

	}
}
