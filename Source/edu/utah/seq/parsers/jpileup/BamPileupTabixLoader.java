package edu.utah.seq.parsers.jpileup;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.query.QueryIndexFileLoader;
import edu.utah.seq.vcf.VCFBkz;
import htsjdk.tribble.readers.TabixReader;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Pulls bpileup lines that overlap each vcf record and calculates AF and z-score stats.*/
public class BamPileupTabixLoader implements Runnable{

	//fields
	private boolean failed = false;
	private VCFBkz vbp;
	private TabixReader tabixReader = null;
	private ArrayList<String> vcfLines = new ArrayList<String>();
	private ArrayList<String> modVcfRecords = new ArrayList<String>();
	private ArrayList<String> tooFewSamples = new ArrayList<String>();
	private int minReadCoverage;
	private int minNumSamples;
	private boolean removeNonZScoredRecords;
	private double minimumZScore;
	private double maxSampleAF;
	private boolean replaceQualScore;
	private int numNotScored = 0;
	private int numFailingZscore = 0;
	private int numWithBKAF = 0;
	private HashSet<Integer> sampleIndexesToExclude = null;
	private boolean excludeBKAFs = false;
	private boolean calcAllBkgAfs = false;
	private boolean includeNsPoorQualBpsInAllBkgAfs = false;
	private boolean debug = false;
	private int bpIndelPad = 0;
	private boolean includeIndelCountsForSnvs = false;
	
	//internal
	public Pattern AF = null;
	private static final NumberFormat fourDecimalMax = NumberFormat.getNumberInstance();

	public BamPileupTabixLoader (File bpileupFile, VCFBkz vbp) throws IOException{
		this.vbp = vbp;
		minReadCoverage = vbp.getMinReadCoverage();
		minNumSamples = vbp.getMinNumSamples();
		removeNonZScoredRecords = vbp.isRemoveNonZScoredRecords();
		minimumZScore = vbp.getMinimumZScore();
		maxSampleAF = vbp.getMaxSampleAF();
		replaceQualScore = vbp.isReplaceQualScore();
		tabixReader = new TabixReader(bpileupFile.getCanonicalPath());
		fourDecimalMax.setMaximumFractionDigits(4);
		sampleIndexesToExclude = vbp.getSampleIndexesToExclude();
		excludeBKAFs = vbp.isExcludeBKAFs();
		calcAllBkgAfs = vbp.isCalcAllBkgAfs();
		AF = Pattern.compile( vbp.getAFInfoName()+"=([\\d\\.]+)");
		bpIndelPad = vbp.getBpIndelPad();
		includeIndelCountsForSnvs = vbp.isIncludeIndelCountsForSnvs();
	}
	
	public void run() {	
		try {
			//get next chunk of work
			while (vbp.loadVcfRecords(vcfLines)){ 
				
				//for each
				for (String record: vcfLines) score(record);
				vbp.save(modVcfRecords);
				
				//cleanup
				vcfLines.clear();
				modVcfRecords.clear();
			}
			//update
			vbp.update(numNotScored, numFailingZscore, tooFewSamples, numWithBKAF);
			//reset
			numNotScored = 0;
			numFailingZscore = 0;
			numWithBKAF = 0;
			tooFewSamples.clear();
		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem fetching bpileup lines\n" );
			e.printStackTrace();
		}
	}
	
	private void printFailingRecord(String[] fields, String record) throws IOException{
		if (removeNonZScoredRecords == false) printScoredVcf(fields, record);
		tooFewSamples.add(record);
		numNotScored++;
	}
	
	private void score(String record) throws Exception {
if (debug) IO.el("\nProc "+record);
		//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
		String[] fields = Misc.TAB.split(record);		

		//fetch AF
		Matcher mat = AF.matcher(fields[7]);
		if (mat.find() == false) {
			System.err.println ("WARNING: Failed to parse AF= number from the INFO field in this variant -> "+record+" skipping.");
			printFailingRecord(fields, record);
			return;
		}
		double altAF = Double.parseDouble(mat.group(1));

		//what kind of variant, returns GATCID
		char allele = fetchAllele(fields);

		//calc z-score for Alt AF
		ScoreAfs scoreAfs = null;
		if (allele == 'I' || allele == 'D') {
			//indel scan
			scoreAfs = scoreIndel(fields, altAF, allele);
		}
		else {
			//snv score, just the exact base
			scoreAfs = scoreSnv(fields, altAF, allele);
		}

		//failed to calculate?
		if (scoreAfs == null) {
			printFailingRecord(fields, record);
			return;
		}

		//append min zscore and save

		//fetch sorted, smallest to largest AFs
		Arrays.sort(scoreAfs.afs);

		//were two or more background AF * 1.111111112 >= tum AF seen?
		int len = scoreAfs.afs.length-1;
		double testAf = 0.9*altAF;
		boolean bkafFound = scoreAfs.afs[len]>=testAf && scoreAfs.afs[len-1]>=testAf;

		if (bkafFound && excludeBKAFs) numWithBKAF++;
		else {
			if (bkafFound) fields[6] = modifyFilterField(fields[6]);

			//fetch AFs
			String bkafString = fetchFormattedAFs(scoreAfs.afs);
			String bkzString = Num.formatNumberNoComma(scoreAfs.zscore, 2);
			fields[7] = "BKZ="+bkzString+";BKAF="+bkafString+";"+fields[7];

			if (replaceQualScore) fields[5] = bkzString;
			String modRecord = Misc.stringArrayToString(fields, "\t");

			//threshold it?
			if (minimumZScore == 0) {
				modVcfRecords.add(modRecord);
			}
			else {
				if (scoreAfs.zscore < minimumZScore) numFailingZscore++;
				else {
					modVcfRecords.add(modRecord);
				}
			}
		}
	}
	
	private ScoreAfs scoreIndel(String[] fields, double altAF, char allele) throws Exception {
		
		if (debug) IO.pl("\tScoring INDEL "+altAF+" "+allele);
		String tabixCoor = null;
		//fetch interbase coordinates for del
		if (allele == 'D') {
			int[] startStop = QueryIndexFileLoader.fetchEffectedBps(fields, true);
			if (startStop == null) return null; 
			if (bpIndelPad !=0) {
				startStop[0] = startStop[0]-bpIndelPad;
				if (startStop[0] <0) startStop[0] = 0;
				startStop[1] = startStop[1]+ bpIndelPad;
			}
			if (debug) IO.pl("\tStartEnd "+startStop[0]+ " to "+startStop[1]);
			//pull bpileup records over deleted bps
			tabixCoor = fields[0]+":"+(startStop[0]+2)+"-"+ (startStop[1]);
		}
		else if (allele ==  'I') {
			//single downstream base
			int start = Integer.parseInt(fields[1])+1;
			int stop = start;
			if (bpIndelPad !=0) {
				start = start-bpIndelPad;
				if (start <0) start = 0;
				stop = stop+bpIndelPad;
			}
			tabixCoor = fields[0]+":"+start+"-"+stop;
		}
		
		
		if (debug) IO.pl("\tTabCoor "+tabixCoor);
		TabixReader.Iterator it = fetchInteratorOnCoordinates(tabixCoor);
		if (it == null) {
			if (debug) IO.pl("No bpileup lines for "+Misc.stringArrayToString(fields, "\t"));		
			return null;
		}

		//for each bpileup record, find lowest zscore
		String bpileupLine = null;
		ScoreAfs scoreAfs = null;

		while ((bpileupLine = it.next()) != null){
			if (debug) IO.pl("\nLine "+bpileupLine);

			//parse bpileup line and pull filtered set of samples
			BpileupLine ml = new BpileupLine(bpileupLine);
			ArrayList<Double> altAFs  = calcAltIndelAFs(allele, ml);
			if (altAFs.size() < minNumSamples) {
				if (debug) IO.pl("Too few passing samples to calc bkz "+bpileupLine);
				continue;
			}

			if (debug) IO.pl("\tNumSamp "+ml.getSamples().length+" "+altAFs.size()+" afs "+altAFs);

			//calculate call zscore
			double[] afs = Num.arrayListOfDoubleToArray(altAFs);
			double[] meanStd = calcMeanStdev(afs);			
			double zscore = (altAF-meanStd[0])/meanStd[1];
			if (Double.isInfinite(zscore)) zscore = VCFBkz.zscoreForInfinity;
			if (zscore < VCFBkz.zscoreLessThanZero) zscore = VCFBkz.zscoreLessThanZero;
			if (debug) IO.pl("\tzScore "+zscore);

			//lowest?
			if (scoreAfs == null) scoreAfs = new ScoreAfs(zscore, afs);
			else if (zscore < scoreAfs.zscore) scoreAfs = new ScoreAfs(zscore, afs); 

		}
		if (debug) IO.pl("Doneeeee "+scoreAfs.zscore);
		
		return scoreAfs;
	}

	private ScoreAfs scoreSnv(String[] fields, double testAF, char allele) throws Exception {
		//pull bpileup record, if none return null;
		String tabixCoor = fields[0]+":"+fields[1]+"-"+fields[1];
		TabixReader.Iterator it = fetchInteratorOnCoordinates(tabixCoor);
		if (it == null) {
			return null;
		}

		String bpileupLine = it.next();
		if (bpileupLine == null) {
			return null;
		}

		//parse bpileup line calculate alt AF for samples not germline (het or homo) for the alt
		BpileupLine ml = new BpileupLine(bpileupLine);
		ArrayList<Double> altAFs = calcAltSnvAFs(allele, ml);
		if (altAFs.size() < minNumSamples) {
			return null;
		}		

		//calculate call zscore
		double[] afs = Num.arrayListOfDoubleToArray(altAFs);
		double[] meanStd = calcMeanStdev(afs);			
		double zscore = (testAF-meanStd[0])/meanStd[1];
		if (Double.isInfinite(zscore)) zscore = VCFBkz.zscoreForInfinity;
		if (zscore < VCFBkz.zscoreLessThanZero) zscore = VCFBkz.zscoreLessThanZero;
		return new ScoreAfs(zscore, afs);
	}
		
	private class ScoreAfs{
			double zscore;
			double[] afs;
			public ScoreAfs(double zscore, double[] afs) {
				this.zscore = zscore;
				this.afs = afs;
			}
		}

	/**Return GATC or ID for indels*/
	public static char fetchAllele(String[] fields) throws IOException {
		if (fields[4].contains(",") || fields[4].startsWith("<")) throw new IOException("Cannot interpret this alt, deconvolute? "+Misc.stringArrayToString(fields,  "\t"));
		//#CHROM	POS	ID	REF	ALT
		int lenRef = fields[3].length();
		int lenAlt = fields[4].length();
		//snv?
		if (lenRef == 1 && lenAlt == 1) return fields[4].charAt(0);
		//indel
		if (lenRef > lenAlt) return 'D';
		return 'I';
	}

	/**Adds a BKAF fail to the FILTER field.*/
	private static String modifyFilterField(String filterField) {
		//nada or pass?
		if (filterField.length() == 0 || filterField.equals(".") || filterField.toUpperCase().equals("PASS")) return "BKAF";
		else return "BKAF;" +filterField;
	}
	
	private void printScoredVcf(String[] fields, String record) throws IOException{
		if (replaceQualScore) {
			fields[5] = "0";
			modVcfRecords.add(Misc.stringArrayToString(fields, "\t"));
		}
		else modVcfRecords.add(record);
	}
	
	
	public static String fetchFormattedAFs(double[] sortedAFs) {
		
		StringBuilder sb = new StringBuilder(Num.formatNumber(sortedAFs[sortedAFs.length-1], 4));
		for (int i=sortedAFs.length-2; i>-1; i--){
			sb.append(",");
			sb.append(Num.formatNumberJustMax(sortedAFs[i], 4));
		}
		return sb.toString();
	}
	
	private static double[] calcMeanStdev(double[] afs) {
		double mean = 0;
		for (int i=0; i< afs.length; i++) mean += afs[i];
		mean = mean/(double)afs.length;
		double stdev = Num.standardDeviation(afs, mean);
		return new double[]{mean, stdev};
	}
	
	private TabixReader.Iterator fetchInteratorOnCoordinates(String coordinates) {
		TabixReader.Iterator it = null;
		//watch out for no retrieved data error from tabix
		try {
			it = tabixReader.query(coordinates);
		} catch (ArrayIndexOutOfBoundsException e){
		}
		return it;
	}
	
	
	private ArrayList<Double> calcAltSnvAFs (char alt, BpileupLine ml) throws Exception{
		//fetch samples 
		ArrayList<Double> al = new ArrayList<Double>();
		BaseCount[] samples = ml.getSamples();
		for (int i=0; i< samples.length; i++){
			//skip sample?
			if (sampleIndexesToExclude != null && sampleIndexesToExclude.contains(i)) continue;	
			
			//check passing read depth, g+a+t+c+ins+del
			if (samples[i].getPassingReadCoverage() < minReadCoverage) continue;
			
			//calc the sample AF, either the all bkg or just for the particular allele
			double af = 0;
			if (calcAllBkgAfs) af = calcAllBkgAf(samples[i]);
			else {
				//check that alt isn't a het in the normal sample
				double snvCounts = samples[i].getSnvCount(alt);
				af = snvCounts / samples[i].getPassingReadCoverageSnv();
				if (af > maxSampleAF) continue;	
				//include indel counts?
				if (includeIndelCountsForSnvs) af = ((double)samples[i].getIndelCount()+snvCounts) / samples[i].getPassingReadCoverage();
			}
			//watch out for infinity or not a number
			if (Double.isInfinite(af) || Double.isNaN(af)) continue;
			al.add(af);
		}
		return al;
	}
	
	private double calcAllBkgAf(BaseCount sample) {
		//test each base obs AF and see if it's less than the maxSampleAF (and potentially germline), if so add it to the sum.
		double sum =0;
		double count = sample.getPassingReadCoverage();
		
		//A
		double af = (double) sample.a / count;
		if (af <= maxSampleAF) sum+= sample.a; 
		//C
		af = (double) sample.c / count;
		if (af <= maxSampleAF) sum+= sample.c;
		//G
		af = (double) sample.g / count;
		if (af <= maxSampleAF) sum+= sample.g;
		//T
		af = (double) sample.t / count;
		if (af <= maxSampleAF) sum+= sample.t;
		//Del
		af = (double) sample.del / count;
		if (af <= maxSampleAF) sum+= sample.del;
		//Ins
		af = (double) sample.ins / count;
		if (af <= maxSampleAF) sum+= sample.ins;
		//include N's and low quality bps? - prob shouldn't do this.
		if (includeNsPoorQualBpsInAllBkgAfs) {
			count = sample.getTotalReadCoverage();
			sum += sample.n;
			sum += sample.failQual;
		}
		return sum/count;
	}
	
	private ArrayList<Double> calcAltIndelAFs (char alt, BpileupLine ml) throws Exception{
		//fetch samples 
		ArrayList<Double> al = new ArrayList<Double>();
		BaseCount[] allSamples = ml.getSamples();
		for (int i=0; i< allSamples.length; i++){
			//skip sample?
			if (sampleIndexesToExclude != null && sampleIndexesToExclude.contains(i)) continue;	

			//check passing read depth, g+a+t+c+ins+del
			if (allSamples[i].getPassingReadCoverage() < minReadCoverage) continue;

			//calc the sample AF, either the all bkg or just for the indel
			double af = 0;
			if (calcAllBkgAfs) af = calcAllBkgAf(allSamples[i]);
			else {
				//check that alt isn't a het in the normal sample
				if (alt == 'D') af = (double)allSamples[i].del / allSamples[i].getPassingReadCoverage();
				else af = (double)allSamples[i].ins / allSamples[i].getPassingReadCoverage();
				if (af >= maxSampleAF) continue;
				//recalc using all indels, real indels have either ins or del, fp indels often have both
				af = (double)allSamples[i].getIndelCount() / allSamples[i].getPassingReadCoverage();
			}
			//watch out for infinity or not a number
			if (Double.isInfinite(af) || Double.isNaN(af)) continue;
			al.add(af);
		}
		return al;
	}

	public boolean isFailed() {
		return failed;
	}

	public TabixReader getTabixReader() {
		return tabixReader;
	}
}
