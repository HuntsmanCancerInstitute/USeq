package edu.utah.seq.parsers.mpileup;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.query.QueryIndexFileLoader;
import edu.utah.seq.vcf.VCFBackgroundChecker;
import edu.utah.seq.vcf.VCFBackgroundCheckerGraphite;
import htsjdk.tribble.readers.TabixReader;
import util.gen.Misc;
import util.gen.Num;

/**Pulls mpileup lines that overlap each vcf record and calculates AF and z-score stats.*/
public class MpileupTabixLoaderGraphite implements Runnable{

	//fields
	private boolean failed = false;
	private VCFBackgroundCheckerGraphite vbc;
	private TabixReader tabixReader = null;
	private ArrayList<String> vcfLines = new ArrayList<String>();
	private ArrayList<String> modVcfRecords = new ArrayList<String>();
	private ArrayList<String> tooFewSamples = new ArrayList<String>();
	private int minBaseQuality;
	private int minReadCoverage;
	private int minNumSamples;
	private boolean removeNonZScoredRecords;
	private double minimumZScore;
	private double maxSampleAF;
	private boolean replaceQualScore;
	private int numNotScored = 0;
	private int numFailingZscore = 0;
	private boolean verbose;
	
	//internal
	public static final Pattern AF = Pattern.compile(".*AF=([\\d\\.]+).*");
	public static final double zscoreForInfinity = 1000;
	private static final NumberFormat fourDecimalMax = NumberFormat.getNumberInstance();

	

	public MpileupTabixLoaderGraphite (File mpileupFile, VCFBackgroundCheckerGraphite vbc) throws IOException{
		this.vbc = vbc;
		minBaseQuality = vbc.getMinBaseQuality();
		minReadCoverage = vbc.getMinReadCoverage();
		minNumSamples = vbc.getMinNumSamples();
		removeNonZScoredRecords = vbc.isRemoveNonZScoredRecords();
		minimumZScore = vbc.getMinimumZScore();
		maxSampleAF = vbc.getMaxSampleAF();
		replaceQualScore = vbc.isReplaceQualScore();
		verbose = vbc.isVerbose();
		tabixReader = new TabixReader(mpileupFile.getCanonicalPath());
		fourDecimalMax.setMaximumFractionDigits(4);
	}
	
	public void run() {	
		try {
			//get next chunk of work
			while (vbc.loadVcfRecords(vcfLines)){ 
				
				//for each
				for (String record: vcfLines) score(record);
				vbc.saveModifiedVcf(modVcfRecords);
				
				//cleanup
				vcfLines.clear();
				modVcfRecords.clear();
				
			}
			//update
			vbc.update(numNotScored, numFailingZscore, tooFewSamples);
			//reset
			numNotScored = 0;
			numFailingZscore = 0;
			tooFewSamples.clear();
		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem fetching mpileup lines\n" );
			e.printStackTrace();
		}
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
			if (removeNonZScoredRecords == false) printScoredVcf(fields, record);
			tooFewSamples.add(record);
			numNotScored++;
			return;
		}
		
		//pull mpileup records, if none then warn and save vcf record
		int start = startStop[0]+1;
		if (start < 1) start = 1;
		String tabixCoor = fields[0]+":"+start+"-"+(startStop[1]);
		TabixReader.Iterator it = fetchInteratorOnCoordinates(tabixCoor);
		if (it == null) {
			if (removeNonZScoredRecords == false) printScoredVcf(fields, record);
			tooFewSamples.add(record);
			numNotScored++;
			return;
		}
		
		StringBuilder sb = null;
		if (verbose) {
			sb = new StringBuilder ("VCF\t");
			sb.append(record);
			sb.append("\n");
		}
		
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
			if (verbose) {
				addPileupInfo(sb, ml, toExamine, zscore);
				System.out.println(sb.toString());
			}
		}
		
		//was a z-score calculated?
		if (minZScore == Double.MAX_VALUE) {
			if (removeNonZScoredRecords == false) printScoredVcf(fields, record);
			tooFewSamples.add(record);
			numNotScored++;
		}
		
		//append min zscore and save
		else {
			//fetch AFs
			String bkafs = fetchFormattedAFs(minZScoreSamples);
			String bkzString = Num.formatNumberNoComma(minZScore, 2);
			fields[7] = "BKZ="+bkzString+";BKAF="+bkafs+";"+fields[7];
			if (replaceQualScore) fields[5] = bkzString;
			String modRecord = Misc.stringArrayToString(fields, "\t");
			//threshold it?
			if (minimumZScore == 0) modVcfRecords.add(modRecord);
			else {
				if (minZScore < minimumZScore) numFailingZscore++;
				else modVcfRecords.add(modRecord);
			}
		}
	}
	
	private void printScoredVcf(String[] fields, String record) throws IOException{
		if (replaceQualScore) {
			fields[5] = "0";
			modVcfRecords.add(Misc.stringArrayToString(fields, "\t"));
		}
		else modVcfRecords.add(record);
	}
	
	public static String fetchFormattedAFs(MpileupSample[] toExamine) {
		StringBuilder sb = new StringBuilder(Num.formatNumber(toExamine[0].getAlleleFreqNonRefPlusIndels(), 4));
		for (int i=1; i< toExamine.length; i++){
			sb.append(",");
			sb.append(Num.formatNumber(toExamine[i].getAlleleFreqNonRefPlusIndels(), 4));
		}
		return sb.toString();
	}
	
	private static double[] calcMeanStdev(MpileupSample[] samples) {
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
	
	private TabixReader.Iterator fetchInteratorOnCoordinates(String coordinates) {
		TabixReader.Iterator it = null;
		//watch out for no retrieved data error from tabix
		try {
			it = tabixReader.query(coordinates);
		} catch (ArrayIndexOutOfBoundsException e){}
		return it;
	}
	
	private void addPileupInfo(StringBuilder sb, MpileupLine ml, MpileupSample[] toExamine, double zscore) {
		sb.append("\tMpileup\t"+ml.getChr()+"\t"+(ml.getZeroPos()+1) +"\t"+ml.getRef()+"\t"+zscore+"\n");
		//generate sample info
		for (MpileupSample ms: toExamine){
				sb.append("\t\tSample\t");
				ms.appendCounts(sb, true);
				sb.append("\t");
				sb.append(ms.getAlleleFreqNonRefPlusIndels());
				sb.append("\n");
		}
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
			al.add(allSamples[i]);
		}
		
		MpileupSample[] toExamine = new MpileupSample[al.size()];
		al.toArray(toExamine);
		return toExamine;
	}

	public boolean isFailed() {
		return failed;
	}

	public TabixReader getTabixReader() {
		return tabixReader;
	}
}
