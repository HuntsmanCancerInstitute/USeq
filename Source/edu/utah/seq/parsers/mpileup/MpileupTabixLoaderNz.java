package edu.utah.seq.parsers.mpileup;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.query.QueryIndexFileLoader;
import edu.utah.seq.vcf.VCFBkzDepreciated;
import edu.utah.seq.vcf.VCFNz;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.tribble.readers.TabixReader;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Pulls mpileup lines that overlap each vcf record and calculates AF and z-score stats.*/
public class MpileupTabixLoaderNz implements Runnable{

	//fields
	private boolean failed = false;
	private VCFNz vnz;
	private TabixReader tabixReader = null;
	private ArrayList<String> vcfLines = new ArrayList<String>();
	private ArrayList<String> modVcfRecords = new ArrayList<String>();
	private ArrayList<String> tooFewSamples = new ArrayList<String>();
	private ArrayList<int[][]> baseCounts = new ArrayList<int[][]>();
	private int minBaseQuality;
	private boolean removeNonZScoredRecords;
	private double maximumZScore;
	private HashMap<String, double[]> seqMeanStd = null;
	private int numNotScored = 0;
	private int numFailingZscore = 0;
	private boolean verbose;
	private int bpPad = 0;
	private IndexedFastaSequenceFile fasta = null;
	private String workingChromosome = "";
	private String workingSequence = null;
	
	//internal
	private static final NumberFormat fourDecimalMax = NumberFormat.getNumberInstance();

	public MpileupTabixLoaderNz (File mpileupFile, VCFNz vnz) throws IOException{
		this.vnz = vnz;
		minBaseQuality = vnz.getMinBaseQuality();
		removeNonZScoredRecords = vnz.isRemoveNonZScoredRecords();
		maximumZScore = vnz.getMaximumZScore();
		tabixReader = new TabixReader(mpileupFile.getCanonicalPath());
		fourDecimalMax.setMaximumFractionDigits(4);
		bpPad = vnz.getBpPad();
		seqMeanStd = vnz.getSeqMeanStd();
		fasta = new IndexedFastaSequenceFile(vnz.getIndexedFasta());
		if (fasta.isIndexed() == false) Misc.printErrAndExit("\nError: cannot find your xxx.fai index or the multi fasta file isn't indexed\n");

	}
	
	public void run() {	
		try {
			//get next chunk of work
			while (vnz.loadVcfRecords(vcfLines)){ 
				
				//for each
				for (String record: vcfLines) score(record);
				vnz.save(modVcfRecords, baseCounts);
				
				//cleanup
				vcfLines.clear();
				modVcfRecords.clear();
				baseCounts.clear();
			}
			//update
			vnz.update(numNotScored, numFailingZscore);
			//reset
			numNotScored = 0;
			numFailingZscore = 0;
			fasta.close();
		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem parsing mpileup lines\n" );
			e.printStackTrace();
		}
	}
	
	private void printFailingRecord(String[] fields, String record) throws IOException{
		if (removeNonZScoredRecords == false) printScoredVcf(fields, record);
		tooFewSamples.add(record);
		numNotScored++;
	}
	
	private void score(String record) throws Exception {
		//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
		String[] fields = Misc.TAB.split(record);		
		
		checkChromosome(fields[0]);
		
		//interbase coor of effected bps
		int[] startStop = QueryIndexFileLoader.fetchEffectedBps(fields, true);
		if (startStop == null) {
			System.err.println("WARNING: Failed to parse the effected bps for : "+record);
			printFailingRecord(fields, record);
			return;
		}
		
		
		//pull mpileup records, if none then warn and save vcf record
		int start = startStop[0] + 1 - bpPad;
		if (start < 1) start = 1;
		String tabixCoor = fields[0]+":"+start+"-"+(startStop[1]+bpPad);
		TabixReader.Iterator it = fetchInteratorOnCoordinates(tabixCoor);
		if (it == null) {
			printFailingRecord(fields, record);
			return;
		}
		
		StringBuilder sb = null;
		if (verbose) {
			sb = new StringBuilder ("VCF\t");
			sb.append(record);
			sb.append("\n");
		}
		
		//for each mpileup record, find highest z-score
		String mpileupLine = null;
		double maxZScore = -1;

		while ((mpileupLine = it.next()) != null){

			//parse mpilup line and pull first sample
			MpileupLine ml = new MpileupLine(mpileupLine, minBaseQuality);
			if (ml.getChr() == null) throw new IOException ("Failed to parse the mpileup line:\n"+mpileupLine);
			
			MpileupSample toExamine = ml.getSamples()[0];
			double afN = toExamine.getAlleleFreqNs();
			
			//fetch five mer
			int pos = ml.getZeroPos();
			String seq = workingSequence.substring(pos-2, pos+3).toUpperCase();
IO.pl("Seq:"+seq+" Ref:"+ml.getRef()+" afN:"+afN);
			
			//fetch meanStd and calc zscore
			double[] meanStd = seqMeanStd.get(seq);			
			double zscore = (afN-meanStd[0])/meanStd[1];
			if (Double.isInfinite(zscore)) zscore = VCFNz.zscoreForInfinity;
			if (zscore < VCFNz.zscoreLessThanZero) zscore = VCFNz.zscoreLessThanZero;
			
			if (zscore > maxZScore) maxZScore = zscore;

		}
		
		//was a z-score calculated?
		if (minZScore == Double.MAX_VALUE) printFailingRecord(fields, record);
		
		//append min zscore and save
		else {
			//fetch sorted, smallest to largest AFs
			double[] bkgAFs = fetchSortedAFs(minZScoreSamples);
			
			//was a background AF >= tum AF seen?
			boolean bkafFound = (freq <= bkgAFs[bkgAFs.length-1]);
			
			if (bkafFound && excludeBKAFs) {
				numWithBKAF++;
			}
			else {
				if (bkafFound) fields[6] = modifyFilterField(fields[6]);
				
				//fetch AFs
				String bkafString = fetchFormattedAFs(bkgAFs);
				String bkzString = Num.formatNumberNoComma(minZScore, 2);
				fields[7] = "BKZ="+bkzString+";BKAF="+bkafString+";"+fields[7];
				
				if (replaceQualScore) fields[5] = bkzString;
				String modRecord = Misc.stringArrayToString(fields, "\t");
				
				//build counts for R, first is the tumor sample, the latter are the normals
				int[][] counts = fetchCountsForR(minZScoreSamples, numNonRef, depth);
				
				//threshold it?
				if (minimumZScore == 0) {
					modVcfRecords.add(modRecord);
					baseCounts.add(counts);
				}
				else {
					if (minZScore < minimumZScore) numFailingZscore++;
					else {
						modVcfRecords.add(modRecord);
						baseCounts.add(counts);
					}
				}
			}
			
			
			
		}	
	}
	
	private void checkChromosome(String chr) throws IOException {
		if (chr.equals(workingChromosome) == false) {
			workingChromosome = chr;
			//find and load sequence
			ReferenceSequence p = fasta.getSequence(workingChromosome);
			if (p == null ) throw new IOException ("\n\nFailed to find or load a fasta sequence for '"+workingChromosome+"', aborting.\n");
			workingSequence = new String(p.getBases());
		}
	}

	/**Adds a BKAF fail to the FILTER field.*/
	private static String modifyFilterField(String filterField) {
		//nada or pass?
		if (filterField.length() == 0 || filterField.equals(".") || filterField.toUpperCase().equals("PASS")) return "BKAF";
		else return "BKAF;" +filterField;
	}

	public static boolean highBkgrndAFsFound(double tAF, MpileupSample[] toExamine) {
		for (int i=0; i< toExamine.length; i++) if (toExamine[i].getAlleleFreqNonRefPlusIndels() >= tAF) return true;
		return false;
	}
	
	private int[][] fetchCountsForR(MpileupSample[] minZScoreSamples, int numNonRef, double depth ){
		int[] nonRef = new int[minZScoreSamples.length+1];
		int[] total = new int[nonRef.length];
		nonRef[0] = numNonRef;
		total[0] = ((int)depth) - numNonRef;
		int index = 1;
		for (MpileupSample ms: minZScoreSamples){
			nonRef[index] = ms.getNonRefBaseCounts() + ms.getInsertions()+ ms.getDeletions();
			total[index++] = ms.getReadCoverageAll();
		}
		return new int[][]{nonRef,total};
	}
	
	private TabixReader.Iterator fetchInteratorOnCoordinates(String coordinates) {
		TabixReader.Iterator it = null;
		//watch out for no retrieved data error from tabix
		try {
			it = tabixReader.query(coordinates);
		} catch (ArrayIndexOutOfBoundsException e){}
		return it;
	}

	public boolean isFailed() {
		return failed;
	}

	public TabixReader getTabixReader() {
		return tabixReader;
	}
}
