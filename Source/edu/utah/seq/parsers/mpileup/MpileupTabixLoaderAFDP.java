package edu.utah.seq.parsers.mpileup;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.query.QueryIndexFileLoader;
import edu.utah.seq.vcf.VCFMpileupAnnotator;
import htsjdk.tribble.readers.TabixReader;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Pulls mpileup lines that overlap each vcf record and calculates AF and z-score stats.*/
public class MpileupTabixLoaderAFDP implements Runnable{

	//fields
	private boolean failed = false;
	private VCFMpileupAnnotator vbc;
	private TabixReader tabixReader = null;
	private ArrayList<String> vcfLines = new ArrayList<String>();
	private ArrayList<String> modVcfRecords = new ArrayList<String>();
	private int minBaseQuality;
	private boolean verbose;
	private int numDecimals = 0;
	private ArrayList<Double> alleleFreqs = new ArrayList<Double>();
	private ArrayList<Double> readDepth = new ArrayList<Double>();

	//internal
	public static final Pattern AF = Pattern.compile("AF=[\\d\\.]+;*");
	public static final Pattern DP = Pattern.compile("DP=[\\d\\.]+;*");

	public MpileupTabixLoaderAFDP (File mpileupFile, VCFMpileupAnnotator vbc) throws IOException{
		this.vbc = vbc;
		minBaseQuality = vbc.getMinBaseQuality();
		verbose = vbc.isVerbose();
		tabixReader = new TabixReader(mpileupFile.getCanonicalPath());
		numDecimals = vbc.getNumDecimals();
	}

	public void run() {	
		try {
			//get next chunk of work
			while (vbc.loadVcfRecords(vcfLines)){ 

				//for each
				for (String record: vcfLines) scoreMultiSample(record);
				vbc.save(modVcfRecords);

				//cleanup
				vcfLines.clear();
				modVcfRecords.clear();
			}

		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem fetching mpileup lines\n" );
			e.printStackTrace();
		}
	}

	
	private void scoreMultiSample(String record) throws Exception {
		if (verbose) IO.pl("\n"+record);

		//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
		String[] fields = Misc.TAB.split(record);

		//look for , in ref or alt
		if (fields[3].contains(",") || fields[4].contains(",")) throw new IOException ("Found multi ref/alt, none allowed, normalize with vt :\n"+record);
		
		//look for . in alt
		if (fields[4].equals(".")) throw new IOException ("Found '.' for an alt, none allowed:\n"+record);

		//snv, ins, del?
		int type = 0;
		if (fields[3].length() == 1 && fields[4].length() == 1) type = 0;
		else if (fields[3].length() > fields[4].length()) type = 1;
		else if (fields[3].length() < fields[4].length())type = 2;
		else throw new IOException ("Found odd var representation? compound snv? Normalize with vt :\n"+record);

		//find interbase coor of effected bps
		int[] startStop = QueryIndexFileLoader.fetchEffectedBps(fields, true);
		if (startStop == null) throw new IOException ("Failed to parse the effected bps for : "+record);

		ArrayList<Double[]> maxAfDpAL = new ArrayList<Double[]>();
		maxAfDpAL.add(new Double[]{-1.0, 0.0});
		int numSamples = 1;

		//pull mpileup record(s)
		int start = startStop[0]+1;
		if (start < 1) start = 1;
		String tabixCoor = fields[0]+":"+start+"-"+(startStop[1]);
		TabixReader.Iterator it = fetchInteratorOnCoordinates(tabixCoor);
		
		//any records? Might be none
		if (it != null) {

			//for each mpileup record, find highest AF, should only be one for snvs
			String mpileupLine = null;

			while ((mpileupLine = it.next()) != null){

				//parse mpilup line and pull filtered set of samples
				MpileupLine ml = new MpileupLine(mpileupLine, minBaseQuality);
				if (ml.getChr() == null) throw new IOException ("Failed to parse the mpileup line:\n"+mpileupLine);

				//if snv them check position and ref
				if (type == 0){
					int vcfPos = Integer.parseInt(fields[1]);
					if (vcfPos != (ml.getZeroPos()+1))  throw new IOException ("Failed to match snv position with the mpileup line:\n"+mpileupLine+"\n"+record);
					//check ref
					if (fields[3].equals(ml.getRef()) == false) throw new IOException ("Failed to match snv ref with the mpileup line:\n"+mpileupLine+"\n"+record);
				}

				//adjust for multisamples
				MpileupSample[] samples = ml.getSamples();
				if (samples.length > 1 && numSamples == 1) {
					numSamples = samples.length;
					//add more defaults
					for (int i=1; i< numSamples; i++) maxAfDpAL.add(new Double[]{-1.0, 0.0});
				}
				
				//for each sample
				for (int i=0; i< numSamples; i++) {
					//fetch Double[]
					Double[] maxAfDp = maxAfDpAL.get(i);
					
					//calculate AF and DP
					double[] afDp = fetchAfDp(samples[i], fields, type);
					if (afDp[0] > maxAfDp[0]) {
						maxAfDp[0] = afDp[0];
						maxAfDp[1] = afDp[1];
					}
					
					if (verbose) {
						IO.pl("\t"+afDp[0]+"\t"+afDp[1]);
						samples[i].debug();
					}
				}
			}
		}
		else if (verbose) IO.pl("\tFound no mpileup records, appending AF=0;DP=0;");
		
		
		//save first af and dp
		if (maxAfDpAL.size() != 0) {
			Double[] maxAfDp = maxAfDpAL.get(0);
			alleleFreqs.add(maxAfDp[0]);
			readDepth.add(maxAfDp[1]);	
		}
		
		//append vcf info
		String modInfo = modifyInfo(maxAfDpAL, fields[7]);
		fields[7] = modInfo;
		modVcfRecords.add(Misc.stringArrayToString(fields, "\t"));
	}	


	private String modifyInfo(double[] afDp, String info) {
		//remove AF and DP if present
		Matcher mat = AF.matcher(info);
		String modInfo = mat.replaceFirst("");
		Matcher pat = DP.matcher(modInfo);
		String mod2Info = pat.replaceFirst("");
		if (afDp[0] == -1) return "AF=0;DP=0;" +mod2Info;
		return "AF="+Num.formatNumber(afDp[0], numDecimals)+";DP="+(int)afDp[1]+";"+mod2Info;
	}
	
	private String modifyInfo(ArrayList<Double[]> maxAfDpAL, String info) {
		//collect AFs and DPs
		String[] afs = new String[maxAfDpAL.size()];
		String[] dps = new String[afs.length];
		for (int i=0; i< afs.length; i++) {
			Double[] afDp = maxAfDpAL.get(i);
			if (afDp[0] == -1.0) {
				afs[i] = "0";
				dps[i] = "0";
			}
			else {
				afs[i] = Num.formatNumber(afDp[0], numDecimals);
				dps[i] = new Integer((int)afDp[1].doubleValue()).toString();
			}
		}
		String af = Misc.stringArrayToString(afs, ",");
		String dp = Misc.stringArrayToString(dps, ",");
		//remove AF and DP if present
		Matcher mat = AF.matcher(info);
		String modInfo = mat.replaceFirst("");
		Matcher pat = DP.matcher(modInfo);
		String mod2Info = pat.replaceFirst("");
		return "AF="+af+";DP="+dp+";"+mod2Info;
	}

	private static double[] fetchAfDp(MpileupSample sample, String[] fields, int type) throws IOException {
		//pull ref counts from curr pos
		double refCounts = sample.getBaseCounts(MpileupSample.getIndex(sample.getRecord().getRef().charAt(0)));

		//snv?
		if (type == 0){
			int baseIndex = MpileupSample.getIndex(fields[4].charAt(0));
			if (baseIndex == -1) throw new IOException("\nFailed to fetch the base index for "+fields[4].charAt(0)+ 
					" from the Alt "+fields[4]+" in "+Misc.stringArrayToString(fields,  "\t"));
			double altCounts = sample.getBaseCounts(baseIndex);
			double total = altCounts+refCounts;
			return new double[]{altCounts/total, total};
		}

		//del  GAC -> G
		if (type == 1){
			double delCounts = sample.getDeletions();
			double total = refCounts + delCounts;
			return new double[]{delCounts/total, total};
		}

		//must be insertion, G-> GAC
		double insCounts = sample.getInsertions();
		double total = refCounts + insCounts;
		return new double[]{insCounts/total, total};
	}

	private TabixReader.Iterator fetchInteratorOnCoordinates(String coordinates) {
		TabixReader.Iterator it = null;
		//watch out for no retrieved data error from tabix
		try {
			it = tabixReader.query(coordinates);
		} catch (ArrayIndexOutOfBoundsException e){}
		return it;
	}

	private MpileupSample fetchMpileupSample (MpileupLine ml) throws Exception{
		MpileupSample[] allSamples = ml.getSamples();
		if (allSamples.length != 1) throw new Exception ("\nMore than one sample found in the mpileup line "+ml.getLine());
		return allSamples[0];
	}

	public boolean isFailed() {
		return failed;
	}

	public TabixReader getTabixReader() {
		return tabixReader;
	}

	public ArrayList<Double> getAlleleFreqs() {
		return alleleFreqs;
	}

	public ArrayList<Double> getReadDepth() {
		return readDepth;
	}
}
