package edu.utah.seq.parsers.mpileup;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.query.QueryIndexFileLoader;
import edu.utah.seq.vcf.VCFBackgroundChecker;
import edu.utah.seq.vcf.VCFBackgroundCheckerGraphite;
import edu.utah.seq.vcf.VCFBamAnnotator;
import edu.utah.seq.vcf.xml.SimpleVcf;
import htsjdk.tribble.readers.TabixReader;
import util.gen.Misc;
import util.gen.Num;

/**Pulls mpileup lines that overlap each vcf record and calculates R1 and R2 coverage stats.*/
public class MpileupTabixLoaderRs implements Runnable{

	//fields
	private boolean failed = false;
	private VCFBamAnnotator vbc;
	private TabixReader tabixReader = null;
	private ArrayList<String> vcfLines = new ArrayList<String>();
	private ArrayList<String> modVcfRecords = new ArrayList<String>();
	private ArrayList<String> notAnnotated = new ArrayList<String>();
	private int minBaseQuality;
	private boolean verbose;
	

	public MpileupTabixLoaderRs (File mpileupFile, VCFBamAnnotator vbc) throws IOException{
		this.vbc = vbc;
		minBaseQuality = vbc.getMinBaseQuality();
		verbose = vbc.isVerbose();
		tabixReader = new TabixReader(mpileupFile.getCanonicalPath());
	}
	
	public void run() {	
		try {
			//get next chunk of work
			while (vbc.loadVcfRecords(vcfLines)){ 
				
				//for each
				for (String record: vcfLines) score(record);
				
				//add back records
vbc.addProcessedVcfs(modVcfRecords, notAnnotated);

				//cleanup
				vcfLines.clear();
				modVcfRecords.clear();
				notAnnotated.clear();
			}
		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem fetching mpileup lines\n" );
			e.printStackTrace();
		}
	}
	
	private void score(String record) throws Exception {
		//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR
		String[] fields = Misc.TAB.split(record);
		SimpleVcf vcf = new SimpleVcf(record, 0);
		
//System.out.println("\nVCF\t"+record);
		
		//interbase coor of effected bps
		int[] startStop = QueryIndexFileLoader.fetchEffectedBps(fields, true);
		if (startStop == null) {
			System.err.println("Failed to parse the effected bps for "+record);
			notAnnotated.add(record);
			return;
		}
		
//System.out.println("SS\t"+startStop[0]+" "+startStop[1]);
		
		//pull mpileup records, if none then warn and save vcf record
		int start = startStop[0]+1;
		if (start < 1) start = 1;
		String tabixCoor = fields[0]+":"+start+"-"+(startStop[1]);
		TabixReader.Iterator it = fetchInteratorOnCoordinates(tabixCoor);
		if (it == null) {
			System.err.println("Failed to fetch mpileup data for "+record);
			notAnnotated.add(record);
			return;
		}
		
		//load mpileup records
		ArrayList<MpileupLine> mLines = new ArrayList<MpileupLine>();
		String mpileupLine = null;
		while ((mpileupLine = it.next()) != null){
			//parse mpilup line and pull filtered set of samples
			MpileupLine ml = new MpileupLine(mpileupLine, minBaseQuality);
			if (ml.getChr() == null) throw new IOException ("Failed to parse the mpileup line:\n"+record);
			mLines.add(ml);
		}
		
		//snp?
		if (vcf.isSnv()){
			//check there is only one line
			if (mLines.size() !=1) throw new IOException ("Multiple mpileup lines found for a snv?!: "+mLines.size()+"\n"+record);
			//check reference
			MpileupLine ml = mLines.get(0);
			if (vcf.getRef().equals(ml.getRef()) == false) throw new IOException ("Reference of vcf and mpileup differ?!:\n"+mpileupLine +"\n"+record);
			//check position
			if (vcf.getPos() != ml.getZeroPos()) throw new IOException ("The position of vcf and mpileup differ?!:\n"+mpileupLine +"\n"+record);
			//calculate AF and DP for each sample
			MpileupSample[] r1r2 = ml.getSamples();
			
			StringBuilder sb = new StringBuilder(vcf.toString());
			sb.append("\tGT:AF:DP:R1:R2:DPR1:DPR2");
			for (int i=0; i< r1r2.length; i++){
				
				int altBaseIndex = MpileupSample.getIndex(vcf.getAlt().toUpperCase().charAt(0));
				//read one
				int R1 = r1r2[i].getForwardGATC()[altBaseIndex] + r1r2[i].getReverseGATC()[altBaseIndex];
				int DPR1 = r1r2[i].getReadCoverageForwardBases() + r1r2[i].getReadCoverageReverseBases();
				i++;
				//read two
				int R2 = r1r2[i].getForwardGATC()[altBaseIndex] + r1r2[i].getReverseGATC()[altBaseIndex];
				int DPR2 = r1r2[i].getReadCoverageForwardBases() + r1r2[i].getReadCoverageReverseBases();
				//combine
				double AF = (double)(R1+R2)/ (DPR1+DPR2);
				int DP = DPR1 + DPR2;
				
				//make string  GT:AF:DP:R1:R2:DPR1:DPR2, GT is required by some apps
				sb.append("\t./.:");
				sb.append(Num.formatNumberJustMax(AF, 4)); sb.append(":");
				sb.append(DP); sb.append(":");
				sb.append(R1); sb.append(":");
				sb.append(R2); sb.append(":");
				sb.append(DPR1); sb.append(":");
				sb.append(DPR2);
			}
			System.out.println(sb);
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

	public boolean isFailed() {
		return failed;
	}

	public TabixReader getTabixReader() {
		return tabixReader;
	}
}
