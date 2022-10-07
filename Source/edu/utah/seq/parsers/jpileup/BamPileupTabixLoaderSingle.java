package edu.utah.seq.parsers.jpileup;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import edu.utah.seq.query.QueryIndexFileLoader;
import htsjdk.tribble.readers.TabixReader;
import util.bio.annotation.Coordinate;
import util.gen.IO;
import util.gen.Misc;

/**Pulls bpileup lines that overlap vcf records.*/
public class BamPileupTabixLoaderSingle {

	//fields
	private boolean failed = false;
	private TabixReader tabixReader = null;
	private static boolean debug = true;
	private int bpPad = 0;

	//constructor
	public BamPileupTabixLoaderSingle (File bpileupFile, int bpPad) throws IOException{
		tabixReader = new TabixReader(bpileupFile.getCanonicalPath());
		this.bpPad = bpPad;
	}
	
	//primary method, provide tab split vcf record
	//#CHROM POS	ID	REF	ALT	QUAL FILTER	INFO...
	public ArrayList<BpileupLine> fetchBpileupRecords(String[] vcfFields) throws Exception {
			
		//what kind of variant, returns GATCID
		char allele = fetchAllele(vcfFields);

		//fetch the coordinates and make the query
		int[] startStop = fetchStartStop(vcfFields, allele);
		if (startStop == null) throw new Exception("ERROR: failed to parse start stop corr for "+Misc.stringArrayToString(vcfFields, "\t"));
		String tabixCoor = vcfFields[0]+":"+ startStop[0]+ "-"+ startStop[1];
		
		//return the lines
		return fetchLines(tabixCoor);
	}
	
	public ArrayList<BpileupLine> fetchBpileupRecords(Coordinate region) throws Exception {	
		String tabixCoor = region.getTabixSearchCoordinates();
		return fetchLines(tabixCoor);
	}
	
	private ArrayList<BpileupLine> fetchLines(String tabixCoor) throws IOException, Exception {
		ArrayList<BpileupLine> al = new ArrayList<BpileupLine>();
		TabixReader.Iterator it = fetchInteratorOnCoordinates(tabixCoor);
		if (it == null) {
			if (debug) IO.pl("No bpileup lines for "+ tabixCoor);
		}
		else {
			String bpileupLine = null;
			while ((bpileupLine = it.next()) != null)al.add(new BpileupLine(bpileupLine));
		}
		return al;
	}

	private int[] fetchStartStop(String[] fields, char allele) {
		int[] startStop = null;
		if (allele == 'D') {
			startStop = QueryIndexFileLoader.fetchEffectedBps(fields, true);
			if (startStop == null) return null;
			startStop[0] = startStop[0]+2;
		}
		else if (allele == 'I') {
			int s = Integer.parseInt(fields[1])+1;
			startStop = new int[] {s,s};
		}
		//snv
		else {
			int s = Integer.parseInt(fields[1]);
			startStop = new int[] {s,s};
		}
		//pad the bps?
		if (bpPad !=0) {
			startStop[0] = startStop[0]-bpPad;
			if (startStop[0] <0) startStop[0] = 0;
			startStop[1] = startStop[1]+ bpPad;
		}
		return startStop;
	}
	

	/**Return GATC or ID for indels*/
	public static char fetchAllele(String[] fields) throws IOException {
		if (fields[4].contains(",") || fields[4].startsWith("<")) throw new IOException("Cannot interpret multi alts or those with <, deconvolute? "+Misc.stringArrayToString(fields,  "\t"));
		//#CHROM	POS	ID	REF	ALT
		int lenRef = fields[3].length();
		int lenAlt = fields[4].length();
		//snv?
		if (lenRef == 1 && lenAlt == 1) return fields[4].charAt(0);
		//indel
		if (lenRef > lenAlt) return 'D';
		return 'I';
	}
	
	private TabixReader.Iterator fetchInteratorOnCoordinates(String coordinates) {
		TabixReader.Iterator it = null;
		//watch out for no retrieved data error from tabix
		try {
			it = tabixReader.query(coordinates);
		} catch (ArrayIndexOutOfBoundsException e){}
		return it;
	}
	
	/*
	public static void main (String[] args) throws Exception {
		String del = "chr5\t68293834\tSSC_262958\tGGT\tG\t1000\t.\tBKZ=1000";
		String ins = "chr7\t140753338\tSSC_345020\tT\tTGTA\t84.78\t.\tBKZ=84.78";
		String snv = "chr5\t68280559\tSSC_261533\tA\tT\t0\tBKAF\tBKZ=0";
		String[] vars = {del,ins,snv};
		BamPileupTabixLoaderSingle bp = new BamPileupTabixLoaderSingle(new File("/Users/u0028003/Downloads/SSCIndelScanning/15352X17.bp.txt.gz"), 0);
		
		for (String s: vars) {
if (debug) IO.pl("\n"+s);
			ArrayList<BpileupLine> bpl = bp.fetchBpileupRecords(Misc.TAB.split(s));
		}
		
		bp.getTabixReader().close();
	}*/

	public boolean isFailed() {
		return failed;
	}

	public TabixReader getTabixReader() {
		return tabixReader;
	}

	public int getBpPad() {
		return bpPad;
	}

	public void setBpPad(int bpPad) {
		this.bpPad = bpPad;
	}
}
