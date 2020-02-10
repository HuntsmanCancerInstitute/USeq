package edu.utah.seq.query;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import util.gen.IO;
import util.gen.Misc;

public class QueryIndexFileLoader implements Runnable {

	//fields
	private boolean failed = false;
	private QueryIndexer queryIndexer;
	private String chrom;
	private TabixReader reader = null;
	private File sourceFile = null;
	private boolean verbose;
	private ArrayList<IndexRegion> toAdd = new ArrayList<IndexRegion>();
	private HashMap<File, Integer> fileId = null;
	private static final int numToLoad = 10000;
	private int start = 0;
	private int stop = 0;
	public static final Pattern END_POSITION = Pattern.compile(".*END=(\\d+).*", Pattern.CASE_INSENSITIVE);

	public QueryIndexFileLoader (QueryIndexer queryIndexer, String chrom, int start, int stop) throws IOException{
		this.queryIndexer = queryIndexer;
		this.fileId = queryIndexer.getFileId();
		this.chrom = chrom;
		this.start = start;
		this.stop = stop;
		verbose = queryIndexer.isVerbose();
	}

	public void run() {	
		try {
			//get next file to parse
			while ((sourceFile = queryIndexer.getFileToParse()) != null){ 
				toAdd.clear();
				//fetch a reader on slice, start must be > 0, end can be past chrom end
				reader = new TabixReader(sourceFile.toString());
				TabixReader.Iterator it = fetchTReader(chrom, start+1, stop);
				if (it == null) {
					reader.close();
					continue;
				}

				parseDataFile(it);
				
				//add last and clear
				if (toAdd.size()!=0) queryIndexer.addRegions(toAdd);
				reader.close();
			}
		} catch (Exception e) {
			failed = true;
			System.err.println("Error: loading "+chrom+":"+start+"-"+stop+" for "+sourceFile );
			e.printStackTrace();
		} finally {
			if (reader != null) reader.close();
		}
	}

	private void parseDataFile(Iterator it) throws NumberFormatException, IOException {
		//is it a vcf?
		boolean vcf = sourceFile.getName().toLowerCase().endsWith(".vcf.gz");
		int[] startStopSubtract = null;
		if (vcf == false) startStopSubtract = queryIndexer.getSSS(sourceFile);

		Integer id = fileId.get(sourceFile);
		String record = null;
		while ((record = it.next()) != null){
			//parse start and stop bp positions
			int[] startStop;
			if (vcf) startStop = parseStartStopBpCoorVcf(record, verbose);
			else startStop = parseStartStopBpCoor(record, startStopSubtract, verbose);
			if (startStop != null) {
				//check against the chrom seg
				if (startStop[0] < start) startStop[0] = start;
				if (startStop[1] >= stop) startStop[1] = stop;
				if ((startStop[1]-startStop[0]) > 0) {
					toAdd.add(new IndexRegion(startStop[0]-start, startStop[1]-start, id));
					if (toAdd.size() > numToLoad) {
						queryIndexer.addRegions(toAdd);
						toAdd.clear();
					}
				}
			}
		}
	}

	private TabixReader.Iterator fetchTReader(String chrom, int start, int stop) throws Exception{
		if (reader.getChromosomes().contains(chrom)) return reader.query(chrom, start, stop);
		String chrChrom = "chr"+chrom;
		if (reader.getChromosomes().contains(chrChrom)) return reader.query(chrChrom, start, stop);
		return null;
	}

	private static int[] parseStartStopBpCoorVcf(String vcfRecord, boolean verbose) {
		int[] startStop = null;
		try {
			String[] t = Misc.TAB.split(vcfRecord);
			startStop = fetchEffectedBps(t, verbose);
		} catch (NumberFormatException e){}
		return startStop;
	}

	private static int[] parseStartStopBpCoor(String record, int[] startStopSub, boolean verbose) throws NumberFormatException {
		int[] startStop = null;
		try {
			String[] t = Misc.TAB.split(record);
			int start = Integer.parseInt(t[startStopSub[0]]);
			start = start- startStopSub[2];
			int stop = Integer.parseInt(t[startStopSub[1]]);
			startStop = new int[]{start, stop};
		} catch (NumberFormatException e){
			if (verbose) System.err.println("\tWARNING: failed to parse start stop, skipping -> "+record);
		}
		return startStop;
	}

	/**Returns the interbase start stop region of effected bps for simple SNV and INDELs. 
	 * SNV=iPos,iPos+LenRef; INS=iPos,iPos+lenRef+1; DEL=iPos+1,iPos+lenRef; iPos=pos-1.
	 * For multi alts, returns min begin and max end of all combinations.
	 * For alts with < indicating a CNV or trans, attempts to parse the END=number from the INFO column. */
	public static int[] fetchEffectedBps(String[] vcf, boolean verbose) throws NumberFormatException{

		/*NOTE, any changes here, please update the web app's QueryRequest too*/

		//CHROM	POS	ID	REF	ALT	
		//  0    1   2   3   4  
		//put into interbase coordinates
		int iPos = Integer.parseInt(vcf[1]) - 1;
		String ref= vcf[3];
		String alt= vcf[4];

		//any commas/ multi alts?
		if (alt.contains(",") == false) return fetchEffectedBpsSingleAlt(iPos, ref, alt, vcf,true, verbose);

		//OK commas present, thus multi alts, these need to be tested for max effect
		//There is complexity with multi alts, best to deconvolute and left justify! 
		String[] alts = Misc.COMMA.split(alt);
		int begin = Integer.MAX_VALUE;
		int end = -1;		
		for (int i=0; i< alts.length; i++){
			int[] ss = fetchEffectedBpsSingleAlt(iPos, ref, alts[i], vcf, false, verbose);
			//skip it?
			if (ss == null) continue;
			if(ss[0] < begin) begin = ss[0];
			if(ss[1]> end) end = ss[1];
		}

		if (begin == Integer.MAX_VALUE) return null;

		return new int[]{begin, end};
	}

	public static int[] fetchEffectedBpsSingleAlt(int iPos, String ref, String alt, String[] vcf, boolean singleAlt, boolean verbose) throws NumberFormatException{

		/*NOTE, any changes here, please update the web app's QueryRequest too*/

		int begin = -1;
		int end = -1;
		int lenRef = ref.length();
		int lenAlt = alt.length();

		//watch out for < in the alt indicative of a CNV, structural var, or gvcf block
		if (alt.contains("<")){
			//CHROM	POS	ID	REF	ALT QUAL FILTER INFO	
			//  0    1   2   3   4   5      6     7
			Matcher mat = END_POSITION.matcher(vcf[7]);
			if (mat.matches()) end = Integer.parseInt(mat.group(1));
			else {
				if (verbose && singleAlt) System.err.println("\tWARNING: found a < containing alt, yet failed to parse END=number position, skipping if this is the only alt -> "+Misc.stringArrayToString(vcf, "\t"));
				return null;
			}
			begin = iPos;
		}

		//single or multi adjacent snp? return just the changed bps,  GC->AT or G->A or G->.
		else if (lenRef == lenAlt) {
			begin = iPos;
			end = iPos+ lenRef;
		}
		//ins? return the bases on either side of the insertion GC->GATTA or G->ATTA
		else if (lenAlt > lenRef) {
			begin = iPos;
			end = iPos+ lenRef +1;
			//C->TTAAT
			if (ref.charAt(0) != alt.charAt(0)) {
				begin--;
				//if (verbose) System.err.println("\tWARNING: Odd INS vcf record, the first base in the ref and alt should be the same, -> "+Misc.stringArrayToString(vcf, "\t")+ " effected bps "+begin+" "+end);
			}

		}
		//del? return the uneffected bp and those that are deleted to match tabix's behaviour AT->A, ATTCG->ACC
		else if (lenRef > lenAlt) {
			begin = iPos;
			end = iPos + lenRef;
			//AT->C
			if (ref.charAt(0) != alt.charAt(0)) {
				begin--;
				//if (verbose) System.err.println("\tWARNING: Odd DEL vcf record, the first base in the ref and alt should be the same -> "+Misc.stringArrayToString(vcf, "\t")+ " effected bps "+begin+" "+end);
			}
		}
		//odd, shouldn't hit this
		else {
			if (verbose) System.err.println("ERROR: Contact admin! Odd vcf record, can't parse effected bps for -> "+Misc.stringArrayToString(vcf, "\t"));
			return null;
		}

		return new int[]{begin, end};
	}


	public boolean isFailed() {
		return failed;
	}

}
