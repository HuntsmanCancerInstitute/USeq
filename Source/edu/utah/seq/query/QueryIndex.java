package edu.utah.seq.query;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.useq.data.RegionScoreText;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;

public class QueryIndex {

	private DataSources dataSources;
	private TQuery tQuery;
	
	public static final Pattern END_POSITION = Pattern.compile(".*END=(\\d+).*", Pattern.CASE_INSENSITIVE);
	//interbase coordinates, 0 is first base, last base in any region is not included.
	private HashMap<String, HashSet<File>[]> chromFileIndex = new HashMap<String, HashSet<File>[]>();
	
	
	//constructor
	public QueryIndex(TQuery tQuery) throws IOException{
		this.tQuery = tQuery;
		
		//make an empty index of chr and hash[] where every bp can have a hashSet of Files with an overlapping data point
		createChromIndex();
		
		//create a container to hold all the File info for per request filtering
		dataSources = new DataSources();
		
		//ok load em up!
		System.err.println("\nParsing data sources...");
		System.err.println("\tRecordsIndexed\tRecordsSkipped\tFile");
		if (tQuery.getVcfDataFiles() != null) loadChromIndexWithVcf();
		if (tQuery.getMafDataFiles() != null) loadChromIndexWithMaf();
		if (tQuery.getBedDataFiles() != null) loadChromIndexWithBed();
		
		//set in TQuery
		tQuery.setDataSources(dataSources);
	}

	//methods
	private void createChromIndex() throws IOException {
		HashMap<String, RegionScoreText[]> chrLen = Bed.parseBedFile(tQuery.getChrLengthFile(), true, false);

		System.err.println("Building empty chrom index...");
		for (String chr: chrLen.keySet()){
			RegionScoreText[] regions = chrLen.get(chr);
			if (regions.length !=1) throw new IOException("\nError: there can be only one bed region for each chromosome, see "+chr);
			HashSet<File>[] indexes = new HashSet[regions[0].getStop()+1];
			chromFileIndex.put(chr, indexes);
			System.err.println("\t"+chr+"\t"+regions[0].getStop());
		}
	}
	
	private void loadChromIndexWithVcf() throws IOException {
		File[] vcfDataFiles = tQuery.getVcfDataFiles();
		boolean printWarnings = tQuery.isPrintWarnings();
		
		for (int i=0; i< vcfDataFiles.length; i++){

			//load filter with file, extension, parent dir name
			dataSources.addFileToFilter(vcfDataFiles[i]);

			BufferedReader in = IO.fetchBufferedReader(vcfDataFiles[i]);

			String[] t;
			String line;
			String currChrom = "";
			HashSet<File>[] currIndex = null;
			long numLoadedRecords = 0;
			long numRecordsSkipped = 0;
			while ((line = in.readLine()) != null){
				if (line.startsWith("#") == false){
					t = Misc.TAB.split(line);

					//diff chrom? pull index
					if (currChrom.equals(t[0]) == false){
						currIndex = chromFileIndex.get(t[0]);
						if (currIndex == null) {
							if (tQuery.isPrintWarnings()) System.err.println("\tWARNING: Failed to find a chromosome for '"+t[0]+ "' skipping -> "+line);
							numRecordsSkipped++;
							continue;
						}
						currChrom = t[0];
					}

					//fetch start and stop for effected bps.
					int[] startStop = fetchEffectedBps(t, printWarnings);
					if (startStop == null) numRecordsSkipped++;
					else {
						//add in references to source file over the covered bases, stop isn't covered.
						for (int j= startStop[0]; j< startStop[1]; j++){
							if (currIndex[j] == null) currIndex[j] = new HashSet<File>(1);
							currIndex[j].add(vcfDataFiles[i]);
						}
						numLoadedRecords++;
					}
				}
			}

			//clean up and stat incrementing
			in.close();
			System.err.println("\t"+numLoadedRecords+"\t"+numRecordsSkipped+ "\t"+vcfDataFiles[i]);
			dataSources.setRecordsLoaded(dataSources.getRecordsLoaded() + numLoadedRecords); 
			dataSources.setRecordsSkipped(dataSources.getRecordsSkipped()+ numRecordsSkipped);
		}
	}

	/**Returns the interbase start stop region of effected bps for simple SNV and INDELs. 
	 * SNV=iPos,iPos+LenRef; INS=iPos,iPos+lenRef+1; DEL=iPos+1,iPos+lenRef; iPos=pos-1.
	 * For multi alts, returns min begin and max end of all combinations.
	 * For alts with < indicating a CNV or trans, attempts to parse the END=number from the INFO column. */
	public static int[] fetchEffectedBps(String[] vcf, boolean printWarnings) throws IOException{
		//CHROM	POS	ID	REF	ALT	
		//  0    1   2   3   4  
		//put into interbase coordinates
		int iPos = Integer.parseInt(vcf[1]) - 1;
		String ref= vcf[3];
		String alt= vcf[4];

		//any commas/ multi alts?
		if (alt.contains(",") == false) return fetchEffectedBpsNoMultiAlt(iPos, ref, alt, vcf, printWarnings);

		//OK commas present, these need to be tested for max effect
		//There is complexity with multi alts, best to deconvolute and left justify! 
		String[] alts = Misc.COMMA.split(alt);
		int begin = Integer.MAX_VALUE;
		int end = -1;
		for (int i=0; i< alts.length; i++){
			int[] ss = fetchEffectedBpsNoMultiAlt(iPos, ref, alts[i], vcf, printWarnings);
			if (ss == null) return null;
			if(ss[0] < begin) begin = ss[0];
			if(ss[1]> end) end = ss[1];
		}

		return new int[]{begin, end};
	}

	public static int[] fetchEffectedBpsNoMultiAlt(int iPos, String ref, String alt, String[] vcf, boolean printWarnings) throws IOException{
		int begin = -1;
		int end = -1;
		int lenRef = ref.length();
		int lenAlt = alt.length();

		//watch out for < in the alt indicative of a CNV or structural var
		if (alt.contains("<")){
			//CHROM	POS	ID	REF	ALT QUAL FILTER INFO	
			//  0    1   2   3   4   5      6     7
			Matcher mat = END_POSITION.matcher(vcf[7]);
			if (mat.matches()) end = Integer.parseInt(mat.group(1));
			else {
				if (printWarnings) System.err.println("\tWARNING: found a < or > containing alt, failed to parse END=number position, skipping -> "+Misc.stringArrayToString(vcf, "\t"));
				return null;
			}
			begin = iPos;
		}

		//single or multi adjacent snp? return just the changed bps,  GC->AT or G->A
		else if (lenRef == lenAlt) {
			begin = iPos;
			end = iPos+ lenRef;
		}
		//ins? return the bases on either side of the insertion GC->GATTA or G->ATTA
		else if (lenAlt > lenRef) {
			begin = iPos;
			end = iPos+ lenRef +1;
			if (ref.charAt(0) != alt.charAt(0)) {
				if (printWarnings) System.err.println("\tWARNING: Odd INS vcf record, the first base in the ref and alt must be the same, use vt to normalize your variants, skipping -> "+Misc.stringArrayToString(vcf, "\t"));
				return null;
			}

		}
		//del? return the bps that are deleted, AT->A, ATTCG->ACC
		else if (lenRef > lenAlt) {
			begin = iPos+1;
			end = iPos + lenRef;
			if (ref.charAt(0) != alt.charAt(0)) {
				if (printWarnings) System.err.println("\tWARNING: Odd DEL vcf record, the first base in the ref and alt must be the same, use vt to normalize your variants, skipping -> "+Misc.stringArrayToString(vcf, "\t"));
				return null;
			}
		}
		//odd, shouldn't hit this
		else throw new IOException("\nError: Contact admin! Odd vcf record, can't parse effected bps for -> "+Misc.stringArrayToString(vcf, "\t"));

		return new int[]{begin, end};
	}

	/**This executes a search loading the results into the QueryRequest's hash maps.*/
	public void queryFileIndex(QueryRequest qr) {
		long startTime = System.currentTimeMillis();
		long numQueries = 0;
		long numSkippedQueries = 0;
		long numQueriesWithHits = 0;
		long numberHits = 0;
		
		HashMap<File, ArrayList<TabixQuery>> fileTabixQueries = qr.getFileTabixQueries();
		HashMap<String, TabixQuery[]> chrTabixQueries = qr.getChrTabixQueries();

		//for each chromosome of regions to query
		for (String chr: chrTabixQueries.keySet()){

			//check that chr exists in index
			if (chromFileIndex.containsKey(chr) == false){
				int numSkipped = chrTabixQueries.get(chr).length;
				if (tQuery.isPrintWarnings()) System.err.println("\nWARNING: chromosome '"+chr+"' not found in index? Skipping "+numSkipped+" query regions.");
				numSkippedQueries += numSkipped; 
			}

			else {
				HashSet<File>[] index = chromFileIndex.get(chr);
				//for each region
				TabixQuery[] regions = chrTabixQueries.get(chr);
				numQueries += regions.length;
				for (TabixQuery tq: regions){
					HashSet<File> fileHits = intersect(index, tq);
					if (fileHits.size() !=0) {
						numberHits += fileHits.size();
						numQueriesWithHits++;
					}
					addHits(fileHits, fileTabixQueries, tq);
				}
			}
		}
		if (tQuery.isPrintStats()){
			long diffTime = System.currentTimeMillis() -startTime;
			System.err.println("\nQuery Stats Pre Filtering:");
			System.err.println(numQueries+ "\tNum user queries");
			System.err.println(numSkippedQueries+ "\tNum skipped user queries");
			System.err.println(numQueriesWithHits+ "\tNum queries with hits");
			System.err.println(numberHits+ "\tNum hits");
			System.err.println(diffTime+"\tMillisec to complete index search\n");
		}
	}

	/**Adds the TabixQuery to an ArrayList associated with a file resource to fetch the data from.*/
	private static void addHits(HashSet<File> hits, HashMap<File, ArrayList<TabixQuery>> toQuery, TabixQuery tq) {
		for (File fHit: hits){
			ArrayList<TabixQuery> al = toQuery.get(fHit);
			if (al == null){
				al = new ArrayList<TabixQuery>();
				toQuery.put(fHit, al);
			}
			al.add(tq);
		}
	}

	/**Checks the begin and end for out of bounds then uses a hash to collapse all the Files found to have region that overlaps the query.*/
	private static HashSet<File> intersect(HashSet<File>[] index, TabixQuery tq) {
		HashSet<File> hits = new HashSet<File>();
		int begin = tq.getStart();
		if (begin < 0) begin = 0;
		int end = tq.getStop();
		if (end >= index.length) end = index.length;
		for (int i=begin; i<end; i++) {
			if (index[i] != null) hits.addAll(index[i]);
		}
		return hits;
	}

	private void loadChromIndexWithBed() throws IOException {
		File[] bedDataFiles = tQuery.getBedDataFiles();
		
		for (int i=0; i< bedDataFiles.length; i++){

			//load filter with file, extension, parent dir name
			dataSources.addFileToFilter(bedDataFiles[i]);

			BufferedReader in = IO.fetchBufferedReader(bedDataFiles[i]);

			String[] t;
			String line;
			String currChrom = "";
			HashSet<File>[] currIndex = null;
			long numLoadedRecords = 0;
			long numSkippedRecords = 0;
			while ((line = in.readLine()) != null){
				if (line.startsWith("#") == false){
					t = Misc.TAB.split(line);
					//diff chrom?
					if (currChrom.equals(t[0]) == false){
						currIndex = chromFileIndex.get(t[0]);
						if (currIndex == null) {
							if (tQuery.isPrintWarnings()) System.err.println("\nError: Failed to find a chromosome for '"+t[0]+ "' from "+bedDataFiles[i]+" skipping line "+line);
							numSkippedRecords++;
							continue;
						}
						currChrom = t[0];
					}

					//parse start and stop, note interbase coordinates
					int start = Integer.parseInt(t[1]);
					int stop = Integer.parseInt(t[2]);

					//add in references to source file over the covered bases, stop isn't covered.
					for (int j= start; j< stop; j++){
						if (currIndex[j] == null) currIndex[j] = new HashSet<File>(1);
						currIndex[j].add(bedDataFiles[i]);
					}
					numLoadedRecords++;
				}
			}
			//clean up and stat incrementing
			in.close();
			System.err.println("\t"+numLoadedRecords+"\t"+numSkippedRecords+"\t"+bedDataFiles[i]);
			dataSources.setRecordsLoaded(dataSources.getRecordsLoaded() + numLoadedRecords); 
			dataSources.setRecordsSkipped(dataSources.getRecordsSkipped()+ numSkippedRecords);
		}
	}

	private void loadChromIndexWithMaf() throws IOException {
		File[] mafDataFiles = tQuery.getMafDataFiles();
		
		for (int i=0; i< mafDataFiles.length; i++){

			//load filter with file, extension, parent dir name
			dataSources.addFileToFilter(mafDataFiles[i]);

			BufferedReader in = IO.fetchBufferedReader(mafDataFiles[i]);

			String[] t;
			String line;
			String currChrom = "";
			int chromIndex = -1;
			int startIndex = -1;
			int endIndex = -1;
			HashSet<File>[] currIndex = null;
			long numLoadedRecords = 0;
			long numSkippedRecords = 0;
			while ((line = in.readLine()) != null){

				//skip blanks and comments
				if (line.trim().length() == 0 || line.startsWith("#")) continue;

				//set indexes, if these don't get set they an array index out of bounds exception will be thrown below
				if (line.startsWith("Hugo_Symbol")){
					//Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_position End_position
					//   0               1          2       3          4            5            6
					HashMap<String, Integer> map = new HashMap<String, Integer>();
					String[] header = Misc.TAB.split(line);
					for (int x=0; x< header.length; x++) map.put(header[x], x);
					chromIndex = map.get("Chromosome");
					startIndex = map.get("Start_position");
					endIndex = map.get("End_position");
					continue;
				}

				//dataline
				t = Misc.TAB.split(line);

				//diff chrom?
				if (currChrom.equals(t[chromIndex]) == false){
					currIndex = chromFileIndex.get(t[chromIndex]);
					if (currIndex == null) {
						if (tQuery.isPrintWarnings()) System.err.println("\tWARNING: Failed to find a chromosome for '"+t[chromIndex]+ "' skipping -> "+line);
						numSkippedRecords++;
						continue;
					}
					currChrom = t[chromIndex];
				}

				//parse start and stop, note sub one from start to convert to interbase coordinates
				int start = Integer.parseInt(t[startIndex]) -1;
				if (start < 0) start = 0;
				int stop = Integer.parseInt(t[endIndex]);
				if (stop > currIndex.length) stop = currIndex.length;

				//add in references to source file over the covered bases, stop isn't covered.
				for (int j= start; j< stop; j++){
					if (currIndex[j] == null) currIndex[j] = new HashSet<File>(1);
					currIndex[j].add(mafDataFiles[i]);
				}
				numLoadedRecords++;
			}

			//clean up and stat incrementing
			in.close();
			System.err.println("\t"+numLoadedRecords+"\t"+numSkippedRecords+"\t"+mafDataFiles[i]);
			dataSources.setRecordsLoaded(dataSources.getRecordsLoaded() + numLoadedRecords); 
			dataSources.setRecordsSkipped(dataSources.getRecordsSkipped()+ numSkippedRecords);
		}
	}

	public HashMap<String, HashSet<File>[]> getChromFileIndex() {
		return chromFileIndex;
	}

}
