package edu.utah.seq.query;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.json.JSONObject;
import edu.utah.seq.its.Interval;
import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
import edu.utah.seq.its.IntervalTree;
import edu.utah.seq.useq.data.RegionScoreText;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.readers.TabixReader;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class QueryIndexer {

	//user params
	private File dataDir;
	private File chrNameLength;
	private File indexDir;
	private boolean verbose = true;
	private File bgzip = null;
	private File tabix = null;
	
	//internal fields
	private File[] dataFilesToParse;
	private String[] fileExtToIndex = null;
	private String[] supportedFileExtensions = {"vcf.gz", "bed.gz", "bedGraph.gz", "maf.txt.gz"};
	private HashMap<String, int[]> extensionsStartStopSub = new HashMap<String, int[]>();
	private HashMap<String, RegionScoreText[]> chrLengths;
	private String[] stringHeaderPatternStarts = {"^#.+", "^browser.+", "^track.+", "^Hugo_Symbol.+"};
	public static final Pattern END_POSITION = Pattern.compile(".*END=(\\d+).*", Pattern.CASE_INSENSITIVE);
	private HashMap<String, ArrayList<File>> chrFiles = new HashMap<String, ArrayList<File>>();
	private HashMap<File, Integer> fileId = new HashMap<File, Integer>();
	private long totalRecordsProcessed = 0;
	
	//per chr fields
	private TreeSet<File>[] workingIndex = null;
	private HashMap<String, TreeSet<File>> fileNames2TheirPointers = null;
	private long workingPassed;
	private long workingFailed;
	private RegionScoreText workingCoor;
	private String workingChr;
	
	//constructor
	@SuppressWarnings("unchecked")
	public QueryIndexer(String[] args) {
		long startTime = System.currentTimeMillis();
		
		processArgs(args);
		
		System.out.println("Creating file id hash...");
		createFileIdHash();
		
		System.out.println("Checking tbi files for chr content...");
		createChrFiles();
		
		//for each chromosome
		if (verbose) System.out.println("Indexing records by chr...\n\tchr\t#pass\t#fail\tmemUsed");
		else System.out.print("Indexing records by chr...\n\t"); 
		for (String chr: chrLengths.keySet()){
			//any files?
			ArrayList<File> filesToParse = chrFiles.get(chr);
			if (filesToParse.size() == 0) continue; 
			
			//prep for new chr
			workingPassed = 0;
			workingFailed = 0;
			fileNames2TheirPointers = new HashMap<String, TreeSet<File>>();
			workingChr = chr;
			workingCoor = chrLengths.get(chr)[0];
			workingIndex = new TreeSet[workingCoor.getStop()+1];
			
			//for each file
			for (File f: filesToParse) parse(f);
			
			//save the index
			saveIndex();
			
			//cleanup
			if (verbose) System.out.println("\t"+chr+"\t"+workingPassed+"\t"+workingFailed+"\t"+IO.memory());
			else System.out.print(chr+" ");
			totalRecordsProcessed+= workingPassed;
			totalRecordsProcessed+= workingFailed;
		}
		if (verbose == false) System.out.println(": "+IO.memory());
		//bgzip and tabix
		System.out.println("Compressing master index...");
		bgzipAndTabixIndex();
		
		System.out.println("Saving file headers...");
		saveFileHeadersAndIds();
		
		System.out.print("Saving interval trees...\n\t");
		saveIntervalTrees();
		
		String diffTime = Num.formatNumberOneFraction(((double)(System.currentTimeMillis() -startTime))/1000/60);
		System.out.println("\n"+ diffTime+" Min to parse "+ totalRecordsProcessed+" records and build the query index trees");
	}
	

	private void createChrFiles() {
		//load hash to hold files with a particular chr
		for (String chr: chrLengths.keySet()) chrFiles.put(chr, new ArrayList<File>());
		
		try {
			//for each file
			for (File f: dataFilesToParse){
				//make an index
				File i = new File (f+".tbi");
				TabixIndex ti = new TabixIndex(i);
				//for each chromosome
				for (String chr: chrLengths.keySet()){
					if (ti.containsChromosome(chr)) chrFiles.get(chr).add(f);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem with testing indexes for particular chroms\n");
		}
	}
	

	private void createFileIdHash() {
		int counter = 0;
		for (File f: dataFilesToParse) fileId.put(f, counter++);
	}
	
	public void saveIntervalTrees(){
		try {
			File[] chrBed = IO.extractFiles(indexDir, ".qi.bed.gz");
			HashMap<String, IntervalST<int[]>> chrTrees = new HashMap<String, IntervalST<int[]>>();
			HashMap<String, int[]> idsFileIndex = new HashMap<String, int[]>();
			//for each file
			for (File f: chrBed){
				BufferedReader in = IO.fetchBufferedReader(f);
				IntervalST<int[]> st = new IntervalST<int[]>();
				String line;
				while ((line = in.readLine())!= null){
					String[] t = Misc.TAB.split(line);
					int start = Integer.parseInt(t[1]);
					int stop = Integer.parseInt(t[2]);
					int[] ids = idsFileIndex.get(t[3]);
					if (ids == null){
						ids = stringToInts(t[3], Misc.COMMA);
						idsFileIndex.put(t[3], ids);
					}
					//the end is included in IntervalST so sub 1 from end
					st.put(new Interval1D(start, stop-1), ids);
				}
				in.close();

				//save tree
				String chr = f.getName().replace(".qi.bed.gz", "");
				
				chrTrees.put(chr, st);
				System.out.print(chr+" ");
			}
			System.out.println(": "+IO.memory());
			
			File t = new File (indexDir, "its.obj");
			IO.saveObject(t, chrTrees);
			
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem saving interval trees, aborting.\n");
		}
	}
	
	/**Given a String of ints delimited by something, will parse or return null.*/
	public static int[] stringToInts(String s, Pattern pat){
		String[] tokens = pat.split(s);
		int[] num = new int[tokens.length];
		try {
			for (int i=0; i< tokens.length; i++){
				num[i] = Integer.parseInt(tokens[i]);
			}
			return num;
		} catch (Exception e){
			return null;
		}
	}


	public void saveIndex(){
		try {
			File queryIndexFile = new File(indexDir, workingChr+".qi.bed");
			PrintWriter out = new PrintWriter(new FileWriter((queryIndexFile)));

			//collect bps with the same set of files
			//start first entry, should always contain one or more files
			int start = Integer.MIN_VALUE;
			String priorIds = "";
			int stop = Integer.MIN_VALUE;
			
			for (int i=0; i< workingIndex.length; i++){
				if (workingIndex[i] == null) {
					//in a block? write out old if not the first
					if (start != Integer.MIN_VALUE){
						out.print(workingChr); out.print("\t");
						out.print(start); out.print("\t");
						out.print(i); out.print("\t");
						out.println(priorIds);

						//reset start 
						priorIds = "";
						start = Integer.MIN_VALUE;
					}
				}

				//nope hash present
				else {
					//save last good position
					stop = i;
					//different ids?
					String currentIds = fetchSortedIds(workingIndex[i]);
					if (priorIds.equals(currentIds) == false){
						//write out old if not the first
						if (start != Integer.MIN_VALUE){
							out.print(workingChr); out.print("\t");
							out.print(start); out.print("\t");
							out.print(stop); out.print("\t");
							out.println(priorIds);
						}
						//start new
						priorIds = currentIds;
						start = i;
					}
				}
			}
			//save last?
			if (start != Integer.MIN_VALUE){
				out.print(workingChr); out.print("\t");
				out.print(start); out.print("\t");
				out.print(stop+1); out.print("\t");
				out.println(priorIds);
			}
			//close it
			out.close();
		} catch (IOException e){
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem saving the index for "+workingChr+", aborting.\n");
		}
	}
	
	

	/**Returns the int id concat rep for all the files in the sorted TreeSet*/
	private String fetchSortedIds(TreeSet<File> files) {
		int[] ids = new int[files.size()];
		int counter = 0;
		for (File f: files) ids[counter++] = fileId.get(f);
		return Misc.intArrayToString(ids, ",");
	}

	private void parse(File f) {
		try {

			//create the master hash to use when no overlaps with other annotations
			TreeSet<File> singleton = new TreeSet<File>();
			singleton.add(f);

			//fetch a reader and iterator on the entire chr
			TabixReader reader = new TabixReader(f.toString());
			TabixReader.Iterator it = reader.query(workingChr);

			//what's the file type
			boolean vcf = f.getName().endsWith(".vcf.gz");
			int[] startStopSubtract = null;
			if (vcf == false) startStopSubtract = getSSS(f);

			String record = null;
			while ((record = it.next()) != null){

				//parse start and stop bp positions
				int[] startStop;
				if (vcf) startStop = parseStartStopBpCoorVcf(record);
				else startStop = parseStartStopBpCoor(record, startStopSubtract);
				if (startStop == null) {
					workingFailed++;
					continue;
				}
				workingPassed++;

				//check against the chrom length
				boolean warn = false;
				if (startStop[0] < 0) {
					startStop[0] = 0;
					warn = true;
				}
				if (startStop[1] > workingIndex.length) {
					startStop[1] = workingIndex.length;
					warn = true;
				}
				if (warn && verbose){
					System.err.println("\tWARNING: This record's covered bps were trimmed ("+startStop[0]+" - "+startStop[1]+
							") to match the chrom length ("+workingIndex.length+") -> "+record);
				}

				//add in references to source file over the covered bases, stop isn't covered.
				for (int j= startStop[0]; j< startStop[1]; j++){
					//anything already there?
					if (workingIndex[j] == null) workingIndex[j] = singleton;

					//does the current hash already contain the file due to overlapping regions?
					else if (workingIndex[j].contains(f) == false){
						//ok something there but it doesn't contain this file
						//attempt to find an existing hash with these files
						String key = fetchCompositeKey(workingIndex[j], f);
						TreeSet<File> ts = fileNames2TheirPointers.get(key);
						if (ts == null){
							ts = new TreeSet<File>();
							ts.addAll(workingIndex[j]);
							ts.add(f);
							fileNames2TheirPointers.put(key, ts);
						}
						workingIndex[j] = ts;
					}
					//yes, so don't do anything it's already there
				}
			}
			reader.close();
		}
		catch (IOException e){
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem parsing "+f+", aborting.\n");
		}
	}
		
		private int[] parseStartStopBpCoorVcf(String vcfRecord) {
			int[] startStop = null;
			try {
				String[] t = Misc.TAB.split(vcfRecord);
				startStop = fetchEffectedBps(t);
			} catch (NumberFormatException e){}
			return startStop;
		}

		private int[] parseStartStopBpCoor(String record, int[] startStopSub) throws NumberFormatException {
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
	public int[] fetchEffectedBps(String[] vcf) throws NumberFormatException{
		
		/*NOTE, any changes here, please update the web app's QueryRequest too*/
		
		//CHROM	POS	ID	REF	ALT	
		//  0    1   2   3   4  
		//put into interbase coordinates
		int iPos = Integer.parseInt(vcf[1]) - 1;
		String ref= vcf[3];
		String alt= vcf[4];

		//any commas/ multi alts?
		if (alt.contains(",") == false) return fetchEffectedBpsNoMultiAlt(iPos, ref, alt, vcf);

		//OK commas present, these need to be tested for max effect
		//There is complexity with multi alts, best to deconvolute and left justify! 
		String[] alts = Misc.COMMA.split(alt);
		int begin = Integer.MAX_VALUE;
		int end = -1;
		for (int i=0; i< alts.length; i++){
			int[] ss = fetchEffectedBpsNoMultiAlt(iPos, ref, alts[i], vcf);
			if (ss == null) return null;
			if(ss[0] < begin) begin = ss[0];
			if(ss[1]> end) end = ss[1];
		}

		return new int[]{begin, end};
	}

	public int[] fetchEffectedBpsNoMultiAlt(int iPos, String ref, String alt, String[] vcf) throws NumberFormatException{
		
		/*NOTE, any changes here, please update the web app's QueryRequest too*/
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
				if (verbose) System.err.println("\tWARNING: found a < or > containing alt, failed to parse END=number position, skipping -> "+Misc.stringArrayToString(vcf, "\t"));
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
				if (verbose) System.err.println("\tWARNING: Odd INS vcf record, the first base in the ref and alt must be the same, use vt to normalize your variants, skipping -> "+Misc.stringArrayToString(vcf, "\t"));
				return null;
			}

		}
		//del? return the bps that are deleted, AT->A, ATTCG->ACC
		else if (lenRef > lenAlt) {
			begin = iPos+1;
			end = iPos + lenRef;
			if (ref.charAt(0) != alt.charAt(0)) {
				if (verbose) System.err.println("\tWARNING: Odd DEL vcf record, the first base in the ref and alt must be the same, use vt to normalize your variants, skipping -> "+Misc.stringArrayToString(vcf, "\t"));
				return null;
			}
		}
		//odd, shouldn't hit this
		else {
			if (verbose) System.err.println("ERROR: Contact admin! Odd vcf record, can't parse effected bps for -> "+Misc.stringArrayToString(vcf, "\t"));
			return null;
		}

		return new int[]{begin, end};
	}

	private String fetchCompositeKey(TreeSet<File> set, File file) {
		StringBuilder sb = new StringBuilder();
		for (File f: set) sb.append(f.toString());
		sb.append(file.toString());
		return sb.toString();
	}
	
	public void bgzipAndTabixIndex() {
		File[] beds = IO.extractFiles(indexDir, ".qi.bed");
		for (File bed: beds){
			//compress with bgzip, this will replace any existing compressed file and delete the uncompressed
			String[] cmd = { bgzip.toString(), "-f", bed.toString()};
			String[] output = IO.executeViaProcessBuilder(cmd, false);
			File compBed = new File (bed+".gz");
			if (output.length != 0 || bed.exists() == true || compBed.exists() == false){
				Misc.printErrAndExit("\nERROR: Failed to bgzip compress "+bed+"\nMessage: "+Misc.stringArrayToString(output, "\n"));
			}

			//tabix
			cmd = new String[]{ tabix.toString(), "-f", "-p", "bed", compBed.toString() };
			output = IO.executeViaProcessBuilder(cmd, false);
			File indexBed = new File (compBed+".tbi");
			if (output.length != 0 || indexBed.exists() == false){
				Misc.printErrAndExit("\nERROR: Failed to tabix index "+compBed+"\nMessage: "+Misc.stringArrayToString(output, "\n"));
			}
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new QueryIndexer (args);
	}	

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		String useqVersion = IO.fetchUSeqVersion();
		File tabixBinDirectory = null;
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': chrNameLength = new File(args[++i]); break;
					case 'd': dataDir = new File(args[++i]); break;
					case 'i': indexDir = new File(args[++i]); break;
					case 'e': fileExtToIndex = Misc.COMMA.split(args[++i]); break;
					case 'q': verbose = false; break;
					case 't': tabixBinDirectory = new File(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//pull tabix and bgzip executables
		if (tabixBinDirectory == null) Misc.printExit("\nError: please point to the dir containing the tabix and bgzip HTSlib executibles (e.g. /Users/Clinton/BioApps/HTSlib/1.3/bin/ )\n");
		bgzip = new File (tabixBinDirectory, "bgzip");
		tabix = new File (tabixBinDirectory, "tabix");
		if (bgzip.canExecute() == false || tabix.canExecute() == false) Misc.printExit("\nCannot find or execute bgzip or tabix executables from "+bgzip+" "+tabix);
		
		if (chrNameLength == null) Misc.printErrAndExit("\nError: please provide a bed file of chromosome and their max lengths to index. e.g. X 0 155270560\n" );
		if (dataDir == null || dataDir.isDirectory() == false) Misc.printErrAndExit("\nERROR: please provide a directory containing gzipped and tabix indexed bed, vcf, maf.txt, and bedGraph files to index." );
		if (indexDir == null ) Misc.printErrAndExit("\nERROR: please provide a directory in which to write the master query index." );
		indexDir.mkdirs();
		
		parseFileExtensions();
		parseDataSources();
		parseChromLengthFile();
	}	
	
	private void parseChromLengthFile() {
		HashMap<String, RegionScoreText[]> chrLen = Bed.parseBedFile(chrNameLength, true, false);
		chrLengths = new HashMap<String, RegionScoreText[]>(); 
		//remove any chr and check there is only one
		for (String chr: chrLen.keySet()){
			RegionScoreText[] regions = chrLen.get(chr);
			if (regions.length !=1) Misc.printErrAndExit("Error: there can be only one bed region for each chromosome, see "+chr);
			if (chr.startsWith("chr")) {
				String t = chr.substring(3);
				if (chrLen.containsKey(t)) Misc.printErrAndExit("Error: looks like there is both the chrXXX and XXX chromosome versions in your chrLen file, please remove one. Note all chromosome names are trimmed of chr.");
				chrLengths.put(t, regions);
			}
			else chrLengths.put(chr, regions);
		}
	}

	/*Methods to modify with hard coded file extensions*/
	
	private void parseFileExtensions() {
		if (fileExtToIndex == null) fileExtToIndex = supportedFileExtensions;
		else {
			HashSet<String> ok = new HashSet<String>();
			for (String e : supportedFileExtensions) ok.add(e);
			for (String t: fileExtToIndex){
				if (ok.contains(t) == false) Misc.printErrAndExit("\nError: one or more of your requested data source extensions isn't supported, choose from these "+ok+"\n");
			}
		}
		//0 index start stop column info for "bed.gz", "bedGraph.gz", "maf.txt.gz"
		//the last is the number to subtract from the start to convert to interbase
		//note vcf.gz is it's own beast and handled separately
		extensionsStartStopSub.put("bed.gz", new int[]{1,2,0});
		extensionsStartStopSub.put("bedGraph.gz", new int[]{1,2,0});
		extensionsStartStopSub.put("maf.txt.gz", new int[]{5,6,1});
	}
	
	private int[] getSSS(File f) {
		String name = f.getName();
		if (name.endsWith(".bed.gz")) return extensionsStartStopSub.get("bed.gz");
		if (name.endsWith(".maf.txt.gz")) return extensionsStartStopSub.get("maf.txt.gz");
		if (name.endsWith(".bedGraph.gz")) return extensionsStartStopSub.get("bedGraph.gz");
		return null;
	}
	
	private void saveFileHeadersAndIds(){
		try {
			//make Patterns to scan
			Pattern[] headerStarts = new Pattern[stringHeaderPatternStarts.length];
			for (int i =0; i< stringHeaderPatternStarts.length; i++) headerStarts[i] = Pattern.compile(stringHeaderPatternStarts[i]); 
				
			//find all the parsed files
			HashSet<File> parsedFiles = new HashSet<File>();
			for (ArrayList<File> al: chrFiles.values()) parsedFiles.addAll(al);
			
			//save ids with truncated paths
			int toSkip = dataDir.getParentFile().toString().length()+1;
			
			HashMap<String, Integer> fileStringId = new HashMap<String, Integer>();
			HashMap<String, String> fileHeaders = new HashMap<String, String>();
			
			for (File f: parsedFiles){
				String trimmedName = f.toString().substring(toSkip);
				
				//save trmmed name : id
				Integer id = fileId.get(f);
				fileStringId.put(trimmedName, id);
				
				//save trimmed name: json header, have to use string rep since JSONObject is not serializable
				ArrayList<String> header = parseHeader(f, headerStarts);
				JSONObject jo = new JSONObject();
				jo.put("source", trimmedName);
				jo.put("header", header);
				fileHeaders.put(trimmedName, jo.toString(1));				
			}
			
			//save them
			File ids = new File(indexDir, "fileIds.obj");
			IO.saveObject(ids, fileStringId);
			File h = new File(indexDir, "fileHeaders.obj");
			IO.saveObject(h, fileHeaders);
			
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: parsing headers, aborting\n");
		}
	}
	
	private ArrayList<String> parseHeader(File f, Pattern[] headerLineStarts) throws IOException {
		ArrayList<String> header = new ArrayList<String>();
		BufferedReader in = IO.fetchBufferedReader(f);
		String line;
		while ((line = in.readLine()) != null){
			//skip blanks
			line = line.trim();
			if (line.length()==0) continue;
			//scan patterns
			boolean found = false;
			for (Pattern p: headerLineStarts) {
				Matcher m = p.matcher(line);
				if (m.matches()) {
					header.add(line);
					found = true;
					break;
				}
			}
			if (found == false) break;
		}
		in.close();
		return header;
	}
	
	private void parseDataSources(){
		System.out.println("Searching for tabix data sources...");

		File[][] toCombine = new File[fileExtToIndex.length][];
		for (int i=0; i< fileExtToIndex.length; i++) toCombine[i] = IO.fetchFilesRecursively(dataDir, fileExtToIndex[i]);

		dataFilesToParse = IO.collapseFileArray(toCombine);
		int numFiles = dataFilesToParse.length;

		dataFilesToParse = returnFilesWithTabix(dataFilesToParse);
		int numFilesWithIndex = dataFilesToParse.length;

		if (numFilesWithIndex == 0) Misc.printErrAndExit("\nERROR: failed to find any data sources with xxx.tbi indexes in your dataDir?!\n");
		System.out.println("\t"+numFiles+" Data sources ("+ Misc.stringArrayToString(fileExtToIndex, ", ")+")");
		System.out.println("\t"+numFilesWithIndex+" Data sources with xxx.tbi indexes");
	}
	
	public static File[] returnFilesWithTabix(File[] tabixFiles) {
		ArrayList<File> goodFiles = new ArrayList<File>();
		for (File tb: tabixFiles){
			File index = new File (tb.toString()+".tbi");
			if (index.exists()) goodFiles.add(tb);
		}
		//any files?
		if (goodFiles.size()==0) return null;
		File[] toReturn = new File[goodFiles.size()];
		goodFiles.toArray(toReturn);
		return toReturn;
	}
	
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Query Indexer: Nov 2016                           **\n" +
				"**************************************************************************************\n" +
				"Builds index files for Query by recursing through a data directory looking for vcf,\n"+
				"maf, bed, and bedGraph files that have been bgzip compressed and tabix indexed.\n"+
				"Interval trees are build containing regions that overlap with one or more data sources.\n"+
				"These are used by the Query REST service to rapidly identify which data files contain\n"+
				"records that overlap a users query. This app needs >30G RAM to run.\n"+
				

				"\nRequired Params:\n"+
				"-c A bed file of chromosomes and their lengths (e.g. 21 0 48129895) to use to \n"+
				"     building the intersection index. For multiple builds and species, remove any 'chr'\n"+
				"     prefixes (it is ignored during searching), and provide only one named chrom with\n"+
				"     the max observed length.\n"+
				"-d A data directory containing gzipped tabixed xxx.vcf.gz, xxx.bed.gz, xxx.bedGraph.gz,\n"+
				"     and xxx.maf.txt.gz files. Recurses through all sub directories. Be sure to\n"+
				"     normalize and decompose_blocksub all vcf records, see \n"+
				"     http://genome.sph.umich.edu/wiki/Vt\n"+
				"-t Full path directory containing the compiled bgzip and tabix executables. See\n" +
				"      https://github.com/samtools/htslib\n"+


				"\nExample: java -Xmx64G -jar pathToUSeq/Apps/QueryIndexer -c b37ChrLen.bed \n"+
				"     -t TabixDataFiles/ | gzip > results.json.gz  \n\n" +

				"**************************************************************************************\n");
	}

}
