
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;

import edu.utah.seq.data.*;
import edu.utah.seq.useq.data.RegionScoreText;
import util.gen.*;
import java.util.*;


import util.bio.annotation.*;

/**Parses a Soap version 1 alignment txt file into point chromPointData, split by chromosome and strand.
 * For each sequence a single hit is assigned to the center position of the read.  Files are saved using the bar format.
 * Final positions are in interbase coordinates (0 start, stop excluded).
 * @author david.nix@hci.utah.edu 
 **/
public class SoapV1Parser {
	//fields
	private File[] dataFiles;
	private File saveDirectory;
	private File workingFile;
	private String versionedGenome;
	private HashMap <String, ArrayList<Point>> chromPointData;
	private HashMap <String, ArrayList<RegionScoreText>> chromNamedStartStop;
	private Pattern tab = Pattern.compile("\\t");
	private ArrayList<File> tempPointDataDirectories = new ArrayList<File>();
	private ArrayList<File> tempBedDirectories = new ArrayList<File>();
	private int readLength;
	private int maxNumberBestMatches = 1;
	private int minimumReadLength = 17;
	private int numberPassingAlignments;
	private Histogram histogram;
	private boolean sumIdenticalPositions = false;
	private boolean makeHistogramOnAllData = true;

	//constructors
	public SoapV1Parser(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		System.out.println("\n"+maxNumberBestMatches+"\tMax number best matches");
		System.out.println(minimumReadLength+"\tMinimum read length");
		System.out.println("\nConverting...");

		//make histogram
		histogram = new Histogram(0,37,37);

		//for each file, parse and split by strand and chromosome and save to disk
		for (int i=0; i< dataFiles.length; i++){
			//set working objects and parse tag file text
			workingFile = dataFiles[i];
			System.out.println("\t"+workingFile);

			//split file to chromosome strand specific temp files
			if (parseWorkingFile() == false) Misc.printExit("\nError: failed to parse, aborting.\n");

			//write temp data to disk
			writePointDataToDisk();
			chromPointData = null;
			writeNamedStartStopToDisk();
			chromNamedStartStop = null;
		}

		//generate a HashMap of strandedChromosomes and PointData files to merge
		HashMap<String,ArrayList<File>> chrFiles = collectStrandedFiles(tempPointDataDirectories);

		//load, sort, make point chromPointData, and save
		System.out.print("Sorting and saving PointData...");
		makePointData(chrFiles);
		System.out.println();
		
		//write bed
		System.out.print("Sorting and saving bed files...");
		chrFiles = collectStrandedFiles(tempBedDirectories);
		makeBedFiles(chrFiles);
		System.out.println();
		
		//cleanup
		for (int i=0; i< tempPointDataDirectories.size(); i++) IO.deleteDirectory(tempPointDataDirectories.get(i));
		for (int i=0; i< tempBedDirectories.size(); i++) IO.deleteDirectory(tempBedDirectories.get(i));

		//stats
		System.out.println("Stats...");
		System.out.println("\t"+numberPassingAlignments+"\tAlignments passing filters");
		//print histogram
		System.out.println("\nRead length distribution:");
		histogram.printScaledHistogram();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	public static HashMap<String,ArrayList<File>> collectStrandedFiles(ArrayList<File> files){
		HashMap<String,ArrayList<File>> chrFiles = new HashMap<String,ArrayList<File>>();
		//for each temp directory get names and load into hash	
		for (int i=0; i< files.size(); i++){
			HashMap<String, File> x = IO.fetchNamesAndFiles(files.get(i));
			//add to merged
			Iterator<String> it = x.keySet().iterator();
			while (it.hasNext()){
				String chr = it.next();
				//fetch ArrayList
				ArrayList<File> al;
				if (chrFiles.containsKey(chr)) al = chrFiles.get(chr);
				else {
					al = new ArrayList<File>();
					chrFiles.put(chr, al);
				}
				//add File
				al.add(x.get(chr));
			}
		}
		return chrFiles;
	}

	/**Makes the PointData hash after sorting, writes to disk.*/
	public void makePointData(HashMap<String,ArrayList<File>> chrFiles){
		ComparatorPointPosition comp = new ComparatorPointPosition();
		File pointDataDir = new File (saveDirectory, "PointData");
		pointDataDir.mkdir();
		//for each stranded chromosome
		Iterator<String> it = chrFiles.keySet().iterator();
		while (it.hasNext()){
			String chromName = it.next();
			//parse strand and chromosome
			int len = chromName.length();
			String strand = chromName.substring(len-1);
			String chromosome = chromName.substring(0, len-1);
			//fetchFiles, load and merge
			ArrayList<File> al = chrFiles.get(chromName);
			ArrayList<Point> pointsAL = new ArrayList<Point>();
			for (int i=0; i< al.size(); i++){
				pointsAL.addAll((ArrayList<Point>)IO.fetchObject(al.get(i)));
			}
			//convert to Point[] and sort
			Point[] points = new Point[pointsAL.size()];
			pointsAL.toArray(points);
			Arrays.sort(points, comp);		
			//make notes
			HashMap <String,String> notes = new HashMap <String,String> ();
			notes.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
			notes.put(BarParser.SOURCE_TAG, IO.concatinateFileFullPathNames(dataFiles, ","));
			notes.put(BarParser.STRAND_TAG, strand);
			notes.put(BarParser.READ_LENGTH_TAG, readLength+"");
			notes.put(BarParser.DESCRIPTION_TAG, "Generated by running the SoapV1Parser on Soap alignment =file(s), the position is assigned to the middle of the read, interbase coordinates.");
			//make an Info object  public Info (String text, String versionedGenome, String chromosome, String strand, int readLength, HashMap<String,String> notes){
			Info info = new Info(chromName, versionedGenome, chromosome, strand, readLength, notes);
			//make pd
			PointData pd = Point.extractPositionScores(points);			
			pd.setInfo(info);
			//merge?
			if (sumIdenticalPositions){
				ArrayList<PointData> pdAL = new ArrayList<PointData>();
				pdAL.add(pd);
				pd = PointData.mergePointData(pdAL, true, true);
			}
			//write to file
			pd.writePointData(pointDataDir);
			//cleanup
			pd = null;
			al = null;
			points = null;
			System.out.print(".");
		}
	}

	/**Makes the composite bed file after sorting, writes to disk.*/
	public void makeBedFiles(HashMap<String,ArrayList<File>> chrFiles){
		try {
			File bedDir = new File (saveDirectory, "Bed_"+versionedGenome);
			bedDir.mkdir();
			//for each stranded chromosome
			Iterator<String> it = chrFiles.keySet().iterator();
			while (it.hasNext()){
				String chromName = it.next();
				//parse strand and chromosome
				int len = chromName.length();
				String strand = chromName.substring(len-1);
				String chromosome = chromName.substring(0, len-1);
				File bedFile = new File (bedDir,chromosome+"_"+strand+"_.bed");
				PrintWriter out = new PrintWriter (new FileWriter (bedFile));
				//fetchFiles, load and merge
				ArrayList<File> al = chrFiles.get(chromName);
				ArrayList<RegionScoreText> ssAL = new ArrayList<RegionScoreText>();
				for (int i=0; i< al.size(); i++){
					ssAL.addAll((ArrayList<RegionScoreText>)IO.fetchObject(al.get(i)));
				}
				//convert to RegionScoreText[] and sort
				RegionScoreText[] ss = new RegionScoreText[ssAL.size()];
				ssAL.toArray(ss);
				Arrays.sort(ss);		
				//write to bed file (chr, start, stop, text, score, strand)
				for (int i=0; i < ss.length; i++){
					out.println(chromosome+ "\t"+ ss[i].getStart()+ "\t"+ ss[i].getStop()+ "\t"+ ss[i].getText()+ "\t"+ ss[i].getScore()+ "\t"+ strand);
				}
				out.close();
				IO.zipAndDelete(bedFile);
				System.out.print(".");
			}
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}


	public void writePointDataToDisk(){
		//make tempDirectory to hold stranded chromosomes
		String rnd = Passwords.createRandowWord(8);
		File dir = new File (saveDirectory, "TempDirPoint_"+rnd);
		dir.mkdir();
		tempPointDataDirectories.add(dir);

		//iterate through hashmap saving serialized objects
		Iterator<String> it = chromPointData.keySet().iterator();
		while (it.hasNext()){
			String name = it.next();
			File file = new File (dir, name);
			IO.saveObject(file, chromPointData.get(name));
		}
	}

	public void writeNamedStartStopToDisk(){ 
		//make tempDirectory to hold stranded chromosomes
		String rnd = Passwords.createRandowWord(8);
		File dir = new File (saveDirectory, "TempDirBed_"+rnd);
		dir.mkdir();
		tempBedDirectories.add(dir);

		//iterate through hashmap saving serialized objects
		Iterator<String> it = chromNamedStartStop.keySet().iterator();
		while (it.hasNext()){
			String name = it.next();
			File file = new File (dir, name);
			IO.saveObject(file, chromNamedStartStop.get(name));
		}
	}


	public boolean parseWorkingFile(){
		try{
			//get reader
			BufferedReader in = IO.fetchBufferedReader(workingFile);
			String line;
			String[] tokens = null;
			chromPointData = new HashMap<String, ArrayList<Point>>();
			chromNamedStartStop = new HashMap <String, ArrayList<RegionScoreText>> ();
			int length =0;
			numberPassingAlignments = 0;

			/*read in lines
			id seq qual #bestHits nada lengthTrimmedRead strand chrom 1basedRefSeq typeOfMatch changes
			0   1   2        3     4            5           6     7         8          9          10
			 */
			while ((line = in.readLine()) !=null){
				tokens = tab.split(line);
				//length
				length = Integer.parseInt(tokens[5]);
				if (makeHistogramOnAllData) histogram.count(length);
				if (length < minimumReadLength) continue;
				//bestHits
				int numberBestHits = Integer.parseInt(tokens[3]);
				if (numberBestHits > maxNumberBestMatches) continue;
				if (makeHistogramOnAllData == false) histogram.count(length);
				//increment counter
				numberPassingAlignments++;
				//make chrStrand
				String chrStrand = tokens[7]+tokens[6];
				//save Point
				ArrayList<Point> al;
				if (chromPointData.containsKey(chrStrand)) al = chromPointData.get(chrStrand);
				else {
					al = new ArrayList<Point>();
					chromPointData.put(chrStrand, al);
				}
				//calc middle
				int start = Integer.parseInt(tokens[8]) -1;
				int middle = (int)Math.round((double)length/2) + start;
				//add Point
				al.add(new Point(middle, 1));
				//save for bed file
				String seq = tokens[1];
				if (tokens[6].equals("-")) seq = util.bio.seq.Seq.reverseComplementDNA(seq);
				ArrayList<RegionScoreText> nssAL;
				if (chromNamedStartStop.containsKey(chrStrand)) nssAL = chromNamedStartStop.get(chrStrand);
				else {
					nssAL = new ArrayList<RegionScoreText>();
					chromNamedStartStop.put(chrStrand, nssAL);
				}
				//start, stop, text, score
				nssAL.add(new RegionScoreText(start, start+length, 1f/(float)numberBestHits, seq));
			}
			//set max read length?
			if (length > readLength) readLength = length;
			in.close();
			return true;
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SoapV1Parser(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': forExtraction = new File(args[i+1]); i++; break;
					case 'v': versionedGenome = args[i+1]; i++; break;
					case 'r': saveDirectory = new File (args[i+1]); i++; break;
					case 'x': maxNumberBestMatches = Integer.parseInt(args[++i]); break;
					case 'm': minimumReadLength = Integer.parseInt(args[++i]); break;
					case 's': sumIdenticalPositions = true; break;
					case 'p': makeHistogramOnAllData = false; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//pull files
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".txt");
		tot[1] = IO.extractFiles(forExtraction,".txt.zip");
		tot[2] = IO.extractFiles(forExtraction,".txt.gz");

		dataFiles = IO.collapseFileArray(tot);
		if (dataFiles == null || dataFiles.length==0) dataFiles = IO.extractFiles(forExtraction);
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.txt(.zip/.gz) file(s)!\n");
		if (versionedGenome == null) Misc.printExit("\nPlease enter a genome version recognized by UCSC, see http://genome.ucsc.edu/FAQ/FAQreleases.\n");
		if (saveDirectory == null)  Misc.printExit("\nPlease provide a directory in which to save the parsed data.\n");
		else if (saveDirectory.exists() == false) saveDirectory.mkdir();
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                SoapV1Parser: Feb 2009                            **\n" +
				"**************************************************************************************\n" +
				"Splits and converts Soap version 1 alignment xxx.txt files into center position binary\n" +
				"PointData xxx.bar files. Interbase coordiantes (zero based, stop excluded).\n" +
				"These can be directly viewed in IGB.\n\n" +

				"-v Versioned Genome (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-f The full path directory/file text of your Soap xxx.txt(.zip or .gz) file(s).\n" +
				"-r Full path directory text for saving the results.\n"+
				"-x Maximum number of best matches, defaults to 1.\n"+
				"-m Miminum read length, defaults to 17.\n"+
				"-s Sum identical PointData positions. This should not be used for any downstream USeq\n" +
				"      applications, only for visualization.\n"+
				"-p Make read length histogram on reads that pass filters, defaults to all.\n"+


				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/SoapV1Parser -f /Soap/Run7/\n" +
				"     -v H_sapiens_Mar_2006 -x 5 -m 20\n\n" +

		"**************************************************************************************\n");

	}	

}
