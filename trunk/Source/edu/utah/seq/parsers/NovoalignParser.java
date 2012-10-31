
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;

import edu.utah.seq.data.*;
import edu.utah.seq.useq.data.Region;
import edu.utah.seq.useq.data.RegionScoreText;
import util.gen.*;
import java.util.*;

import util.bio.annotation.*;
import util.bio.seq.Seq;


/**Parses a Novoalign alignment txt file into point chromPointData, split by chromosome and strand.
 * For each sequence a single hit is assigned to the center position of the read.  Files are saved using the bar format.
 * Final positions are in interbase coordinates (0 start, stop excluded).
 * Recommend: novoalign -r0.015 -q5 -d database -f pbrFiles"_prb.txt.gz" | grep chr >> results
 * If a chrAdapter.fasta was included will remove any reads that have hits to it.
 * @author david.nix@hci.utah.edu 
 **/
public class NovoalignParser{
	//fields
	private File[] dataFiles;
	private File saveDirectory;
	private File bedDirectory;
	private File pointDataDirectory;
	private File tempDirectory;
	private File workingFile;
	private String versionedGenome;
	private HashMap <String, DataOutputStream> chromOut = new HashMap <String, DataOutputStream>();
	private Pattern tab = Pattern.compile("\\t");
	private Pattern spliceJunction = Pattern.compile("(.+)_(\\d+)_(\\d+)");
	private Pattern underscore = Pattern.compile("_");
	private Pattern colorspace = Pattern.compile("\\w[\\.\\d]+");
	private Pattern whiteSpace = Pattern.compile("\\s");
	private Pattern leadingNumber = Pattern.compile("^(\\d+).+");
	private HashMap <String, File> tempChrData = new HashMap <String, File>();
	private HashMap <String, File> tempChrSpliceJunctionData = new HashMap <String, File>();
	private ComparatorPointPosition comparator = new ComparatorPointPosition();
	private float probThreshold = 13;
	private float alignmentThreshold = 60;
	private int numberPassingAlignments = 0;
	private int numberAlignments = 0;
	private long numberAlignedBasePairs = 0;
	private int minimumBasePairQuality = 30;
	private String chromosomePrefix = ">chr";
	private boolean skipUnderScoredChromosomes = false;
	private String adapterName = "chrAdapt";
	private File geneRegionFile;
	private HashMap<String,Region[]> geneRegions;
	private boolean ignoreStrandForSpliceJunctions = false;
	private boolean warningCalled = false;
	private boolean printBadLines = true;
	private boolean justPrintStats = false;
	private boolean verbose = true;

	//for filtering novoalign files for repeats near each other
	private boolean filterForPairedRepeats = false;
	private int maxRepeats = 10;
	private int maxDistance = 1000;
	private int chromIndex = -1;
	private int positionIndex = -1;
	private int oriIndex = -1;
	private int errorIndex = -1;
	
	//histogram
	private Histogram histogram = new Histogram(15,110, 110-15);
	private double[] bpErrors = new double[150];
	private double[] bpCounts = new double[150];


	//constructors
	/**For incorporation into RNA-Seq and ChIP-Seq apps.*/
	public NovoalignParser (File saveDirectory, File[] dataFiles, float minimumMappingQuality, float maximumAlignmentScore, String versionedGenome, File geneRegionFile){
		this.dataFiles = dataFiles;
		this.saveDirectory = saveDirectory;
		this.probThreshold = minimumMappingQuality;
		this.alignmentThreshold = maximumAlignmentScore; 
		this.versionedGenome = versionedGenome;
		this.geneRegionFile = geneRegionFile;
		verbose = false;

		tempDirectory = new File (saveDirectory,"TempFilesDelme");
		tempDirectory.mkdir();
		bedDirectory = new File (saveDirectory, "Bed_"+versionedGenome);
		bedDirectory.mkdir();
		pointDataDirectory = new File (saveDirectory, "PointData");
		pointDataDirectory.mkdir();

		//flattened gene regions for scaling? Not stranded
		if (geneRegionFile != null){
			geneRegions = Region.parseStartStops(geneRegionFile, 0, 0, 0);
		}

		//for each file, parse and save to disk	
		for (int i=0; i< dataFiles.length; i++){
			//set working objects and parse tag file text
			workingFile = dataFiles[i];
			//split file to chromosome strand specific temp files
			boolean  parsed = parseWorkingFile();
			if (parsed == false) Misc.printExit("\n\tError: failed to parse, aborting.\n");
		}
		//close the writers
		closeWriters();
		//load, sort, make point chromPointData, and save
		makePointDataAndBedFiles();
		if (skipUnderScoredChromosomes == false) saveSpliceJunctionBedFile();
		//cleanup
		IO.deleteDirectory(tempDirectory);
	}


	/**Stand alone.*/
	public NovoalignParser(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		System.out.println("\n"+probThreshold+"\tPosterior probability threshold");
		System.out.println(alignmentThreshold+"\tAlignment score threshold");
		System.out.println("\nParsing and filtering...");

		//for each file, parse and save to disk	
		for (int i=0; i< dataFiles.length; i++){
			//set working objects and parse tag file text
			workingFile = dataFiles[i];
			System.out.print("\t"+workingFile);

			//split file to chromosome strand specific temp files
			boolean  parsed;
			if (filterForPairedRepeats) {
				parsed = parseWorkingFileForRepeats(); 
				System.out.println("\t"+numberAlignments+"\tAlignments");
				System.out.println("\t"+numberPassingAlignments+"\tAlignments passing filters");
			}
			parsed = parseWorkingFile(); 
			if (parsed == false) Misc.printExit("\n\tError: failed to parse, aborting.\n");
			System.out.println();
		}

		//stop?
		if (filterForPairedRepeats) System.exit(0);

		//close the writers
		closeWriters();

		//stats
		System.out.println("Stats:");
		System.out.println("\t"+numberAlignments+"\tAlignments");
		System.out.println("\t"+numberPassingAlignments+"\tAlignments passing filters");
		System.out.println("\t"+numberAlignedBasePairs+"\tAligned Q"+ minimumBasePairQuality+ " bps passing filters");
		System.out.println("\nRead length histogram:\nLength\tCounts\tScaledStars");
		histogram.setTrimLabelsToSingleInteger(true);
		histogram.printScaledHistogram();
		
		//print error rates
		System.out.println("\nPer base error rates:\nPosition\tErrors\tObservations\tErrorRate");
		for (int i=0; i< bpCounts.length; i++){
			if (bpCounts[i] == 0) break;
			double rate = bpErrors[i]/bpCounts[i];
			System.out.println(i+"\t" + (int)bpErrors[i] + "\t" +(int)bpCounts[i]+ "\t"+Num.formatNumber(rate, 3));
		}

		//load, sort, make point chromPointData, and save?
		if (justPrintStats == false){
			System.out.print("\nLoading, sorting and saving PointData and bed files for ");
			makePointDataAndBedFiles();
			if (skipUnderScoredChromosomes == false) saveSpliceJunctionBedFile();
			System.out.println();
		}
		//cleanup
		IO.deleteDirectory(tempDirectory);

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	/**Closes writers.*/
	public void closeWriters(){
		try{
			Iterator<DataOutputStream> it = chromOut.values().iterator();
			while (it.hasNext()) it.next().close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void makePointDataAndBedFiles(){
		//for each composite data file in non splice junction data
		Iterator<String> it = tempChrData.keySet().iterator();
		while (it.hasNext()){
			String chromStrand = it.next();
			System.out.print(chromStrand+" ");
			File chromDataFile = tempChrData.get(chromStrand);
			//parse strand and chrom
			String strand = chromStrand.substring(chromStrand.length()-1);
			String chrom = chromStrand.substring(0, chromStrand.length()-1);
			//get StartStopScores and sort
			RegionScoreText[] sss = RegionScoreText.oldLoadBinary_DEPRECIATED(chromDataFile, true);
			//make files
			saveBedFile(sss,chrom,strand);
			savePointData(sss,chrom,strand);
		}
		System.out.println();
	}

	/**Returns Point[] where splice read assigned to last bp in first exon.*/
	public Point[] fetchSpliceJunctionPoints(String chrom, String strand){
		File chromDataFile = tempChrSpliceJunctionData.get(chrom+strand);
		if (chromDataFile == null) return null;
		RegionScoreText[] sss = RegionScoreText.oldLoadBinary_DEPRECIATED(chromDataFile, false);
		//make Point[]
		Point[] points = new Point[sss.length];
		for (int i=0; i<sss.length; i++) points[i] = new Point(sss[i].getStart()-1, sss[i].getScore());
		return points;
	}

	/**Makes the PointData, sorts and saves.*/
	public void savePointData(RegionScoreText[] sss, String chrom, String strand){
		//make Point[]
		Point[] points = new Point[sss.length];
		long totalReadLengthSize = 0;
		for (int i=0; i<sss.length; i++){
			//add to length
			totalReadLengthSize += sss[i].getLength();
			//make point
			points[i] = new Point(sss[i].getMiddle(), sss[i].getScore());
		}
		//average read length
		int averageReadLength = (int)(totalReadLengthSize/points.length);

		//any splice junction points to add?
		Point[] spliceJunctionPoints = fetchSpliceJunctionPoints(chrom, strand);
		if (spliceJunctionPoints !=null){
			Point[] merge = new Point[points.length+ spliceJunctionPoints.length];
			System.arraycopy(points, 0, merge, 0, points.length);
			System.arraycopy(spliceJunctionPoints, 0, merge, points.length, spliceJunctionPoints.length);
			points = merge;
		}

		Arrays.sort(points, comparator);

		//make notes
		HashMap <String,String> notes = new HashMap <String,String> ();
		notes.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
		notes.put(BarParser.SOURCE_TAG, IO.concatinateFileFullPathNames(dataFiles, ","));
		notes.put(BarParser.STRAND_TAG, strand);
		notes.put(BarParser.READ_LENGTH_TAG, +averageReadLength+"");
		notes.put(BarParser.UNIT_TAG, "Probability of originating from the given location.");
		notes.put(BarParser.DESCRIPTION_TAG, "Generated by running the NovoalignParser on Novoalign alignment file(s), the position is assigned to the middle of the read, interbase coordinates. The posterior probability thresholds was set to "+probThreshold);
		//make an Info object  public Info (String text, String versionedGenome, String chromosome, String strand, int readLength, HashMap<String,String> notes){
		Info info = new Info(chrom+strand, versionedGenome, chrom, strand, averageReadLength, notes);
		info.setNumberObservations(points.length);
		//make pd
		PointData pd = Point.extractPositionScores(points);			
		pd.setInfo(info);
		//write to file
		pd.writePointData(pointDataDirectory);
		//cleanup
		points = null;
	}

	/**Makes and saves the data to a bed file.*/
	public void saveBedFile(RegionScoreText[] sss, String chrom, String strand){
		try {
			File bedFile = new File (bedDirectory, chrom+strand+".bed");
			PrintWriter out = new PrintWriter (new FileWriter (bedFile));
			for (int i=0; i< sss.length; i++){
				out.println(chrom+ "\t"+ sss[i].getStart()+ "\t"+ sss[i].getStop()+ "\t"+ sss[i].getText()+ "\t"+ sss[i].getScore()+ "\t"+ strand);
			}
			out.close();
			IO.zipAndDelete(bedFile);
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	/**Makes and saves the data to a gff file.*/
	public void saveSpliceJunctionGFFFile(){
		try {
			//make file
			File file = new File (saveDirectory, "spliceJunctions.gff");
			PrintWriter out = new PrintWriter (new FileWriter (file));
			//for each composite data file in splice junction data
			Iterator<String> it = tempChrSpliceJunctionData.keySet().iterator();
			while (it.hasNext()){
				String chromStrand = it.next();
				File chromDataFile = tempChrSpliceJunctionData.get(chromStrand);
				//parse strand and chrom
				String strand = chromStrand.substring(chromStrand.length()-1);
				String chrom = chromStrand.substring(0, chromStrand.length()-1);
				//get StartStopScores and sort
				RegionScoreText[] sss = RegionScoreText.oldLoadBinary_DEPRECIATED(chromDataFile, true);
				//count number with same start stop
				LinkedHashMap<String, Integer> coorNum = new LinkedHashMap<String, Integer>();
				for (int i=0; i< sss.length; i++){
					String coor = sss[i].getStart()+"_"+sss[i].getStop();
					if (coorNum.containsKey(coor)) {
						Integer num = coorNum.get(coor);
						num = new Integer(num.intValue() + 1);
						coorNum.put(coor, num);
					}
					else coorNum.put(coor, new Integer(1));
				}
				//print
				Iterator<String> it2 = coorNum.keySet().iterator();
				while (it2.hasNext()){
					String startStop = it2.next();
					int counts = coorNum.get(startStop).intValue();
					String[] ss = underscore.split(startStop);
					//note gff uses 1bp coord, stop included so add 1bp to start, set start=stop
					int exonAEnd = Integer.parseInt(ss[0]);
					int exonAStart = exonAEnd;
					int exonBStart = Integer.parseInt(ss[1])+1;
					int exonBEnd = exonBStart;
					//print line
					String left = chrom+"\tUSeq\tSpliceJunction\t"; 
					String right = "\t"+counts+"\t"+strand+"\t.\t"+counts+"_"+startStop;
					out.println(left +exonAStart+"\t"+exonAEnd+right);
					out.println(left +exonBStart+"\t"+exonBEnd+right);
				}
			}
			out.close();
			IO.zipAndDelete(file);
			if (verbose) System.out.print(".");
		} catch (Exception e){
			e.printStackTrace();
		}
	}


	/**Makes and saves the data to a bed file.*/
	public void saveSpliceJunctionBedFile(){
		if (verbose) System.out.println();
		try {
			//make file
			File file = new File (saveDirectory, "spliceJunctions.bed");
			PrintWriter out = new PrintWriter (new FileWriter (file));
			out.println("track text=SpliceJunctions description=\"SpliceJunctions\" useScore=1 color=255,255,0");
			//for each composite data file in splice junction data
			Iterator<String> it = tempChrSpliceJunctionData.keySet().iterator();
			HashSet<String> testedChroms = new HashSet<String>();
			while (it.hasNext()){
				String chromStrand = it.next();
				File chromDataFile = tempChrSpliceJunctionData.get(chromStrand);				
				//parse strand and chrom
				String strand = chromStrand.substring(chromStrand.length()-1);
				String chrom = chromStrand.substring(0, chromStrand.length()-1);
				//already tested?
				if (ignoreStrandForSpliceJunctions && testedChroms.contains(chrom)) continue;
				testedChroms.add(chrom);
				//get StartStopScores and sort
				RegionScoreText[] sss = RegionScoreText.oldLoadBinary_DEPRECIATED(chromDataFile, false);
				//add in other strand?
				if (ignoreStrandForSpliceJunctions){
					String oppStrand = "-";
					if (strand.equals("-"))oppStrand = "+";
					File oppChromData = tempChrSpliceJunctionData.get(chrom+oppStrand);
					if (oppChromData != null){
						RegionScoreText[] oppSSS = RegionScoreText.oldLoadBinary_DEPRECIATED(oppChromData,false);
						RegionScoreText[] merge = new RegionScoreText[sss.length+ oppSSS.length];
						System.arraycopy(sss, 0, merge, 0, sss.length);
						System.arraycopy(oppSSS, 0, merge, sss.length, oppSSS.length);
						sss = merge;
					}
					strand = ".";
				}
				Arrays.sort(sss);
				//count number with same start stop
				LinkedHashMap<String, Integer> coorNum = new LinkedHashMap<String, Integer>();
				for (int i=0; i< sss.length; i++){
					String coor = sss[i].getStart()+"_"+sss[i].getStop();
					if (coorNum.containsKey(coor)) {
						Integer num = coorNum.get(coor);
						num = new Integer(num.intValue() + 1);
						coorNum.put(coor, num);
					}
					else coorNum.put(coor, new Integer(1));
				}
				//print
				Iterator<String> it2 = coorNum.keySet().iterator();
				ArrayList<Bed> chrmBedLines = new ArrayList<Bed>();
				while (it2.hasNext()){
					String startStop = it2.next();
					int counts = coorNum.get(startStop).intValue();
					String[] ss = underscore.split(startStop);
					//note gff uses 1bp coord, stop included so add 1bp to start, set start=stop
					int start = Integer.parseInt(ss[0])-1;
					int end = Integer.parseInt(ss[1])+1;
					//save line
					//out.println(chrom+ "\t"+start+"\t"+stop+"\t"+startStop+"\t"+counts+"\t"+strand);
					chrmBedLines.add(new Bed(chrom, start, end, counts+":"+startStop, (counts* -1), strand.charAt(0)));
				}
				//scale?
				Bed[] bedLines;
				if (geneRegions!= null) bedLines = scaleBedLines(chrmBedLines);
				else {
					bedLines =  new Bed[chrmBedLines.size()];
					chrmBedLines.toArray(bedLines);
					//invert scores
					for (int z=0; z< bedLines.length; z++) {
						bedLines[z].setScore(bedLines[z].getScore() * -1);
					}
				}
				//print
				if (bedLines != null) for (int i=0; i< bedLines.length; i++) out.println(bedLines[i]);
			}
			out.close();
			IO.zipAndDelete(file);
			if (verbose) System.out.print(".");
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public Bed[] scaleBedLines(ArrayList<Bed> chrmBedLines){
		//get chrom
		String chrom = chrmBedLines.get(0).getChromosome();
		//get gene regions, not stranded
		Region[] chromGeneRegions = geneRegions.get(chrom);
		if (chromGeneRegions == null){
			System.out.println("\nWarning! No gene regions found for "+chrom);
			return null;
		}
		//covert AL to Bed[]
		Bed[] bed = new Bed[chrmBedLines.size()];
		chrmBedLines.toArray(bed);
		//for each gene region find intersecting Beds
		ArrayList<Bed> intersectingBedRegions = new ArrayList<Bed>();
		for (int i=0; i< chromGeneRegions.length; i++){
			intersectingBedRegions.clear();
			int start = chromGeneRegions[i].getStart();
			int stop = chromGeneRegions[i].getStop() -1;
			double maxScore = -1;
			for (int j=0; j< bed.length; j++){
				if (bed[j].intersects(start, stop)) {
					//check and set score
					double score = bed[j].getScore();
					if (score >0); //System.out.println("Error, your gene regions file contains overlapping regions, run the MergeRegions app on this file and restart.");
					else score = -1 * score;
					bed[j].setScore(score);
					if (maxScore < score) maxScore = score;
					intersectingBedRegions.add(bed[j]);
				}
			}
			//get scores
			int number = intersectingBedRegions.size();
			if (number == 0) continue;
			if (number == 1){
				intersectingBedRegions.get(0).setScore(1000);
			}
			else {
				double[] scores = new double[intersectingBedRegions.size()];
				for (int k=0; k<intersectingBedRegions.size(); k++) scores[k] = intersectingBedRegions.get(k).getScore();
				//convert
				int[] scaledScores = Num.scale100To1000(scores);
				//set
				for (int k=0; k<intersectingBedRegions.size(); k++) intersectingBedRegions.get(k).setScore(scaledScores[k]);
			}
		}
		//convert remainder to negatives
		boolean foundNegs = false;
		for (int i=0; i< bed.length; i++){
			if (bed[i].getScore()<0){
				foundNegs = true;
				bed[i].setScore(1000);
			}
		}
		if (foundNegs && warningCalled == false) {
			System.out.println("\nWarning! Some of your splice junctions did not fall within one of your gene regions?!");
			warningCalled = true;
		}
		return bed;
	}

	public boolean parseWorkingFileForRepeats(){
		try{
			String line;
			String[] tokens = null;
			int counter = 0;
			/*read in lines to find out if the miRNA scoring mode was turned on
			text type read qual matchClass alignScore prob (miRNAScore) chrom position ori   .   .   .   error
			 0    1     2   3     4            5        6                 7       8     9    10  11  12  13
			                                                     7        8       9     10   11  12  13  14
			 */
			BufferedReader in = IO.fetchBufferedReader(workingFile);
			while ((line = in.readLine()) !=null){
				//does it contain 'chr'?
				if (line.startsWith("#") == false){
					if (line.contains(chromosomePrefix)){
						tokens = tab.split(line);
						if (tokens[7].startsWith(chromosomePrefix)){
							chromIndex = 7;
							positionIndex = 8;
							oriIndex = 9;
							errorIndex = 13;
							break;
						}
						else if (tokens[8].startsWith(chromosomePrefix)){
							chromIndex = 8;
							positionIndex = 9;
							oriIndex = 10;
							errorIndex = 14;
							break;
						}
						else Misc.printExit("\nProblem identifying chromosome column? This should start with '"+chromosomePrefix+"' and reside in column 7 or 8. See -> "+line);
					}
					if (counter++ == 10000) {
						System.err.println("\nProblem identifing chromosome column? No '"+chromosomePrefix+"' found in 1st 10000 lines? Assuming chr column is 7th.\n");
						chromIndex = 7;
						positionIndex = 8;
						oriIndex = 9;
						errorIndex = 13;
						break;
					}
				}
			}


			//reset BufferedReader
			in = IO.fetchBufferedReader(workingFile);
			numberAlignments = 0;
			numberPassingAlignments = 0;
			counter = 0;
			String name = Misc.removeExtension(workingFile.getName())+"_repFilt"+maxRepeats+".txt";
			String workingReadName = "";
			ArrayList<String[]> lines = new ArrayList<String[]>();
			PrintWriter out = new PrintWriter (new FileWriter( new File(saveDirectory,name)));
			while ((line = in.readLine()) !=null){
				try {
					if (++counter == 25000){
						if (verbose) System.out.print(".");
						counter = 0;
					}
					if (line.startsWith("#")) {
						out.println(line);
						continue;
					}
					tokens = tab.split(line);
					numberAlignments++;
					//unique?
					if (tokens[4].equals("U")){
						out.println(line);
						numberPassingAlignments++;
						continue;
					}
					//if not repeat then skip
					if (tokens[4].equals("R") == false)	continue;

					//check read name
					String testReadName = tokens[0];
					//if same add to ArrayList
					if (testReadName.equals(workingReadName)) lines.add(tokens);
					//not same so process
					else {
						//print if less than max otherwise delete
						processLines(lines, out);
						lines.clear();
						lines.add(tokens);
						workingReadName = testReadName;
					}
				} catch (NumberFormatException ne){
					System.err.println("\nBad line skipping -> "+line);
				}
			}
			//print last?
			processLines(lines, out);
			in.close();
			out.close();
			return true;
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
	}


	private void processLines(ArrayList<String[]> lines, PrintWriter out) throws NumberFormatException{
		int numLines = lines.size();
		//meets limit?
		if (numLines <= maxRepeats && numLines!=0){
			//all same chromosome and within maxDistance
			String[] first = lines.get(0);
			String chrom = first[chromIndex];
			int position = Integer.parseInt(first[positionIndex]);
			for (int i=1; i< numLines; i++) {
				String[] test = lines.get(i);
				//chrom
				if (chrom.equals(test[chromIndex]) == false) return;
				//dist
				int testPos = Integer.parseInt(test[positionIndex]);
				if (Math.abs(position-testPos) > maxDistance) return;
			}
			//print
			for (int i=0; i< numLines; i++) {
				String[] test = lines.get(i);
				out.print(test[0]);
				for (int x =1; x< test.length; x++){
					out.print("\t");
					out.print(test[x]);
				}
				out.println();
				numberPassingAlignments++;
			}
		}
	}

	public boolean parseWorkingFile(){
		try{
			//get reader
			BufferedReader in = IO.fetchBufferedReader(workingFile);
			String line;
			String[] tokens = null;

			//find out if the miRNA scoring mode was turned on
			int counter = 0;
			int chromIndex = -1;
			int positionIndex = -1;
			int oriIndex = -1;
			while ((line = in.readLine()) !=null){
				//does it contain 'chr'?
				if (line.startsWith("#") == false){
					if (line.contains(chromosomePrefix)){
						tokens = tab.split(line);
						//is this a colorspace file?
						if (colorspace.matcher(tokens[2]).matches()) {
							return parseWorkingColorSpaceFile();
						}
						else if (tokens[7].startsWith(chromosomePrefix)){
							chromIndex = 7;
							positionIndex = 8;
							oriIndex = 9;
							errorIndex = 13;
							break;
						}
						else if (tokens[8].startsWith(chromosomePrefix)){
							chromIndex = 8;
							positionIndex = 9;
							oriIndex = 10;
							errorIndex = 14;
							break;
						}
						else Misc.printExit("\nProblem identifying chromosome column? This should start with '"+chromosomePrefix+"' and reside in column 7 or 8. See -> "+line);
					}
					if (counter++ == 10000) {
						System.err.println("\nProblem identifing chromosome column? No '"+chromosomePrefix+"' found in 1st 10000 lines? Assuming chr column is 7th.\n");
						chromIndex = 7;
						positionIndex = 8;
						oriIndex = 9;
						errorIndex = 13;
						break;
					}
				}
			}

			/*read in lines
			text type read qual matchClass alignScore prob (miRNAScore) chrom position ori
			 0    1     2   3     4            5        6                 7       8     9  
			                                                     7        8       9     10
			 */
			//reset BufferedReader
			in = IO.fetchBufferedReader(workingFile);
			counter = 0;
			while ((line = in.readLine()) !=null){
				try {
					if (++counter == 25000){
						if (verbose) System.out.print(".");
						counter = 0;
					}
					if (line.startsWith("#")) continue;
					tokens = tab.split(line);
					numberAlignments++;

					if (tokens.length < 13 || tokens.length > 18 || tokens[chromIndex].contains(adapterName)) {
						if (printBadLines) System.err.println("Bad line, tokens, skipping -> "+line);
						continue;
					}

					//check probability threshold
					float probScore = Float.parseFloat(tokens[6]);
					if (probScore < probThreshold) continue;

					//check alignment threshold
					float alignmentScore = Float.parseFloat(tokens[5]);
					if (alignmentScore > alignmentThreshold) continue;

					//check that the seq and quality score string are of the same size!  Novoaligner stinks! 3hrs to find their bug
					int seqLength = tokens[2].length();
					if (seqLength != tokens[3].length() && tokens[3].length()!=0) {
						if (printBadLines) System.err.println("\nBad line, read/ qual length skipping -> "+line);
						continue;
					}

					//increment seq length histogram
					histogram.count(seqLength);
					
					//increment bp scored 
					for (int i=0; i<seqLength; i++) bpCounts[i]++;
					
					//process errors?
					if (tokens.length > errorIndex && tokens[errorIndex].length() !=0){
						String[] errors = whiteSpace.split(tokens[errorIndex]);
						for (int i=0; i< errors.length; i++){
							Matcher mat = leadingNumber.matcher(errors[i]);
							if (mat.matches() == false) {
								//throw new Exception ("Failed to match error pattern?! "+tokens[errorIndex]);
								//bisulfite notation should be skipped
								break;
							}
							int num = Integer.parseInt(mat.group(1)) -1;
							bpErrors[num]++;
						}
					}
					

					//parse chromosome
					String chromosome = tokens[chromIndex].substring(1);

					//increment counter
					numberPassingAlignments++;

					//count bases that pass minimum quality threshold
					int[] qscores = Seq.convertScores(tokens[3]);
					for (int i=0; i< qscores.length; i++) if (qscores[i]>= minimumBasePairQuality) numberAlignedBasePairs++;

					//continue?
					if (justPrintStats == false){
						//parse strand
						String strand = "+";
						if (tokens[oriIndex].equals("R")) strand = "-";
						int size = tokens[2].length();

						//calc start stop and chrStrand
						int start;
						int stop;
						String chrStrand;
						boolean spliceJunctionPresent;

						Matcher mat = spliceJunction.matcher(chromosome);
						if (mat.matches()){
							if (skipUnderScoredChromosomes) continue;
							spliceJunctionPresent = true;
							chromosome = mat.group(1);
							start = Integer.parseInt(mat.group(2));
							stop = Integer.parseInt(mat.group(3));
							chrStrand = chromosome+strand+"SpliceJunction";
						}
						else {
							spliceJunctionPresent = false;
							start = Integer.parseInt(tokens[positionIndex]) -1;
							stop = start + size;
							chrStrand = chromosome+strand;
						}

						float convertedScore = (float)(1- Num.antiNeg10log10(probScore));

						//get PrintWriter
						DataOutputStream dos;
						if (chromOut.containsKey(chrStrand)) dos = chromOut.get(chrStrand);
						else {
							//make and set file
							File f = new File(tempDirectory, chrStrand);
							if (spliceJunctionPresent) tempChrSpliceJunctionData.put(chromosome+strand, f);
							else tempChrData.put(chromosome+strand, f);
							dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(f)));
							chromOut.put(chrStrand, dos);
						}
						//write out
						dos.writeInt(start);
						dos.writeInt(stop);
						dos.writeFloat(convertedScore);
						dos.writeInt(size+ tokens[3].length()+1);
						//write concatinate of the sequence ; quality scores
						dos.writeBytes(tokens[2]+";"+tokens[3]);
					}
				} catch (NumberFormatException ne){
					System.err.println("\nBad line skipping -> "+line);
				}
			}
			in.close();
			return true;
		} catch (Exception e){
			System.err.println("\nError parsing Novoalign file or writing split binary chromosome files.\nToo many open files? Too many chromosomes? " +
			"If so then login as root and set the default higher using the ulimit command (e.g. ulimit -n 10000)\n");
			e.printStackTrace();
			return false;
		}
	}

	public boolean parseWorkingColorSpaceFile(){
		try{
			//get reader
			BufferedReader in = IO.fetchBufferedReader(workingFile);
			String line;
			String[] tokens = null;

			//find out if the miRNA scoring mode was turned on
			int counter = 0;
			int chromIndex = 9;
			int positionIndex = 10;
			int oriIndex = 11;
			int probabilityIndex =8;
			int alignmentScoreIndex =7;
			int readSequenceIndex =4;
			int readQualityIndex =5;


			/*read in lines
			text type csRead csQuality read qual matchClass alignScore prob chrom position ori pairChrom pairPosition pairStrand mismatches
			 0    1     2       3        4   5       6           7       8     9    10      11    12          13          14        15
			 */
			//reset BufferedReader
			in = IO.fetchBufferedReader(workingFile);
			counter = 0;
			while ((line = in.readLine()) !=null){
				try {
					if (++counter == 25000){
						if (verbose) System.out.print(".");
						counter = 0;
					}
					if (line.startsWith("#")) continue;
					tokens = tab.split(line);
					numberAlignments++;

					if (tokens.length < 15 || tokens.length > 18 || tokens[chromIndex].contains(adapterName)) {
						if (printBadLines) System.err.println("Bad line, tokens, skipping -> "+line);
						continue;
					}

					//check probability threshold
					float probScore = Float.parseFloat(tokens[probabilityIndex]);
					if (probScore < probThreshold) continue;

					//check alignment threshold
					float alignmentScore = Float.parseFloat(tokens[alignmentScoreIndex]);
					if (alignmentScore > alignmentThreshold) continue;

					//check that the seq and quality score string are the same size!
					int size = tokens[readSequenceIndex].length();
					if (size != tokens[readQualityIndex].length()) {
						if (printBadLines) System.err.println("\nBad line, read/ qual length, skipping -> "+line);
						continue;
					}

					//increment histogram
					histogram.count(size);

					//parse chromosome
					String chromosome = tokens[chromIndex].substring(1);

					//increment counter
					numberPassingAlignments++;

					//parse strand
					String strand = "+";
					if (tokens[oriIndex].equals("R")) strand = "-";


					//calc start stop and chrStrand
					int start;
					int stop;
					String chrStrand;
					boolean spliceJunctionPresent;

					Matcher mat = spliceJunction.matcher(chromosome);
					if (mat.matches()){
						if (skipUnderScoredChromosomes) continue;
						spliceJunctionPresent = true;
						chromosome = mat.group(1);
						start = Integer.parseInt(mat.group(2));
						stop = Integer.parseInt(mat.group(3));
						chrStrand = chromosome+strand+"SpliceJunction";
					}
					else {
						spliceJunctionPresent = false;
						start = Integer.parseInt(tokens[positionIndex]) -1;
						stop = start + size;
						chrStrand = chromosome+strand;
					}

					float convertedScore = (float)(1- Num.antiNeg10log10(probScore));

					//get PrintWriter
					DataOutputStream dos;
					if (chromOut.containsKey(chrStrand)) dos = chromOut.get(chrStrand);
					else {
						//make and set file
						File f = new File(tempDirectory, chrStrand);
						if (spliceJunctionPresent) tempChrSpliceJunctionData.put(chromosome+strand, f);
						else tempChrData.put(chromosome+strand, f);
						dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(f)));
						chromOut.put(chrStrand, dos);
					}

					//save data
					dos.writeInt(start);
					dos.writeInt(stop);
					dos.writeFloat(convertedScore);
					dos.writeInt(size+ tokens[3].length()+1);
					//write concatinate of the sequence ; quality scores
					dos.writeBytes(tokens[readSequenceIndex]+";"+tokens[readQualityIndex]);
				} catch (NumberFormatException ne){
					System.err.println("\nBad line skipping -> "+line);
				}
			}
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
		new NovoalignParser(args);
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
					case 'p': probThreshold = Float.parseFloat(args[++i]); break;
					case 'q': alignmentThreshold = Float.parseFloat(args[++i]); break;
					case 'c': chromosomePrefix = args[++i]; break;
					case 'n': adapterName = args[++i]; break;
					case 'i': ignoreStrandForSpliceJunctions = true; break;
					case 'g': geneRegionFile = new File(args[++i]); break;
					case 'h': printBadLines = true; break;
					case 's': justPrintStats = true; break;
					case 'z': filterForPairedRepeats = true; break;
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
		tempDirectory = new File (saveDirectory,"TempFilesDelme");
		tempDirectory.mkdir();
		bedDirectory = new File (saveDirectory, "Bed_"+versionedGenome);
		bedDirectory.mkdir();
		pointDataDirectory = new File (saveDirectory, "PointData");
		pointDataDirectory.mkdir();

		//flattened gene regions for scaling? Not stranded
		if (geneRegionFile != null){
			geneRegions = Region.parseStartStops(geneRegionFile, 0, 0, 0);
		}
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Novoalign Parser: Jan 2011                            **\n" +
				"**************************************************************************************\n" +
				"Parses Novoalign xxx.txt(.zip/.gz) files into center position binary PointData xxx.bar\n" +
				"files, xxx.bed files, and if appropriate, a splice junction bed file. For the later,\n" +
				"create a gene regions bed file and run it through the MergeRegions application to\n" +
				"collapse overlapping transcripts. We recommend using the following settings while\n" +
				"running Novoalign 'novoalign -r0.2 -q5 -d yourDataBase -f your_prb.txt | grep '>chr' >\n" +
				"yourResultsFile.txt'. NP works with native, colorspace, and miRNA novoalignments. \n" +

				"\nOptions:\n"+
				"-v Versioned Genome (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-f The full path directory/file text of your Novoalign xxx.txt(.zip or .gz) file(s).\n" +
				"-r Full path directory text for saving the results.\n"+
				"-p Posterior probability threshold (-10Log10(prob)) of being incorrect, defaults to 13\n"+
				"      (0.05). Larger numbers are more stringent. The parsed scores are delogged and\n" +
				"      converted to 1-prob.\n"+
				"-q Alignment score threshold, smaller numbers are more stringent, defaults to 60\n"+
				"-c Chromosome prefix, defaults to '>chr'.\n"+
				"-i Ignore strand when making splice junctions.\n"+
				"-g (Optional) Full path gene region bed file (chr start stop...) containing gene\n" +
				"      regions to use in scaling intersecting splice junctions.\n"+
				"-s Just print alignment stats, don't save any data.\n"+


				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/NovoalignParser -f /Novo/Run7/\n" +
				"     -v H_sapiens_Mar_2006 -p 20 -q 30 -r /Novo/Run7/mRNASeq/ -i -g\n" +
				"     /Anno/Hg18/mergedUCSCKnownGenes.bed \n\n" +

		"**************************************************************************************\n");

	}



	public int getNumberPassingAlignments() {
		return numberPassingAlignments;
	}	

}
