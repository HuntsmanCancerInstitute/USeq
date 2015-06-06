package edu.utah.seq.data;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.parsers.BarParser;
import edu.utah.seq.useq.ArchiveInfo;
import edu.utah.seq.useq.SliceInfo;
import edu.utah.seq.useq.USeqUtilities;
import edu.utah.seq.useq.data.*;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import util.bio.annotation.Bed;
import util.bio.annotation.ExportIntergenicRegions;
import util.gen.*;

/**Per base read coverage
 * @author Nix
 * */
public class Sam2USeq {

	//user fields
	private File[] samFiles;
	private File tempDirectory;
	private File logFile;
	private File perRegionCoverageStats = null;
	private String versionedGenome;
	private boolean makeRelativeTracks = true;
	private boolean scaleRepeats = false;
	private float minimumMappingQuality = 0;
	private float maximumAlignmentScore = 1000;
	private boolean uniquesOnly = false;
	private boolean stranded = false;
	private boolean makeAveReadLengthGraph = false;
	private float minimumCounts = 0;
	private int minimumLength = 0;
	private double maximumCoverageCalculated = 101;

	//internal
	private ChromData chromData;
	private String adapter = "chrAdapt";
	private String phiX = "chrPhiX";
	private HashMap <String, ChromData> chromDataHash = new HashMap <String, ChromData>();
	private HashMap<String,RegionScoreText[]> regions = null;
	private ArrayList<String> countedChromosomes = new ArrayList<String>();
	private Histogram histogram;
	private ArrayList<File> files2Zip = new ArrayList<File>();
	private float scalar = 0;
	private Pattern cigarSub = Pattern.compile("(\\d+)([MSDHN])");
	private Pattern supportedCigars = Pattern.compile(".*[^\\dMIDSHN].*");
	private Pattern cigarWithRepeats = Pattern.compile("(.+\\D)(\\d+)$");
	private long numberPassingAlignments = 0;
	private double numberPassingAlignmentsForScaling = 0;
	private long numberAlignments = 0;
	private int rowChunkSize = 10000;
	private PrintWriter bedOut = null;
	private boolean verbose = true;
	private File useqOutputFile;
	private Gzipper perRegionsGzipper = null;
	private HashMap<Long,Long> baseCoverageHist = new HashMap<Long,Long>();
	private File barDirectory;
	private File chromDataFile = null;
	

	/**For stand alone app.*/
	public Sam2USeq(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		//process args
		processArgs(args);

		doWork();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");
	}

	/**For integration with the RNASeq app*/
	public Sam2USeq (File[] samFiles, File useqOutputFile, String versionedGenome, boolean makeRelativeTracks, boolean stranded, boolean scaleRepeats, boolean verbose){
		this.samFiles = samFiles;
		this.useqOutputFile = useqOutputFile;
		this.versionedGenome = versionedGenome;
		this.makeRelativeTracks = makeRelativeTracks;
		this.stranded = stranded;
		this.scaleRepeats = scaleRepeats;
		this.verbose = verbose;

		doWork();
	}

	/**For integration with the MergePairedAlignments application.*/
	public Sam2USeq (){

	}

	public void doWork(){
		//make tempDir
		if (chromDataFile == null) tempDirectory = new File (samFiles[0].getParentFile(), "TempDir_"+Passwords.createRandowWord(8));
		else tempDirectory = new File (chromDataFile.getParentFile(), "TempDir_"+Passwords.createRandowWord(8));
		tempDirectory.mkdir();

		//break up sam files? or have these already been loaded from binaries
		if (chromDataFile == null){
			//split sam files by chromosome
			if (verbose) System.out.println("\nSplitting sam files by chromsome");
			splitSamBamFiles();
			double percent = (double)numberPassingAlignments/(double)numberAlignments;
			if (verbose) {
				System.out.println(numberAlignments+" Alignments");
				System.out.println(numberPassingAlignments+" Alignments passing filters ("+Num.formatPercentOneFraction(percent)+")");
			}
		}

		//make coverage track
		if (verbose) {
			if (makeAveReadLengthGraph) System.out.print("\nMaking average alignment lenght coverage tracks");
			else System.out.print("\nMaking depth coverage tracks");
		}
		makeCoverageTracks();

		if (minimumCounts !=0) bedOut.close();

		//finish read depth coverage stats
		finishReadDepthStats();

		if (perRegionsGzipper != null) perRegionsGzipper.closeNoException();

	}

	public void finishReadDepthStats(){
		if (regions != null){
			//increment chroms not scanned
			for (String chromStrand: countedChromosomes) regions.remove(chromStrand);
			for (RegionScoreText[] chromRegions: regions.values()){
				for (RegionScoreText r: chromRegions){
					int length = r.getLength();					
					for (int i=0; i < length; i++) {
						histogram.count(0);
						if (!baseCoverageHist.containsKey((long)0)) {
							baseCoverageHist.put((long)0,(long) 0);
						}
						baseCoverageHist.put((long)0,baseCoverageHist.get((long)0) + 1);
					}
				}
			}

			//print histogram and summary stats
			printStats();

		}
	}


	public void splitSamBamFiles(){
		SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
		
		for (File samFile: samFiles){
			if (verbose) System.out.print("\t"+samFile.getName());

			int counter =0;
			String currentChromStrand = "";
			ChromData data = null;
			try {
				SamReader samReader = factory.open(samFile);
				SAMRecordIterator it = samReader.iterator();

				while (it.hasNext()) {
					SAMRecord sam = it.next();
					//print status blip
					if (++counter == 2000000){
						if (verbose) System.out.print(".");
						counter = 0;
					}

					//is it aligned?
					if (sam.getReadUnmappedFlag()) continue;

					numberAlignments++;

					//does it pass the vendor qc?
					if (sam.getReadFailsVendorQualityCheckFlag()) continue;

					//chromosome, converting to chr
					String chromosome = sam.getReferenceName();
					
					/*Don't do this its causing issues with b37
					if (chromosome.startsWith("chr") == false) chromosome = "chr"+chromosome;
					*/

					//skip phiX and adapter
					if (chromosome.startsWith(phiX) || chromosome.startsWith(adapter)) continue;

					//does it pass the score thresholds?
					List<SAMTagAndValue> attributes = sam.getAttributes();
					int alignmentScore = Integer.MIN_VALUE;
					for (SAMTagAndValue tagVal : attributes){
						String tag = tagVal.tag;
						if (tag.equals("AS")){
							alignmentScore = (Integer)tagVal.value;
							break;
						}
					}
					if (alignmentScore != Integer.MIN_VALUE){
						if (alignmentScore > maximumAlignmentScore) continue;
					}
					int mappingQuality = sam.getMappingQuality();
					if (mappingQuality < minimumMappingQuality) continue;

					//check for unique alignments? Not sure this works unless dups have been marked
					if (uniquesOnly && sam.getDuplicateReadFlag()) continue;

					//increment counter
					numberPassingAlignments++;

					//make chromosome strand
					String strand = null;
					if (stranded){
						if (sam.getReadNegativeStrandFlag()) strand = "-";
						else strand = "+";
					}
					else strand = ".";
					String chromosomeStrand = chromosome+strand;

					//check cigar
					String cigar = sam.getCigarString();
					checkCigar(cigar, sam);

					//add number of repeats to cigar string for subsequent repeat scaling?
					double forScaling = 1;
					if (scaleRepeats){
						Object o = sam.getAttribute("NH");
						if (o != null)  {
							int numRepeats = (Integer)o;
							cigar = cigar + numRepeats;
							forScaling = 1/(double)numRepeats;
						}
					}
					numberPassingAlignmentsForScaling+= forScaling;

					//get start and end
					int start = sam.getUnclippedStart() -1; 
					int end = sam.getAlignmentEnd();

					//get ChromData
					if (currentChromStrand.equals(chromosomeStrand) == false){
						currentChromStrand = chromosomeStrand;
						//already present?
						if (chromDataHash.containsKey(currentChromStrand)) {
							//yes so fetch
							data = chromDataHash.get(currentChromStrand);
						}
						else {
							//no thus make new and set
							File f = new File(tempDirectory, currentChromStrand+"ChromData");
							DataOutputStream dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(f)));
							data = new ChromData (start, end, chromosome, strand, f, dos);
							chromDataHash.put(currentChromStrand, data);
						}
					}

					//set first and last?
					if (start < data.firstBase) data.firstBase = start;
					if (end > data.lastBase) data.lastBase = end;

					//save data start and cigar
					data.out.writeInt(start);
					data.out.writeUTF(cigar);
				}
				samReader.close();
				if (verbose) System.out.println();
			} catch (Exception e){
				System.err.println("\nError parsing sam file or writing split binary chromosome files.\n\nToo many open files exception? Too many chromosomes? " +
						"If so then login as root and set the default higher using the ulimit command (e.g. ulimit -n 10000)\n");
				e.printStackTrace();
				System.exit(1);
			}
			if (verbose) System.out.println();
		}

		closeWriters();

		//set scalar
		scalar = (float)(numberPassingAlignmentsForScaling/ 1000000.0);
	}

	/**Closes writers.*/
	public void closeWriters(){
		try{
			Iterator<ChromData> it = chromDataHash.values().iterator();
			while (it.hasNext()) it.next().out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void checkCigar(String cigar, SamAlignment sam){
		Matcher mat = supportedCigars.matcher(cigar);
		if (mat.matches()) Misc.printErrAndExit("\nUnsupported cigar string in -> \n"+sam.toString());
	}

	public void makeCoverageTracks(){	

		//for each ChromData, these might be stranded
		for (ChromData cd : chromDataHash.values()){
			chromData = cd;
			String s = "";
			if (chromData.strand.equals(".") == false) s= chromData.strand;
			if (verbose) System.out.print(" "+chromData.chromosome+ s);
			makeCoverageTrack();
		}
		if (verbose) System.out.println();

		//make useq archive
		writeReadMeTxt();

		if (useqOutputFile == null) {
			String zipName = null;
			File dir = null;
			if (chromDataFile != null) {
				dir = chromDataFile.getParentFile();
				zipName = USeqUtilities.capitalizeRemoveExtension(dir) +USeqUtilities.USEQ_EXTENSION_WITH_PERIOD;
			}
			else {
				if (samFiles.length == 1) {
					zipName= USeqUtilities.capitalizeRemoveExtension(samFiles[0]) +USeqUtilities.USEQ_EXTENSION_WITH_PERIOD;
					dir = samFiles[0].getParentFile();
				}
				else {
					dir = samFiles[0].getParentFile();
					zipName = USeqUtilities.capitalizeRemoveExtension(dir) +USeqUtilities.USEQ_EXTENSION_WITH_PERIOD;
				}
			}
			useqOutputFile = new File (dir, zipName);
		}

		File[] files = new File[files2Zip.size()];
		files2Zip.toArray(files);
		if (USeqUtilities.zip(files, useqOutputFile)== false) {
			useqOutputFile.delete();
			USeqUtilities.deleteDirectory(tempDirectory);
			USeqUtilities.printErrAndExit("\nProblem zipping data for "+samFiles[0]);
		}
		USeqUtilities.deleteDirectory(tempDirectory);
	}

	/**Sums M bases in cigar.*/
	public int alignmentLength (String cigar){
		//for each cigar block
		Matcher mat = cigarSub.matcher(cigar);
		int alignmentLength = 0;
		while (mat.find()){
			String call = mat.group(2);
			//a match
			if (call.equals("M")) {
				int numberBases = Integer.parseInt(mat.group(1));
				alignmentLength+= numberBases;
			}
		}
		return alignmentLength;
	}

	public void makeCoverageTrack(){	
		//make array to hold counts
		int firstBase = chromData.firstBase;
		int lastBase = chromData.lastBase;
		float[] baseCounts = new float[lastBase - firstBase];
		float[] alignmnentLengths = new float[baseCounts.length];

		//fetch dis
		DataInputStream dis = null;
		try {
			dis = new DataInputStream(new BufferedInputStream(new FileInputStream(chromData.binaryFile)));

			while (true){
				//read binary
				int start =  dis.readInt() - firstBase;
				String cigar = dis.readUTF();

				//scaling repeats?
				float numRepeats = 1;
				if (scaleRepeats){
					Matcher c = cigarWithRepeats.matcher(cigar);
					if (c.matches()){
						numRepeats = 1.0f/Float.parseFloat(c.group(2));
						cigar = c.group(1);
					}
				}

				//for each cigar block
				Matcher mat = cigarSub.matcher(cigar);
				int alignmentLength = alignmentLength(cigar);
				while (mat.find()){
					String call = mat.group(2);
					int numberBases = Integer.parseInt(mat.group(1));
					//a match
					if (call.equals("M")) {
						for (int i = 0; i< numberBases; i++) {
							baseCounts[start] += numRepeats;
							alignmnentLengths[start] += alignmentLength;
							start++;
						}
					}
					//just advance for all but insertions which should be skipped via the failure to match
					else start += numberBases;
				}
			}
		} catch (EOFException eof){	
			PositionScore[] positions = null;
			try {
				//calculate read coverage over interrogated regions
				if (regions != null){
					String chromStrand = chromData.chromosome+chromData.strand;
					if (chromData.strand.equals(".")) countedChromosomes.add(chromData.chromosome);
					else countedChromosomes.add(chromStrand);
					//get regions
					RegionScoreText[] chrRegions = regions.get(chromStrand);
					if (chrRegions == null) chrRegions = regions.get(chromData.chromosome);
					if (chrRegions != null) {
						//for each region
						for (RegionScoreText r: chrRegions){
							int start = r.getStart() - firstBase;
							int stop = r.getStop() - firstBase;
							float[] counts = new float[r.getLength()];
							double num10 = 0;
							double num20 = 0;
							int index = 0;
							//before counted bases? past end? Add zeros to all bases
							for (int i=start; i< stop; i++){
								//before or after scored bases
								if (i < 0 || i >= baseCounts.length) {
									histogram.count(0);
									if (!baseCoverageHist.containsKey((long)0)) {
										baseCoverageHist.put((long)0,(long) 0);
									}
									baseCoverageHist.put((long)0,baseCoverageHist.get((long)0) + 1);
									//baseCoverage.add(zeroShort);
									counts[index++] = 0.0f;
								}
								//nope inside
								else {
									histogram.count(baseCounts[i]);
									long bc = (long)baseCounts[i];
									if (!baseCoverageHist.containsKey(bc)) {
										baseCoverageHist.put(bc, (long)0);
									}
									baseCoverageHist.put(bc,baseCoverageHist.get(bc) + 1);
									//baseCoverage.add((short)baseCounts[i]);
									counts[index++] = baseCounts[i];
									if (baseCounts[i] >= 20){
										num10++;
										num20++;
									}
									else if (baseCounts[i] >= 10) num10++;
								}
							}
							double total = (double) counts.length;
							String fraction10 = Num.formatNumber(num10/total, 2);
							String fraction20 = Num.formatNumber(num20/total, 2);
							Arrays.sort(counts);
							perRegionsGzipper.println(chromData.chromosome+"\t"+chromData.strand+"\t"+r.getStart()+"\t"+r.getStop()+"\t"+
									r.getText()+"\t"+fraction10+"\t"+fraction20+"\t"+Num.statFloatArrayWithSizeChecks(counts)); 
						}
					}
				}
			} catch (IOException e){
				e.printStackTrace();
				Misc.printErrAndExit("\nError writing data to gzipper for per region coverage stats.\n");
			}
			//write out bar point data?
			if (barDirectory != null) saveBarPointData(firstBase, baseCounts);

			//make ave read length graph?
			if (makeAveReadLengthGraph){
				//make averages
				for (int i=0; i< baseCounts.length; i++) {
					if (baseCounts[i] != 0 ) alignmnentLengths[i] = alignmnentLengths[i] / baseCounts[i];
				}
				//make PositionScore[]
				positions = makeStairStepGraph(firstBase, alignmnentLengths, false);
				//write data to disk
				SliceInfo sliceInfo = new SliceInfo(chromData.chromosome, chromData.strand,0,0,0,null);
				PositionScoreData psd = new PositionScoreData (positions, sliceInfo);
				psd.sliceWritePositionScoreData(rowChunkSize, tempDirectory, files2Zip);

			}
			//do they want relative read coverage graphs and good block counts?
			else if (makeRelativeTracks && minimumCounts !=0){
				//first make positions without scaling
				positions = makeStairStepGraph(firstBase, baseCounts, false);
				//make blocks 
				makeGoodBlocks(firstBase, baseCounts, chromData.chromosome, chromData.strand);
				//now make positions with scaling, this will alter the baseCounts
				positions = makeStairStepGraph(firstBase, baseCounts, true);
				//write data to disk
				SliceInfo sliceInfo = new SliceInfo(chromData.chromosome, chromData.strand,0,0,0,null);
				PositionScoreData psd = new PositionScoreData (positions, sliceInfo);
				psd.sliceWritePositionScoreData(rowChunkSize, tempDirectory, files2Zip);
			}

			else {
				//make PositionScore[]
				positions = makeStairStepGraph(firstBase, baseCounts, makeRelativeTracks);

				//write data to disk
				SliceInfo sliceInfo = new SliceInfo(chromData.chromosome, chromData.strand,0,0,0,null);
				PositionScoreData psd = new PositionScoreData (positions, sliceInfo);
				psd.sliceWritePositionScoreData(rowChunkSize, tempDirectory, files2Zip);

				if (minimumCounts !=0) makeGoodBlocks(firstBase, baseCounts, chromData.chromosome, chromData.strand);
			}

			//delete file if this was generated by splitting sams
			if (samFiles.length != 0) chromData.binaryFile.delete();
		}
		catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nCan't make Read Coverage track for "+chromData.binaryFile.getName());
		} finally {
			if (dis != null) {
				try {
					dis.close();
				} catch (IOException ignored) {}
			}
		}
	}

	/**Saves bar formatted point data for use in apps like AggregatePlotter.*/
	private void saveBarPointData(int firstBase, float[] baseCounts) {
		//create non zero Point data
		ArrayList<Point> al = new ArrayList<Point>();
		double scoreTotal = 0;
		for (int i=0; i< baseCounts.length; i++){
			if (baseCounts[i] != 0) {
				al.add(new Point (i+firstBase, baseCounts[i]));
				scoreTotal += baseCounts[i];
			}
		}
		Point[] pts = new Point[al.size()];
		al.toArray(pts);

		PointData hitTrack = Point.extractPositionScores(pts);
		Info info = new Info();
		info.setChromosome(chromData.chromosome);
		info.setStrand(chromData.strand);
		info.setNumberObservations(pts.length);
		info.setScoreTotal(scoreTotal);
		info.setVersionedGenome(versionedGenome);
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = new HashMap<String,String>();		
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
		//save in info
		info.setNotes(map);
		hitTrack.setInfo(info);
		hitTrack.writePointData(barDirectory);	
	}

	public void makeGoodBlocks(int firstBase, float[] baseCount, String chromosome, String strand){
		String nameScoreStrand = "\t.\t0\t"+strand;
		boolean[] falseMask = new boolean[baseCount.length];
		Arrays.fill(falseMask, true);
		for (int i=0; i< baseCount.length; i++){
			if (baseCount[i] >= minimumCounts) falseMask[i] = false;
		}
		int[][] blocks = ExportIntergenicRegions.fetchFalseBlocks(falseMask, 0, 0);
		for (int j=0; j< blocks.length; j++){
			int length = 1+ blocks[j][1] - blocks[j][0];
			if (length >= minimumLength) bedOut.println(chromosome+"\t"+(blocks[j][0]+ firstBase)+"\t"+(blocks[j][1]+ firstBase +1) + nameScoreStrand);
		}
	}

	/**Makes a stairstep graph from base count data.  The firstBase is added to the index to create a bp position.  Note if you scaleCounts then the original float[] is modified!*/
	public PositionScore[] makeStairStepGraph(int firstBase, float[] baseCount, boolean scaleCounts){
		//scale it?
		if (scaleCounts){
			for (int i=0; i< baseCount.length; i++){
				if (baseCount[i] !=0) baseCount[i] = baseCount[i]/scalar;
			}
		}

		ArrayList<PositionScore> psAL = new ArrayList<PositionScore>();
		PositionScore ps = null;

		//find first non zero score
		float hits = 0;
		int basePosition = 0;
		int i = 0;
		for (; i< baseCount.length; i++){
			if (baseCount[i] !=0){
				hits = baseCount[i];
				basePosition = firstBase +i;
				break;
			}
		}

		//set beginning zero point?
		if (basePosition !=0) {
			int basePositionMinusOne = basePosition-1;
			int priorPosition = -1;
			if (ps != null) priorPosition = ps.getPosition();
			if (priorPosition != basePositionMinusOne) {
				ps = new PositionScore(basePositionMinusOne, 0);
				psAL.add(ps);
			}
		}
		//set first data point
		ps = new PositionScore(basePosition,hits);
		psAL.add(ps);

		//for each subsequent
		for (int m=i+1; m<baseCount.length; m++){
			//add Point only if hits differs
			if (baseCount[m] != hits){
				//add old if not the same position as prior
				if (ps.getPosition() != basePosition){
					ps = new PositionScore(basePosition,hits);
					psAL.add(ps);
				}			
				//set new
				basePosition = m + firstBase;
				hits = baseCount[m];
				ps = new PositionScore(basePosition,hits);
				psAL.add(ps);	
			}
			else basePosition = m + firstBase;
		}

		//add last
		if (ps.getPosition() != basePosition) {
			ps = new PositionScore(basePosition,hits);
			psAL.add(ps);
		}

		//set ending point
		ps = new PositionScore(basePosition+1,0);
		psAL.add(ps);

		//return array
		PositionScore[] p = new PositionScore[psAL.size()];
		psAL.toArray(p);

		return p;
	}

	public static void printSam(SAMRecord sam){
		System.out.println(sam+" "+sam.getCigarString());
		System.out.println(sam.getAlignmentStart()+" - "+sam.getAlignmentEnd());
		System.out.println(sam.getUnclippedStart()+" - "+sam.getUnclippedEnd());

	}

	public void checkCigar(String cigar, SAMRecord sam){
		Matcher mat = supportedCigars.matcher(cigar);
		if (mat.matches()) Misc.printErrAndExit("\nUnsupported cigar string in -> \n"+sam.toString());
	}


	private void writeReadMeTxt(){
		try {
			ArchiveInfo ai = new ArchiveInfo(versionedGenome, null, verbose);
			//set data type, graph or region
			ai.setDataType(ArchiveInfo.DATA_TYPE_VALUE_GRAPH);
			ai.setInitialGraphStyle(ArchiveInfo.GRAPH_STYLE_VALUE_STAIRSTEP);
			//set text file source
			String source = null;
			if (samFiles.length == 1) source = samFiles[0].toString();
			else if (samFiles.length == 0) source = chromDataFile.getParent();
			else source = samFiles[0].getParent();
			ai.setOriginatingDataSource(source);
			//set color
			ai.setInitialColor("#0066FF");
			//write
			File readme = ai.writeReadMeFile(tempDirectory);
			files2Zip.add(0, readme);
		} catch (IOException e){
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new Sam2USeq(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		if (verbose) System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		File regionFile = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': forExtraction = new File(args[i+1]); i++; break;
					case 'p': perRegionCoverageStats = new File(args[i+1]); i++; break;
					case 'v': versionedGenome = args[i+1]; i++; break;
					case 'r': makeRelativeTracks = false; break;
					case 'd': barDirectory = new File(args[i+1]); i++; break;
					case 's': stranded = true; break;
					case 'e': scaleRepeats = true; break;
					case 'k': makeAveReadLengthGraph = true; break;
					case 'm': minimumMappingQuality = Float.parseFloat(args[++i]); break;
					case 'a': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'c': minimumCounts  = Float.parseFloat(args[++i]); break;
					case 'l': minimumLength  = Integer.parseInt(args[++i]); break;
					case 'x': maximumCoverageCalculated  = (Double.parseDouble(args[++i]) + 1.0); break;
					case 'b': regionFile = new File(args[i+1]); i++; break;
					case 'o': logFile = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		File[][] tot = new File[4][];
		tot[0] = IO.extractFiles(forExtraction,".sam");
		tot[1] = IO.extractFiles(forExtraction,".sam.gz");
		tot[2] = IO.extractFiles(forExtraction,".sam.zip");
		tot[3] = IO.extractFiles(forExtraction,".bam");
		samFiles = IO.collapseFileArray(tot);
		if (samFiles == null || samFiles.length ==0 || samFiles[0].canRead() == false) {
			//look for chromData.obj
			chromDataFile = new File (forExtraction, "chromData.obj");
			if (chromDataFile.exists()) {
				if (stranded || scaleRepeats) Misc.printErrAndExit("\nError: cannot perform a stranded or scaled repeat analysis with ChromData from MergePairedAlignments.  Aborting.\n");
				loadChromData();
			}
			else Misc.printExit("\nError: cannot find your xxx.sam(.zip/.gz) file(s)!\n");
		}

		//barDirectory for saving point data?
		if (barDirectory!= null) barDirectory.mkdirs();

		//stranded and regionFile?
		if (regionFile !=null){
			regions = Bed.parseBedFile(regionFile, stranded == false, false);
			
			/*Don't do this, it breaks b37 stuff!
			 * HashMap<String,RegionScoreText[]> tempRegions = Bed.parseBedFile(regionFile, stranded == false, false);
			regions = new HashMap<String,RegionScoreText[]>();
			for (String chrom: tempRegions.keySet()) {
				if (!chrom.startsWith("chr")) {
					String modChrom = "chr" + chrom;
					regions.put(modChrom, tempRegions.get(chrom));
				} else {
					regions.put(chrom,tempRegions.get(chrom));
				}
			}*/
			
			histogram = new Histogram(0, maximumCoverageCalculated, (int)maximumCoverageCalculated);
			//watch out for stranded analysis
			if (stranded) {
				for (String s : regions.keySet()){ 
					if (s.endsWith(".")) Misc.printErrAndExit("\nError: cannot perform a stranded analysis with non stranded interrogated regions.  Aborting.\n");
				}
			}
			//start up gzipper for saving the individual region read coverage stats
			if (perRegionCoverageStats == null){
				String name = Misc.removeExtension(regionFile.getName()) +"_RegionCoverageStats.txt.gz";
				perRegionCoverageStats = new File (regionFile.getParentFile(), name);
			}
			try {
				perRegionsGzipper = new Gzipper(perRegionCoverageStats);
				perRegionsGzipper.println("Chr\tStrand\tStart\tStop\tInfo\tFracBPs>=10\tFracBPs>=20\tMean\tMedian\tStdDev\tMin\tMax\t10th\t90th");
			} catch (Exception e) {
				e.printStackTrace();
				Misc.printErrAndExit("Failed instantiating a gziper for saving the individual read coverage stats.");
			} 

		}

		//genome version?
		if (versionedGenome == null) Misc.printErrAndExit("\nPlease provide a versioned genome (e.g. H_sapiens_Mar_2006).\n");

		if (minimumCounts !=0){
			String name = "regions";
			if (samFiles.length == 1) name = Misc.removeExtension(samFiles[0].getName())+"_Regions";
			name = name+minimumCounts+"C"+minimumLength+"L.bed";
			File bed = new File (samFiles[0].getParentFile(), name);
			try {
				bedOut = new PrintWriter (new FileWriter (bed));
			} catch (IOException e) {
				e.printStackTrace();
				Misc.printErrAndExit("\nProblem making PrintWriter!\n");
			}
		}

		//settings
		if (verbose) {
			System.out.println("Settings:");
			if (makeAveReadLengthGraph){
				System.out.println(makeAveReadLengthGraph+"\tMake average alignment length graph.");
			}
			else {
				System.out.println(makeRelativeTracks+"\tMake relative covererage tracks.");
				System.out.println(scaleRepeats+"\tScale repeat alignments by the number of matches.");
			}
			System.out.println(stranded+"\tMake stranded covererage tracks.");
			System.out.println(minimumMappingQuality+"\tMinimum mapping quality score.");
			System.out.println(maximumAlignmentScore+"\tMaximum alignment score.");
			if (minimumCounts !=0){
				System.out.println(minimumCounts+"\tMinimum counts.");
				System.out.println(minimumLength+"\tMinimum length.");
			}
			if (regionFile!= null) System.out.println(regionFile.getName()+"\tInterrogated region file.");
		}


	}	

	private void loadChromData() {
		ChromDataSave[] cds = (ChromDataSave[]) IO.fetchObject(chromDataFile);
		ChromData[] cd = new ChromData[cds.length];

		//load hash
		for (int i=0; i< cds.length; i++){
			cd[i] = cds[i].getChromData();
			//check if data file is present
			File dataFile = new File (chromDataFile.getParentFile(), cd[i].binaryFile.getName());
			if (dataFile.exists() == false) Misc.printErrAndExit("Error: cannot find the binary data file for ChromData "+dataFile);
			//if (dataFile.exists()){
				cd[i].binaryFile = dataFile;
				chromDataHash.put(cd[i].chromosome+".", cd[i]);
			//}
			
		}

		//load counts
		File countsFile = new File (chromDataFile.getParentFile(), "count.obj");
		if (countsFile.exists() == false) Misc.printErrAndExit("Error: cannot find the binary count.obj file "+countsFile);
		Long count = (Long) IO.fetchObject(countsFile);
		numberPassingAlignmentsForScaling = count.longValue();
		
		//set scalar
		scalar = (float)(((double)numberPassingAlignmentsForScaling)/ 1000000.0);
	}

	public void printStats(){
		try {
			PrintStream oldOut = System.out;
			if (logFile != null) {
				System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream(logFile))));		
			}
			System.out.println("\nInterrogated region read depth coverage statistics");
			float[] counts = histogram.getBinCountsFloat();
			float total = histogram.getTotalBinCountsFloat();

			System.out.println("BaseCoverage\tObservedBasesWithGivenCoverage\tFractionObserved\tFractionObservedWithGivenOrMoreCoverage");
			double numCounts = 0;
			for (int i=0; i< counts.length; i++){
				double fract = (double)counts[i] / total;
				double cumFrac = (total-numCounts)/ total;
				System.out.println(i+"\t"+counts[i]+"\t"+ Num.formatNumber(fract, 3) +"\t"+Num.formatNumber(cumFrac, 3));
				numCounts += counts[i];
				if (numCounts == total) break;
			}
			System.out.println("\nTotal interrogated bases\t"+ total);

			//print summary stats
			System.out.println("Mean Coverage\t"+this.calcHistMean()+"\nMedian Coverage\t"+this.calcHistMedian()+"\nMinimum\t"+calcHistMin()+"\nMaximum\t"+calcHistMax());

			if (logFile != null){
				System.out.close();
				System.setOut(oldOut);
			}
		} catch (FileNotFoundException e) {
			System.out.println("Could not open log file");
			e.printStackTrace();
			System.exit(1);
		}

	}

	private long calcHistMin() {
		return Collections.min(this.baseCoverageHist.keySet());
	}

	private long calcHistMax() {
		return Collections.max(this.baseCoverageHist.keySet());
	}

	private double calcHistMean() {
		long totalPositions = 0;

		long totalCoverage = 0;
		for (long k: this.baseCoverageHist.keySet()) {
			long count = this.baseCoverageHist.get(k);
			totalPositions += count;
			totalCoverage += (count * k);
		}
		System.out.println("Total postitions: " + totalPositions);
		double mean = (double)totalCoverage / totalPositions; 
		return mean;

	}

	private double calcHistMedian() {
		int totalPositions = 0;
		for (long i: this.baseCoverageHist.values()) {
			totalPositions += i;
		}
		int midPosition = totalPositions / 2; 

		ArrayList<Long> keySet = new ArrayList<Long>();
		keySet.addAll(this.baseCoverageHist.keySet());
		Collections.sort(keySet);

		double median = 0;

		if (totalPositions % 2 == 1) {
			long count = 0;
			for (Long key: keySet) {
				long elements = this.baseCoverageHist.get(key);
				if (midPosition > count && midPosition < count+elements) {
					median = key;
					break;
				}
				count += elements;
			}
		} else {
			long count = 0;
			long e1 = 0;
			long e2 = 0;
			boolean e1found = false;
			boolean e2found = false;
			for (long key: keySet) {
				long elements = this.baseCoverageHist.get(key);
				if (!e1found && (midPosition > count) && (midPosition < (count+elements))) {
					e1 = key;
					e1found = true;
				} 
				if (!e2found && ((midPosition+1) > count) && ((midPosition+1) < (count + elements))) {
					e2 = key;
					e2found = true;
				}

				if (e1found && e2found) {
					break;
				}
				count += elements;
			}

			median = (e1 + e2) / 2.0;
		}
		return median;
	}


	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Sam 2 USeq : June 2015                            **\n" +
				"**************************************************************************************\n" +
				"Generates per base read depth stair-step graph files for genome browser visualization.\n" +
				"By default, values are scaled per million mapped reads with no score thresholding. Can\n" +
				"also generate a list of regions that pass a minimum coverage depth.\n\n" +

				"Required Options:\n"+
				"-f Full path to a bam or a sam file (xxx.sam(.gz/.zip OK) or xxx.bam) or directory\n" +
				"      containing such. Multiple files are merged. Also works with a directory of\n"+
				"      ChromData from MergePairedAlignments.\n"+
				"-v Versioned Genome (ie H_sapiens_Mar_2006, D_rerio_Jul_2010), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +

				"\nDefault Options:\n"+
				"-s Generate strand specific coverage graphs.\n"+
				"-m Minimum mapping quality score. Defaults to 0, bigger numbers are more stringent.\n" +
				"      This is a phred-scaled posterior probability that the mapping position of read\n" +
				"      is incorrect.\n" +
				"-a Maximum alignment score. Defaults to 1000, smaller numbers are more stringent.\n"+
				"-r Don't scale graph values. Leave as actual read counts. \n"+
				"-e Scale repeat alignments by dividing the alignment count at a given base by the\n" +
				"      total number of genome wide alignments for that read.  Repeat alignments are\n" +
				"      thus given fractional count values at a given location. Requires that the IH\n" +
				"      tag was set.\n"+
				"-b Path to a region bed file (tab delim: chr start stop ...) to use in calculating\n" +
				"      read coverage statistics.  Be sure these do not overlap! Run the MergeRegions app\n" +
				"      if in doubt.\n"+ 
				"-x Maximum read coverage stats calculated, defaults to 100, for use with -b.\n"+
				"-p Path to a file for saving per region coverage stats. Defaults to variant of -b.\n"+
				"-c Print regions that meet a minimum # counts, defaults to 0, don't print.\n"+
				"-l Print regions that also meet a minimum length, defaults to 0.\n"+
				"-o Path to log file.  Write coverage statistics to a log file instead of stdout.\n" +
				"-k Make average alignment length graph instead of read depth.\n"+
				"-d Full path to a directory for saving bar binary PointData, defaults to not saving.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/Sam2USeq -f /Data/SamFiles/ -r\n"+
				"     -v H_sapiens_Feb_2009 -b ccdsExons.bed.gz \n\n"+

				"**************************************************************************************\n");

	}

	public long getNumberPassingAlignments() {
		return numberPassingAlignments;
	}		

}
