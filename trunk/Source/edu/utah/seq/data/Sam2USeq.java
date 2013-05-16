package edu.utah.seq.data;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import edu.utah.seq.data.sam.SamAlignment;
import edu.utah.seq.useq.ArchiveInfo;
import edu.utah.seq.useq.SliceInfo;
import edu.utah.seq.useq.USeqUtilities;
import edu.utah.seq.useq.data.PositionScore;
import edu.utah.seq.useq.data.PositionScoreData;
import edu.utah.seq.useq.data.Region;
import edu.utah.seq.useq.data.RegionScoreText;
import net.sf.samtools.*;
import net.sf.samtools.SAMRecord.SAMTagAndValue;
import util.bio.annotation.Bed;
import util.bio.annotation.Coordinate;
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
	private String versionedGenome;
	private boolean makeRelativeTracks = true;
	private boolean scaleRepeats = false;
	private float minimumMappingQuality = 0;
	private float maximumAlignmentScore = 1000;
	private boolean uniquesOnly = false;
	private boolean stranded = false;
	private float minimumCounts = 0;
	private int minimumLength = 0;

	//internal
	private ChromData chromData;
	private String adapter = "chrAdapt";
	private String phiX = "chrPhiX";
	private HashMap <String, ChromData> chromDataHash = new HashMap <String, ChromData>();
	private HashMap<String,Region[]> regions = null;
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
	private ArrayList<Short> baseCoverage = new ArrayList<Short>();
	private Short zeroShort = new Short((short)0);

	private int maxNumberBases;


	/**For stand alone app.*/
	public Sam2USeq(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();
		//process args
		processArgs(args);

		doWork();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		if (verbose) System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
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

	public void doWork(){
		//make tempDir
		tempDirectory = new File (samFiles[0].getParentFile(), "TempDir_"+Passwords.createRandowWord(8));
		tempDirectory.mkdir();

		//split sam files by chromosome
		if (verbose) System.out.println("\nSplitting sam files by chromsome");
		splitSamBamFiles();
		double percent = (double)numberPassingAlignments/(double)numberAlignments;

		if (verbose) {
			System.out.println(numberAlignments+" Alignments");
			System.out.println(numberPassingAlignments+" Alignments passing filters ("+Num.formatPercentOneFraction(percent)+")\n");
		}

		//make coverage track
		if (verbose) System.out.print("Making depth coverage tracks");
		makeCoverageTracks();

		if (minimumCounts !=0) bedOut.close();
		
		
		//finish read depth coverage stats
		finishReadDepthStats();
	}
	
	public void finishReadDepthStats(){
		if (regions != null){
			
			//increment chroms not scanned
			for (String chromStrand: countedChromosomes) regions.remove(chromStrand);
			for (Region[] chromRegions: regions.values()){
				for (Region r: chromRegions){
					int length = r.getLength();					
					for (int i=0; i < length; i++) {
						histogram.count(0);
						baseCoverage.add(zeroShort);
					}
				}
			}
			
			//print histogram and summary stats
			printStats();

		}
	}


	public void splitSamBamFiles(){
		for (File samFile: samFiles){
			if (verbose) System.out.print("\t"+samFile.getName());

			SAMFileReader samReader = null;
			int counter =0;
			String currentChromStrand = "";
			ChromData data = null;
			try {
				samReader = new SAMFileReader(samFile);
				samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
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

					//skip phiX and adapter
					if (sam.getReferenceName().startsWith(phiX) || sam.getReferenceName().startsWith(adapter)) continue;

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
					String chromosome = sam.getReferenceName();
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
						Object o = sam.getAttribute("IH");
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
							File f = new File(tempDirectory, currentChromStrand);
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

	private class ChromData{
		//fields
		int firstBase;
		int lastBase;
		DataOutputStream out;
		File binaryFile;
		String chromosome;
		String strand;

		//constructors
		private ChromData (int firstBase, int lastBase, String chromosome, String strand, File binaryFile, DataOutputStream out){
			this.firstBase = firstBase;
			this.lastBase = lastBase;
			this.chromosome = chromosome;
			this.strand = strand;
			this.binaryFile = binaryFile;
			this.out = out;
		}
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
			if (samFiles.length == 1) zipName= USeqUtilities.capitalizeRemoveExtension(samFiles[0], "sam") +USeqUtilities.USEQ_EXTENSION_WITH_PERIOD;
			else zipName = USeqUtilities.capitalizeRemoveExtension(samFiles[0].getParentFile(), "sam") +USeqUtilities.USEQ_EXTENSION_WITH_PERIOD;
			useqOutputFile = new File (samFiles[0].getParentFile(), zipName);
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

	public void makeCoverageTrack(){	
		//make array to hold counts
		int firstBase = chromData.firstBase;
		int lastBase = chromData.lastBase;
		float[] baseCounts = new float[lastBase - firstBase];

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
				while (mat.find()){
					String call = mat.group(2);
					int numberBases = Integer.parseInt(mat.group(1));
					//a match
					if (call.equals("M")) {
						for (int i = 0; i< numberBases; i++) baseCounts[start++] += numRepeats;
						if (numberBases > maxNumberBases) maxNumberBases = numberBases;
					}
					//just advance for all but insertions which should be skipped via the failure to match
					else start += numberBases;

				}


			}
		} catch (EOFException eof){	
			PositionScore[] positions = null;
			
			//calculate read coverage over interrogated regions
			if (regions != null){
				String chromStrand = chromData.chromosome+chromData.strand;
				countedChromosomes.add(chromStrand);
				//get regions
				Region[] chrRegions = regions.get(chromStrand);
				if (chrRegions != null) {
					//for each region
					for (Region r: chrRegions){
						int start = r.getStart() - firstBase;
						int stop = r.getStop() - firstBase;
						//before counted bases? past end? Add zeros to all bases
						for (int i=start; i< stop; i++){
							//before or after scored bases
							if (i < 0 || i >= baseCounts.length) {
								histogram.count(0);
								baseCoverage.add(zeroShort);
							}
							//nope inside
							else {
								histogram.count(baseCounts[i]);
								baseCoverage.add((short)baseCounts[i]);
							}
						}
					}
				}
			}

			//do they want relative read coverage graphs and good block counts?
			if (makeRelativeTracks && minimumCounts !=0){
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

			//delete file
			chromData.binaryFile.delete();
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
					case 'v': versionedGenome = args[i+1]; i++; break;
					case 'r': makeRelativeTracks = false; break;
					case 's': stranded = true; break;
					case 'e': scaleRepeats = true; break;
					case 'm': minimumMappingQuality = Float.parseFloat(args[++i]); break;
					case 'a': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'c': minimumCounts  = Float.parseFloat(args[++i]); break;
					case 'l': minimumLength  = Integer.parseInt(args[++i]); break;
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
		if (samFiles == null || samFiles.length ==0 || samFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.sam(.zip/.gz) file(s)!\n");

		//stranded and regionFile?
		if (regionFile !=null){
			regions = Bed.parseRegions(regionFile, stranded == false);
			histogram = new Histogram(0, 101, 101);
			//watch out for stranded analysis
			if (stranded) {
				for (String s : regions.keySet()){ 
					if (s.endsWith(".")) Misc.printErrAndExit("\nError: cannot perform a stranded analysis with non stranded interrogated regions.  Aborting.\n");
				}
			}
		}
		
		//genome version?
		if (versionedGenome == null) Misc.printErrAndExit("\nPlease provide a versioned genome (e.g. H_sapiens_Mar_2006).\n");

		//look for point directories
		if (samFiles == null || samFiles.length == 0) Misc.printExit("\nError: cannot find your sam files?\n");
		

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
			System.out.println(makeRelativeTracks+"\tMake relative covererage tracks.");
			System.out.println(scaleRepeats+"\tScale repeat alignments by the number of matches.");
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
	
	public void printStats(){
		try {
			PrintStream oldOut = System.out;
			if (logFile != null) {
				System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream(logFile))));		
			}
			System.out.println("\nInterrogated region read depth coverage statistics");
			int[] counts = histogram.getBinCounts();
			double total = histogram.getTotalBinCounts();

			System.out.println("BaseCoverage\tObservedBasesWithGivenCoverage\tFractionObserved\tFractionObservedWithGivenOrMoreCoverage");
			double numCounts = 0;
			for (int i=0; i< counts.length; i++){
				double fract = (double)counts[i] / total;
				double cumFrac = (total-numCounts)/ total;
				System.out.println(i+"\t"+counts[i]+"\t"+ Num.formatNumber(fract, 3) +"\t"+Num.formatNumber(cumFrac, 3));
				numCounts += counts[i];
				if (numCounts == total) break;
			}
			System.out.println("\nTotal interrogated bases\t"+(int)total);
			
			//print summary stats

			 short[] c = Num.arrayListOfShortToArray(baseCoverage);
			 Arrays.sort(c);
			//calc mean
			 double mean = Num.mean(c);
			 short min = c[0];
			 short max = c[c.length-1];
			//calc median
			double median = Num.median(c);
			System.out.println("Mean Coverage\t"+mean+"\nMedian Coverage\t"+median+"\nMinimum\t"+min+"\nMaximum\t"+max);
			
			
			System.out.close();
			System.setOut(oldOut);
		} catch (FileNotFoundException e) {
			System.out.println("Could not open log file");
			e.printStackTrace();
			System.exit(1);
		}

	}
	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Sam 2 USeq : April 2013                            **\n" +
				"**************************************************************************************\n" +
				"Generates per base read depth stair-step graph files for genome browser visualization.\n" +
				"By default, values are scaled per million mapped reads with no score thresholding. Can\n" +
				"also generate a list of regions that pass a minimum coverage depth.\n\n" +

				"Required Options:\n"+
				"-f Full path to a bam or a sam file (xxx.sam(.gz/.zip OK) or xxx.bam) or directory\n" +
				"      containing such. Multiple files are merged.\n"+
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
				"-c Print regions that meet a minimum # counts, defaults to 0, don't print.\n"+
				"-l Print regions that also meet a minimum length, defaults to 0.\n"+
				"-o Path to log file.  Write coverage statistics to a log file instead of stdout.\n" +

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/Sam2USeq -f /Data/SamFiles/ -r\n"+
				"     -v H_sapiens_Feb_2009 -b ccdsExons.bed.gz \n\n"+

		"**************************************************************************************\n");

	}

	public long getNumberPassingAlignments() {
		return numberPassingAlignments;
	}		

}
