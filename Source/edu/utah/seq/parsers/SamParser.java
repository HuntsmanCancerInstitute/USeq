
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import edu.utah.seq.data.*;
import edu.utah.seq.data.sam.MalformedSamAlignmentException;
import edu.utah.seq.data.sam.SamAlignment;
import util.gen.*;
import java.util.*;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord.SAMTagAndValue;


/**Parses a SAM or BAM files into point chromPointData, split by chromosome and strand.
 * For each sequence a single hit is assigned to the center position of the read.  Files are saved using the bar format.
 * The bar score value is set to one.
 * Final positions are in interbase coordinates (0 start, stop excluded).
 * @author david.nix@hci.utah.edu 
 **/
public class SamParser{
	//fields
	private File[] dataFiles;
	private File saveDirectory;
	private File pointDataDirectory;
	private File tempDirectory;
	private File workingFile;
	private String versionedGenome;
	private HashMap <String, DataOutputStream> chromOut = new HashMap <String, DataOutputStream>();
	private HashMap <String, Integer> chromLength = new HashMap <String, Integer>();
	private HashMap <String, Integer> chromStart = new HashMap <String, Integer>();
	private HashMap <String, File> tempChrData = new HashMap <String, File>();
	private float minimumMappingQuality = 13;
	private float maximumAlignmentScore = 60;
	private int numberAlignments = 0;
	private int numberUnmapped = 0;
	private int numberFailingVendorQC = 0;
	private int numberAdapter = 0;
	private int numberPhiX = 0;
	private int numberFailingAlignmentScore = 0;
	private int numberFailingMappingQualityScore = 0;
	private int numberPassingAlignments = 0;
	private int maxReadLength = 0;
	private boolean verbose = true;

	//constructors
	//for ChIPSeq app
	public SamParser (File saveDirectory, File[] dataFiles, float minimumMappingQuality, float maximumAlignmentScore, String versionedGenome){
		this.dataFiles = dataFiles;
		this.pointDataDirectory = saveDirectory;
		this.minimumMappingQuality = minimumMappingQuality;
		this.maximumAlignmentScore = maximumAlignmentScore; 
		this.versionedGenome = versionedGenome;
		verbose = false;

		tempDirectory = new File (saveDirectory,"TempFilesDelme");
		tempDirectory.mkdir();
		pointDataDirectory = new File (saveDirectory, "PointData");
		pointDataDirectory.mkdir();

		//for each file, parse and save to disk	
		for (int i=0; i< dataFiles.length; i++){
			//set working objects and parse tag file text
			workingFile = dataFiles[i];
			//watch out for bam index files
			if (workingFile.getName().endsWith(".bai")) continue;
			//split file to chromosome strand specific temp files
			boolean parsed;
			if (isSamFormat(workingFile)) parsed = parseWorkingSAMFile();
			else parsed = parseWorkingBAMFile();
			if (parsed == false) Misc.printExit("\n\tError: failed to parse, aborting.\n");
		}
		//close the writers
		closeWriters();
		//load and save
		makePointData();
		//cleanup
		IO.deleteDirectory(tempDirectory);
	}
	
	//for RNASeq app
	public SamParser (File pointDataDirectory, File samFile, float minimumMappingQuality, float maximumAlignmentScore, String versionedGenome, boolean verbose){
		dataFiles = new File[]{samFile};
		
		this.pointDataDirectory = pointDataDirectory;
		pointDataDirectory.mkdirs();
		this.minimumMappingQuality = minimumMappingQuality;
		this.maximumAlignmentScore = maximumAlignmentScore; 
		this.versionedGenome = versionedGenome;
		this.verbose = verbose;

		tempDirectory = new File (pointDataDirectory,"TempFilesDelme");
		tempDirectory.mkdir();

		workingFile = samFile;

		//split file to chromosome strand specific temp files
		boolean parsed;
		if (isSamFormat(workingFile)) parsed = parseWorkingSAMFile();
		else parsed = parseWorkingBAMFile();
		if (parsed == false) Misc.printExit("\n\tError: failed to parse, aborting.\n");

		//close the writers
		closeWriters();
		//load and save
		makePointData();
		//cleanup
		IO.deleteDirectory(tempDirectory);
	}

	public SamParser(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		System.out.println("Minimum mapping quality threshold = "+minimumMappingQuality);
		System.out.println("Maximum alignment score = "+maximumAlignmentScore);
		
		System.out.println("\nParsing and filtering:");

		//for each file, parse and save to disk	
		for (int i=0; i< dataFiles.length; i++){
			//set working objects and parse tag file text
			workingFile = dataFiles[i];
			System.out.print("\t"+workingFile.getName()+" ");

			//split file to chromosome strand specific temp files
			boolean parsed;
			if (isSamFormat(workingFile)) parsed = parseWorkingSAMFile();
			else parsed = parseWorkingBAMFile();
			if (parsed == false) Misc.printExit("\n\tError: failed to parse, aborting.\n");
			System.out.println();

		}

		//close the writers
		closeWriters();
		
		System.out.println("\n\t"+maxReadLength+" Maximum read length");

		//load, sort, make point chromPointData, and save
		System.out.print("\nLoading and saving PointData to "+ saveDirectory);
		makePointData();
		System.out.println();

		//cleanup
		IO.deleteDirectory(tempDirectory);

		//stats
		double fractionPassing = ((double)numberPassingAlignments)/((double)numberAlignments);
		System.out.println("Stats (some flags aren't set so be suspicious of zero read catagories):");
		System.out.println("\t"+numberAlignments+"\tTotal # Alignments");
		System.out.println("\t"+numberPassingAlignments+"\tAlignments passing filters ("+Num.formatPercentOneFraction(fractionPassing)+")");
		System.out.println("\t\t"+numberUnmapped+"\t# Unmapped");
		System.out.println("\t\t"+numberFailingVendorQC+"\t# Failing vendor/ platform QC");
		System.out.println("\t\t"+numberAdapter+"\t# Aligning to the adapters");
		System.out.println("\t\t"+numberPhiX+"\t# Aligning to phiX");
		System.out.println("\t\t"+numberFailingAlignmentScore+"\t# Failing to pass the alignment score.");
		System.out.println("\t\t"+numberFailingMappingQualityScore+"\t# Failing to pass the mapping quality score.");

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

	public void makePointData(){
		//for each composite data file in non splice junction data
		Iterator<String> it = tempChrData.keySet().iterator();
		while (it.hasNext()){
			String chromStrand = it.next();
			File chromDataFile = tempChrData.get(chromStrand);
			int maxLength = chromLength.get(chromStrand);
			int minStart = chromStart.get(chromStrand);
			//parse strand and chrom
			String strand = chromStrand.substring(chromStrand.length()-1);
			String chrom = chromStrand.substring(0, chromStrand.length()-1);
			//get int[]
			int[] hits = loadBinary(chromDataFile, minStart, maxLength);
			//make files
			savePointData(hits,chrom,strand, minStart);
		}
	}

	/**Loads a binary file containing int
	 * @return an array! or null if something bad happened.*/
	public static int[] loadBinary(File binaryFile, int minStart, int maxLength){
		int[] hits = new int[maxLength+1 - minStart];
		DataInputStream dis = null;
		int pos = 0;
		try {
			dis = new DataInputStream(new BufferedInputStream(new FileInputStream(binaryFile)));
			while (true) {
				pos = dis.readInt() - minStart;
				hits[pos]++;
			}
		} catch (EOFException eof){
			return hits;
		}
		catch (Exception e){
			System.out.println("\nBad binary file "+binaryFile);
			System.out.println("\nPos "+pos);
			e.printStackTrace();
			return null;
		} finally {
			if (dis != null) {
				try {
					dis.close();
				} catch (IOException ignored) {
				}
			}
		}
	}

	/**Makes the PointData, sorts and saves.*/
	public void savePointData(int[] hits, String chrom, String strand, int start){
		//calc total number of points
		int sum = Num.sumIntArray(hits);
		
		//make Point[], already sorted!
		Point[] points = new Point[sum];
		int index =0;
		for (int pos =0; pos < hits.length; pos++){
			if (hits[pos] !=0){
				int numToMake = hits[pos];
				Point p = new Point(pos + start, 1.0f);
				for (int i=0; i< numToMake; i++) points[index++] = p;
			}
		}

		//make notes
		HashMap <String,String> notes = new HashMap <String,String> ();
		notes.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
		notes.put(BarParser.SOURCE_TAG, IO.concatinateFileFullPathNames(dataFiles, ","));
		notes.put(BarParser.STRAND_TAG, strand);
		notes.put(BarParser.READ_LENGTH_TAG, +maxReadLength+"");
		notes.put(BarParser.UNIT_TAG, "Mapping quality score.");
		notes.put(BarParser.DESCRIPTION_TAG, "Generated by running the SamBamParser on Sam/Bam alignment file(s), the position is assigned to the middle of the read, interbase coordinates.");
		//make an Info object  public Info (String text, String versionedGenome, String chromosome, String strand, int readLength, HashMap<String,String> notes){
		Info info = new Info(chrom+strand, versionedGenome, chrom, strand, maxReadLength, notes);
		info.setNumberObservations(points.length);
		//make pd
		PointData pd = Point.extractPositionScores(points);			
		pd.setInfo(info);
		//write to file
		pd.writePointData(pointDataDirectory);
		//cleanup
		points = null;
	}

	/**Uses picard's SAMFileReader to parse BAM file.  This requires a header.*/
	public boolean parseWorkingBAMFile(){
		SAMFileReader samReader = null;
		try {
			samReader = new SAMFileReader(workingFile);
			samReader.setValidationStringency(ValidationStringency.SILENT);
			SAMRecordIterator it = samReader.iterator();
			while (it.hasNext()) {
				SAMRecord sam = it.next();
				numberAlignments++;

				//is it aligned?
				if (sam.getReadUnmappedFlag()){
					numberUnmapped++;
					continue;
				}
				//does it pass the vendor qc?
				if (sam.getReadFailsVendorQualityCheckFlag()){
					numberFailingVendorQC++;
					continue;
				}

				//skip phiX and adapter
				if (sam.getReferenceName().startsWith("chrPhiX")){
					numberPhiX++;
					continue;
				}
				if (sam.getReferenceName().startsWith("chrAdapt")){
					numberAdapter++;
					continue;
				}

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
				boolean failedScore = false;
				if (alignmentScore != Integer.MIN_VALUE){
					if (alignmentScore > maximumAlignmentScore){
						numberFailingAlignmentScore++;
						failedScore = true;
					}
				}
				int mappingQuality = sam.getMappingQuality();
				if (mappingQuality < minimumMappingQuality){
					numberFailingMappingQualityScore++;
					failedScore = true;
				}
				if (failedScore) continue;
				
				//increment counter
				numberPassingAlignments++;

				//parse strand
				String strand = "+";
				if (sam.getReadNegativeStrandFlag()) strand = "-";

				//calc start stop and chrStrand
				int start = sam.getAlignmentStart();
				
				//int stop = sam.getAlignmentEnd();
				//double alignmentLength = stop - start;
				//switching to use length of the read, not the Alignment End to avoid huge lengths with splice junction reads.
				
				double alignmentLength = sam.getReadLength();
				int mid = (int)Math.round((alignmentLength/2.0) + start);
				
				//reset max seq length?
				if (alignmentLength > maxReadLength) maxReadLength = (int)alignmentLength;
				
				String chrStrand = sam.getReferenceName() + strand;
				
				//check max value
				int maxLength = 0;
				if (chromLength.containsKey(chrStrand)) maxLength = chromLength.get(chrStrand);
				if (mid > maxLength) chromLength.put(chrStrand, new Integer(mid));
				//check min value
				int minStart = Integer.MAX_VALUE;
				if (chromStart.containsKey(chrStrand)) minStart = chromStart.get(chrStrand);
				if (mid < minStart) chromStart.put(chrStrand, new Integer(mid));

				//get PrintWriter
				DataOutputStream dos;
				if (chromOut.containsKey(chrStrand)) dos = chromOut.get(chrStrand);
				else {
					//make and set file
					File f = new File(tempDirectory, chrStrand);
					tempChrData.put(chrStrand, f);
					dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(f)));
					chromOut.put(chrStrand, dos);
				}

				//save data
				dos.writeInt(mid);

			}
		} catch (Exception e) {
			System.err.println("\nError parsing BAM file or writing split binary chromosome files.\n");
			e.printStackTrace();
			return false;
		} 
		return true;
	}


	public boolean parseWorkingSAMFile(){
		BufferedReader in = null;
		PrintWriter out = null;
		int numBadLines = 0;
		try {
			in = IO.fetchBufferedReader(workingFile);
			String line;
			while ((line=in.readLine())!= null) {
				
				line = line.trim();
				if (line.length() == 0 || line.startsWith("@")) continue;
				numberAlignments++;
				SamAlignment sa;
				try {
					sa = new SamAlignment(line, true);
				} catch (MalformedSamAlignmentException e) {
					System.err.println("Skipping malformed sam alignment -> "+e.getMessage());
					if (numBadLines++ > 1000) Misc.printErrAndExit("Aboring: too many malformed SAM alignments");
					continue;
				}

				//is it aligned?
				if (sa.isUnmapped()) {
					numberUnmapped++;
					continue;
				}
				//does it pass the vendor qc?
				if (sa.failedQC()) {
					numberFailingVendorQC++;
					continue;
				}

				//skip phiX and adapter
				if (sa.getReferenceSequence().startsWith("chrPhiX")){
					numberPhiX++;
					continue;
				}
				if (sa.getReferenceSequence().startsWith("chrAdapt")){
					numberAdapter++;
					continue;
				}

				//does it pass the scores threshold?
				boolean failedScore = false;
				int alignmentScore = sa.getAlignmentScore();
				if (alignmentScore != Integer.MIN_VALUE){
					if (alignmentScore > maximumAlignmentScore){
						numberFailingAlignmentScore++;
						failedScore = true;
					}
				}

				//was likelihood score of being derived from the given map location set?
				float mappingQuality = sa.getMappingQuality();
				if (mappingQuality < minimumMappingQuality){
					numberFailingMappingQualityScore++;
					failedScore = true;
				}
				
				if (failedScore) continue;

				//increment counter
				numberPassingAlignments++;

				//parse strand
				String strand = "+";
				if (sa.isReverseStrand()) strand = "-";

				//reset max seq length?
				int seqLength = sa.getSequence().length();
				if (seqLength > maxReadLength) maxReadLength = seqLength;

				//calc start stop and chrStrand
				double start = sa.getPosition();
				double stop = start + sa.countLengthOfAlignment();
				int mid = (int)Math.round(   ((stop-start)/2.0) + start     );
				
				String chrStrand = sa.getReferenceSequence()+ strand;
				
				//check max value
				int maxLength = 0;
				if (chromLength.containsKey(chrStrand)) maxLength = chromLength.get(chrStrand);
				if (mid > maxLength) chromLength.put(chrStrand, new Integer(mid));
				//check min value
				int minStart = Integer.MAX_VALUE;
				if (chromStart.containsKey(chrStrand)) minStart = chromStart.get(chrStrand);
				if (mid < minStart) chromStart.put(chrStrand, new Integer(mid));

				//get PrintWriter
				DataOutputStream dos;
				if (chromOut.containsKey(chrStrand)) dos = chromOut.get(chrStrand);
				else {
					//make and set file
					File f = new File(tempDirectory, chrStrand);
					tempChrData.put(chrStrand, f);
					dos = new DataOutputStream(new BufferedOutputStream (new FileOutputStream(f)));
					chromOut.put(chrStrand, dos);
				}

				//save data
				dos.writeInt(mid);

			}
		} catch (Exception e) {
			System.err.println("\nError parsing SAM file or writing split binary chromosome files.\n");
			e.printStackTrace();
			return false;
		} finally {
			try {
				if (in != null) in.close();
				if (out != null) out.close();
			} catch (IOException e) {}
		}
		return true;
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SamParser(args);
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
					case 'm': minimumMappingQuality = Float.parseFloat(args[++i]); break;
					case 'a': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		if (versionedGenome == null) Misc.printErrAndExit("\nPlease provide a versioned genome (e.g. H_sapiens_Mar_2006).\n");

		//pull files
		if (forExtraction == null ) Misc.printExit("\nError: cannot find your xxx.sam(.zip/.gz) or xxx.bam file(s)!\n");
		File[][] tot = new File[4][];
		tot[0] = IO.extractFiles(forExtraction,".sam");
		tot[1] = IO.extractFiles(forExtraction,".sam.gz");
		tot[2] = IO.extractFiles(forExtraction,".sam.zip");
		tot[3] = IO.extractFiles(forExtraction,".bam");

		dataFiles = IO.collapseFileArray(tot);
		if (dataFiles == null || dataFiles.length==0) dataFiles = IO.extractFiles(forExtraction);
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.sam(.zip/.gz) or xxx.bam file(s)!\n");
		if (versionedGenome == null) Misc.printExit("\nPlease enter a genome version recognized by UCSC, see http://genome.ucsc.edu/FAQ/FAQreleases.\n");
		if (saveDirectory == null)  {
			if (dataFiles.length == 1) saveDirectory = IO.makeDirectory(dataFiles[0], "_SBP");
			if (saveDirectory == null) {
				saveDirectory = new File (dataFiles[0].getParentFile(), "SBParser_"+Passwords.createRandowWord(5));
			}
		}
		saveDirectory.mkdir();
		tempDirectory = new File (saveDirectory,"TempFilesDelme");
		tempDirectory.mkdir();
		pointDataDirectory = new File (saveDirectory, "PointData");
		pointDataDirectory.mkdir();

	}	
	
	public static boolean isSamFormat(File f){
		String name = f.getName().toLowerCase();
		if (name.endsWith(".bam")) return false;
		if (name.endsWith(".sam")) return true;
		if (name.endsWith(".sam.gz")) return true;
		if (name.endsWith(".sam.zip")) return true;
		return false;
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                Sam Parser: June 2013                             **\n" +
				"**************************************************************************************\n" +
				"Parses SAM and BAM files into alignment center position PointData xxx.bar files.\n" +
				"For RNASeq data, first run the SamTranscriptomeParser to convert splice junction\n" +
				"coordinates to genomic coordinates and set -m to 0 below.\n"+

				"\nOptions:\n"+
				"-v Versioned Genome (ie H_sapiens_Mar_2006), see UCSC Browser,\n"+
				"      http://genome.ucsc.edu/FAQ/FAQreleases.\n" +
				"-f The full path file or directory containing xxx.sam(.gz/.zip OK) or xxx.bam file(s).\n" +
				"      Multiple files will be merged.\n" +
				"-r Full path directory for saving the results.\n"+
				"-m Minimum mapping quality score. Defaults to 13, bigger numbers are more stringent.\n" +
				"      This is a phred-scaled posterior probability that the mapping position of read\n" +
				"      is incorrect. For RNA-Seq data from the SamTranscriptomeParser, set this to 0.\n" +
				"-a Maximum alignment score. Defaults to 60, smaller numbers are more stringent.\n"+


				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/SamParser -f /Novo/Run7/\n" +
				"     -v C_elegans_May_2008 -m 0 -a 120  \n\n" +

		"**************************************************************************************\n");

	}

	public int getNumberPassingAlignments() {
		return numberPassingAlignments;
	}	

}
