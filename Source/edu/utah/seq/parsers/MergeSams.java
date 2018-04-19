
package edu.utah.seq.parsers;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import java.io.*;
import java.util.regex.*;
import util.gen.*;
import java.util.*;
import edu.utah.seq.data.sam.*;

/**
 * @author david.nix@hci.utah.edu 
 **/
public class MergeSams{
	//fields
	private File[] dataFiles;
	private File saveFile;
	private File outputFile;
	private String samHeader = null;
	private boolean deleteTempSamFile = true;
	private float maximumAlignmentScore = 300;
	private float minimumMappingQualityScore = 0;
	private int numberAlignments = 0;
	private int numberUnmapped = 0;
	private int numberFailingVendorQC = 0;
	private int numberPassingAlignments = 0;
	private int numberFailingAlignmentScore = 0;
	private int numberFailingMappingQualityScore = 0;
	private int numberAdapter = 0;
	private int numberPhiX = 0;
	private boolean verbose = true;
	private boolean addProgramArgs = false;
	private boolean appendFileIndex = true;
	
	private Gzipper samOut;
	private boolean saveBadReads = false;
	private HashMap <String, Integer> chromLength = new HashMap <String, Integer>();
	private HashSet<String> readGroups = new HashSet<String>();
	private HashSet<String> programParams = new HashSet<String>();
	private String programArguments;
	private static final Pattern TAB = Pattern.compile("\\t");
	private String genomeVersion = null;
	private SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	
	
	//constructors
	public MergeSams(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			if (verbose) System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public MergeSams(File[] samBam, File mergedBam, boolean verbose, boolean appendFileIndex){
		try {
			dataFiles = samBam;
			saveFile = mergedBam;
			this.verbose = verbose;
			this.appendFileIndex = appendFileIndex;
			doWork();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void doWork() throws IOException{
		//make print writer
		outputFile = new File(saveFile+".temp.sam.gz");
		if (deleteTempSamFile) outputFile.deleteOnExit();
		samOut = new Gzipper(outputFile);
		
		//add header
		if (samHeader == null) {
			for (int i=0; i< dataFiles.length; i++){
				SamReader sr = factory.open(dataFiles[i]);
				parseHeader(sr.getFileHeader().getTextHeader().trim());
				sr.close();
			}
			samHeader = fetchSamHeader();
		}
		samOut.println(samHeader);

		//for each file, parse and save to disk	
		if (verbose) System.out.println("\nParsing, filtering, and merging SAM files...");
		for (int i=0; i< dataFiles.length; i++){
			if (verbose) System.out.print("\t"+dataFiles[i].getName());
			if (appendFileIndex) parseSam(dataFiles[i],i+"_"); 
			else parseSam(dataFiles[i],""); 
			if (verbose) System.out.println();
		}

		//close the writers
		samOut.close();
		
		if (verbose) System.out.println("\nSorting and writing bam output with Picard's SortSam...");
		new PicardSortSam (outputFile, saveFile);


		//stats
		double fractionPassing = ((double)numberPassingAlignments)/((double)numberAlignments);
		if (verbose) { 
			System.out.println("\nStats:\n");
			System.out.println("\t"+numberAlignments+"\tTotal # Alignments from raw sam file");
			System.out.println("\t"+numberPassingAlignments+"\tAlignments passing filters ("+Num.formatPercentOneFraction(fractionPassing)+")");
			System.out.println("\t\t"+numberUnmapped+"\t# Unmapped Reads");
			System.out.println("\t\t"+numberFailingVendorQC+"\t# Alignments failing vendor/ platform QC and or malformed");
			System.out.println("\t\t"+numberFailingAlignmentScore+"\t# Alignments failing alignment score");
			System.out.println("\t\t"+numberFailingMappingQualityScore+"\t# Alignments failing mapping quality score");
			System.out.println("\t\t"+numberAdapter+"\t# Adapter alignments");
			System.out.println("\t\t"+numberPhiX+"\t# PhiX alignments");
			System.out.println();
		}

	}

	public void parseHeader(String textHeader){
		String[] lines = textHeader.split("\\n");
		for (int i=0; i< lines.length; i++) {
			lines[i] = lines[i].trim();
			if (lines[i].length() == 0) continue;
			
			if (lines[i].startsWith("@SQ")) {
				//parse genome version, chrom and length
				String chrom = null;
				int length = 0;
				String[] tokens = TAB.split(lines[i]);
				for (String t : tokens){
					if (t.startsWith("SN:")) chrom = t.substring(3);
					else if (t.startsWith("LN:")) length = Integer.parseInt(t.substring(3));
					else if (genomeVersion == null && t.startsWith("AS:")) genomeVersion = t;
				}
				//set chrom and length?
				if (chrom != null){
					int currLength = 0;
					if (chromLength.containsKey(chrom)) currLength = chromLength.get(chrom);
					if (length > currLength) chromLength.put(chrom, new Integer(length));
				}
			}
			else if (lines[i].startsWith("@RG")) readGroups.add(lines[i]);
			else if (lines[i].startsWith("@PG")) programParams.add(lines[i]);
		}	
	}

	public boolean parseSam(File samFile, String nameAppender){
		SAMRecord sam = null;
		try {
			SamReader sr = factory.open(samFile);
			SAMRecordIterator it = sr.iterator();
			int dotCounter = 0;
			while (it.hasNext()) {
				if (verbose && ++dotCounter > 1000000){
					System.out.print(".");
					dotCounter = 0;
				}
				try {
					sam = it.next();
					numberAlignments++;

					//is it aligned?
					if (sam.getReadUnmappedFlag()){
						numberUnmapped++;
						if (saveBadReads) samOut.println(nameAppender+ sam.getSAMString().trim());
						continue;
					}
					//does it pass the vendor qc?
					if (sam.getReadFailsVendorQualityCheckFlag()){
						numberFailingVendorQC++;
						if (saveBadReads) samOut.println(nameAppender+ sam.getSAMString().trim());
						continue;
					}

					//skip phiX and adapter
					if (sam.getReferenceName().startsWith("chrPhiX")){
						numberPhiX++;
						if (saveBadReads) samOut.println(nameAppender+ sam.getSAMString().trim());
						continue;
					}
					if (sam.getReferenceName().startsWith("chrAdapt")){
						numberAdapter++;
						if (saveBadReads) samOut.println(nameAppender+ sam.getSAMString().trim());
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
					if (alignmentScore != Integer.MIN_VALUE){
						if (alignmentScore > maximumAlignmentScore){
							numberFailingAlignmentScore++;
							if (saveBadReads) samOut.println(nameAppender+ sam.getSAMString().trim());
							continue;
						}
					}
					int mappingQuality = sam.getMappingQuality();
					if (mappingQuality < minimumMappingQualityScore){
						numberFailingMappingQualityScore++;
						if (saveBadReads) samOut.println(nameAppender+ sam.getSAMString().trim());
						continue;
					}

					//OK, it passes, increment counter
					numberPassingAlignments++;
					samOut.println(nameAppender+ sam.getSAMString().trim());

				} catch (IllegalArgumentException i){
					numberFailingVendorQC++;
				}

			}
			sr.close();

		} catch (Exception e) {
			if (sam != null) System.err.println("Problem processing -> "+sam.getSAMString());
			e.printStackTrace();

			return false;
		}
		return true;
	}

	public String fetchSamHeader() {
		ArrayList<String> al = new ArrayList<String>();
		//add unsorted
		al.add("@HD\tVN:1.0\tSO:unsorted");

		//as sq lines for each chromosome @SQ	SN:chr10	AS:mm9	LN:129993255
		String gv = "";
		if (genomeVersion != null) gv = "\t" +genomeVersion;
		for (String chromosome: chromLength.keySet()){
			int length = chromLength.get(chromosome);
			al.add("@SQ\tSN:"+chromosome+ gv+ "\tLN:"+length);
		}
		//add programs
		al.add("@PG\tID:MergeSams\tCL: args "+programArguments);
		if (addProgramArgs) for (String pg: programParams) al.add(pg);
		
		//add read group
		for (String rg: readGroups) al.add(rg);
		
		return Misc.stringArrayListToString(al, "\n");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MergeSams(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		File replacementHeaderFile = null;
		String useqVersion = IO.fetchUSeqVersion();
		programArguments = useqVersion+" "+Misc.stringArrayToString(args, " ");
		System.out.println("\n"+useqVersion+" Arguments: "+ Misc.stringArrayToString(args, " ") +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'd': forExtraction = new File(args[++i]); break;
					case 'f': saveBadReads = true; break;
					case 'v': verbose = false; break;
					case 's': saveFile = new File(args[++i]); break;
					case 't': deleteTempSamFile = false; break;
					case 'h': replacementHeaderFile = new File(args[++i]); break;
					case 'a': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'm': minimumMappingQualityScore = Float.parseFloat(args[++i]); break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}


		//pull files
		File[][] tot = new File[2][];
		tot[0] = IO.extractFiles(forExtraction,".sam.gz");
		tot[1] = IO.extractFiles(forExtraction,".bam");

		dataFiles = IO.collapseFileArray(tot);
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.bam or xxx.sam.gz file(s)!\n");

		//check save file
		if (saveFile == null) {
			saveFile = new File (forExtraction , "merged.bam");
		}
		else if (saveFile.getName().endsWith(".bam") == false) Misc.printErrAndExit("\nError: your indicated save file doesn't end in xxx.bam?\n");

		//load header?
		if (replacementHeaderFile != null){
			StringBuilder sb = new StringBuilder();
			String[] lines = IO.loadFile(replacementHeaderFile);
			if (lines.length ==0) Misc.printErrAndExit("\nError: replacement header contains no comment lines?\n");
			for (String l: lines) {
				sb.append(l);
				sb.append("\n");
			}
			//add program
			sb.append("@PG\tID:MergeSams\tCL: args "+programArguments);
			samHeader = sb.toString();
		}


		//print info
		System.out.println(maximumAlignmentScore+ "\tMaximum alignment score.");
		System.out.println(minimumMappingQualityScore+ "\tMinimum mapping quality score.");
		System.out.println(saveBadReads +"\tSave alignments failing filters.");

	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                 Merge Sams: May 2017                             **\n" +
				"**************************************************************************************\n" +
				"Merges sam and bam files. Adds a consensus header if one is not provided. These may\n"+
				"not work with GATK or Picard downstream apps, good for USeq.\n"+

				"\nOptions:\n"+
				"-d The full path to a directory containing xxx.bam or xxx.sam.gz files to merged.\n" +

				"\nDefault Options:\n"+
				"-s Save file, must end in xxx.bam, defaults merge.bam in -d.\n"+
				"-a Maximum alignment score. Defaults to 300, smaller numbers are more stringent.\n" +
				"      Approx 30pts per mismatch.\n"+
				"-m Minimum mapping quality score, defaults to 0 (no filtering), larger numbers are\n" +
				"      more stringent. Set to 13 or more to require near unique alignments. DO NOT set\n"+
				"      for alignments parsed by the SamTranscriptomeParser!\n"+
				"-f Save reads failing filters, defaults to tossing them.\n"+
				"-h Full path to a txt file containing a sam header, defaults to autogenerating the\n"+
				"      header from the sam/bam headers.\n"+
				"-t Don't delete temp xxx.sam.gz file.\n"+
				"-p Add program arguments to header, defaults to deleting, note duplicate cause Picard\n"+
				"      apps to fail.\n"+
				"-q Quiet, print only errors.\n"+

				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/MergeSams -f /Novo/Run7/\n" +
				"     -m 20 -a 120  \n\n" +

				"**************************************************************************************\n");

	}

	public int getNumberPassingAlignments() {
		return numberPassingAlignments;
	}	
}
