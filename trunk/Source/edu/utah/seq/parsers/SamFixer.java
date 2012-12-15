
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import edu.utah.seq.data.*;
import edu.utah.seq.useq.data.Region;
import edu.utah.seq.useq.data.RegionScoreText;
import util.gen.*;
import java.util.*;
import util.bio.annotation.*;
import edu.utah.seq.data.sam.*;

/**Parses, filters, merges, and fixes Sam alignment txt files
 * @author david.nix@hci.utah.edu 
 **/
public class SamFixer{
	//fields
	private File[] dataFiles;
	private File samFile;
	private Pattern spliceJunction = Pattern.compile("(.+)_(\\d+)_(\\d+)");
	private Pattern spliceJunctionLine = Pattern.compile("^@SQ\\s+SN:\\w+_\\d+_\\d+[GATCgatc]*\\s.+");
	private Pattern toSkip = Pattern.compile("^[@#]");
	private float minimumMappingQuality = 0;
	private float maximumAlignmentScore = 1000;
	private int numberAlignments = 0;
	private int numberUnmapped = 0;
	private int numberFailingVendorQC = 0;
	private int numberAdapter = 0;
	private int numberFailingScore = 0;
	private int numberSpliceJunctions = 0;
	private int numberFailedSpliceJunctions = 0;
	private int numberPassingAlignments = 0;
	private String adapterName = "chrAdapt";
	private boolean rnaSeqData = false;
	private int spliceJunctionRadius = 0;
	private PrintWriter out = null;
	private String argString = "";
	private boolean stripMDField = true;
	private boolean removeUnMapped = false;
	private boolean removePoorQualityReads = true;
	private BufferedReader workingBufferedReader;

	//constructors
	public SamFixer(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			System.out.println("Minimum mapping quality threshold = "+minimumMappingQuality);
			System.out.println("Maximum alignment score = "+maximumAlignmentScore);
			System.out.println("Remove MD: fields = "+stripMDField);
			System.out.println("Remove poor quality reads = "+removePoorQualityReads);
			System.out.println("Remove unmapped reads = "+removeUnMapped);
			System.out.println("RNA-Seq data present = "+rnaSeqData);
			if (rnaSeqData && spliceJunctionRadius !=0) System.out.println("Splice junction radius = "+spliceJunctionRadius);

			System.out.println("\nParsing and filtering...");

			//make writer
			out = new PrintWriter ( new FileWriter (samFile));

			//collect headers and print
			printHeader();

			//for each file, parse and write to disk	
			for (int i=0; i< dataFiles.length; i++){
				System.out.println("\t"+dataFiles[i]);
				parseSamFile(dataFiles[i]); 
			}

			//close writer
			out.close();

			//stats
			double fractionPassing = ((double)numberPassingAlignments)/((double)numberAlignments);
			System.out.println("\nStats:");
			System.out.println("\t"+numberAlignments+"\tTotal # Alignments");
			System.out.println("\t"+numberPassingAlignments+"\tAlignments passing filters ("+Num.formatPercentOneFraction(fractionPassing)+")");
			System.out.println("\t\t"+numberUnmapped+"\t# Unmapped");
			System.out.println("\t\t"+numberFailingVendorQC+"\t# Failing vendor/ platform QC");
			System.out.println("\t\t"+numberAdapter+"\t# Aligning to the adapters");
			System.out.println("\t\t"+numberFailingScore+"\t# Failing to pass the mapping quality score ("+minimumMappingQuality+") or the alignment score ("+maximumAlignmentScore+").");
			if (rnaSeqData) {
				System.out.println("\t\t"+numberSpliceJunctions+"\t# converted splice junctions");
				System.out.println("\t\t"+numberFailedSpliceJunctions+"\t# non-converted splice junctions");
			}

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
			System.out.println("\nDone! "+Math.round(diffTime)+" seconds");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void printHeader(){
		Pattern numberChrom = Pattern.compile(".+SN:chr.+\\s.+");
		try {
			LinkedHashSet<String> header = new LinkedHashSet<String>();
			String pg = null;
			String rg = null;
			String hd = null;
			//for each file, parse header and load into hash	
			for (int i=0; i< dataFiles.length; i++){
				BufferedReader in  = IO.fetchBufferedReader(dataFiles[i]);
				String line;
				//for each line in the start of the file
				while ((line=in.readLine())!= null) {
					line = line.trim();
					if (line.length() == 0 || line.startsWith("#")) continue;
					if (line.startsWith("@")) {
						//skip splice junction chromosomes
						if (spliceJunctionRadius !=0 && spliceJunctionLine.matcher(line).matches()) continue;
						//skip multiple occurrences of @PG
						if (line.startsWith("@PG")) {
							pg = line;
							continue;
						}
						//read group? skip those that are the same
						if (line.startsWith("@RG")) {
							if (rg == null) rg = line;
							else if (rg.equals(line) == false) rg = rg + "\n"+line;
							continue;
						}
						//sort line?
						if (line.startsWith("@HD")){
							hd = line;
							continue;
						}
						//look for chromosomes with just numbers
						if (line.startsWith("@SQ")){
							Matcher mat = numberChrom.matcher(line);
							if (mat.matches()==false){
								System.err.println("\t\tWARNING: chromosome name does not start with 'chr', this is incompatible with DAS/2 data distribution. Correct before uploading into GenoPub -> "+line);
							}
						}
						//add line
						header.add(line);
					}
					else break;
				}
				//close in
				in.close();
			}
			
			//write out header
			
			//missing sorted?
			if (hd == null) out.println("@HD\tVN:1.0\tSO:unsorted");
			else out.println(hd);
			
			Iterator<String> it = header.iterator();
			while (it.hasNext()) out.println(it.next());

			//add PGs
			if (pg != null) out.println(pg);
			//out.println("@PG\tID:SamFixer\tVN:"+version+"\tCL:java -jar SamFixer "+argString);

			//add dummy read group if not present
			if (rg == null) out.println("@RG\tID:unknownReadGroup\tPG:SamFixer\tSM:unknownSample");
			else out.println(rg);

		} catch (Exception e) {
			System.err.println("\nError parsing SAM/BAM file header\n");
			e.printStackTrace();
		} 
	}


	public void parseSamFile(File file){
		try {
			BufferedReader in = IO.fetchBufferedReader(file);
			String line;
			//for each line in the file
			SamAlignment prior = null;
			while ((line=in.readLine())!= null) {
				line = line.trim();
				if (line.length() == 0 || toSkip.matcher(line).find()) continue;

				numberAlignments++;
				SamAlignment sa;
				try {
					sa = new SamAlignment(line, false);
				} catch (MalformedSamAlignmentException e) {
					System.err.println("Skipping, malformed-> "+line);
					continue;
				}

				//does it pass the vendor qc?
				if (removePoorQualityReads && sa.failedQC()) {
					numberFailingVendorQC++;
					continue;
				}

				//is it hitting the adapter chromosome?
				String chromosome = sa.getReferenceSequence();
				if (chromosome.contains(adapterName)) {
					numberAdapter++;
					continue;
				}

				//does it pass the scores threshold?
				int alignmentScore = sa.getAlignmentScore();
				if (alignmentScore != Integer.MIN_VALUE){
					if (alignmentScore > maximumAlignmentScore){
						numberFailingScore++;
						continue;
					}
				}

				//was likelihood score of being derived from the given map location set?
				float mappingQuality = sa.getMappingQuality();
				if (mappingQuality < minimumMappingQuality){
					numberFailingScore++;
					continue;
				}

				//is it a splice junction?
				boolean printIt = true;
				if (rnaSeqData){
					Matcher mat = spliceJunction.matcher(chromosome);
					boolean spliceJunctionPresent = mat.matches();
					if (spliceJunctionPresent){
						//paired? if so then must abort! 
						if (sa.isPartOfAPairedAlignment()) Misc.printErrAndExit("\nError: cannot fix paired RNA-Seq data containing splices. \n");
						//attempt to correct the splice junction to match genomic coordinates
						if (spliceJunctionRadius != 0) {
							printIt = sa.checkAndConvertSpliceJunction(spliceJunctionRadius);
							if (printIt){
								numberSpliceJunctions++;
							}
							else {
								numberFailedSpliceJunctions++;
								//conversion failed so make it unmapped if not part of a pair
								if (sa.isPartOfAPairedAlignment() == false ) {
									sa.convertToUnmapped();
									printIt = true;
								}
								//part of a pair so skip alignment
								else {
									printIt = false;
								}
							}
						}
					}
				}

				//write out sam alignment
				if (printIt) {
					//is it aligned?
					if (sa.isUnmapped()) {
						numberUnmapped++;
						if (removeUnMapped == false) sa.fixUnMappedNovo();
						else continue;
					}
					//fix mate info
					if (sa.isMateUnMapped()) sa.setUnMappedMate();
					if (stripMDField) out.println(sa.toStringNoMDField());
					else out.println(sa);
					//increment counter
					numberPassingAlignments++;
				}
			}
			in.close();
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing SAM/BAM file\n");
		} 

	}
	
	//For block parsing
	/*
	public SamAlignment[] fetchSamAlignmentBlock(){
		String line;
		SamAlignment sa = null;
		try {
			while ((line=workingBufferedReader.readLine())!= null) {
				line = line.trim();
				if (line.length() == 0 || toSkip.matcher(line).find()) continue;
				sa = null;
				try {
					sa = new SamAlignment(line);
				} catch (MalformedSamAlignmentException e) {
					System.err.println("Skipping, malformed-> "+line);
					continue;
				}
				return sa;
			}
		} catch (IOException e){
			e.printStackTrace();
		}
		return null;
	}

	
	public boolean processSamAlignment (SamAlignment sa){
		try {
			numberAlignments++;
			//does it pass the vendor qc?
			if (removePoorQualityReads && sa.failedQC()) {
				numberFailingVendorQC++;
				return false;
			}

			//is it hitting the adapter chromosome?
			String chromosome = sa.getReferenceSequence();
			if (chromosome.contains(adapterName)) {
				numberAdapter++;
				return false;
			}

			//does it pass the scores threshold?
			int alignmentScore = sa.getAlignmentScore();
			if (alignmentScore != Integer.MIN_VALUE){
				if (alignmentScore > maximumAlignmentScore){
					numberFailingScore++;
					return false;
				}
			}

			//was likelihood score of being derived from the given map location set?
			float mappingQuality = sa.getMappingQuality();
			if (mappingQuality < minimumMappingQuality){
				numberFailingScore++;
				return false;
			}

			//is it a splice junction?
			if (rnaSeqData){
				Matcher mat = spliceJunction.matcher(chromosome);
				boolean spliceJunctionPresent = mat.matches();
				if (spliceJunctionPresent){
					//attempt to correct the splice junction to match genomic coordinates
					if (spliceJunctionRadius != 0) {
						if (sa.checkAndConvertSpliceJunction(spliceJunctionRadius)) numberSpliceJunctions++;
						else {
							//conversion failed so make it unmapped
							sa.convertToUnmapped();
							numberFailedSpliceJunctions++;
						}
					}
				}
			}

			//is it aligned?
			if (sa.isUnmapped()) {
				numberUnmapped++;
				sa.fixUnMappedNovo();
				here review
				if (removeUnMapped == false) out.println(sa);
				return false;
			}
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
		return true;
	}

	public void parseSamFile2(File file){
		try {
			workingBufferedReader = IO.fetchBufferedReader(file);

			while (true) {

				SamAlignment saOne = fetchNextSamAlignment();
				if (saOne == null) return;


				//is it part of a pair?
				if (sa.isPartOfAPairedAlignment()){
					//read in second alignment
					SamAlignment saTwo = fetchNextSamAlignment();
				}


				//write out sam alignment
				if (printIt) {
					//fix mate info
					if (sa.isMateUnMapped()) sa.fixUnMappedMateNovo();
					if (stripMDField) out.println(sa.toStringNoMDField());
					else out.println(sa);
					//increment counter
					numberPassingAlignments++;
				}
			}
			workingBufferedReader.close();
		} catch (Exception e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nError parsing SAM/BAM file\n");
		} 

	}*/

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SamFixer(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		argString = Misc.stringArrayToString(args, " ");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+ argString +"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': forExtraction = new File(args[i+1]); i++; break;
					case 's': samFile = new File (args[i+1]); i++; break;
					case 'm': minimumMappingQuality = Float.parseFloat(args[++i]); break;
					case 'a': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'c': spliceJunctionRadius = Integer.parseInt(args[++i]); rnaSeqData = true; break;
					case 'u': removeUnMapped = true; break;
					case 'q': removePoorQualityReads = false; break;
					case 'd': stripMDField = false; break;
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//pull files
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".sam");
		tot[1] = IO.extractFiles(forExtraction,".sam.gz");
		tot[2] = IO.extractFiles(forExtraction,".sam.zip");

		dataFiles = IO.collapseFileArray(tot);
		if (dataFiles == null || dataFiles.length==0) dataFiles = IO.extractFiles(forExtraction);
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.sam(.zip/.gz) file(s)!\n");
		if (samFile == null)  Misc.printExit("\nError: please enter a full path sam output file (e.g. -s /save/sam/here/myFixedSamFile.sam \n");

		//rna seq data?
		if (rnaSeqData){
			if (spliceJunctionRadius < 1) Misc.printErrAndExit("\nError: please enter a splice junction radius used in making the splice " +
			"juction fasta file. This might be present in the sam header. Otherwise contact the person who ran the alignments.\n");
		}

	}	

	public static final String version = "August 2011";

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                 Sam Fixer: August 2011                           **\n" +
				"**************************************************************************************\n" +
				"Parses, filters, merges, and fixes xxx.sam files.\n"+

				"\nOptions:\n"+
				"-f The full path file or directory containing xxx.sam(.gz/.zip OK) file(s). Multiple \n" +
				"      files will be merged.\n" +
				"-s Full path file name for saving the fixed sam file.\n"+

				"\nDefault Options:\n"+
				"-m Minimum mapping quality score. Defaults to 0, bigger numbers are more stringent.\n" +
				"      This is a phred-scaled posterior probability that the mapping position of read\n" +
				"      is incorrect.\n" +
				"-a Maximum alignment score. Defaults to 1000, smaller numbers are more stringent.\n"+
				"-d Don't strip optional MD fields from alignments, defaults to removing these.\n"+
				"-u Remove unmapped reads.\n"+
				"-q Don't remove poor quality reads.\n"+
				"-c Convert splice-junctions to genomic coordinates, by providing a splice junction\n" +
				"      radius. Only works for single read RNA-Seq data where a splice junction fasta\n" +
				"      file was included in the alignments from the USeq MakeSpliceJunctionFasta app.\n" +
				"      This does NOT work for paired RNA-Seq data.\n"+


				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/SamParser -f /Novo/Run7/\n" +
				"     -m 20 -a 120 -s /Novo/Run7/mergedFixed.sam  -c 46 -u\n\n" +

		"**************************************************************************************\n");

	}

	public int getNumberPassingAlignments() {
		return numberPassingAlignments;
	}	

}
