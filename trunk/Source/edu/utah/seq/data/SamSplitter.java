
package edu.utah.seq.data;
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


/**Splits and sorts a sam file randomly into two parts.
 * @author david.nix@hci.utah.edu 
 **/
public class SamSplitter{
	//fields
	private File samFile;
	private String samHeader = null;
	private float maximumAlignmentScore = 240;
	private float minimumMappingQualityScore = 0;
	private int numberAlignments = 0;
	private int numberUnmapped = 0;
	private int numberExceedingEnd = 0;
	private int numberFailingVendorQC = 0;
	private int numberPassingAlignments = 0;
	private int numberFailingAlignmentScore = 0;
	private int numberFailingMappingQualityScore = 0;
	private int numberAdapter = 0;
	private int numberPhiX = 0;

	private Gzipper samOut1;
	private Gzipper samOut2;
	private String programArguments;
	private SamReaderFactory factory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX).validationStringency(ValidationStringency.SILENT);
	private Random random = new Random();
	private int numberOne = 0;
	private int numberTwo = 0;
	
	//constructors
	public SamSplitter(String[] args){
		try {
			long startTime = System.currentTimeMillis();
			processArgs(args);

			doWork();

			//finish and calc run time
			double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
			System.out.println("\nDone! "+Math.round(diffTime)+" Min\n");

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void doWork() throws IOException{
		File saveDir = samFile.getParentFile();
		String samName = Misc.removeExtension(samFile.getName());
		//make print writer
		File temp1 = new File(saveDir, samName+"_1.sam.gz");
		samOut1 = new Gzipper(temp1);
		File temp2 = new File(saveDir, samName+"_2.sam.gz");
		samOut2 = new Gzipper(temp2);
		
		//add header
		SamReader sr = factory.open(samFile);
		samHeader = sr.getFileHeader().getTextHeader().trim();
		samHeader = samHeader+"\n"+ "@PG\tID:SplitSam\tCL: args "+programArguments;
		samOut1.println(samHeader);
		samOut2.println(samHeader);
		
		parseSam(sr);
		
		sr.close();
		samOut1.close();
		samOut2.close();

		//stats
		double fractionPassing = ((double)numberPassingAlignments)/((double)numberAlignments);
		System.out.println("\nStats:\n");
		System.out.println("\t"+numberAlignments+"\tTotal # Alignments from raw sam file");
		System.out.println("\t"+numberPassingAlignments+"\tAlignments passing filters ("+Num.formatPercentOneFraction(fractionPassing)+")");
		System.out.println("\t\t"+numberUnmapped+"\t# Unmapped Reads");
		System.out.println("\t\t"+numberFailingVendorQC+"\t# Alignments failing vendor/ platform QC and or malformed");
		System.out.println("\t\t"+numberFailingAlignmentScore+"\t# Alignments failing alignment score");
		System.out.println("\t\t"+numberFailingMappingQualityScore+"\t# Alignments failing mapping quality score");
		System.out.println("\t\t"+numberExceedingEnd+"\t# Alignments exceeding the length of the chromosome in the header");
		System.out.println("\t\t"+numberAdapter+"\t# Adapter alignments");
		System.out.println("\t\t"+numberPhiX+"\t# PhiX alignments");
		System.out.println("\t\t"+numberOne+"\t# Split1 Alignments");
		System.out.println("\t\t"+numberTwo+"\t# Split2 alignments");
		System.out.println();

	}

	public boolean parseSam(SamReader sr){
		SAMRecord sam = null;
		try {
			SAMRecordIterator it = sr.iterator();
			int dotCounter = 0;
			ArrayList<SAMRecord> al = new ArrayList<SAMRecord>();
			String lastName = null;
			while (it.hasNext()) {
				if (++dotCounter > 1000000){
					System.out.print(".");
					dotCounter = 0;
				}
					sam = it.next();
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
					if (alignmentScore != Integer.MIN_VALUE){
						if (alignmentScore > maximumAlignmentScore){
							numberFailingAlignmentScore++;
							continue;
						}
					}
					int mappingQuality = sam.getMappingQuality();
					if (mappingQuality < minimumMappingQualityScore){
						numberFailingMappingQualityScore++;
						continue;
					}

					//OK, it passes, increment counter
					numberPassingAlignments++;
					
					//first passing record?
					if (lastName == null) lastName = sam.getReadName();
					
					//same as last name
					String currName = sam.getReadName();
					if (currName.equals(lastName)) al.add(sam);
					else {
						//print old
						print(al);
						//reset for new
						lastName = currName;
						al.clear();
						al.add(sam);
					}
			}
			//print last
			print(al);
			
			sr.close();

		} catch (Exception e) {
			if (sam != null) System.err.println("Problem processing -> "+sam.getSAMString());
			e.printStackTrace();

			return false;
		}
		return true;
	}

	

	private void print(ArrayList<SAMRecord> al) throws IOException{
		if (random.nextBoolean()) {
			for (SAMRecord sam : al) samOut1.println(sam.getSAMString().trim());
			numberOne+= al.size();
		}
		else {
			for (SAMRecord sam : al) samOut2.println(sam.getSAMString().trim());
			numberTwo+= al.size();
		}
		
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new SamSplitter(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
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
					case 's': samFile = new File(args[++i]); break;
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

		//pull file
		if (samFile == null) Misc.printErrAndExit("\nError: please enter a sam or bam file to split.\n");
		String name = samFile.getName();
		if (name.endsWith("bam") == false && name.endsWith("sam.gz") == false && name.endsWith("sam") == false) {
			Misc.printErrAndExit("\nError: please enter a sam or bam file to split.\n");
		}


	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                SamSplitter: Oct 2014                             **\n" +
				"**************************************************************************************\n" +
				"Randomly splits a sam or bam file in ~1/2. To maintain pairs, sort by queryname.\n"+

				"\nOptions:\n"+
				"-s The full path to a queryname sorted xxx.bam or xxx.sam (.gz OK) file.\n" +

				"\nDefault Options:\n"+
				"-a Maximum alignment score. Defaults to 240, smaller numbers are more stringent.\n" +
				"      Approx 30pts per mismatch.\n"+
				"-m Minimum mapping quality score, defaults to 0 (no filtering), larger numbers are\n" +
				"      more stringent. Set to 13 or more to require near unique alignments. DO NOT set\n"+
				"      for alignments parsed by the SamTranscriptomeParser!\n"+

				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/SamSplitter -f /Novo/Run7/exome.bam\n" +
				"     -m 20 -a 120  \n\n" +

				"**************************************************************************************\n");

	}	
}
