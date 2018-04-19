package edu.utah.seq.vcf.sim;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.data.SamSubsampler;
import edu.utah.seq.parsers.MergeSams;
import edu.utah.seq.parsers.SamAlignmentExtractor;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

/**Generates mixes of alignments for tumor normal simulations with BamBlaster
 * @author Nix
 * */
public class BamMixer {

	//fields
	private File pairedVariantBamFile;
	private File unPairedVariantBamFile;
	private File unModifiedMatchingBamFile;
	private File unModifiedNoMatchBamFile;
	private File saveDirectory;
	private double[] fractions = {0.025, 0.05, 0.1, 0.2};
	private boolean verbose = false;
	
	//internal fields
	private SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
	private MergeSams variantMS;
	private int[] countForVar = new int[fractions.length];
	private int[] countForUnMod = new int[fractions.length];
	private File mergedVariantBam;
	private File modifiedPairedVariantBamFile;
	private File unmodifiedPairedVariantBamFile; 
	private static final Pattern BB = Pattern.compile(":BB");
	private Pattern numUnder = Pattern.compile("^\\d[\\d_]+");
	private Pattern trailingBB = Pattern.compile(":BB-");
	
	public BamMixer (String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);
		
		doWork();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" min\n");
	}
	
	public void doWork(){
		
		System.out.println("Splitting paired alignment file into modified and unmodified records...");
		int[] modUnMod = splitPairedVariantByMod();
		System.out.println("\t"+modUnMod[0]+"\t# modified alignments");
		System.out.println("\t"+modUnMod[1]+"\t# unmodified mate alignments");
		System.out.println("\t"+modUnMod[2]+"\t# not aligned or secondary");
		
		System.out.println("\nFiltering, merging, and sorting paired and unpaired variant aligments...");
		mergePairedAndUnPaired();
		
		System.out.println("\nCounting and subsampling bams...\n\tFracVar\t#Var\t#NonVar");
		
		//find lower number of passing alignments
		int numGoodAlign = countPassingAlignments(unModifiedMatchingBamFile);
		if (numGoodAlign > variantMS.getNumberPassingAlignments()) numGoodAlign = variantMS.getNumberPassingAlignments();
		
		//calc actual # of alignments needed
		calculateAlignmentNumbersForFractions(numGoodAlign);
		
		//subsample
		File[] subVar = subsampleBam(countForVar, mergedVariantBam);
		File[] subNorm = subsampleBam(countForUnMod, unModifiedMatchingBamFile);
		
		//merge
		System.out.println("\nMerging and sorting subsampled bams with filtered bam...");
		mergeFinals(subVar, subNorm);
	}

	private File[] subsampleBam(int[] alignmentCounts, File bam) {
		SamSubsampler ss = new SamSubsampler(bam, verbose);
		File[] sub = new File[alignmentCounts.length];
		String name = Misc.removeExtension(bam.getName());
		for (int i=0; i< alignmentCounts.length; i++){
			sub[i] = new File(saveDirectory, alignmentCounts[i]+"_"+name+".temp.sam.gz");
			sub[i].deleteOnExit();
			ss.printRandomAlignments(alignmentCounts[i], sub[i]);
		}
		return sub;
	}

	private void calculateAlignmentNumbersForFractions(double numGoodAlign) {
		countForVar = new int[fractions.length];
		countForUnMod = new int[fractions.length];
		for (int i=0; i< fractions.length; i++){
			if (fractions[i] > 1) Misc.printErrAndExit("\nError: one of your fractions is greater than 1.0, aborting.\n");
			countForVar[i] = (int)Math.round(fractions[i] * numGoodAlign);
			countForUnMod[i] = (int)Math.round((1-fractions[i]) * numGoodAlign);
			System.out.println("\t"+fractions[i]+ "\t" +countForVar[i]+ "\t" +countForUnMod[i]);
		}
	}
	
	/**Removes leading and trailing info from BamMixer
	 * e.g. 0_HWI-D00294:322:CATY4ANXX:6:1101:1166:34209:BB-INS_178536299_C_CT_2  ->  HWI-D00294:322:CATY4ANXX:6:1101:1166:34209*/
	public void stripNameNumber(SAMRecord sam){
		String name = sam.getReadName();
		int start = 0; 
		int stop = name.length();
		
		Matcher mat = numUnder.matcher(name);
		if (mat.find()) start = mat.end();

		Matcher trailing = trailingBB.matcher(name);
		if (trailing.find()) stop = trailing.start();
		
		String newName = name.substring(start, stop);
		sam.setReadName(newName);
	}


	public int countPassingAlignments(File bamSam) {
		SamReader reader = readerFactory.open(bamSam);
		int numPassingAlignments = 0;
		SAMRecordIterator i = reader.iterator();
		while (i.hasNext()) {
			SAMRecord sam = i.next();
			if (passThresholds(sam) == false) continue;
			numPassingAlignments++;
		}
		i.close();
		return numPassingAlignments;
	}
	
	public static boolean passThresholds(SAMRecord sam){
		//any problematic flags?
		if (sam.getReadUnmappedFlag() || sam.getNotPrimaryAlignmentFlag() || sam.getReadFailsVendorQualityCheckFlag()) return false;
		return true;
	}

	/*The paired alignment file contains alignments with the variant and a matched mate that may or may not also be modified so need to split these.
	 * returns numMod and numUnMod*/
	private int[] splitPairedVariantByMod(){
		int numMod = 0;
		int numUnMod = 0;
		int numNotAligned = 0;
		SamReader reader = readerFactory.open(pairedVariantBamFile);
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		writerFactory.setTempDirectory(saveDirectory);
		modifiedPairedVariantBamFile = new File (saveDirectory, "modPair.temp.bam");
		modifiedPairedVariantBamFile.deleteOnExit();
		SAMFileWriter modWriter = writerFactory.makeBAMWriter(reader.getFileHeader(), false, modifiedPairedVariantBamFile);
		
		unmodifiedPairedVariantBamFile = new File (saveDirectory, "unModPair.temp.bam");
		unmodifiedPairedVariantBamFile.deleteOnExit();
		SAMFileWriter unmodWriter = writerFactory.makeBAMWriter(reader.getFileHeader(), false, unmodifiedPairedVariantBamFile);

		//Iterate through the sam file 
		SAMRecordIterator i = reader.iterator();
		String first = "_1";
		String second = "_2";
		String both = "_12";
		String toSearch = null;
		while (i.hasNext()) {
			SAMRecord sam = i.next();
			//check it
			if (sam.isSecondaryOrSupplementary() || sam.getReadUnmappedFlag() || sam.getProperPairFlag() == false) {
				numNotAligned++;
				continue;
			}
			
			String name = sam.getReadName();
			String[] tokens = BB.split(name);
			if (tokens.length == 1) Misc.printErrAndExit("\nERROR: a ':BB' delimiter was not found in the read name for this modified alignment:\n"+sam.getSAMString());
			
			//is this the first or second read
			if (sam.getFirstOfPairFlag()) toSearch = first;
			else toSearch = second;
			boolean isModified = false;
			//scan for trailing _12, _1, _2
			for (int x=1; x< tokens.length; x++){
				if (tokens[x].endsWith(both) || tokens[x].endsWith(toSearch)){
					isModified = true;
					break;
				}
			}
			
			//strip out leading and trailing header info
			stripNameNumber(sam);
			
			//write out
			if (isModified) {
				modWriter.addAlignment(sam);
				numMod++;
			}
			else {
				unmodWriter.addAlignment(sam);
				numUnMod++;
			}
		}
		
		//close the IO
		try {
			reader.close();
			i.close();
			modWriter.close();
			unmodWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
			Misc.printErrAndExit("\nERROR: problem closing readers or writers when spliting the paired alignments\n");
		}
		
		return new int[] {numMod, numUnMod, numNotAligned};
	}
	
	private void mergePairedAndUnPaired() {
		File[] toMerge = {modifiedPairedVariantBamFile, unPairedVariantBamFile};
		mergedVariantBam = new File(saveDirectory, "mergedVariants.bam");
		variantMS = new MergeSams(toMerge, mergedVariantBam, verbose, false);
		System.out.println("\t"+variantMS.getNumberPassingAlignments()+"\ttotal # modified alignments");
	}
	
	private void mergeFinals(File[] subVar, File[] subNorm){
		//for each fraction
		for (int i=0; i< fractions.length; i++){
			File f = new File (saveDirectory, fractions[i]+"Var.bam");
			File[] toMerge = {unModifiedNoMatchBamFile, unmodifiedPairedVariantBamFile, subVar[i], subNorm[i]};
			System.out.println("\t"+subVar[i].getName()+"\t"+subNorm[i].getName()+"\t"+unModifiedNoMatchBamFile.getName()+"\t"+unPairedVariantBamFile.getName());
			new MergeSams(toMerge, f, verbose, false);
		}
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BamMixer(args);
	}

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		String frac = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'p': pairedVariantBamFile = new File(args[++i]); break;
					case 's': unPairedVariantBamFile = new File(args[++i]); break;
					case 'u': unModifiedMatchingBamFile = new File(args[++i]); break;
					case 'f': unModifiedNoMatchBamFile = new File(args[++i]); break;
					case 'r': saveDirectory = new File(args[++i]); break;
					case 'v': verbose = true; break;
					case 'm': frac = args[++i]; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//diff fractions than defaults?
		if (frac != null) fractions = Num.parseDoubles(Misc.COMMA.split(frac));
		
		//check bams
		if (pairedVariantBamFile == null || pairedVariantBamFile.canRead() == false) Misc.printErrAndExit("Error: please proved a sam/bam alignment file containing paired alignments from aligning the paired fastq.gz files from your BamBlaster run.\n");
		if (unPairedVariantBamFile == null || unPairedVariantBamFile.canRead() == false) Misc.printErrAndExit("Error: please proved a sam/bam alignment file containing unpaired alignments from aligning the unpaired fastq.gz file from your BamBlaster run.\n");
		if (unModifiedMatchingBamFile == null || unModifiedMatchingBamFile.canRead() == false) Misc.printErrAndExit("Error: please proved a path to the xxx_unmodified.bam from your BamBlaster run.\n");
		if (unModifiedNoMatchBamFile == null || unModifiedNoMatchBamFile.canRead() == false) Misc.printErrAndExit("Error: please proved a path to the xxx_filtered.bam from your BamBlaster run.\n");
		SamAlignmentExtractor.lookForBaiIndexes(new File[]{pairedVariantBamFile,unPairedVariantBamFile,unModifiedMatchingBamFile,unModifiedNoMatchBamFile}, true);
		
		//check save directory
		if (saveDirectory == null) Misc.printErrAndExit("Error: please provide a path to a directory for saving the results.");
		saveDirectory.mkdirs();
		
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                               Bam Mixer : June 2014                              **\n" +
				"**************************************************************************************\n" +
				"Combines bam alignment files in different fractions to simulate multiple variant\n"+
				"frequencies. Run BamBlaster first.\n\n"+

				"Required:\n"+
				"-r Path to a directory to save the results\n" +
				"-u Path to the xxx_unmodified.bam from your BamBlaster run\n"+
				"-f Path to the xxx_filtered.bam from your BamBlaster run\n"+
				"-p Path to your realigned paired end bam\n"+
				"-s Path to your realigned single end bam\n"+
				
				"\nOptional:\n"+
				"-m Fractions to mix in the variant alignments, comma delimited, no spaces, defaults to\n"+
				"     0.025,0.05,0.1,0.2\n"+
				"-v Verbose output.\n"+

				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/BamMixer -r ~/TumorSim/\n"+
				"    -u ~/bb_unmodified.bam -f ~/bb_filtered.bam -p ~/bb_paired.bam -s ~/bb_single.bam \n\n" +

				"**************************************************************************************\n");
	}
}
