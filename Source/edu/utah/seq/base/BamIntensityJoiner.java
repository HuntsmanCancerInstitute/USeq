package edu.utah.seq.base;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord.SAMTagAndValue;

import util.gen.*;

/**Takes a sam/ bam alignment file and the bam intensity file from running the modified IlluminaBasecallsToSam app and joins them.
 * @author Nix
 * */
public class BamIntensityJoiner {

	//user fields
	private File intensityFile;
	private File alignmentFile;
	private File resultsFile;
	private float minimumMappingQuality = 20;
	private float maximumAlignmentScore = 240;
	private String mdField = null;
	private int subSample = 0;
	
	//internal
	private int numberAlignments = 0;
	private int numberUnmapped = 0;
	private int numberFailingVendorQC = 0;
	private int numberAdapter = 0;
	private int numberPhiX = 0;
	private int numberFailingAlignmentScore = 0;
	private int numberFailingMappingQualityScore = 0;
	private int numberPassingAlignments = 0 ;
	private int numberFailingMDField = 0;

	private SAMRecordIterator alignmentInterator;
	private SAMRecordIterator intensityInterator;
	private SAMFileWriter outputSam;


	//constructor
	public BamIntensityJoiner(String[] args){
		long startTime = System.currentTimeMillis();
		
		processArgs(args);
		
		makeIO();

		System.out.print("Parsing");
		parseFiles();
		
		closeIO();

		//stats
		double fractionPassing = ((double)numberPassingAlignments)/((double)numberAlignments);
		System.out.println("\nStats (some flags aren't set so be suspicious of zero read catagories):");
		System.out.println("\t"+numberAlignments+"\tTotal # Alignments");
		System.out.println("\t"+numberPassingAlignments+"\tAlignments passing filters ("+Num.formatPercentOneFraction(fractionPassing)+")");
		System.out.println("\t\t"+numberUnmapped+"\t# Unmapped");
		System.out.println("\t\t"+numberFailingVendorQC+"\t# Failing vendor/ platform QC");
		System.out.println("\t\t"+numberAdapter+"\t# Aligning to the adapters");
		System.out.println("\t\t"+numberPhiX+"\t# Aligning to phiX");
		System.out.println("\t\t"+numberFailingAlignmentScore+"\t# Failing to pass the alignment score -> "+(int)maximumAlignmentScore);
		System.out.println("\t\t"+numberFailingMappingQualityScore+"\t# Failing to pass the mapping quality score -> "+(int)minimumMappingQuality);
		if (mdField != null)System.out.println("\t\t"+numberFailingMDField+"\t# Failing MD field -> "+mdField);

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	public void makeIO(){
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader alignmentReader = new SAMFileReader(alignmentFile);
		alignmentInterator = alignmentReader.iterator();
		SAMFileReader intensityReader = new SAMFileReader(intensityFile);
		intensityInterator = intensityReader.iterator();
		outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(alignmentReader.getFileHeader(), true, resultsFile);
	}
	
	public void closeIO(){
		outputSam.close();
		alignmentInterator.close();
		intensityInterator.close();
	}
	
	public boolean parseFiles(){
		try {
			//for each alignment
			int dotCounter = 0;
			int counter = 0;
			while (alignmentInterator.hasNext()) {
				if (++dotCounter > 1000000){
					System.out.print(".");
					dotCounter = 0;
				}
				
				SAMRecord samAlignment = alignmentInterator.next();
				numberAlignments++;

				//is it aligned?
				if (samAlignment.getReadUnmappedFlag()){
					numberUnmapped++;
					continue;
				}
				//does it pass the vendor qc?
				if (samAlignment.getReadFailsVendorQualityCheckFlag()){
					numberFailingVendorQC++;
					continue;
				}

				//skip phiX and adapter
				if (samAlignment.getReferenceName().startsWith("chrPhiX")){
					numberPhiX++;
					continue;
				}
				if (samAlignment.getReferenceName().startsWith("chrAdapt")){
					numberAdapter++;
					continue;
				}

				//does it pass the score thresholds?
				List<SAMTagAndValue> attributes = samAlignment.getAttributes();
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
				int mappingQuality = samAlignment.getMappingQuality();
				if (mappingQuality < minimumMappingQuality){
					numberFailingMappingQualityScore++;
					failedScore = true;
				}
				
				//check MD field?
				if (mdField != null){
					String md = (String) samAlignment.getAttribute("MD");
					if (md != null){
						if (md.equals(mdField) == false) {
							numberFailingMDField++;
							continue;
						}
					}
				}
				
				if (failedScore) continue;
				
				//increment counter
				numberPassingAlignments++;
				
				//print a subsampling of the data?
				if (subSample !=0){
					if (counter == subSample){
						counter = 0;
					}
					else {
						counter++;
						continue;
					}
				}
				
				//fetch intensities sam object
				String alignmentReadName = samAlignment.getReadName();
				boolean alignmentFirstPair = samAlignment.getFirstOfPairFlag();
				SAMRecord samIntensity = null;
				while (intensityInterator.hasNext()){
					SAMRecord test = intensityInterator.next();
					//same name?
					if (test.getReadName().equals(alignmentReadName)){
						//same read?
						if (alignmentFirstPair == test.getFirstOfPairFlag()){
							samIntensity = test;
							break;
						}
					}
				}
				
				if (samIntensity == null) Misc.printErrAndExit("\nError: no intensity sam record found matching alignment sam record.  Did you synchronized the alignment output to match the fastq input?\n");
				
				//set intensities
				samAlignment.setAttribute("IA", samIntensity.getSignedShortArrayAttribute("IA"));
				samAlignment.setAttribute("IC", samIntensity.getSignedShortArrayAttribute("IC"));
				samAlignment.setAttribute("IG", samIntensity.getSignedShortArrayAttribute("IG"));
				samAlignment.setAttribute("IT", samIntensity.getSignedShortArrayAttribute("IT"));
				
				//write it
				outputSam.addAlignment(samAlignment);

			}
			System.out.println();
		} catch (Exception e) {
			System.err.println("\nError parsing sam files.\n");
			e.printStackTrace();
			return false;
		} 
		return true;
	}








	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BamIntensityJoiner(args);
	}		


	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'i': intensityFile = new File(args[++i]); break;
					case 'a': alignmentFile = new File(args[++i]); break;
					case 'r': resultsFile = new File(args[++i]); break;
					case 'q': minimumMappingQuality = Float.parseFloat(args[++i]); break;
					case 's': maximumAlignmentScore = Float.parseFloat(args[++i]); break;
					case 'm': mdField = args[++i]; break;
					case 'u': subSample = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}


	}	
	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Bam Intensity Joiner : July 2013                        **\n" +
				"**************************************************************************************\n" +
				"Extracts base level intensity information from the output of modified Picard\n" +
				"IlluminaBaseCallsToSam app and inserts this into an alignment file. Be sure to \n" +
				"syncronize the alignment output (e.g. -oSync in novoalign) so it is in the same order\n" +
				"as the intensity data.\n\n" +

				"Options:\n"+
				"-a Full path to sam/bam alignment file with header.\n"+
				"-i Full path to bam intensity file from running the modified IlluminaBasecallsToSam.\n"+
				"-r Full path bam file for saving the merged results.\n"+
				"-q Minimum mapping quality score. Defaults to 20, bigger numbers are more stringent.\n" +
				"      This is a phred-scaled posterior probability that the mapping position of read\n" +
				"      is incorrect. For RNA-Seq data from the SamTranscriptomeParser, set this to 0.\n" +
				"-s Maximum alignment score. Defaults to 240, smaller numbers are more stringent.\n"+
				"-m Filter for particular MD fields.\n"+
				"-u Sub sample data, printing only every XXX alignment.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/BamIntensityJoiner -u 10000 -m 101 -a\n" +
				"      /Alignments/8341X.sam.gz -i /Ints/8341X.bam -r /Merged/8341.bam\n\n" +

		"**************************************************************************************\n");

	}


}
