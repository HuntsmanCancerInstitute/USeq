package edu.utah.seq.vcf;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.Gzipper;
import util.gen.IO;
import util.gen.Misc;

public class MultiSampleVCFFilter {

	//user defined fields
	private File[] vcfFiles;
	private boolean filterAnySample = false;
	private boolean controlHomozygousFilter = false;
	private boolean oneOrMorePassingCohortFilter = false;
	private boolean filterRecordQuality = false;
	private boolean requireOneObservationInCases = false;
	private int sampleMinimumReadDepthDP = 0;
	private int sampleMinimumGenotypeQualityGQ = 0;
	private float recordMinimumQUAL = 0;
	private boolean printSampleNames = false;
	private String[] controlSampleNames = null;
	private String[] cohortSampleNames = null;
	

	
	public MultiSampleVCFFilter(String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);

		printOptions();

		//for each file
		System.out.println("\nFile\tFilterType\tStarting#\tEnding#");
		for (int i=0; i< vcfFiles.length; i++){
			VCFParser parser = new VCFParser(vcfFiles[i], true, true);

			//set everything to pass (note this won't change the original when you print because printing grabs the original record line)
			parser.setFilterFieldOnAllRecords(VCFRecord.PASS);
			
			if (filterRecordQuality){
				int[] startEndCounts = filterRecordQuality(parser);
				System.out.println(vcfFiles[i].getName()+ "\tRecordQuality\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
			}

			if (filterAnySample) {
				int[] startEndCounts = filterAnySample(parser);
				System.out.println(vcfFiles[i].getName()+ "\tAnySample\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
			}

			if (controlHomozygousFilter){
				int[] startEndCounts = controlHomozygousFilter(parser);
				System.out.println(vcfFiles[i].getName()+ "\tAnyControlHomozygous\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
			}

			if (oneOrMorePassingCohortFilter){
				int[] startEndCounts = filterAnySampleCohort(parser);
				System.out.println(vcfFiles[i].getName()+ "\tAnyCohortSample\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
			}
			
			if (requireOneObservationInCases) {
				int[] startEndCounts = requireOneObservationFilter(parser);
				System.out.println(vcfFiles[i].getName()+ "\tAtLeastOneObservationAboveThresholds\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
			}
			
			

			//print good and bad records 
			parser.printRecords(VCFRecord.PASS);
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private void printSampleNames() {
		System.out.println("File\tSampleNames");
		for (int i=0; i< vcfFiles.length; i++){
			VCFParser parser = new VCFParser(vcfFiles[i], false, false);
			System.out.println(vcfFiles[i].getName()+ "\t"+ Misc.stringArrayToString(parser.getSampleNames(), ","));
		}
	}

	/**Sets passing records to fail if any of the control samples that pass the read depth and genotype quality and are also homozygous for the non reference allele. 
	 * Returns int[startingNumPassing, endingNumPassing]*/
	private int[] controlHomozygousFilter(VCFParser parser) {
		//fetch indexes for controls
		int[] controlSampleIndex = fetchControlIndexes(parser);
		//fetch records
		VCFRecord[] records = parser.getVcfRecords();
		int startingRecordNumber = parser.countMatchingVCFRecords(VCFRecord.PASS);

		//filter
		for (VCFRecord test : records){
			//is it a passing record?
			if (test.getFilter().equals(VCFRecord.FAIL)) continue;

			boolean passes = true;
			VCFSample[] samples = test.getSample();
			for (int i=0; i< controlSampleIndex.length; i++){
				VCFSample s = samples[controlSampleIndex[i]];
				if (s.isNoCall()== true ||  
						s.getReadDepthDP() < sampleMinimumReadDepthDP || 
						s.getGenotypeQualityGQ() < sampleMinimumGenotypeQualityGQ) continue;
				//check if homozygous non reference
				if (s.getGenotypeGT().equals("1/1")){
					passes = false;
					break;
				}
			}
			if (passes == false) test.setFilter(VCFRecord.FAIL);
		}

		int numStillPassing = parser.countMatchingVCFRecords(VCFRecord.PASS);
		return new int[]{startingRecordNumber, numStillPassing};
	}
	
	
	private int[] requireOneObservationFilter(VCFParser parser) {
		//get indexes for cases
		int[] cohortSampleIndex = this.fetchCohortIndexes(parser);
		
		VCFRecord[] records = parser.getVcfRecords();
		int startingRecordNumber = parser.countMatchingVCFRecords(VCFRecord.PASS);
		
		//Run filtering
		for (VCFRecord test: records) {
			//Make sure it's a passing record
			if (test.getFilter().equals(VCFRecord.FAIL)) {
				continue;
			}
			
			boolean passes = false;
			VCFSample[] samples = test.getSample();
			for (int i=0; i<cohortSampleIndex.length; i++) {
				VCFSample s = samples[cohortSampleIndex[i]];
				if ((s.isNoCall())) {
					continue;
				}
				if ((s.getGenotypeGT().equals("0/1") || s.getGenotypeGT().equals("1/1")) && 
						s.getReadDepthDP() >= sampleMinimumReadDepthDP &&
						s.getGenotypeQualityGQ() >= sampleMinimumGenotypeQualityGQ) {
					passes = true;
					break;
				}
			}
			if (passes == false) {
				test.setFilter(VCFRecord.FAIL);
			}
		}
		
		int numStillPassing = parser.countMatchingVCFRecords(VCFRecord.PASS);
		return new int[] {startingRecordNumber, numStillPassing};
	}


	private int[] fetchControlIndexes(VCFParser parser) {
		int[] indexes = new int[controlSampleNames.length];
		String[] sampleNames = parser.getSampleNames();
		for (int i=0; i< indexes.length; i++){
			boolean found = false;
			for (int j=0; j< sampleNames.length; j++){
				if (controlSampleNames[i].equals(sampleNames[j])){
					indexes[i] = j;
					found = true;
					break;
				}
			}
			if (found == false) Misc.printErrAndExit("\nCannot find a matching control sample name for '"+
					controlSampleNames[i]+"' in the list of sample names : "+Misc.stringArrayToString(sampleNames, ",") +"\n");
		}
		return indexes;
	}


	/**Sets passing records to fail if no cohort/ affected sample makes the read depth and genotype quality thresholds. 
	 * Returns int[startingNumPassing, endingNumPassing]*/
	private int[] filterAnySampleCohort(VCFParser parser) {
		//fetch indexes for controls
		int[] cohortSampleIndex = fetchCohortIndexes(parser);
		//fetch records
		VCFRecord[] records = parser.getVcfRecords();
		int startingRecordNumber = parser.countMatchingVCFRecords(VCFRecord.PASS);

		//filter
		for (VCFRecord test : records){
			//is it a passing record?
			if (test.getFilter().equals(VCFRecord.FAIL)) continue;

			boolean passes = false;
			VCFSample[] samples = test.getSample();
			for (int i=0; i< cohortSampleIndex.length; i++){
				VCFSample s = samples[i];
				if (s.isNoCall()== false &&  
						s.getReadDepthDP() >= sampleMinimumReadDepthDP && 
						s.getGenotypeQualityGQ() >= sampleMinimumGenotypeQualityGQ) {
					passes = true;
					break;
				}
			}
			if (passes == false) test.setFilter(VCFRecord.FAIL);
		}

		int numStillPassing = parser.countMatchingVCFRecords(VCFRecord.PASS);
		return new int[]{startingRecordNumber, numStillPassing};
	}


	private int[] fetchCohortIndexes(VCFParser parser) {
		int[] indexes = new int[cohortSampleNames.length];
		String[] sampleNames = parser.getSampleNames();
		for (int i=0; i< indexes.length; i++){
			boolean found = false;
			for (int j=0; j< sampleNames.length; j++){
				if (cohortSampleNames[i].equals(sampleNames[j])){
					indexes[i] = j;
					found = true;
					break;
				}
			}
			if (found == false) Misc.printErrAndExit("\nCannot find a matching cohort sample name for '"+
					cohortSampleNames[i]+"' in the list of sample names : "+Misc.stringArrayToString(sampleNames, ",") +"\n");
		}
		return indexes;
	}



	/**Sets passing records to fail if no sample passes the sample read depth and genotype quality. Returns int[startingNumPassing, endingNumPassing]*/
	private int[] filterAnySample(VCFParser parser) {
		//fetch records
		VCFRecord[] records = parser.getVcfRecords();
		int startingRecordNumber = parser.countMatchingVCFRecords(VCFRecord.PASS);
		//filter
		for (VCFRecord test : records){
			//is it a passing record?
			if (test.getFilter().equals(VCFRecord.FAIL)) continue;
			boolean passes = false;
			for (VCFSample sample: test.getSample()){
				//is it a passing record?
				if (sample.isNoCall() == false && 
						sample.getReadDepthDP() >= sampleMinimumReadDepthDP && 
						sample.getGenotypeQualityGQ() >= sampleMinimumGenotypeQualityGQ) {
					passes = true;
					break;
				}
			}
			if (passes == false) test.setFilter(VCFRecord.FAIL);
		}
		int numStillPassing = parser.countMatchingVCFRecords(VCFRecord.PASS);
		return new int[]{startingRecordNumber, numStillPassing};
	}
	
	/**Sets passing records to fail if no sample passes the sample read depth and genotype quality. Returns int[startingNumPassing, endingNumPassing]*/
	private int[] filterRecordQuality(VCFParser parser) {
		//fetch records
		VCFRecord[] records = parser.getVcfRecords();
		int startingRecordNumber = parser.countMatchingVCFRecords(VCFRecord.PASS);
		//filter
		for (VCFRecord test : records){
			//is it a passing record?
			if (test.getFilter().equals(VCFRecord.FAIL)) continue;
			if (test.getQuality() < recordMinimumQUAL) test.setFilter(VCFRecord.FAIL);
		}
		int numStillPassing = parser.countMatchingVCFRecords(VCFRecord.PASS);
		return new int[]{startingRecordNumber, numStillPassing};
	}



	private void printOptions() {
		System.out.println("Options:");
		System.out.println(filterRecordQuality+"\tFail records with QUAL scores < "+ recordMinimumQUAL);
		System.out.println(filterAnySample + "\tPass records where any sample passses thresholds");
		System.out.println(controlHomozygousFilter + "\tFail records where any control is homozygous non reference and also passes the sample thresholds");
		System.out.println(oneOrMorePassingCohortFilter + "\tFail records where no cohort/ affected sample passes the sample thresholds");
		if (requireOneObservationInCases) {
			System.out.println(requireOneObservationInCases+"\tFail records where no cohort sample is homozygous or heterozgous for alt allele above the sample thresholds");
		}
		if (controlHomozygousFilter) System.out.println(Misc.stringArrayToString(controlSampleNames, ",")+ "\tControl sample names");
		if (oneOrMorePassingCohortFilter) System.out.println(Misc.stringArrayToString(cohortSampleNames, ",")+ "\tCohort/ affected sample names");
		System.out.println(sampleMinimumReadDepthDP + "\tMinimum sample read depth DP");
		System.out.println(sampleMinimumGenotypeQualityGQ + "\tMinimum sample genotype quality GQ");
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MultiSampleVCFFilter(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': forExtraction = new File(args[++i]); break;
					case 'a': filterAnySample = true; break;
					case 'b': controlHomozygousFilter = true; break;
					case 'c': oneOrMorePassingCohortFilter = true; break;
					case 'e': requireOneObservationInCases = true; break;
					case 'g': sampleMinimumGenotypeQualityGQ = Integer.parseInt(args[++i]); break;
					case 'r': sampleMinimumReadDepthDP = Integer.parseInt(args[++i]); break;
					case 'd': recordMinimumQUAL = Float.parseFloat(args[++i]); filterRecordQuality = true; break;
					case 's': printSampleNames = true; break;
					case 'o': cohortSampleNames = args[++i].split(","); break;
					case 'n': controlSampleNames = args[++i].split(","); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//pull files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please indicate a vcf file to filter.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		vcfFiles = IO.collapseFileArray(tot);
		if (vcfFiles == null || vcfFiles.length ==0 || vcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz) file(s)!\n");

		if (printSampleNames) {
			printSampleNames();
			System.out.println();
			System.exit(0);
		}

		if (controlHomozygousFilter){
			if (controlSampleNames == null || controlSampleNames.length ==0) Misc.printExit("\nError: please enter a comma delimited list (no spaces) of the control sample names.\n");
		}

		if (oneOrMorePassingCohortFilter){
			if (cohortSampleNames == null || cohortSampleNames.length ==0) Misc.printExit("\nError: please enter a comma delimited list (no spaces) of the cohort/ affected sample names.\n");
		}
		
		
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Multi Sample VCF Filter  : March 2013                   **\n" +
				"**************************************************************************************\n" +
				"Splits vcf file(s) containing multiple sample records into those that pass and fail\n" +
				"the tests below. This works with VCFv4.1 files created by the GATK package. Note, the\n" +
				"records are not modified. There is an incompatibility with the Tabix gzip function and\n" +
				"java gzip reader on Linux.  If you find premature termination of your vcf file try\n" +
				"uncompressing the file.\n\n" +

				"Options:\n"+
				"-v Full path to a sorted multi sample vcf file or directory containing such\n" +
				"      (xxx.vcf(.gz/.zip OK)). \n"+
				"-a Fail records where no sample passes the sample thresholds.\n"+
				"-b Fail records where any of the control samples that pass the sample thresholds also\n" +
				"      contain the homozygous non reference allele. Requires setting -n .\n"+
				"-c Fail records where none of the cohort/ affected samples pass the sample thresholds.\n" +
				"      Requires setting -o .\n"+
				"-e Fail records where none of the cohort/ affected samples have alt observations\n" +
				"      above the specified GQ and DP threshholds.\n" +
				"-d Minimum record QUAL score, defaults to 0, recommend >=20 .\n"+
				"-g Minimum sample genotype quality GQ, defaults to 0, recommend >= 20 .\n"+
				"-r Minimum sample read depth DP, defaults to 0, recommend >=10 .\n"+
				"-o Comma delimited (no spaces) list of cohort/ affected sample names.\n"+
				"-n Comma delimited (no spaces) list of control sample names.\n"+
				"-s Print sample names and exit.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/MultiSampleVCFFilter -b -c -g 20 -r 10\n" +
				"     -d 20 -n norm5,norm6,norm7  -o cancer1,cancer2,cancer3,cancer4\n\n"+

		"**************************************************************************************\n");

	}

}
