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
	private File vcfInFile;
	private File vcfOutFile;
	private String pathToTabix = "/tomato/app/tabix/";
	private boolean filterAnySample = false;
	private boolean passing = true;
	private boolean compressOutput = true;
	private boolean controlHomozygousFilter = false;
	private boolean oneOrMorePassingCohortFilter = false;
	private boolean filterRecordQuality = false;
	private boolean filterRecordVQSLOD = false;
	private boolean failNonPassRecords = false;
	private boolean requireOneObservationInCases = false;
	private int sampleMinimumReadDepthDP = 0;
	private int sampleMinimumGenotypeQualityGQ = 0;
	private float recordMinimumQUAL = 0;
	private float recordMinimumVQSLOD = 0;
	private boolean printSampleNames = false;
	private String[] controlSampleNames = null;
	private String[] cohortSampleNames = null;




	public MultiSampleVCFFilter(String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);

		printOptions();

		//for each file
		System.out.println("\nFile\tFilterType\tStarting#\tEnding#");

		VCFParser parser = new VCFParser(vcfInFile, true, true, true);
		
		//set everything to pass (note this won't change the original when you print because printing grabs the original record line)?
		if (failNonPassRecords) {
			parser.setFilterFieldPeriodToTextOnAllRecords(VCFRecord.PASS);
			int pass  = parser.countMatchingVCFRecords(VCFRecord.PASS);
			int total = parser.getVcfRecords().length;
			System.out.println(vcfInFile.getName()+ "\tPASSRecordFILTER\t"+total+"\t"+pass);
		}
		else parser.setFilterFieldOnAllRecords(VCFRecord.PASS);

		if (filterRecordQuality){
			int[] startEndCounts = filterRecordQuality(parser);
			System.out.println(vcfInFile.getName()+ "\tRecordQuality\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
		}

		if (filterRecordVQSLOD){
			int[] startEndCounts = filterRecordVQSLOD(parser);
			System.out.println(vcfInFile.getName()+ "\tRecordVQSLOD\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
		}

		if (filterAnySample) {
			int[] startEndCounts = filterAnySample(parser);
			System.out.println(vcfInFile.getName()+ "\tAnySample\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
		}

		if (controlHomozygousFilter){
			int[] startEndCounts = controlHomozygousFilter(parser);
			System.out.println(vcfInFile.getName()+ "\tAnyControlHomozygous\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
		}

		if (oneOrMorePassingCohortFilter){
			int[] startEndCounts = filterAnySampleCohort(parser);
			System.out.println(vcfInFile.getName()+ "\tAnyCohortSample\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
		}

		if (requireOneObservationInCases) {
			int[] startEndCounts = requireOneObservationFilter(parser);
			System.out.println(vcfInFile.getName()+ "\tAtLeastOneObservationAboveThresholds\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
		}



		//print good or bad records 
		if (this.passing) {
			parser.printFilteredRecords(this.vcfOutFile,VCFRecord.PASS);
		} else {
			parser.printFilteredRecords(this.vcfOutFile, VCFRecord.FAIL);
		}

		if (this.compressOutput) {
			VCFUtilities.createTabix(this.vcfOutFile, this.pathToTabix);
		}

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private void printSampleNames() {
		System.out.println("File\tSampleNames");

		VCFParser parser = new VCFParser(vcfOutFile, false, false, false);
		System.out.println(vcfOutFile.getName()+ "\t"+ Misc.stringArrayToString(parser.getSampleNames(), ","));
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

	/**Sets passing records to fail if the don't pass the record quality. Returns int[startingNumPassing, endingNumPassing]*/
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

	/**Sets passing records to fail if no sample passes the sample read depth and genotype quality. Returns int[startingNumPassing, endingNumPassing]*/
	private int[] filterRecordVQSLOD(VCFParser parser) {
		try {
			//fetch records
			VCFRecord[] records = parser.getVcfRecords();
			int startingRecordNumber = parser.countMatchingVCFRecords(VCFRecord.PASS);
			//filter
			for (VCFRecord test : records){
				//is it a passing record?
				if (test.getFilter().equals(VCFRecord.FAIL)) continue;
				float score = test.getInfoObject().getInfoFloat("VQSLOD");
				if (score < recordMinimumVQSLOD) test.setFilter(VCFRecord.FAIL);
			}
			int numStillPassing = parser.countMatchingVCFRecords(VCFRecord.PASS);
			return new int[]{startingRecordNumber, numStillPassing};
			
		} catch (Exception e) {
			System.err.println("\nProblem parsing VQSLOD from INFO? Was the GATK ApplyRecalibration run on your vcf file?\n");
			e.printStackTrace();
			System.exit(0);
		}
		return null;
	}



	private void printOptions() {
		System.out.println("Options:");
		System.out.println(failNonPassRecords + "\tFail records where the original FILTER field is not 'PASS' or '.'");
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

		File inputFile = null;
		String outputFile = null;


		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': inputFile = new File(args[++i]); break;
					case 'p': outputFile = args[++i]; break;
					case 't': this.pathToTabix = args[++i]; break;
					case 'f': passing = false; break;
					case 'a': filterAnySample = true; break;
					case 'b': controlHomozygousFilter = true; break;
					case 'c': oneOrMorePassingCohortFilter = true; break;
					case 'e': requireOneObservationInCases = true; break;
					case 'g': sampleMinimumGenotypeQualityGQ = Integer.parseInt(args[++i]); break;
					case 'i': failNonPassRecords = true; break;
					case 'r': sampleMinimumReadDepthDP = Integer.parseInt(args[++i]); break;
					case 'd': recordMinimumQUAL = Float.parseFloat(args[++i]); filterRecordQuality = true; break;
					case 'q': recordMinimumVQSLOD = Float.parseFloat(args[++i]); filterRecordVQSLOD = true; break;
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

		if (inputFile == null) {
			System.out.println("Input file was not specified, exiting.\n");
			System.exit(1);
		} else if (!inputFile.exists()) {
			System.out.println("Input file does not exist, exiting.\n");
			System.exit(1);
		} else if (Pattern.matches(".+?.vcf",inputFile.getName())) {
			this.vcfInFile = inputFile;
		} else if (Pattern.matches(".+?.vcf.gz",inputFile.getName())) {
			this.vcfInFile = VCFUtilities.unzipTabix(inputFile,this.pathToTabix);
		} else {
			System.out.println("Input file does not appear to be a XXX.vcf/XXX.vcf.gz file.\n");
			System.exit(1);
		}


		if (outputFile == null) {
			System.out.println("Output file was no specified, exiting.\n");
			System.exit(1);
		} else if (Pattern.matches(".+?.vcf",outputFile)) {
			this.compressOutput = false;
			this.vcfOutFile = new File(outputFile);
		} else if (Pattern.matches(".+?.vcf.gz",outputFile)) {
			File vcfOutComp = new File(outputFile);
			if (vcfOutComp.exists()) {
				System.out.println("Tabix won't overwrite an existing file, rename the output or delete exisiting file.\n");
				System.exit(1);
			}

			this.vcfOutFile = new File(outputFile.substring(0,outputFile.length()-3));

			this.compressOutput = true;
		} else {
			System.out.println("Output file does not appear to be a XXX.vcf/XXX.vcf.gz file.\n");
			System.exit(1);
		}

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
				"**                          Multi Sample VCF Filter  : April 2013                   **\n" +
				"**************************************************************************************\n" +
				"Filters a vcf file containing multiple sample records into those that pass or fail\n" +
				"the tests below. This works with VCFv4.1 files created by the GATK package. Note, the\n" +
				"records are not modified. \n\n" +

				"Required:\n"+
				"-v Full path to a sorted multi sample vcf file (xxx.vcf/xxx.vcf.gz)). \n"+
				"-p Full path to the output VCF (xxx.vcf/xxx.vcf.gz).  Specifying xxx.vcf.gz will\n" +
				"      compress and index the VCF using tabix (set -t too).\n\n" +
				
				"Optional:\n" +
				"-f Print out failing records, defaults to printing those passing the filters.\n" +
				"-a Fail records where no sample passes the sample thresholds.\n"+
				"-b Fail records where any of the control samples that pass the sample thresholds also\n" +
				"      contain the homozygous non reference allele. Requires setting -n .\n"+
				"-c Fail records where none of the cohort/ affected samples pass the sample thresholds.\n" +
				"      Requires setting -o .\n"+
				"-e Fail records where none of the cohort/ affected samples have alt observations\n" +
				"      above the specified GQ and DP threshholds.\n" +
				"-i Fail records where the original FILTER field is not 'PASS' or '.'\n"+
				"-d Minimum record QUAL score, defaults to 0, recommend >=20 .\n"+
				"-g Minimum sample genotype quality GQ, defaults to 0, recommend >= 20 .\n"+
				"-r Minimum sample read depth DP, defaults to 0, recommend >=10 .\n"+
				"-q Minimum VQSLOD, defaults to 0, recommend >=4 . Requires that the GATK\n" +
				"      ApplyRecalibration was run on the vcf file.\n"+
				"-o Comma delimited (no spaces) list of cohort/ affected sample names.\n"+
				"-n Comma delimited (no spaces) list of control sample names.\n"+
				"-s Print sample names and exit.\n"+
				"-t Path to tabix\n" +

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/MultiSampleVCFFilter -v 9901R.vcf -b -c\n" +
				"     -i -g 20 -r 10 -o 9901R_filtered.vcf -d 20 -n norm5,norm6,norm7 \n" +
				"     -o cancer1,cancer2,cancer3,cancer4\n\n"+

		"**************************************************************************************\n");

	}

}
