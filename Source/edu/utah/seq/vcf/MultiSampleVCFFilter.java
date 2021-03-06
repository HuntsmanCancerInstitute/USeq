package edu.utah.seq.vcf;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utah.seq.useq.data.RegionScoreText;
import util.bio.annotation.Bed;
import util.gen.IO;
import util.gen.Misc;
import util.gen.Num;

public class MultiSampleVCFFilter {

	//user defined fields
	private File vcfInFile;
	private File vcfOutFile;
	private String pathToTabix = null;
	private boolean filterAnySample = false;
	private boolean passing = true;
	private boolean compressOutput = false;
	private boolean filterRecordQuality = false;
	private boolean failNonPassRecords = false;
	private int sampleMinimumReadDepthDP = 0;
	private int sampleMaximumReadDepthDP = Integer.MAX_VALUE;
	private int sampleMinimumReadDepthAD = 0;
	private float sampleMinimumGenotypeQualityGQ = 0;
	private float recordMinimumQUAL = 0;
	private boolean printSampleNames = false;
	private boolean filterByGenotype = false;
	private boolean failNoCall = false;
	private HashMap<String,RegionScoreText[]> regions = null;
	private String[][] samplesByGroup = null;
	private String[] flagsByGroup = null;
	

	public MultiSampleVCFFilter(String[] args){
		long startTime = System.currentTimeMillis();

		processArgs(args);

		printOptions();
		
		//Count vcf records
		ArrayList<File> tempVcfFiles = new ArrayList<File>();
		int recordCount = VCFUtilities.countReads(vcfInFile);
		if (recordCount == 0) return;
		int chunks = recordCount / VCFUtilities.readsToChunk + 1;
		
		File tempDir = null;
		if (chunks > 1) {
			System.out.println("\nFile too big to intersect all in one go, splitting into " + chunks + " chunks\n");
			//make tempDirectory
			tempDir = new File("tempVCF");
			if (tempDir.exists()) {
				IO.deleteDirectory(tempDir);
			}
			tempDir.mkdir();
		}
				

		//for each file or file chunk
		System.out.println("File\tFilterType\tStarting#\tEnding#");
		for (int i=0;i<chunks;i++) {
			
			//create VCF parser
			VCFParser parser = null;
			File tempVcf = null;
			if (tempDir != null){
				tempVcf = new File(tempDir,"tempVcf_" + i + ".vcf");
				tempVcfFiles.add(tempVcf);
				parser =  new VCFParser(vcfInFile, true, true, true, i,VCFUtilities.readsToChunk);
			}
			else parser = new VCFParser(vcfInFile,true, true, true);
			parser.appendChrFixMT(); //this is ok since only the original record is printed out, not the modified.  
			
			String fileName = vcfInFile.getName();
			if (chunks !=1) fileName = i+"_"+fileName;
			
			//set all to pass? or keep incoming status?
			if (failNonPassRecords) {
				parser.setFilterFieldPeriodToTextOnAllRecords(VCFRecord.PASS);
				int pass  = parser.countMatchingVCFRecords(VCFRecord.PASS);
				int total = parser.getVcfRecords().length;
				System.out.println(fileName+ "\tPASSRecordFILTER\t"+total+"\t"+pass);
			}
			else parser.setFilterFieldOnAllRecords(VCFRecord.PASS);

			//filter on regions?
			if (regions != null){			
				int[] startEndCounts = parser.filterVCFRecords(regions);
				System.out.println(fileName+ "\tRegionFilter\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
			}
			
			if (filterRecordQuality){
				int[] startEndCounts = filterRecordQuality(parser);
				System.out.println(fileName+ "\tRecordQuality\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
			}
	
			if (filterAnySample) {
				int[] startEndCounts = filterAnySample(parser);
				System.out.println(fileName+ "\tAnySample\t"+startEndCounts[0]+"\t"+startEndCounts[1]);
			}
			
			if (filterByGenotype) {
				int[] startEndCounts = filterByGenotype(parser);
				System.out.println(fileName + "\tGenotypeFilter\t" + startEndCounts[0]+"\t"+startEndCounts[1]);
			}
	
			//print good or bad records 
			File out = vcfOutFile;
			if (tempDir != null) out = tempVcf;
			if (passing) parser.printFilteredRecords(out,VCFRecord.PASS);
			else parser.printFilteredRecords(out, VCFRecord.FAIL);
			
			parser = null;
		}
		
		//chunking?
		if (tempDir != null){
			//Merge vcf file
			VCFUtilities.mergeVcf(tempVcfFiles, vcfOutFile);	
			//delete temp files
			IO.deleteDirectory(tempDir);
		}

		if (compressOutput) VCFUtilities.createTabix(vcfOutFile, pathToTabix);

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private void printSampleNames() {
		System.out.println("File\tSampleNames");
		VCFParser parser = new VCFParser(vcfOutFile, false, false, false);
		System.out.println(vcfOutFile.getName()+ "\t"+ Misc.stringArrayToString(parser.getSampleNames(), ","));
	}

	/**Fails records where groups don't match their flags */
	private int[] filterByGenotype(VCFParser parser) {
		//fetch records
		VCFRecord[] records = parser.getVcfRecords();
		int startingRecordNumber = parser.countMatchingVCFRecords(VCFRecord.PASS);
		
		//filter
		for (VCFRecord test : records){
			//is it a passing record?
			if (test.getFilter().equals(VCFRecord.FAIL)) continue;
			
			boolean passFlag;
			boolean globalPass = true;
			
			for (int i=0; i<flagsByGroup.length; i++) {
				passFlag = true;
				
				int[] groupSamples = this.fetchGroupIndexes(parser, i);
				VCFSample[] samples = test.getSample();
				
				for (int j=0; j<groupSamples.length; j++) {
					VCFSample s = samples[groupSamples[j]];
					
					if (s.isNoCall()== true || s.getReadDepthDP() < sampleMinimumReadDepthDP || s.getReadDepthDP() > sampleMaximumReadDepthDP || s.getGenotypeQualityGQ() < sampleMinimumGenotypeQualityGQ) {
						if (failNoCall) {
							passFlag = false;
						} 
						continue;
					} 
					
					if (!checkFlag(s.getGenotypeGT(),this.flagsByGroup[i])) {
						passFlag = false;
					} 
				}
				
				if (passFlag != true) {
					globalPass = false;
				}
			}
			
			if (globalPass != true) {
				test.setFilter(VCFRecord.FAIL);
			}
			
		}

		int numStillPassing = parser.countMatchingVCFRecords(VCFRecord.PASS);
		return new int[]{startingRecordNumber, numStillPassing};
	}
	
	private boolean checkFlag(String genotype,String flag) {
		boolean retVal =  false;
		if (flag.equals("W") && genotype.equals("0/0")) {
			retVal = true;
		} else if (flag.equals("H") && genotype.equals("0/1")) {
			retVal = true;
		} else if (flag.equals("M") && genotype.equals("1/1")) {
			retVal = true;
		} else if (flag.equals("-W") && !genotype.equals("0/0")) {
			retVal = true;
		} else if (flag.equals("-H") && !genotype.equals("0/1")) {
			retVal = true;
		} else if (flag.equals("-M") && !genotype.equals("1/1")) {
			retVal = true;
		}
		
		return retVal;
	}

	
	private int[] fetchGroupIndexes(VCFParser parser, int group) {
		int[] indexes = new int[this.samplesByGroup[group].length];
		String[] sampleNames = parser.getSampleNames();
		for (int i=0; i < indexes.length; i++) {
			boolean found = false;
			for (int j=0; j<sampleNames.length; j++) {
				if (this.samplesByGroup[group][i].equals(sampleNames[j])) {
					indexes[i] = j;
					found = true;
					break;
				}
			}
			if (found == false) Misc.printErrAndExit("\nCannot find a matching sample for '" + this.samplesByGroup[group][i] + 
					"' in the list of sample names "+Misc.stringArrayToString(sampleNames, ",") +"\n");
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
				int dp = sample.getReadDepthDP();				
				int ad = 0;
				String adString = sample.getAlleleCount();
				if (adString != null) ad = Integer.parseInt(sample.getAlternateCounts());
				else {
					//look for info DP4, from SamTools, 
					//DP4=1,4,5,11 high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases
					String dp4 = test.getInfoObject().getInfo("DP4");
					if (dp4 != null){
						int[] baseCounts = Num.parseInts(dp4, Misc.COMMA);
						if (baseCounts.length == 4) {
							ad = baseCounts[2]+baseCounts[3];
							if (dp == -1) dp = ad + baseCounts[0]+baseCounts[1];
						}
					}
				}

				if (sample.isNoCall() == false && 
						dp >= sampleMinimumReadDepthDP && 
						dp <= sampleMaximumReadDepthDP &&
						sample.getGenotypeQualityGQ() >= sampleMinimumGenotypeQualityGQ &&
						ad >= sampleMinimumReadDepthAD) {
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

	private void printOptions() {
		System.out.println("Options:");
		System.out.println(failNonPassRecords + "\tFail records where the original FILTER field is not 'PASS' or '.'");
		System.out.println(filterRecordQuality+"\tFail records with QUAL scores < "+ recordMinimumQUAL);
		System.out.println(filterAnySample + "\tPass records where any sample passses thresholds");
		boolean regionFilter = (regions != null);
		System.out.println(regionFilter+"\tFail records that don't intersect the user provided regions");
		
		System.out.println(sampleMinimumReadDepthDP + "\tMinimum sample read depth DP");
		if (sampleMaximumReadDepthDP == Integer.MAX_VALUE) System.out.println("unlim\tMaximum sample read depth DP");
		else System.out.println(sampleMaximumReadDepthDP + "\tMaximum sample read depth DP");
		System.out.println(sampleMinimumReadDepthAD + "\tMinimum sample allele read depth AD or DP4");
		System.out.println(sampleMinimumGenotypeQualityGQ + "\tMinimum sample genotype quality GQ");
		
		if (filterByGenotype) {
			for (int i=0; i<flagsByGroup.length;i++) {
				System.out.println("\nGroup: " + i);
				System.out.println("Flag: " + flagsByGroup[i]);
				for (int j=0; j<samplesByGroup[i].length; j++) {
					System.out.println("\tSample: " + samplesByGroup[i][j]);
				}
			}
		}
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
		File regionFile = null;
		String[] sampleNames = null;
		String[] groups = null;

		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'v': inputFile = new File(args[++i]); break;
					case 'c': regionFile = new File(args[++i]); break;
					case 'p': outputFile = args[++i]; break;
					case 't': pathToTabix = args[++i]; break;
					case 'f': passing = false; break;
					case 'a': filterAnySample = true; break;
					case 'g': sampleMinimumGenotypeQualityGQ = Float.parseFloat(args[++i]); filterAnySample = true; break;
					case 'i': failNonPassRecords = true; break;
					case 'b': filterByGenotype = true; break;
					case 'r': sampleMinimumReadDepthDP = Integer.parseInt(args[++i]); filterAnySample = true; break;
					case 'x': sampleMaximumReadDepthDP = Integer.parseInt(args[++i]); filterAnySample = true; break;
					case 'y': sampleMinimumReadDepthAD = Integer.parseInt(args[++i]); filterAnySample = true;break;
					case 'd': recordMinimumQUAL = Float.parseFloat(args[++i]); filterRecordQuality = true; break;
					case 's': printSampleNames = true; break;
					case 'n': sampleNames = args[++i].split(","); break;
					case 'u': groups = args[++i].split(","); break;
					case 'l': flagsByGroup = args[++i].split(","); break;
					case 'e': failNoCall = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		
		//region file
		if (regionFile != null){
			regions = Bed.parseBedFile(regionFile, true, true);
		}
		
		//Check sample name / group / flag information
		//Make sure sample number == group numbers
		if (filterByGenotype) {
			if (sampleNames == null) {
				System.out.println("If filtering by genotype, you must specify sample names (-n)\n");
				System.exit(1);
			}
			if (groups == null) {
				System.out.println("If filtering by genotype, you must specify groups (-u)\n");
				System.exit(1);
			}
			if (this.flagsByGroup == null) {
				System.out.println("If filtering by genotype, you must specify flags by group (-l)\n");
				System.exit(1);
			}
			
			int total = 0;
			for (String group: groups) {
				total += Integer.parseInt(group);
			}
			if (total != sampleNames.length) {
				System.out.println("Error: The number of samples listed does not match the number of group assignments. Samples: " + sampleNames.length + " Groups: " + total + "\n");
				System.exit(1);
			}
			
			if (groups.length != flagsByGroup.length) {
				System.out.println("Error: The number of flags does not match the number of groups. Groups: " + groups.length + " Flags: " + flagsByGroup.length + "\n");
				System.exit(1);
			}
			
			
			//Check flag
			HashMap<String,Integer> flagChoice = new HashMap<String,Integer>();
			flagChoice.put("W", 0);
			flagChoice.put("H", 0);
			flagChoice.put("M", 0);
			flagChoice.put("-W", 0);
			flagChoice.put("-H", 0);
			flagChoice.put("-M", 0);
			
			for (String flag: flagsByGroup) {
				if (!flagChoice.containsKey(flag)) {
					System.out.println("Error: Don't recognize the specified flag: " + flag + "\n");
					System.exit(1);
				} else {
					int count = flagChoice.get(flag);
					flagChoice.put(flag, count++);
				}
			}
			
			for (Integer val: flagChoice.values()) {
				if (val > 1) {
					System.out.println("Error: Flag specified more than once\n");
					System.exit(1);
				}
			}
			
			int nameIndex = 0;
			int groupIndex = 0;
			this.samplesByGroup = new String[groups.length][];
			for (String group: groups) {
				int count = Integer.parseInt(group);
				this.samplesByGroup[groupIndex] = new String[count];
				int sampleIndex = 0;
				for (int i = nameIndex; i< nameIndex+count; i++) {
					this.samplesByGroup[groupIndex][sampleIndex] = sampleNames[i];
					sampleIndex += 1;
				}
				nameIndex += count;
				groupIndex += 1;
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
			if (pathToTabix != null) vcfInFile = VCFUtilities.unzipTabix(inputFile,pathToTabix);
			else vcfInFile = inputFile;
		} else {
			System.out.println("Input file does not appear to be a XXX.vcf/XXX.vcf.gz file.\n");
			System.exit(1);
		}


		if (outputFile == null) {
			String name = Misc.removeExtension(inputFile.getName());
			vcfOutFile = new File (inputFile.getParentFile(), name+"_Filt.vcf");
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
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                            Multi Sample VCF Filter  : July 2015                  **\n" +
				"**************************************************************************************\n" +
				"Filters a vcf file containing multiple sample records into those that pass or fail the\n" +
				"tests below. This works with VCFv4.1 files created by the GATK package. Note, the \n" +
				"records are not modified. If the number of records in the VCF file is greater than \n" +
				"500000, the VCF file is intersected in chunks.\n\n" +

				"Required:\n"+
				"-v Full path to a sorted single or multi sample vcf file (xxx.vcf/xxx.vcf.gz)). Note,\n"+
				"       Java often fails to parse tabix compressed vcf files.  Best to uncompress.\n\n"+
								
				"Optional:\n" +
				"-p Full path to an output VCF (xxx.vcf or xxx.vcf.gz).  Specifying xxx.vcf.gz will \n" + 
				"       compress and index the VCF using tabix (set -t too). Defaults to input_Filt.vcf\n" +
				"-f Print out failing records, defaults to printing those passing the filters.\n" +
				"-a Fail records where no sample passes the sample thresholds.\n"+
				"-i Fail records where the original FILTER field is not 'PASS' or '.'\n"+
				"-c Fail records that don't intersect the regions in this bed file, full path.\n"+
				"-b Filter by genotype flags.  -n, -u and -l must be set.\n" +
				"-n Sample names ordered by category.  \n" +
				"-u Number of samples in each category.  \n" +
				"-l Requirement flags for each category. All samples that pass the specfied filters \n" + 
				"       must meet the flag requirements, or the variant isn't reported.  At least one \n" +
				"       sample in each group must pass the specified filters, or the variant isn't \n" +
				"       reported.\n" +
				"		   a) 'W' : homozygous common \n" +
				"		   b) 'H' : heterozygous \n" +
				"		   c) 'M' : homozygous rare \n" +
				"		   d) '-W' : not homozygous common \n" +
				"	 	   e) '-H' : not heterozygous \n" +
				"		   f) '-M' : not homozygous rare\n" +
				"-e Strict genotype matching. If this is selected, records with no-call samples \n" +
				"       or samples falling below either minimum sample genotype quality (-g) or \n" +
				"       minimum sample read depth (-r) won't be reported. Only samples listed in (-n)\n" +
				"       will be checked \n"+
				"-d Minimum record QUAL score, defaults to 0, recommend >=20 \n"+
				"-g Minimum sample genotype quality GQ, defaults to 0, recommend >= 20 \n"+
				"-r Minimum sample read depth DP, defaults to 0, recommend >=10 \n"+
				"-x Maximum sample read depth DP, defaults to unlimited\n"+
				"-y Minimum sample allele count read depth AD or DP4, defaults to 0\n"+
				"-s Print sample names and exit.\n"+
				"-t Path to tabix.\n" +

				"\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/MultiSampleVCFFilter \n" +
				"       -v DEMO.passing.vcf -p DEMO.intersection.vcf -c exomeV4.bed -b \n" +
				"       -n SRR504516,SRR776598,SRR504515,SRR504517,SRR504483 -u 2,2,1 -l M,H,-M \n\n" +

				"**************************************************************************************\n");

	}

}
