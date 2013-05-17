package edu.utah.seq.vcf;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.analysis.OverdispersedRegionScanSeqs;
import edu.utah.seq.useq.data.RegionScoreText;

import util.bio.annotation.Bed;
import util.bio.annotation.ExportIntergenicRegions;
import util.bio.seq.Seq;
import util.gen.*;

/**Compares variant lists, uses the vcf QUAL score to filter.
 * @author Nix
 * */
public class VCFComparator {

	//user fields
	private File vcfKey;
	private File bedKey;
	private File[] testVcfFiles;
	private File vcfTest;
	private File bedTest;
	private File saveDirectory = null;
	private boolean requireGenotypeMatch = false;
	private boolean removeSNPs = false;
	private boolean removeNonSNPs = false;
	private boolean removeNonPass = false;
	private boolean useVQSLOD = false;

	private HashMap<String,RegionScoreText[]> keyRegions = null;
	private HashMap<String,RegionScoreText[]> testRegions = null;
	private HashMap<String,RegionScoreText[]> commonRegions = null;
	private VCFParser keyParser = null;
	private VCFParser testParser;
	private VCFMatch[] testMatchingVCF;
	private VCFRecord[] testNonMatchingVCF;
	private VCFRecord[] keyNonMatchingVCF;
	private StringBuilder results;
	private String options;
	private long keyBps;
	private long testBps;
	private long commonBps;
	private int numberUnfilteredKeyVariants;
	private ArrayList<Float> tprAL = new ArrayList<Float>();
	private ArrayList<Float> fdrAL = new ArrayList<Float>();
	private ArrayList<ScoredCalls> scoredCallsAL = new ArrayList<ScoredCalls>();

	//constructor
	public VCFComparator(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);

		for (int i=0; i< testVcfFiles.length; i++){
			vcfTest = testVcfFiles[i];
			results = new StringBuilder();

			//add options to results?
			printOptions();
			if (saveDirectory != null) results.append(options);

			//parse vcf file
			System.out.println("Parsing and filtering variant data for common interrogated regions...");
			parseFilterFiles();

			//if useVQSLOD
			if (useVQSLOD){
				try {
					testParser.setRecordVQSLODAsScore();
				} catch (Exception e) {
					System.err.println("\nProblem parsing VQSLOD from INFO? Was the GATK ApplyRecalibration run on your vcf file?\n");
					e.printStackTrace();
					System.exit(0);
				}
			}
			else {
				//set record QUAL score as thresholding score for roc curve data
				testParser.setRecordQUALAsScore();
			}

			//compare calls in common interrogated regions
			System.out.println("Comparing calls...");
			thresholdAndCompareCalls();
			
			//make scoredCalls
			scoredCallsAL.add(new ScoredCalls(Misc.removeExtension(vcfTest.getName())));

			//call after comparing!
			if (saveDirectory != null) printParsedDatasets();

		}
		
		//print out composite scoredCalls?
		if (saveDirectory !=null && scoredCallsAL.size() >1) printScoredCalls();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private void printScoredCalls() {
		ScoredCalls[] sc = new ScoredCalls[scoredCallsAL.size()];
		scoredCallsAL.toArray(sc);
		
		//find max row length
		int maxLength = 0;
		for (ScoredCalls s: sc){
			if (s.fdr.length > maxLength) maxLength = s.fdr.length;
		}
		
		//make writer
		String filter = "All.xls";
		if (removeSNPs) filter = "NonSNP.xls";
		else if (removeNonSNPs) filter = "SNP.xls";
		File f = new File (saveDirectory, "fdrTprSummary"+filter);
		PrintWriter out;
		
		try {
			out = new PrintWriter( new FileWriter(f));
			
			//print header
			out.print("dFDR\t");
			out.print(sc[0].name);
			for (int i=1; i< sc.length; i++){
				out.print("\tdFDR\t");
				out.print(sc[i].name);
			}
			out.println();
			
			//print first row with 1 for fdr
			out.print("1.0\t");
			out.print(sc[0].tpr[0]);
			for (int i=1; i< sc.length; i++){
				out.print("\t1.0\t");
				out.print(sc[i].tpr[0]);  
			}
			out.println();

			//print data
			//for each row
			for (int i=0; i< maxLength; i++){
				
				//print first sample
				if (i>= sc[0].fdr.length) out.print("\t");
				else {
					out.print(sc[0].fdr[i]);
					out.print("\t");
					out.print(sc[0].tpr[i]);
				}
				//for each sample
				for (int j=1; j< sc.length; j++){
					//past length?
					if (i>= sc[j].fdr.length) out.print("\t\t");
					else {
						out.print("\t");
						out.print(sc[j].fdr[i]);
						out.print("\t");
						out.print(sc[j].tpr[i]);
					}
				}
				out.println();
			}

			
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
	}

	public void thresholdAndCompareCalls(){
		//starting totals
		float totalKey = keyParser.getVcfRecords().length;		
		
		//clear old results
		tprAL.clear();
		fdrAL.clear();

		//intersect and split test into matching and non matching
		intersectVCF();

		//sort by score smallest to largest
		Arrays.sort(testNonMatchingVCF, new ComparatorVCFRecordScore());
		Arrays.sort(testMatchingVCF);

		results.append("QUALThreshold\tNumMatchTest\tNumNonMatchTest\tFDR=nonMatchTest/(matchTest+nonMatchTest)\tdecreasingFDR\tTPR=matchTest/totalKey\tFPR=nonMatchTest/totalKey\tPPV=matchTest/(matchTest+nonMatchTest)\n");

		//first do without thresholds
		float oldScore = testNonMatchingVCF[0].getScore();
		int numNonMatches = testNonMatchingVCF.length;
		int numMatches = testMatchingVCF.length;
		float oldFDR = ((float)numNonMatches)/((float)(numMatches + numNonMatches));
		String res = formatResults(Float.MIN_NORMAL, totalKey, oldFDR, numMatches, numNonMatches);
		results.append(res.toString());
		results.append("\n");

		//for each score in the nonMatching
		for (int i=0; i< testNonMatchingVCF.length; i++){
			float score = testNonMatchingVCF[i].getScore();
			if (score == oldScore) continue;
			numNonMatches = testNonMatchingVCF.length - i;
			numMatches = countNumberMatches(score);
			float fdr = ((float)numNonMatches)/((float)(numMatches + numNonMatches));
			if (fdr < oldFDR) oldFDR = fdr;
			res = formatResults(score, totalKey, oldFDR, numMatches, numNonMatches);
			results.append(res.toString());
			results.append("\n");
			if (numNonMatches == 0 || numMatches == 0) break;
			oldScore = score;
		}

		if (saveDirectory == null) System.out.println("\n"+results);
	}

	private int countNumberMatches(float score) {
		for (int i=0; i< testMatchingVCF.length; i++){
			if (testMatchingVCF[i].getScore() >= score){
				return testMatchingVCF.length - i;
			}
		}
		return 0;
	}
	
	private class ScoredCalls{
		String name;
		float[] tpr;
		float[] fdr;
		
		public ScoredCalls(String name){
			this.name = name;
			tpr = Num.arrayListOfFloatToArray(tprAL);
			fdr = Num.arrayListOfFloatToArray(fdrAL);
		}
	}

	public String formatResults(float threshold, float totalKey, float ratchetFDR, float intTest, float nonIntTest ){
		StringBuilder sb = new StringBuilder();
		//threshold
		if (threshold == Float.MIN_NORMAL) sb.append("none");
		else sb.append(threshold);
		sb.append("\t");
		sb.append((int)intTest); sb.append("\t");
		sb.append((int)nonIntTest); sb.append("\t");
		//fdr nonIntTest/totalTest
		sb.append(nonIntTest/(nonIntTest + intTest)); sb.append("\t");
		//ratchet fdr (always decreasing or the prior FDR)
		sb.append(ratchetFDR); sb.append("\t");
		fdrAL.add(ratchetFDR);
		//tpr intTest/totalKey
		float tpr = intTest/totalKey;
		tprAL.add(tpr);
		sb.append(intTest/totalKey); sb.append("\t");
		//fpr nonIntTest/totalKey
		sb.append(nonIntTest/totalKey); sb.append("\t");
		//ppv intTest/totalTest
		sb.append(intTest/(nonIntTest + intTest));
		return sb.toString();
	}

	public void intersectVCF(){
		//set all records to fail
		keyParser.setFilterFieldOnAllRecords(VCFRecord.FAIL);
		testParser.setFilterFieldOnAllRecords(VCFRecord.FAIL);

		ArrayList<VCFMatch> matches = new ArrayList<VCFMatch>();
		ArrayList<VCFRecord> testNonMatches = new ArrayList<VCFRecord>();

		//for each test record
		for (String chr: testParser.getChromosomeVCFRecords().keySet()){			
			VCFLookUp key = keyParser.getChromosomeVCFRecords().get(chr);
			VCFLookUp test = testParser.getChromosomeVCFRecords().get(chr);
			if (key == null) {
				//add all test to nonMatch
				for (VCFRecord r: test.getVcfRecord()) testNonMatches.add(r);
				continue;
			}
			countMatches(key, test, matches, testNonMatches);
		} 
		if (matches.size() == 0) Misc.printExit("\tNo matching vcf records?! Aborting.\n");

		//set arrays
		testMatchingVCF = new VCFMatch[matches.size()];
		testNonMatchingVCF = new VCFRecord[testNonMatches.size()];
		matches.toArray(testMatchingVCF);
		testNonMatches.toArray(testNonMatchingVCF);

		//Add failing key records to nonMatches
		ArrayList<VCFRecord> keyNonMatches = new ArrayList<VCFRecord>();
		for (VCFRecord v : keyParser.getVcfRecords()){
			if (v.getFilter().equals(VCFRecord.FAIL)) keyNonMatches.add(v);
		}
		keyNonMatchingVCF = new VCFRecord[keyNonMatches.size()];
		keyNonMatches.toArray(keyNonMatchingVCF);

	}

	public void countMatches(VCFLookUp key, VCFLookUp test, ArrayList<VCFMatch> matches, ArrayList<VCFRecord> nonMatches){
		//for each record in the test
		int[] posTest = test.getBasePosition();
		VCFRecord[] vcfTest = test.getVcfRecord();
		for (int i=0; i< vcfTest.length; i++){
			//fetch records from key
			VCFRecord[] matchingKey = key.fetchVCFRecords(posTest[i], posTest[i]+1);
			//no records found
			if (matchingKey == null) {
				nonMatches.add(vcfTest[i]);
			}
			else {
//System.out.println(vcfTest[i]);				
				//for each match
				boolean matchFound = false;
				for (int x=0; x< matchingKey.length; x++){
					//check to see if it matches
					if (vcfTest[i].matchesAlternateAlleleGenotype(matchingKey[x], requireGenotypeMatch)) {
//System.out.println("MatchFound");
						matchingKey[x].setFilter(VCFRecord.PASS);
						vcfTest[i].setFilter(VCFRecord.PASS);
						matches.add(new VCFMatch(matchingKey[x], vcfTest[i]));
						matchFound = true;
						break;
					}
				}

				if (matchFound == false) {
//System.out.println("NoMatch");					
					nonMatches.add(vcfTest[i]);
				}
			}
		}
	}

	public long overlapRegions(){
		commonRegions = new HashMap<String,RegionScoreText[]> ();
		long numberCommonBases = 0;
		for (String chr: testRegions.keySet()){
			//fetch arrays
			RegionScoreText[] key = keyRegions.get(chr);
			if (key == null) continue;
			RegionScoreText[] test = testRegions.get(chr);

			//find last base
			int lastBase = RegionScoreText.findLastBase(test);
			int lastBaseKey = RegionScoreText.findLastBase(key);
			if (lastBaseKey > lastBase) lastBase = lastBaseKey;

			//make booleans, initially all false, flip to true if covered;
			boolean[] keyBases = new boolean[lastBase+1];
			boolean[] testBases = new boolean[lastBase+1];
			for (int i=0; i< key.length; i++){
				for (int j=key[i].getStart(); j < key[i].getStop(); j++) keyBases[j] = true;
			}
			for (int i=0; i< test.length; i++){
				for (int j=test[i].getStart(); j < test[i].getStop(); j++) testBases[j] = true;
			}
			boolean[] good = new boolean[lastBase];
			Arrays.fill(good, true);
			for (int i=0; i<lastBase; i++){
				//both true, then set false, otherwise set true
				if (keyBases[i] == true && testBases[i] == true) good[i] = false;
			}
			//not intergenic, so add one to end.
			int[][] blocks = ExportIntergenicRegions.fetchFalseBlocks(good, 0, 0);
			RegionScoreText[] common = new RegionScoreText[blocks.length];
			for (int i=0; i< blocks.length; i++){
				common[i] = new RegionScoreText(blocks[i][0], blocks[i][1]+1, 0.0f, null);
				numberCommonBases+= common[i].getLength();
			}
			//add to hash
			commonRegions.put(chr, common);
		}
		return numberCommonBases;
	}

	public static HashMap<String, RegionScoreText[]> fixRegionChromosomeNames(HashMap<String, RegionScoreText[]> hash){
		HashMap<String, RegionScoreText[]> fixed = new HashMap<String, RegionScoreText[]>();
		for (String chr : hash.keySet()){
			RegionScoreText[] regions = hash.get(chr);
			if (chr.startsWith("chr")) fixed.put(chr, regions);
			else fixed.put("chr"+chr, regions);
		}
		return fixed;
	}

	public void parseFilterFiles(){

		//key regions
		if (keyRegions== null){
			keyRegions = Bed.parseBedFile(bedKey, true);
			keyRegions = fixRegionChromosomeNames(keyRegions);
			keyBps = RegionScoreText.countBases(keyRegions);
		}

		String res = keyBps +"\tInterrogated bps in key\n";
		results.append(res);

		//test regions
		if (testRegions == null){
			testRegions = Bed.parseBedFile(bedTest, true);
			testRegions = fixRegionChromosomeNames (testRegions);
			testBps = RegionScoreText.countBases(testRegions);
		}

		res = testBps +"\tInterrogated bps in test\n";
		results.append(res);

		//find common intersected regions common
		if (commonRegions == null) commonBps = overlapRegions();

		res = commonBps +"\tInterrogated bps in common\n";
		results.append(res);

		if (keyParser == null){			
			keyParser = new VCFParser(vcfKey, true, true, false);
			if (removeNonPass){				
				keyParser.setFilterFieldPeriodToTextOnAllRecords(VCFRecord.PASS);
				keyParser.filterVCFRecords(VCFRecord.PASS);
			}			
			keyParser.appendChr();			
			if (removeSNPs) keyParser.removeSNPs();
			if (removeNonSNPs) keyParser.removeNonSNPs();
			numberUnfilteredKeyVariants = keyParser.getVcfRecords().length;		
			keyParser.filterVCFRecords(commonRegions);
		}
		res = numberUnfilteredKeyVariants +"\tKey variants\n";
		results.append(res);
		
		res = keyParser.getVcfRecords().length +"\tKey variants in shared regions\n";
		results.append(res);
		
		res = keyParser.calculateTiTvRatio() +"\tShared key variants Ti/Tv\n";
		results.append(res);
		
		if (keyParser.getVcfRecords().length == 0) {
			System.out.println(results);
			Misc.printErrAndExit("\nNo key variants in shared regions? Aboring.\n");
		}

		testParser = new VCFParser(vcfTest, true, true, useVQSLOD);
		if (removeNonPass){
			testParser.setFilterFieldPeriodToTextOnAllRecords(VCFRecord.PASS);
			testParser.filterVCFRecords(VCFRecord.PASS);
		}
		testParser.appendChr();
		if (removeSNPs) testParser.removeSNPs();
		if (removeNonSNPs) testParser.removeNonSNPs();
		res = testParser.getVcfRecords().length +"\tTest variants\n";
		results.append(res);

		testParser.filterVCFRecords(commonRegions);
		res = testParser.getVcfRecords().length +"\tTest variants in shared regions\n";
		results.append(res);
		res = testParser.calculateTiTvRatio() +"\tShared test variants Ti/Tv\n";
		results.append(res);
		
		results.append("\n");
		
		if (testParser.getVcfRecords().length == 0) {
			System.out.println(results);
			Misc.printErrAndExit("\nNo test variants in shared regions? Aboring.\n");
		}
		
	}

	public void printParsedDatasets(){
		try {
			String filter = "All_";
			if (removeSNPs) filter = "NonSNP_";
			else if (removeNonSNPs) filter = "SNP_";

			//print common regions
			String bedKeyName = Misc.removeExtension(bedKey.getName());
			String bedTestName = Misc.removeExtension(bedTest.getName());
			File commonBed = new File (saveDirectory, "shared_"+bedKeyName+"_"+bedTestName+".bed.gz");
			Gzipper bed = new Gzipper(commonBed);
			for (String chr: commonRegions.keySet()){
				RegionScoreText[] r = commonRegions.get(chr);
				for (RegionScoreText x : r) bed.println(x.getBedLineJustCoordinates(chr));
			}
			bed.close();

			//sort arrays by chromosome and position
			VCFRecord[][] keyTestMatches = VCFMatch.split(testMatchingVCF);
			Arrays.sort(keyTestMatches[0]);
			Arrays.sort(keyTestMatches[1]);
			Arrays.sort(testNonMatchingVCF);
			Arrays.sort(keyNonMatchingVCF);

			//print matching and non matching vcf files
			String keyName = Misc.removeExtension(vcfKey.getName()); 
			String testName = Misc.removeExtension(vcfTest.getName());
			File matchingKey = new File (saveDirectory, "match_"+filter+keyName+"_"+testName+".vcf.gz");
			File noMatchingKey = new File (saveDirectory, "noMatch_"+filter+keyName+"_"+testName+".vcf.gz");
			keyParser.printRecords(keyTestMatches[0],matchingKey);
			keyParser.printRecords(keyNonMatchingVCF,noMatchingKey);

			File matchingTest = new File (saveDirectory, "match_"+filter+testName+"_"+keyName+".vcf.gz");
			File noMatchingTest = new File (saveDirectory, "noMatch_"+filter+testName+"_"+keyName+".vcf.gz");
			testParser.printRecords(keyTestMatches[1],matchingTest);
			testParser.printRecords(testNonMatchingVCF,noMatchingTest);

			//print results 
			File intersection = new File (saveDirectory, "comparison_"+filter+keyName+"_"+testName+".xls");
			IO.writeString(results.toString(), intersection);

		} catch (Exception e){

		}
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new VCFComparator(args);
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
					case 'a': vcfKey = new File(args[++i]); break;
					case 'b': bedKey = new File(args[++i]); break;
					case 'c': forExtraction = new File(args[++i]); break;
					case 'd': bedTest = new File(args[++i]); break;
					case 'p': saveDirectory = new File(args[++i]); break;
					case 'g': requireGenotypeMatch = true; break;
					case 's': removeNonSNPs = true; break;
					case 'v': useVQSLOD = true; break;
					case 'n': removeSNPs = true; break;
					case 'e': removeNonPass = true; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//checkfiles
		if (vcfKey == null || vcfKey.canRead() == false || 
				bedKey == null || bedKey.canRead() == false ||
				bedTest == null || bedTest.canRead() == false ) Misc.printErrAndExit("\nError: looks like you are missing or cannot read one of the four required files!\n");

		//pull files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please provide a test vcf file or directory containing such to compare against the key.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".vcf");
		tot[1] = IO.extractFiles(forExtraction,".vcf.gz");
		tot[2] = IO.extractFiles(forExtraction,".vcf.zip");
		testVcfFiles = IO.collapseFileArray(tot);
		if (testVcfFiles == null || testVcfFiles.length ==0 || testVcfFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx.vcf(.zip/.gz) file(s)!\n");


		if (saveDirectory != null){
			saveDirectory.mkdirs();
			if (saveDirectory.isDirectory() == false || saveDirectory.exists() == false) Misc.printErrAndExit("\nCannot find or make your save directory?! "+saveDirectory);
		}

		if (removeNonSNPs == true && removeSNPs == true) Misc.printErrAndExit("\nError: looks like you are throwing out all of your data by removing SNPs and non SNPs?! One or the other, not both.\n");
	}	

	private void printOptions() {
		StringBuilder res = new StringBuilder();
		res.append("VCF Comparator Settings:\n\n");
		res.append(vcfKey+"\tKey vcf file\n");
		res.append(bedKey+"\tKey interrogated regions file\n");
		res.append(vcfTest+"\tTest vcf file\n");
		res.append(bedTest+"\tTest interrogated regions file\n");
		res.append(saveDirectory+"\tSave directory for parsed datasets\n");
		res.append(requireGenotypeMatch+"\tRequire matching genotypes\n");
		res.append(useVQSLOD+"\tUse record VQSLOD score as ranking statistic\n");
		res.append(removeNonPass+ "\tExclude non PASS or . records\n");
		boolean all = removeSNPs == false && removeNonSNPs == false;
		if (all) res.append(all+"\tCompare all variant\n");
		else if (removeSNPs) res.append(removeSNPs+"\tCompare non-SNP variants, not SNPs\n");
		else if (removeNonSNPs) res.append(removeNonSNPs+"\tCompare SNPs, not non-SNP variants\n");
		res.append("\n");
		options = res.toString();
		System.out.print(options);
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             VCF Comparator : April 2013                          **\n" +
				"**************************************************************************************\n" +
				"Compares test vcf file(s) against a gold standard key of trusted vcf calls. Only calls\n" +
				"that fall in the common interrogated regions are compared. WARNING tabix gzipped files\n" +
				"often fail to parse correctly with java. Seeing odd error messages? Try uncompressing.\n\n" +

				"Required Options:\n"+
				"-a VCF file for the key dataset (xxx.vcf(.gz/.zip OK)).\n"+
				"-b Bed file of interrogated regions for the key dataset (xxx.bed(.gz/.zip OK)).\n"+
				"-c VCF file for the test dataset (xxx.vcf(.gz/.zip OK)). May also provide a directory\n" +
				"       containing xxx.vcf(.gz/.zip OK) files to compare.\n"+
				"-d Bed file of interrogated regions for the test dataset (xxx.bed(.gz/.zip OK)).\n"+

				"\nOptional Options:\n"+
				"-g Require the genotype to match, defaults to scoring a match when the alternate\n" +
				"       allele is present.\n"+
				"-v Use VQSLOD score as ranking statistic in place of the QUAL score.\n"+
				"-s Only compare SNPs, defaults to all.\n"+
				"-n Only compare non SNPs, defaults to all.\n"+
				"-p Provide a full path directory for saving the parsed data. Defaults to not saving.\n"+
				"-e Exclude test and key records whose FILTER field is not . or PASS. Defaults to\n" +
				"       scoring all.\n"+

				"\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/VCFComparator -a /NIST/NA12878/key.vcf\n" +
				"       -b /NIST/NA12878/regions.bed.gz -c /EdgeBio/Exome/testHaploCaller.vcf.zip\n" +
				"       -d /EdgeBio/Exome/NimbleGenExomeV3.bed -g -v -s -e -p /CompRes/ \n\n"+

		"**************************************************************************************\n");

	}
}
