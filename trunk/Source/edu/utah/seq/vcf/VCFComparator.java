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
	private File vcfTest;
	private File bedTest;
	private File saveDirectory = null;
	private boolean requireGenotypeMatch = false;
	private boolean removeSNPs = false;
	private boolean removeNonSNPs = false;
	
	private HashMap<String,RegionScoreText[]> keyRegions;
	private HashMap<String,RegionScoreText[]> testRegions;
	private HashMap<String,RegionScoreText[]> commonRegions;
	private VCFParser keyParser;
	private VCFParser testParser;
	private VCFMatch[] testMatchingVCF;
	private VCFRecord[] testNonMatchingVCF;
	private VCFRecord[] keyNonMatchingVCF;
	private StringBuilder results = new StringBuilder();
	private String options;
	
	private float[] minMaxScoreThresholds;

	//constructor
	public VCFComparator(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		//parse vcf file
		System.out.println("Parsing and filtering variant data for common interrogated regions...\n");
		parseFilterFiles();
		
		//set record QUAL score as thresholding score for roc curve data
		minMaxScoreThresholds = testParser.setRecordQUALAsScore();
		
		//compare calls in common interrogated regions
		System.out.println("\nComparing calls...");
		thresholdAndCompareCalls();
		
		if (saveDirectory != null) printParsedDatasets();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	public void thresholdAndCompareCalls(){
			//starting totals
			float totalKey = keyParser.getVcfRecords().length;
			
			//intersect and split test into matching and non matching
			intersectVCF();
			
			//sort by score smallest to largest
			Arrays.sort(testNonMatchingVCF, new ComparatorVCFRecordScore());
			Arrays.sort(testMatchingVCF);
			
			results.append("QUALThreshold\tNumMatchTest\tNumNonMatchTest\tTPR=matchTest/totalKey\tFPR=nonMatchTest/totalKey\tFDR=nonMatchTest/(matchTest+nonMatchTest)\tPPV=matchTest/(matchTest+nonMatchTest)\n");
			String r = formatResults(-1, totalKey, testMatchingVCF.length, testNonMatchingVCF.length);
			results.append(r.toString());
			results.append("\n\n");
			
			//for each score in the nonMatching
			float oldScore = Float.MIN_NORMAL;
			for (int i=0; i< testNonMatchingVCF.length; i++){
				float score = testNonMatchingVCF[i].getScore();
				if (score == oldScore) continue;
				int numNonMatches = testNonMatchingVCF.length - i;
				int numMatches = countNumberMatches(score);
				String res = formatResults(score, totalKey, numMatches, numNonMatches);
				results.append(res.toString());
				results.append("\n");
				if (numNonMatches == 0 || numMatches == 0) break;
				oldScore = score;
			}
			
			if (saveDirectory == null) System.out.println("\n"+results);
	}
	
	private int countNumberMatches(float score) {
		for (int i=0; i< testMatchingVCF.length; i++){
			if (testMatchingVCF[i].score >= score){
				return testMatchingVCF.length - i;
			}
		}
		return 0;
	}

	public String formatResults(float threshold, float totalKey, float intTest, float nonIntTest ){
		StringBuilder sb = new StringBuilder();
		sb.append(threshold); sb.append("\t");
		sb.append((int)intTest); sb.append("\t");
		sb.append((int)nonIntTest); sb.append("\t");
		//tpr intTest/totalKey
		sb.append(Num.formatNumber(intTest/totalKey, 3)); sb.append("\t");
		//fpr nonIntTest/totalKey
		sb.append(Num.formatNumber(nonIntTest/totalKey, 3)); sb.append("\t");
		//fdr nonIntTest/totalTest
		sb.append(Num.formatNumber(nonIntTest/(nonIntTest + intTest), 3)); sb.append("\t");
		//ppv intTest/totalTest
		sb.append(Num.formatNumber(intTest/(nonIntTest + intTest), 3));
		return sb.toString();
	}
	
	public void intersectVCF(){
		//set all records to fail
		keyParser.setFilterFieldOnAllRecords(VCFRecord.FAIL);
		testParser.setFilterFieldOnAllRecords(VCFRecord.FAIL);
		
		ArrayList<VCFMatch> matches = new ArrayList<VCFMatch>();
		ArrayList<VCFRecord> nonMatches = new ArrayList<VCFRecord>();
		//for each test record
		for (String chr: testParser.getChromosomeVCFRecords().keySet()){
			VCFLookUp key = keyParser.getChromosomeVCFRecords().get(chr);
			VCFLookUp test = testParser.getChromosomeVCFRecords().get(chr);
			if (key == null) {
				//add all test to nonMatch
				for (VCFRecord r: test.getVcfRecord()) nonMatches.add(r);
				continue;
			}
			countMatches(key, test, matches, nonMatches);
		} 
		if (matches.size() == 0) Misc.printExit("\tNo matching vcf records?! Aborting.\n");
		//set arrays
		testMatchingVCF = new VCFMatch[matches.size()];
		testNonMatchingVCF = new VCFRecord[nonMatches.size()];
		matches.toArray(testMatchingVCF);
		nonMatches.toArray(testNonMatchingVCF);
		
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
				//check to see only one of the key was found
				if (matchingKey.length !=1) Misc.printErrAndExit("\nMore than one vcf record in the key was found to match \ntest\t"+vcfTest[i]+"\nkey[0]\t"+matchingKey[0]);
				//check to see if it matches
				if (vcfTest[i].matchesAlternateAlleleGenotype(matchingKey[0], requireGenotypeMatch)) {
					matchingKey[0].setFilter(VCFRecord.PASS);
					vcfTest[i].setFilter(VCFRecord.PASS);
					matches.add(new VCFMatch(matchingKey[0], vcfTest[i]));
				}
				else nonMatches.add(vcfTest[i]);
			}
		}
	}
	
	private class VCFMatch implements Comparable<VCFMatch>{
		float score;
		VCFRecord key;
		VCFRecord test;
		
		private VCFMatch (VCFRecord key, VCFRecord test){
			this.key = key;
			this.test = test;
			score = test.getScore();
		}
		/**Sorts by score, smallest to largest*/
		public int compareTo(VCFMatch second) {
			if (this.score < second.score) return -1;
			if (this.score > second.score) return 1;
			return 0;
		}
	}
	
	public int overlapRegions(){
		commonRegions = new HashMap<String,RegionScoreText[]> ();
		int numberCommonBases = 0;
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
			for (int i=0; i<lastBase; i++){
				//both true, then set false, otherwise set true
				if (keyBases[i] == true && testBases[i] == true) testBases[i] = false;
				else testBases[i] = true;
			}
			//not intergenic, so add one to end.
			int[][] blocks = ExportIntergenicRegions.fetchFalseBlocks(testBases, 0, 0);
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
	
	public void parseFilterFiles(){
		keyRegions = Bed.parseBedFile(bedKey, true);
		int numberBasesInKey = RegionScoreText.countBases(keyRegions);
		String res = numberBasesInKey +"\tInterrogated bps in key\n";
		System.out.print(res);
		results.append(res);
		testRegions = Bed.parseBedFile(bedTest, true);
		int numberBasesInTest = RegionScoreText.countBases(testRegions);
		res = numberBasesInTest +"\tInterrogated bps in test\n";
		System.out.print(res);
		results.append(res);
		
		//find common intersected regions common 
		int num = overlapRegions();
		res = num +"\tInterrogated bps in common\n";
		System.out.print(res);
		results.append(res);
		
		keyParser = new VCFParser(vcfKey, true, true);
		if (removeSNPs) keyParser.removeSNPs();
		if (removeNonSNPs) keyParser.removeNonSNPs();
		res = keyParser.getVcfRecords().length+"\tKey variants\n";
		System.out.print(res);
		results.append(res);
		
		keyParser.filterVCFRecords(commonRegions);
		res = keyParser.getVcfRecords().length +"\tKey variants in shared regions\n";
		System.out.print(res);
		results.append(res);
		
		testParser = new VCFParser(vcfTest, true, true);
		if (removeSNPs) testParser.removeSNPs();
		if (removeNonSNPs) testParser.removeNonSNPs();
		res = testParser.getVcfRecords().length +"\tTest variants\n";
		System.out.print(res);
		results.append(res);
		
		testParser.filterVCFRecords(commonRegions);
		res = testParser.getVcfRecords().length +"\tTest variants in shared regions\n";
		System.out.print(res);
		results.append(res);
		results.append("\n");
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

			//print matching and non matching vcf files
			String keyName = Misc.removeExtension(vcfKey.getName()); 
			File matchingKey = new File (saveDirectory, "match_"+filter+keyName+".vcf.gz");
			File noMatchingKey = new File (saveDirectory, "noMatch_"+filter+keyName+".vcf.gz");
			//keyParser.printRecords(VCFRecord.PASS, matchingKey, noMatchingKey);

			String testName = Misc.removeExtension(vcfTest.getName()); 
			File matchingTest = new File (saveDirectory, "match_"+filter+testName+".vcf.gz");
			File noMatchingTest = new File (saveDirectory, "noMatch_"+filter+testName+".vcf.gz");
			//testParser.printRecords(VCFRecord.PASS, matchingTest, noMatchingTest);

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
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'a': vcfKey = new File(args[++i]); break;
					case 'b': bedKey = new File(args[++i]); break;
					case 'c': vcfTest = new File(args[++i]); break;
					case 'd': bedTest = new File(args[++i]); break;
					case 'p': saveDirectory = new File(args[++i]); break;
					case 'g': requireGenotypeMatch = true; break;
					case 's': removeNonSNPs = true; break;
					case 'n': removeSNPs = true; break;
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
				vcfTest == null || vcfTest.canRead() == false ||
				bedTest == null || bedTest.canRead() == false ) Misc.printErrAndExit("\nError: looks like you are missing or cannot read one of the four required files!\n");
		if (saveDirectory != null){
			saveDirectory.mkdirs();
			if (saveDirectory.isDirectory() == false || saveDirectory.exists() == false) Misc.printErrAndExit("\nCannot find or make your save directory?! "+saveDirectory);
		}
		
		if (removeNonSNPs == true && removeSNPs == true) Misc.printErrAndExit("\nError: looks like you are throwing out all of your data by removing SNPs and non SNPs?! One or the other, not both.\n");
		
		printOptions();
		
		//add options to results?
		if (saveDirectory != null) results.append(options);
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
		boolean all = removeSNPs == false && removeNonSNPs == false;
		res.append(all+"\tCompare all variant\n");
		if (all == false){
			res.append(removeSNPs+"\tCompare non-SNP variants, not SNPs\n");
			res.append(removeNonSNPs+"\tCompare SNPs, not non-SNP variants\n");
		}
		res.append("\n");
		options = res.toString();
		System.out.print(options);
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             VCF Comparator : March 2013                          **\n" +
				"**************************************************************************************\n" +
				"Compares a test vcf file against a gold standard key of trusted vcf calls. Only calls\n" +
				"that fall in the common interrogated regions are compared. \n\n" +

				"Required Options:\n"+
				"-a VCF file for the key dataset (xxx.vcf(.gz/.zip OK)).\n"+
				"-b Bed file of interrogated regions for the key dataset (xxx.bed(.gz/.zip OK)).\n"+
				"-c VCF file for the test dataset (xxx.vcf(.gz/.zip OK)).\n"+
				"-d Bed file of interrogated regions for the test dataset (xxx.bed(.gz/.zip OK)).\n"+
				
				"\nOptional Options:\n"+
				"-g Require the genotype to match, defaults to scoring a match when then alternate\n" +
				"       allele is present.\n"+
				"-s Only compare SNPs, defaults to all.\n"+
				"-n Only compare non SNPs, defaults to all.\n"+
				"-p Provide a full path directory for saving the parsed data.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/VCFComparator\n\n"+

		"**************************************************************************************\n");

	}
}
