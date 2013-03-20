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
	private boolean requireGenotypeMatch = false;
	
	private HashMap<String,RegionScoreText[]> keyRegions;
	private HashMap<String,RegionScoreText[]> testRegions;
	private HashMap<String,RegionScoreText[]> commonRegions;
	private VCFParser keyParser;
	private VCFParser testParser;
	private VCFMatch[] testMatchingVCF;
	private VCFRecord[] testNonMatchingVCF;
	
	private float[] minMaxScoreThresholds;

	
	
	//constructor
	public VCFComparator(String[] args){
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);
		
		//parse vcf file
		System.out.println("Parsing and filtering variant data for common interrogated regions...");
		parseFilterFiles();
		
		//set genotypeQualityGQ from sampleIndexForScore as VCFRecord score for thresholding
		minMaxScoreThresholds = testParser.setRecordQUALAsScore();
		System.out.println(minMaxScoreThresholds[0]+"\tMinimum test score");
		System.out.println(minMaxScoreThresholds[1]+"\tMaximum test score");
		
		
		//compare calls in common interrogated regions
		System.out.println("\nComparing calls...");
		thresholdAndCompareCalls();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	public void printStats(){
		System.out.println();
		
	}
	
	public void thresholdAndCompareCalls(){
			//starting totals
			float totalKey = keyParser.getVcfRecords().length;
			
			//intersect and split test into matching and non matching
			intersectVCF();
			
			//sort by score smallest to largest
			Arrays.sort(testNonMatchingVCF, new ComparatorVCFRecordScore());
			Arrays.sort(testMatchingVCF);
			
			System.out.println("QUALThreshold\tNumMatchTest\tNumNonMatchTest\tTPR=matchTest/totalKey\tFPR=nonMatchTest/totalKey\tFDR=nonMatchTest/(matchTest+nonMatchTest)\tPPV=matchTest/(matchTest+nonMatchTest)");
			String r = formatResults(-1, totalKey, testMatchingVCF.length, testNonMatchingVCF.length);
			System.out.println(r);
			System.out.println();
			
			//for each score in the nonMatching
			float oldScore = Float.MIN_NORMAL;
			for (int i=0; i< testNonMatchingVCF.length; i++){
				float score = testNonMatchingVCF[i].getScore();
				if (score == oldScore) continue;
				int numNonMatches = testNonMatchingVCF.length - i;
				int numMatches = countNumberMatches(score);
				String res = formatResults(score, totalKey, numMatches, numNonMatches);
				System.out.println(res);
				if (numNonMatches == 0 || numMatches == 0) break;
				oldScore = score;
			}
			
		
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
		sb.append(Num.formatNumber(nonIntTest/(nonIntTest + intTest), 3));
		//ppv intTest/totalTest
		sb.append(Num.formatNumber(intTest/(nonIntTest + intTest), 3));
		return sb.toString();
	}
	
	public void intersectVCF(){
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
				if (vcfTest[i].matchesAlternateAlleleGenotype(matchingKey[0], requireGenotypeMatch)) matches.add(new VCFMatch(matchingKey[0], vcfTest[i]));
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
				for (int j=key[i].getStart(); j <= key[i].getStop(); j++) keyBases[j] = true;
			}
			for (int i=0; i< test.length; i++){
				for (int j=test[i].getStart(); j <= test[i].getStop(); j++) testBases[j] = true;
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
		System.out.println(numberBasesInKey +"\tInterrogated bps in key");
		testRegions = Bed.parseBedFile(bedTest, true);
		int numberBasesInTest = RegionScoreText.countBases(testRegions);
		System.out.println(numberBasesInTest +"\tInterrogated bps in test");
		
		//find common intersected regions common 
		int num = overlapRegions();
		System.out.println(num +"\tInterrogated bps in common ");
		
		keyParser = new VCFParser(vcfKey, true, true);
		System.out.println(keyParser.getVcfRecords().length+"\tKey variants");
		//remove all non intersecting records
		keyParser.filterVCFRecords(commonRegions);
		System.out.println(keyParser.getVcfRecords().length +"\tKey variants in common regions");
		
		testParser = new VCFParser(vcfTest, true, true);
		System.out.println(testParser.getVcfRecords().length +"\tTest variants");
		testParser.filterVCFRecords(commonRegions);
		System.out.println(testParser.getVcfRecords().length +"\tTest variants in common regions");
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
					//case 'n': requireGenotypeMatch = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//checkfiles
		if (vcfKey == null || vcfKey.canRead() == false || 
				bedKey == null || bedKey.canRead() == false ||
				vcfTest == null || vcfTest.canRead() == false ||
				bedTest == null || bedTest.canRead() == false ) Misc.printExit("\nError: looks like you are missing or cannot read one of the four required files!\n");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                             VCF Comparator : March 2013                          **\n" +
				"**************************************************************************************\n" +
				"Compares a test vcf file against a gold standard key vcf file. \n\n" +

				"Options:\n"+
				"-a VCF file for the key dataset (xxx.vcf(.gz/.zip OK)).\n"+
				"-b Bed file of interrogated regions for the key dataset (xxx.bed(.gz/.zip OK)).\n"+
				"-c VCF file for the test dataset (xxx.vcf(.gz/.zip OK)).\n"+
				"-d Bed file of interrogated regions for the test dataset (xxx.bed(.gz/.zip OK)).\n"+
				"-n Don't require the genotype to match, just the presence of the alternate allele.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/\n\n"+

		"**************************************************************************************\n");

	}
}
