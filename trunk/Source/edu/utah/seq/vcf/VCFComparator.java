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

/**Compares variant lists
 * @author Nix
 * */
public class VCFComparator {

	//user fields
	private File vcfKey;
	private File bedKey;
	private File vcfTest;
	private File bedTest;
	private boolean requireGenotypeMatch = true;
	
	private HashMap<String,RegionScoreText[]> keyRegions;
	private HashMap<String,RegionScoreText[]> testRegions;
	private HashMap<String,RegionScoreText[]> commonRegions;
	private VCFParser keyParser;
	private VCFParser testParser;
	private int sampleIndexForScore = 0;
	private float[] minMaxScoreThresholds;
	
	private ArrayList<String> matches = new ArrayList<String>();
	private ArrayList<String> misMatches = new ArrayList<String>();
	private int totalNumberArraySnps = 0;
	private int numberArraySnpsFailingMinimumScore = 0;
	private int numberArraySnpsWithNoVcf = 0;
	private int numberVCFsFailingMinimumScore = 0;
	private int numberSNPMatches = 0;
	private int numberSNPNoMatches;
	private static final Pattern UNDERSCORE = Pattern.compile("_");
	private String url;
	
	
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
		minMaxScoreThresholds = testParser.setGenotypeQualityGQScore(sampleIndexForScore);
		
		//compare calls in common interrogated regions
		System.out.print("\nComparing calls...");
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
			float totalTest = testParser.getVcfRecords().length;
			//intersect and remove non intersecting
			compareCalls();
			float subKey = keyParser.getVcfRecords().length;
			float subTest = testParser.getVcfRecords().length;
			
			String r = formatResults(minMaxScoreThresholds[0], totalKey, totalTest, subKey, subTest);
			System.out.println(r);
			
		
	}
	
	public String formatResults(float threshold, float totalKey, float totalTest, float intKey, float intTest){
		StringBuilder sb = new StringBuilder();
		//threshold
		sb.append(threshold); sb.append("\t");
		sb.append((int)totalKey); sb.append("\t");
		sb.append((int)totalTest); sb.append("\t");
		sb.append((int)intKey); sb.append("\t");
		sb.append((int)intTest);
		return sb.toString();
	}
	
	public void compareCalls(){
		//set fail
		keyParser.setFilterFieldOnAllRecords(VCFRecord.FAIL);
		testParser.setFilterFieldOnAllRecords(VCFRecord.FAIL);
		
		//for each test record
		for (String chr: testParser.getChromosomeVCFRecords().keySet()){
			VCFLookUp key = keyParser.getChromosomeVCFRecords().get(chr);
			if (key == null) continue;
			VCFLookUp test = testParser.getChromosomeVCFRecords().get(chr);
			countMatches(key, test);
		}
		
		//reduce to those that pass
		keyParser.filterVCFRecords(VCFRecord.PASS);
		testParser.filterVCFRecords(VCFRecord.PASS);
		
	}
	
	public void countMatches(VCFLookUp key, VCFLookUp test){
		
		//for each record in the key
		int numKey = key.getBasePosition().length;
		int numTest = test.getBasePosition().length;
		VCFLookUp walk = key;
		VCFLookUp lookup = test;
		if (numKey > numTest){
			walk = test;
			lookup = key;
		}
		int[] pos = walk.getBasePosition();
		VCFRecord[] rec = walk.getVcfRecord();
		for (int i=0; i< pos.length; i++){
			//fetch records from lookup, should only be one
			VCFRecord[] matchingLookup = lookup.fetchVCFRecords(pos[i], pos[i]+1);
			if (matchingLookup == null) continue;
			if (matchingLookup.length !=1) Misc.printErrAndExit("\nMore than one vcf record found to match ->\n\t"+rec[i]+"\n\t"+matchingLookup[0]);
			if (rec[i].matchesAlternateAlleleGenotype(matchingLookup[0], requireGenotypeMatch)){
				rec[i].setFilter(VCFRecord.PASS);
				matchingLookup[0].setFilter(VCFRecord.PASS);
			}
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
		
		testParser = new VCFParser(vcfTest, false, true);
		System.out.println(keyParser.getVcfRecords().length +"\tTest variants");
		testParser.filterVCFRecords(commonRegions);
		System.out.println(keyParser.getVcfRecords().length +"\tTest variants in common regions");
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
				"Compares a test vcf file against a gold standard key vcf file. Be sure each of your\n" +
				"region files don't contain overlaping regions, run the MergeRegions application if\n" +
				"needed. BETAAAAAAAAAA \n\n" +

				"Options:\n"+
				"-a VCF file for the key dataset (xxx.vcf(.gz/.zip OK)).\n"+
				"-b Bed file of interrogated regions for the key dataset (xxx.bed(.gz/.zip OK)).\n"+
				"-c VCF file for the test dataset (xxx.vcf(.gz/.zip OK)).\n"+
				"-d Bed file of interrogated regions for the test dataset (xxx.bed(.gz/.zip OK)).\n"+
				//"-n Don't require the genotype to match, just the presence of the alternate allele.\n"+

				"\n"+

				"Example: java -Xmx1500M -jar pathTo/USeq/Apps/\n\n"+

		"**************************************************************************************\n");

	}
}
