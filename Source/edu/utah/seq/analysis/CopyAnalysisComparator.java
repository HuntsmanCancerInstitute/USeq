package edu.utah.seq.analysis;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import edu.utah.seq.its.Interval;
import edu.utah.seq.its.IntervalTree;
import edu.utah.seq.useq.data.RegionScoreText;
import edu.utah.seq.vcf.ComparatorVCFRecordScore;
import edu.utah.seq.vcf.VCFLookUp;
import edu.utah.seq.vcf.VCFMatch;
import edu.utah.seq.vcf.VCFParser;
import edu.utah.seq.vcf.VCFRecord;
import edu.utah.seq.vcf.VCFSample;
import util.apps.MergeRegions;
import util.bio.annotation.Bed;
import util.bio.annotation.Coordinate;
import util.bio.annotation.ExportIntergenicRegions;
import util.bio.parsers.MergeAdjacentRegions;
import util.gen.*;

/**Compares call lists of regions against a key calculating pseudo ROC data for graphing.
 * @author Nix
 * */
public class CopyAnalysisComparator {

	//user fields
	private File keyRegionsFile;
	private File[] testRegionFiles;
	private File testRegionFile;
	private File saveDirectory = null;
	private int minimumBpOverlapKeyTest = 1;
	private int mergeBpGap = 1000;
	private int minimumMergedRegionSize = 2500;
	private boolean takeAbsOfTestRegionScores = false;

	private HashMap<String,IntervalTree<RegionScoreText>> keyIntervalTrees = null;
	private HashMap<String,RegionScoreText[]> keyRegions = null;
	private HashMap<String,RegionScoreText[]> testRegions = null;
	private long keyBps;
	private int keyRegionCount;

	private long testBps;
	private int testRegionCount;

	private ArrayList<Float> tprAL = null;
	private ArrayList<Float> fdrAL = null;
	private ArrayList<ScoredCalls> scoredCallsAL = new ArrayList<ScoredCalls>();
	private String[] fixedFdrLines = null;
	double[] fixedFdr = new double[] {0.2, 0.15, 0.10};
	private String headerLine = "Threshold\tKeyMatches\tKeyNonMatches\tTestMatches\tTestNonMatches\t"
			+ "FDR=TestNonMatches/TotalTest\t"
			+ "DecreasingFDR\t"
			+ "Recall TPR=KeyMatches/TotalKey\t"
			+ "FPR=KeyNonMatchs/TotalKey\t"
			+ "Precision PPV=TestMatches/TotalTest\t"
			+ "F-score=HarmonicMean(Precision,Recall)";

	//constructor
	public CopyAnalysisComparator(String[] args) throws Exception{
		//start clock
		long startTime = System.currentTimeMillis();

		//process args
		processArgs(args);

		//parse the key
		IO.pl("Parsing the key file...");
		keyRegions = Bed.parseBedFile(keyRegionsFile, true, true);
		keyBps = RegionScoreText.countBases(keyRegions);
		keyRegionCount = countRegions(keyRegions);
		keyIntervalTrees = createIntervalTreesForBedCalls(keyRegions);
		IO.pl("\t"+ keyRegionCount+" # Keys, "+keyBps+" # Key BPs");


		//for each test file
		for (int i=0; i< testRegionFiles.length; i++){
			testRegionFile = testRegionFiles[i];

			//parse test region file
			IO.pl("\nParsing and test region data for "+testRegionFile.getName());
			testRegions = Bed.parseBedFile(testRegionFile, true, true);
			if (takeAbsOfTestRegionScores) setAbsValOfRegionScores(testRegions);
			testBps = RegionScoreText.countBases(testRegions);
			testRegionCount = countRegions(testRegions);
			IO.pl("\tInput\t"+ testRegionCount+" # Test Regions, "+testBps+" # Test BPs");
			
			//merge the regions and assign max score
			testRegions = makeMergeAssignHighestScore(testRegions);
			testBps = RegionScoreText.countBases(testRegions);
			testRegionCount = countRegions(testRegions);
			IO.pl("\tMerged\t"+ testRegionCount+" # Test Regions, "+testBps+" # Test BPs");
			saveRegions(testRegions, new File(saveDirectory, Misc.removeExtension(testRegionFile.getName())+"_CACMerged.bed.gz"));

			//compare calls in common interrogated regions
			IO.pl("Intersecting key with merged regions...");
			tprAL = new ArrayList<Float>();
			fdrAL = new ArrayList<Float>();
			thresholdAndCompareCalls();

			//make scoredCalls
			scoredCallsAL.add(new ScoredCalls(Misc.removeExtension(testRegionFile.getName())));

		}

		//print out composite scoredCalls?
		if (scoredCallsAL.size() >1) printScoredCalls();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		IO.pl("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	private void saveRegions(HashMap<String, RegionScoreText[]> regions, File file) throws Exception {
		Gzipper out = new Gzipper(file);
		for (String chr: regions.keySet()){
			for (RegionScoreText r: regions.get(chr)) {
				out.println(r.getBedLine(chr));
			}
		}
		out.close();
	}

	private HashMap<String, RegionScoreText[]> makeMergeAssignHighestScore( HashMap<String, RegionScoreText[]> toMerge) throws IOException {
		IO.pl("Merging test regions...");
		String name = Misc.removeExtension(testRegionFile.getName());

		//MergeRegions
		//public MergeRegions (File[] regionFiles, File mergedFile, boolean verbose, boolean printHeader){
		File mr = new File (saveDirectory, name+".mr.tmp.bed");
		if (mr.exists() == false) {
			IO.pl("\tMergingRegions...");
			new MergeRegions (new File[] {testRegionFile}, mr, false, false);
		}

		//MergeAdjacentRegions
		File mar = new File (saveDirectory, name+".mar.tmp.bed.gz");
		if (mar.exists() == false) {
			IO.pl("\tMergingAdjacentRegions...");
			String[] args = {"-b", mr.toString(), "-r", mar.toString(), "-m", new Integer(mergeBpGap).toString(), "-q"};
			new MergeAdjacentRegions (args);
		}

		//load the merge and size select
		IO.pl("\tSize selecting merged regions...");
		HashMap<String, RegionScoreText[]> mergedRegions = Bed.parseBedFile(mar, true, false);
		sizeSelect(mergedRegions);
		
		IO.pl("\tAssigning max score to each merged test region...");
		HashMap<String, IntervalTree<RegionScoreText>> trees = createIntervalTreesForBedCalls(toMerge);

		//for each merged region chr
		for (String chr: mergedRegions.keySet()){
			//fetch big regions
			RegionScoreText[] bigRs = mergedRegions.get(chr);

			//fetch the interval tree from the key
			IntervalTree<RegionScoreText> key = trees.get(chr);

			//for each mergedRegion look for intersection with the tree, interbase coordinates
			for (RegionScoreText r: bigRs) {
				ArrayList<RegionScoreText> matchingKey = key.search(r.getStart(), r.getStop());
				if (matchingKey.size() ==0) throw new IOException("\nFAILED to find matching small regions that intersect the merged region "+r.toString());
				//for each key match add it
				float maxScore = 0;
				for (RegionScoreText k: matchingKey) {
					if (k.getScore()> maxScore) maxScore = k.getScore();

				}
				r.setScore(maxScore);
				//IO.pl(r.getBedLine(chr));
			}
		}
		return mergedRegions;

	}

	private void sizeSelect(HashMap<String, RegionScoreText[]> mergedRegions) {
		ArrayList<RegionScoreText> toKeep = new ArrayList<RegionScoreText>();
		int unfilteredCount = 0;
		int postFilteredCount = 0;
		for (String chr: mergedRegions.keySet()) {
			RegionScoreText[] regions = mergedRegions.get(chr);
			unfilteredCount+= regions.length;
			toKeep.clear();
			for (RegionScoreText r: regions) if (r.getLength()>= minimumMergedRegionSize) toKeep.add(r);
			if (toKeep.size() == 0) mergedRegions.remove(chr);
			else {
				regions = new RegionScoreText[toKeep.size()];
				toKeep.toArray(regions);
				mergedRegions.put(chr, regions);
				postFilteredCount+= regions.length;
			}
		}
		IO.pl("\t\t"+unfilteredCount+" -> "+postFilteredCount);
	}

	private void setAbsValOfRegionScores(HashMap<String, RegionScoreText[]> regions) {
		for (RegionScoreText[] r : regions.values()){
			for (int i=0; i< r.length; i++) r[i].setScore(Math.abs(r[i].getScore())) ;
		}
		
	}

	private int countRegions(HashMap<String, RegionScoreText[]> r) {
		int count = 0;
		for (String chr: r.keySet()) {
			count += r.get(chr).length;
		}
		return count;
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
		File f = new File (saveDirectory, "fdrTprSummary.xls");
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
			e.printStackTrace();
			System.exit(1);
		}
	}

	public void thresholdAndCompareCalls(){	
		IO.pl(headerLine);
		//find all thresholds, sort smallest to largest
		TreeSet<Float> thresholds = fetchTestScores();

		//calculate counts for each threshold
		float priorFdr = 1;
		for (Float threshold: thresholds) {
			//#KeyMatches, #KeyNonMatches, #TestMatches, #TestNonMatches
			int[] counts = intersect(threshold);
			priorFdr = printResultsCnv(threshold, counts, priorFdr);
		}
	}
	
	public float printResultsCnv(float threshold, int[] counts, float priorFdr){
		StringBuilder sb = new StringBuilder();
		//threshold
		sb.append(threshold);
		sb.append("\t");
		
		//#KeyMatches, #KeyNonMatches, #TestMatches, #TestNonMatches
		//     0             1              2              3
		for (int i=0; i<counts.length; i++) {
			sb.append(counts[i]);
			sb.append("\t");
		}

		//fdr nonIntTest/totalTest
		float fdr = (float)counts[3]/(float)(counts[2]+counts[3]);
		sb.append(fdr); sb.append("\t");
		
		//ratchet fdr (always decreasing or the prior FDR)
		if (priorFdr < fdr) fdr = priorFdr;
		sb.append(fdr); sb.append("\t");
		fdrAL.add(fdr);
		
		//tpr intKey/totalKey
		float tpr = (float)counts[0]/(float)keyRegionCount;
		sb.append(tpr); sb.append("\t");
		tprAL.add(tpr);
		
		//FPR=NonMatchKey/TotalKey
		float fpr = (float)counts[1]/(float)keyRegionCount;
		sb.append(fpr); sb.append("\t");
		
		//ppv intTest/totalTest
		float ppv = (float)counts[2]/(float)(counts[2]+counts[3]);
		sb.append(ppv); sb.append("\t");
		
		//Recall, TRUTH.TP / (TRUTH.TP + TRUTH.FN), same as tpr
		//Precision, QUERY.TP / (QUERY.TP + QUERY.FP)
		//F-score harmonic mean of tpr (recall) and ppv (precision)
		double hm = Num.harmonicMean(new double[] {tpr, ppv});
		sb.append((float)hm);
		IO.pl(sb);
		
		return fdr;
	}

	private TreeSet<Float> fetchTestScores() {
		TreeSet<Float> uni = new TreeSet<Float>();
		uni.add(new Float(0.0f));
		for (RegionScoreText[] r : testRegions.values()){
			for (int i=0; i< r.length; i++) uni.add(r[i].getScore());
		}
		return uni;
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

	public int[] intersect(float scoreThreshold){

		//could just use counters
		ArrayList<Bed> testMatches = new ArrayList<Bed>();
		ArrayList<Bed> testNonMatches = new ArrayList<Bed>();
		HashMap<String, Bed> keyMatches = new HashMap<String, Bed>();

		//for each chr of test records
		for (String chr: testRegions.keySet()){
			RegionScoreText[] t = testRegions.get(chr);
			ArrayList<RegionScoreText> filteredTests = filterRegions(t, scoreThreshold);
			if (filteredTests.size()==0) continue;

			//fetch the interval tree from the key
			IntervalTree<RegionScoreText> key = keyIntervalTrees.get(chr);
			if (key == null) {
				//no match for any of the tests so add all to nonMatch
				for (RegionScoreText r: filteredTests) testNonMatches.add(new Bed(chr, '.', r));
			}
			else {
				//for each test Region look for intersection with the tree, interbase coordinates
				for (RegionScoreText r: filteredTests) {
					ArrayList<RegionScoreText> matchingKey = key.search(r.getStart(), r.getStop());
					//no potential matches?
					if (matchingKey.size() == 0) testNonMatches.add(new Bed(chr, '.', r));
					//potential matches found, test each one for required overlap
					else {
						boolean matchFound = false;
						//for each key match add it
						for (RegionScoreText k: matchingKey) {
							int overlap = r.bpsIntersection(k);
							if (overlap >= minimumBpOverlapKeyTest || overlap == 0) {
								matchFound = true;
								Bed kb = new Bed(chr, '.', k);
								String kbString = kb.toStringNoStrand();
								keyMatches.put(kbString, kb);
							}
						}
						if (matchFound )testMatches.add(new Bed(chr, '.', r));
						else testNonMatches.add(new Bed(chr, '.', r));
					}
				}
			}
		} 
		//return results
		//#KeyMatches, #KeyNonMatches, #TestMatches, #TestNonMatches 
		//IO.pl("\nTestMatches\n"+testMatches);
		//IO.pl("\nTestNonMatches\n"+testNonMatches);
		
		return new int[] {keyMatches.size(), keyRegionCount-keyMatches.size(), testMatches.size(), testNonMatches.size()}; 

	}

	private ArrayList<RegionScoreText> filterRegions(RegionScoreText[] tests, float scoreThreshold) {
		ArrayList<RegionScoreText> toReturn = new ArrayList<RegionScoreText>();
		for (RegionScoreText r: tests) if (r.getScore()>= scoreThreshold) toReturn.add(r);
		return toReturn;
	}

	private HashMap<String,IntervalTree<RegionScoreText>> createIntervalTreesForBedCalls(HashMap<String,RegionScoreText[]> regionsToTree) {
		//make HashMap of trees
		HashMap<String,IntervalTree<RegionScoreText>> trees = new HashMap<String,IntervalTree<RegionScoreText>>();
		for (String chr : regionsToTree.keySet()){
			RegionScoreText[] regions = regionsToTree.get(chr);
			ArrayList<Interval<RegionScoreText>> ints = new ArrayList();
			for (int i =0; i< regions.length; i++) {				
				ints.add(new Interval<RegionScoreText>(regions[i].getStart(), regions[i].getStop(), regions[i]));
			}
			IntervalTree<RegionScoreText> tree = new IntervalTree(ints, false);
			trees.put(chr, tree);
		}
		return trees;
	}

	public static void main(String[] args) throws Exception {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new CopyAnalysisComparator(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		IO.pl("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File forExtraction = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'k': keyRegionsFile = new File(args[++i]); break;
					case 't': forExtraction = new File(args[++i]); break;
					case 'g': mergeBpGap = Integer.parseInt(args[++i]); break;
					case 'o': minimumBpOverlapKeyTest = Integer.parseInt(args[++i]); break;
					case 'm': minimumMergedRegionSize = Integer.parseInt(args[++i]); break;
					case 'a': takeAbsOfTestRegionScores = true; break;
					case 's': saveDirectory = new File(args[++i]); break;
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
		if (keyRegionsFile == null) Misc.printErrAndExit("\nError: please provide a bed file of key regions.\n");

		//pull test region files
		if (forExtraction == null || forExtraction.canRead() == false) Misc.printExit("\nError: please provide a test region file xxx.bed(.zip/.gz) or directory containing such to compare against the key.\n");
		File[][] tot = new File[3][];
		tot[0] = IO.extractFiles(forExtraction,".bed");
		tot[1] = IO.extractFiles(forExtraction,".bed.gz");
		tot[2] = IO.extractFiles(forExtraction,".bed.zip");
		testRegionFiles = IO.collapseFileArray(tot);
		if (testRegionFiles == null || testRegionFiles.length ==0 || testRegionFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your test xxx.bed(.zip/.gz) file(s)!\n");
		if (saveDirectory == null) Misc.printErrAndExit("\nCannot find your save directory?! "+saveDirectory);
		saveDirectory.mkdirs();
		printOptions();
	}	

	private void printOptions() {
		StringBuilder res = new StringBuilder();
		res.append("Copy Analysis Comparator Settings:\n");
		res.append("  "+ keyRegionsFile.getName()+"\tKey bed file\n");
		for (File f: testRegionFiles) res.append("  "+ f.getName()+"\tTest regions file\n");
		res.append("  "+ saveDirectory.getName()+"\tSave directory for parsed datasets\n");
		res.append("  "+ takeAbsOfTestRegionScores+"\tTake the absolute value of each region score\n");
		res.append("  "+ minimumBpOverlapKeyTest+"\tMin BP key vs test overlap to score a match\n");
		res.append("  "+ mergeBpGap+"\tMaximim BP gap for merging adjacent regions\n");
		res.append("  "+ minimumMergedRegionSize+"\tMin BP merged region size\n");
		IO.pl(res.toString());
	}

	public static void printDocs(){
		IO.pl("\n" +
				"**************************************************************************************\n" +
				"**                             Copy Analysis Comparator : April 2023                **\n" +
				"**************************************************************************************\n" +
				"Compares test vcf file(s) against a gold standard key of trusted vcf calls. Only calls\n" +
				"that fall in the common interrogated regions are compared. WARNING tabix gzipped files\n" +
				"often fail to parse correctly with java. Seeing odd error messages? Try uncompressing.\n"+
				"Be sure a score is provided in the QUAL field.\n\n" +
				
				"Note, since many test regions can intersect one key region, the calculation of summary scores have been adjusted from standard.\n"+

				"Required Options:\n"+
				"-k VCF file for the key dataset (xxx.vcf(.gz/.zip OK)).\n"+
				"-b Bed file of interrogated regions for the key dataset (xxx.bed(.gz/.zip OK)).\n"+
				"-t VCF file for the test dataset (xxx.vcf(.gz/.zip OK)). May also provide a directory\n" +
				"       containing xxx.vcf(.gz/.zip OK) files to compare.\n"+
				"-c Bed file of interrogated regions for the test dataset (xxx.bed(.gz/.zip OK)).\n"+
				"-p Path to a directory for saving the results.\n"+
				"-m Minimum bp length of overlap to score an intersection as a match, defaults to 1bp\n"+

				"\nExample: java -Xmx10G -jar pathTo/USeq/Apps/VCFComparator -a /NIST/NA12878/key.vcf\n" +
				"       -b /NIST/NA12878/regions.bed.gz -c /EdgeBio/Exome/testHaploCaller.vcf.zip\n" +
				"       -d /EdgeBio/Exome/NimbleGenExomeV3.bed -g -v -s -e -p /CompRes/ \n\n"+

				"**************************************************************************************\n");

	}
}
