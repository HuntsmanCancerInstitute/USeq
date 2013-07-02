package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import trans.main.WilcoxonRankSumTest;
import trans.main.WilcoxonSignedRankTest;
import util.gen.*;
import edu.utah.seq.data.ComparatorPointAscendingScore;
import edu.utah.seq.data.ComparatorPointPosition;
import edu.utah.seq.data.MethylationArraySamplePair;
import edu.utah.seq.data.Point;
import edu.utah.seq.data.PointData;
import edu.utah.seq.data.SmoothingWindow;
import edu.utah.seq.useq.data.Region;

/** Application for looking at user defined regions for changes in methylation.  I've hacked in a non paired analysis on top of the original app.
 * @author Nix
 * */
public class MethylationArrayDefinedRegionScanner {

	//fields
	private LinkedHashMap<String,SmoothingWindow[]> chromosomeRegion;
	private int minimumObservations = 3;
	private boolean skipNoDataRegions = false;
	private boolean performPairedAnalysis = true;
	private File bedFile;
	
	private HashMap<String, MethylationArraySamplePair[]> chromSamplePairs;
	private HashMap<String, int[]> chromPositions = new HashMap<String, int[]>();
	private String chromosome;
	private int[] positions;
	private MethylationArraySamplePair[] pairs;
	private SmoothingWindow[] smoothingWindow;
	private ArrayList<SmoothingWindow> smALAll = new ArrayList<SmoothingWindow>();
	private int pseIndex = 0;
	private int pvalueBHFDRIndex = 1;
	private int numObsIndex = 2;
	private int pseTreatementIndex = 3;
	private int pseControlIndex = 4;
	private File dataDirectory;
	private String treatmentPairedSamples;
	private String controlPairedSamples;
	private String genomeVersion;

	//constructor

	public MethylationArrayDefinedRegionScanner(String[] args){	
		long startTime = System.currentTimeMillis();
		
		//set fields
		processArgs(args);

		//launch
		run();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("\nDone! "+Math.round(diffTime)+" minutes\n");
	}

	public void run(){
		//load data hashes
		System.out.println("Loading samples...");
		loadSamplePairs();

		//for each chromosome
		System.out.println("\nRegion scanning...");
		for (String chrom: chromosomeRegion.keySet()){
			//fetch data
			chromosome = chrom;
			positions = chromPositions.get(chromosome);
			if (positions == null) {
				System.err.println("\t"+chromosome+" Skipping! No data.");
				continue;
			}
			pairs = chromSamplePairs.get(chromosome);
			smoothingWindow = chromosomeRegion.get(chrom);
			System.out.println("\t"+chromosome);

			//scan windows
			scanRegions();
		}
		
		//multiple test correct the pvalues with B&H
		System.out.println("\nConverting wilcoxon pvalues to FDRs with B&H...");
		convertPValuesToFDRs();
		
		//print spreadsheet
		printSpreadSheet();

	}
	
	/**Writes out an excel compatible tab delimited spreadsheet with hyperlinks for IGB.*/
	public void printSpreadSheet(){
		try{
			File file = new File(bedFile.getParentFile(), Misc.removeExtension(bedFile.getName())+"_MADRS.xls");
			PrintWriter out = new PrintWriter (new FileWriter (file));
			//print header line
			out.println("#"+genomeVersion+"_IGBHyperLinks\tChr\tStart\tStop\t#Obs\tPseMedianLog2Ratio\tFDR\tTreatmentPse\tControlPse");
			String url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";
			String tab = "\t";
			//for each region
			int index = 0;
			for (String chr: chromosomeRegion.keySet()){
				for (SmoothingWindow sw : chromosomeRegion.get(chr)){
					float[] scores = sw.getScores();
					if (skipNoDataRegions == true) {
						if (scores == null || (scores[numObsIndex] == 0.0f)) {
							continue;
						}
					}
					//url
					int winStart = sw.getStart() - 10000;
					if (winStart < 0) winStart = 0;
					int winEnd = sw.getStop() + 10000;
					out.print(url+chr+"&start="+winStart+"&end="+winEnd+"\",\""+(++index)+"\")\t");
					//chrom
					out.print(chr); out.print(tab);
					//start
					out.print(sw.getStart()); out.print(tab);
					//stop
					out.print(sw.getStop()); out.print(tab);
					//scores = pse, pvalFDR, # obs
					if (scores == null) {
						out.print("0\t0\t0\t0\t0");
					}
					else {
						//obs
						out.print((int)scores[numObsIndex]); out.print(tab);
						//pse
						out.print(scores[pseIndex]); out.print(tab);
						//FDR
						out.print(scores[pvalueBHFDRIndex]); out.print(tab);
						//treatment pse
						out.print(scores[pseTreatementIndex]); out.print(tab);
						//control pse
						out.print(scores[pseControlIndex]); out.print(tab);
					}
					out.println();
				}
			}
			out.close();
		}catch (Exception e){
			System.out.println("\nError: problem printing spreadsheet report");
			e.printStackTrace();
		}
	}

	private void convertPValuesToFDRs() {
		//only convert if enough observations
		if (smALAll.size() < 10) return;
		
		smoothingWindow = new SmoothingWindow[smALAll.size()];
		smALAll.toArray(smoothingWindow);
		
		//collect pvals for B&H corr
		Point[] pvals = new Point[smoothingWindow.length];
		for (int i=0; i< smoothingWindow.length; i++){
			pvals[i] = new Point(i, smoothingWindow[i].getScores()[pvalueBHFDRIndex]);
		}
		//sort by score
		Arrays.sort(pvals, new ComparatorPointAscendingScore());
		//correct
		Point.benjaminiHochbergCorrect(pvals, 0);
		//sort back to original position
		Arrays.sort(pvals, new ComparatorPointPosition());
		//assign FDRs to pVals
		for (int i=0; i< smoothingWindow.length; i++){
			float scores[] = smoothingWindow[i].getScores();
			scores[pvalueBHFDRIndex] = pvals[i].getScore();
			if (scores[pvalueBHFDRIndex] < 0) scores[pvalueBHFDRIndex] = 0;
		}
	}


	private void scanRegions() {
		ArrayList<Float> treatmentAL = new ArrayList<Float>();
		ArrayList<Float> controlAL = new ArrayList<Float>();
		
		//for each region 
		for (int i=0; i< smoothingWindow.length; i++){
			treatmentAL.clear();
			controlAL.clear();
			
			//fetch indexes
			int startIndex = Num.findClosestStartIndex(positions, smoothingWindow[i].getStart());
			int stopIndex = Num.findClosestEndIndex(positions, smoothingWindow[i].getStop());
			if (stopIndex > positions.length) stopIndex = positions.length;
			
			//add scores from each sample pair, stopIndex is included
			for (int x=0; x< pairs.length; x++){
				pairs[x].fetchScoresByIndex(startIndex, stopIndex, treatmentAL, controlAL);
			}
			
			//any obs?
			if (treatmentAL.size() < minimumObservations) {
				smoothingWindow[i].setScores(null);
				continue;
			}
			
			float[] treatment = Num.arrayListOfFloatToArray(treatmentAL);
			float[] control = Num.arrayListOfFloatToArray(controlAL);
			float pseT = pseudoMedianWithChecker(treatment);
			float pseC = pseudoMedianWithChecker(control);
			float pse = 0;
			float pvalue = 0;
			
			//paired?
			if (performPairedAnalysis){
				double[] fraction = Num.ratio(treatment, control);
				Num.log2(fraction);
				pse = pseudoMedianWithChecker(fraction);
				
				//calculate pval if > 4 pairs
				if (treatment.length > 4){
					WilcoxonSignedRankTest w = new WilcoxonSignedRankTest(treatment, control);
					pvalue = (float)w.getTransformedPValue();
				}
			}
			//non-paired
			else{
				pse = Num.log2(pseT/pseC);
				//calc pval if more than 9 obs
				if (treatment.length > 9 && control.length > 9){
					WilcoxonRankSumTest w = new WilcoxonRankSumTest();
					pvalue = (float)w.test(treatment, control);
				}
			}
			
			//scores = pse, pvalFDR, # obs, treatPse, controlPse            
			float[] scores = new float[5];
			scores[pseIndex] = pse;
			scores[pvalueBHFDRIndex] = pvalue;
			scores[numObsIndex] = treatment.length;
			scores[pseTreatementIndex] = pseT;
			scores[pseControlIndex] = pseC;
			smoothingWindow[i].setScores(scores);
			smALAll.add(smoothingWindow[i]);
			
			/*if (smoothingWindow[i].getStart() == 179456){
				System.out.println(smoothingWindow[i]);
				Misc.printArray(treatment);
				Misc.printArray(control);
				System.exit(0);
			}*/
		}
	}
	
	public static float pseudoMedianWithChecker(float[] vals){
		int len = vals.length;
		if (len == 1) return vals[0];
		if (len == 2) return (vals[0] + vals[1])/2.0f;
		return (float)Num.pseudoMedian(vals);
	}
	
	public static float pseudoMedianWithChecker(double[] vals){
		int len = vals.length;
		if (len == 1) return (float)vals[0];
		if (len == 2) return (float)((vals[0] + vals[1])/2.0);
		return (float)Num.pseudoMedian(vals);
	}

	/*Also hacked this to deal with missing pairs. Ugly!*/
	public void loadSamplePairs(){
		//parse names
		String[] treatmentNames = treatmentPairedSamples.split(",");
		String[] controlNames = controlPairedSamples.split(",");
		if (performPairedAnalysis && (treatmentNames.length != controlNames.length)) Misc.printErrAndExit("\nError: the number of treatment and control samples differ?");
		int numSamples = treatmentNames.length;
		if (controlNames.length > numSamples) numSamples = controlNames.length;

		HashMap<String, ArrayList<MethylationArraySamplePair>> pairsHM = new HashMap<String, ArrayList<MethylationArraySamplePair>>();

		for (int i=0; i< numSamples; i++){
			System.out.print("\t");
			HashMap<String, PointData> pdT = null;
			if (i < treatmentNames.length && treatmentNames[i] !=null) {
				System.out.print(treatmentNames[i]);
				pdT = loadPointData(treatmentNames[i]);
			}
			System.out.print("\t");
			HashMap<String, PointData> pdC = null;
			if (i< controlNames.length && controlNames[i] !=null) {
				System.out.print(controlNames[i]);
				pdC = loadPointData(controlNames[i]);
			}
			System.out.println();
			//set version
			if (genomeVersion == null && pdT !=null) genomeVersion = pdT.values().iterator().next().getInfo().getVersionedGenome();

			HashSet<String> allChroms = new HashSet<String>();
			if (pdT!=null) allChroms.addAll(pdT.keySet());
			if (pdC!=null) allChroms.addAll(pdC.keySet());
			for (String chrom: allChroms){
				float[] tScores = null;
				float[] cScores = null;
				if (pdT!=null) tScores = pdT.get(chrom).getScores();
				if (pdC!=null) cScores = pdC.get(chrom).getScores();
				MethylationArraySamplePair sp = new MethylationArraySamplePair(tScores, cScores, performPairedAnalysis);
				ArrayList<MethylationArraySamplePair> al = pairsHM.get(chrom);
				if (al == null){
					al = new ArrayList<MethylationArraySamplePair>();
					pairsHM.put(chrom, al);
				}
				al.add(sp);
			}
		}

		//convert to []s
		chromSamplePairs = new HashMap<String, MethylationArraySamplePair[]>();
		for (String chrom: pairsHM.keySet()){
			ArrayList<MethylationArraySamplePair> al = pairsHM.get(chrom);
			if (performPairedAnalysis && al.size() != treatmentNames.length) Misc.printErrAndExit("\nError: One or more samples are missing a chromosome PointData set.");
			MethylationArraySamplePair[] sp = new MethylationArraySamplePair[numSamples];
			al.toArray(sp);
			chromSamplePairs.put(chrom, sp);
		}
	}

	

	/**Does a variety of checks, also sets the positions for each chrom.*/
	public HashMap<String, PointData> loadPointData(String name){
		File f = new File (dataDirectory, name);
		if (f.exists() == false) Misc.printErrAndExit("\nError: failed to find sample -> "+f);
		if (f.isDirectory() == false) Misc.printErrAndExit("\nError: sample doesn't appear to be a directory -> "+f);
		HashMap<String, PointData> pd = PointData.fetchPointData (f, null, false);
		if (pd == null || pd.keySet().size() == 0) Misc.printErrAndExit("\nError: failed to load PointData from -> "+f);
		//check and or set positions
		for (String chrom: pd.keySet()){
			int[] existingPos = chromPositions.get(chrom);
			int[] pdPos = pd.get(chrom).getPositions();
			if (existingPos == null) chromPositions.put(chrom, pdPos);
			else {
				if (existingPos.length != pdPos.length) Misc.printErrAndExit("\nError: lengths of PointData don't match.");
				for (int i=0; i< existingPos.length; i++){
					if (existingPos[i] != pdPos[i]) Misc.printErrAndExit("\nError: positions in PointData don't match.");
				}
			}
		}
		return pd;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new MethylationArrayDefinedRegionScanner(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		bedFile = null;
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bedFile = new File(args[++i]); break;
					case 'd': dataDirectory = new File (args[++i]); break;
					case 't': treatmentPairedSamples = args[++i]; break;
					case 'c': controlPairedSamples = args[++i]; break;
					case 'z': skipNoDataRegions = true; break;
					case 'o': minimumObservations = Integer.parseInt(args[++i]); break;
					case 'n': performPairedAnalysis = false; break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printErrAndExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printErrAndExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}	
		}
		//look for data dir
		if (dataDirectory == null || dataDirectory.isDirectory() == false) Misc.printErrAndExit("\nError: cannot find your data directory containing sample bar file directories.\n");

		//samples
		if (treatmentPairedSamples == null || controlPairedSamples == null) Misc.printErrAndExit("\nError: please enter at least one paired treatment control sample set to contrast.\n");
		
		//load regions
		if (bedFile == null) Misc.printExit("\nError: please provide a text bed file (tab delimited: chr start stop) of regions to score for differential methylation.\n");
		HashMap<String, Region[]>cr = Region.parseStartStops(bedFile, 0, 0, 0);
		chromosomeRegion = new LinkedHashMap<String, SmoothingWindow[]>();
		for (String chr: cr.keySet()){
			Region[] regions = cr.get(chr);
			SmoothingWindow[] sw = new SmoothingWindow[regions.length];
			for (int i=0; i< regions.length; i++){
				sw[i] = new SmoothingWindow (regions[i].getStart(), regions[i].getStop(), new float[]{0,0,0,0,0,0});
			}
			chromosomeRegion.put(chr, sw);
		}
	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                    Methylation Array Defined Region Scanner: July 2013           **\n" +
				"**************************************************************************************\n" +
				"MADRS takes paired sample PointData representing beta values (0-1) from arrays and\n" +
				"a list of regions to score for differential methylation using a B&H corrected Wilcoxon\n" +
				"signed rank test and pseudo median of the paired log2(treat/control) ratios. Pairs\n" +
				"containing a zero value are ignored. It generates a spreadsheet of statistics for each\n" +
				"region. If a non-paired analysis is selected, a Wilcoxon rank sum test and\n" +
				"log2(pseT/pseC) are calculated on each region. Note this is a very underpowered test\n" +
				"requiring >30 observations to see any significant FDRs.\n" +

				"\nRequired Options:\n"+
				"-b A bed file of regions to score (tab delimited: chr start stop ...)\n"+
				"-d Path to a directory containing individual sample PointData directories, each of\n"+
				"      which should contain chromosome split bar files (e.g. chr1.bar, chr2.bar, ...)\n"+
				"-t Names of the treatment sample directories in -d, comma delimited, no spaces.\n"+
				"-c Ditto but for the control samples, the ordering is critical and describes how to\n"+
				"      pair the samples for a paired analysis.\n"+
				"-o Minimum number paired observations in window, defaults to 3.\n" +
				"-z Skip printing regions with less than minimum observations.\n"+
				"-n Run a non-paired analysis where t and c are treated as groups and pooled. Uneven\n" +
				"      numbers of t and c are allowed.\n"+

				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/MethylationArrayDefinedRegionScanner \n" +
				"     -v H_sapiens_Feb_2009 -d ~/MASS/Bar/ -t Early1,Early2,Early3\n" +
				"     -c Late1,Late2,Late3 \n\n" +

		"**************************************************************************************\n");
	}
}
