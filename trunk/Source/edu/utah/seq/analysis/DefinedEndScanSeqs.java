package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.bio.annotation.Bed;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.*;
import util.gen.*;
import edu.utah.seq.data.*;


/**Takes PointData and generates defined window scan statistics. See the help screen info below.
 * @author Nix
 * */
public class DefinedEndScanSeqs {

	//user defined fields
	private File[] treatmentPointDirs;
	private File[] controlPointDirs;
	private File[] treatmentSpliceFiles;
	private File[] controlSpliceFiles;
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private File bedFile;
	private File refSeqFile;
	private int peakShift = 0;
	private boolean useWeightedReads = false;
	private boolean findReducedRegions = true;
	private int minNumObs2ScoreSpliceJunctions = 100;
	private float minSplicePVal = 10;

	//internal fields
	private HashMap<String,UCSCGeneLine[]> geneModels;
	private UCSCGeneLine[] allGeneLines;
	private ArrayList<UCSCGeneLine> geneLinesWithReadsAL = new ArrayList<UCSCGeneLine>();
	private UCSCGeneLine[] geneLinesWithReads;
	private HashMap<String, File> splitSpliceTreatmentFiles;
	private HashMap<String, File> splitSpliceControlFiles;
	private HashMap<String, Integer> spliceTreatmentCounts;
	private HashMap<String, Integer> spliceControlCounts;
	private int totalNumberJunctionsScored = 0;
	private int halfPeakShift = 0;
	private double numberTreatmentObservations;
	private double numberControlObservations;
	private double expectedFractionUp;
	private double expectedFractionDown;
	private double scalarTC;
	private double scalarCT;
	private double millionMappedTreatmentReads;
	private double millionMappedControlReads;
	private HashMap<String,PointData[]> treatmentPlusPointData;
	private HashMap<String,PointData[]> treatmentMinusPointData;
	private HashMap<String,PointData[]> controlPlusPointData;
	private HashMap<String,PointData[]> controlMinusPointData;
	private int numberCachedBinPVals = 500;
	private float[][] pValsDataUp;
	private float[][] pValsDataDown;
	private float[][] pValsSkew;
	private float[] controlThresholds;
	private int[] controlNumberGenes;
	private String genomeVersion;

	//by chromosome
	private String chromosome;
	private PointData treatmentChromPlus = null;
	private PointData treatmentChromMinus = null;	
	private PointData controlChromPlus = null;
	private PointData controlChromMinus = null;	


	//constructor
	public DefinedEndScanSeqs(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//load gene models
		loadGeneModels();

		//split splice data?	
		if (treatmentSpliceFiles != null && controlSpliceFiles != null){
			System.out.println("Splitting splice junction bed files...");
			splitSpliceTreatmentFiles = splitByChrom(treatmentSpliceFiles);
			splitSpliceControlFiles = splitByChrom(controlSpliceFiles);
		}

		//fetch counts
		System.out.println("Calculating read count stats and binomial p-value look up tables...");
		calculateReadCountStatistics();

		//for each chromosome of gene models
		System.out.print("Scanning regions by chromosome");
		Iterator<String> it = geneModels.keySet().iterator();
		while (it.hasNext()){
			chromosome = it.next();
			scanChromosome();
		}
		System.exit(0);
		System.out.println();

		//control data? run FDR estimations
		if (controlPointDirs != null){
			//load geneLines with those that have observations
			geneLinesWithReads = new UCSCGeneLine[geneLinesWithReadsAL.size()];
			geneLinesWithReadsAL.toArray(geneLinesWithReads);

			//convert binomial pvalues to FDRs
			System.out.println("Converting fuzzy binomial p-values to q-value FDRs in R...");

			//calculate reduce regions and EmpFDRs?
			if (findReducedRegions == false){
				convertPValuesToQValuesFuzzedUp();

				System.out.println("Estimating empirical FDRs...");
				//process controls for empFDRs
				processControlsForEmpFDR();

				//first for regions with enrichment, 
				//scores = pVal,upPVal,downPVal,skew,chiSqrPVal,bps,tSumPlus,tSumMinus,cSumPlus,cSumMinus, realND,mockND,upEmpFDR
				int scoreIndex = 10;
				int fdrIndex = 12;			
				estimateEmpiricalFDRs(geneLinesWithReads, scoreIndex, fdrIndex);

			}
			else convertPValuesToQValuesFuzzedUpDown();

			//calculate Fisher's exact on distribution of reads between the exons
			System.out.println("Estimating read distribution, bonferroni corrected, chi-square p-values in R (very slow) ...");
			estimateDifferencesInReadDistributions();

			//                    0     1      2      3       4        5     6         7         8       9         10     11      12     
			//starting scores = pVal,upPVal,downPVal,skew,chiSqrPVal, bps,tSumPlus,tSumMinus,cSumPlus,cSumMinus, realND,mockND,upEmpFDR

			//                    0        1      2       3         4        5             6                 7         8         9       10        11        12
			//collapse scores to pVal, qValFDR, eFDR, pValSkew, chiSqrPVal, bps, Log2((sumT+1)/(sumC+1)), tSumPlus, tSumMinus, tRPKM, cSumPlus, cSumMinus, cRPKM
			for (int i=0; i< allGeneLines.length; i++){
				float[] scores = allGeneLines[i].getScores();
				float[] col = null;
				//any reads?
				float total = scores[6]+scores[7]+scores[8]+scores[9];
				if (total !=0){
					col = new float[13];
					//pval
					col[0] = scores[0];
					//qValueFDR
					if (scores[1]>=scores[2]) col[1] = scores[1];
					else col[1] = scores[2] * -1;
					//eFDR
					if (findReducedRegions) col[2] = 0;
					else col[2] = scores[11];
					//skew
					col[3] = scores[3];
					//chi sqr
					col[4] = scores[4];
					//bps of exonic seq
					col[5] = scores[5];
					//log2(ratio)
					col[6] = calculateLog2Ratio(scores[6]+scores[7], scores[8]+scores[9], scalarTC, scalarCT);
					//sums
					col[7] = scores[6];
					col[8] = scores[7];
					col[9] = calculateRPKM(millionMappedTreatmentReads, scores[5], scores[6]+scores[7]);

					col[10] = scores[8];
					col[11] = scores[9];
					col[12] = calculateRPKM(millionMappedControlReads, scores[5], scores[8]+scores[9]);
				}
				else col = new float[]{0};
				allGeneLines[i].setScores(col);
			}
			//calculate binomial p-values for splice junctions?
			if (totalNumberJunctionsScored != 0) {
				System.out.println("Scoring splice junctions...");
				scoreSpliceJunctions();
			}
		}
		//sort by score[0] pValue
		System.out.println("Printing results...");
		Arrays.sort(allGeneLines, new UCSCGeneLineComparatorScoreBigToSmall(0));

		//print report
		printGeneModels();

		//print gff 
		//printGeneModelsAsBed();

		//delete temp files?
		if (splitSpliceTreatmentFiles != null){
			it = splitSpliceControlFiles.keySet().iterator();
			while (it.hasNext()) splitSpliceControlFiles.get(it.next()).deleteOnExit();
			it = splitSpliceTreatmentFiles.keySet().iterator();
			while (it.hasNext()) splitSpliceTreatmentFiles.get(it.next()).deleteOnExit();
		}


		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");

	}

	/**Calculates the reads per kb per million mapped reads 
	 * # Observed reads in the region/ bp size of the region / 1000/ total number reads/ 1000000 */
	public float calculateRPKM(double millionTotalMappedReads, double interrogatedRegionBPSize, double numberObservedReadsInRegion){
		double exonicBasesPerKB = interrogatedRegionBPSize/1000;
		double rpkm = numberObservedReadsInRegion/exonicBasesPerKB/millionTotalMappedReads;
		return new Double(rpkm).floatValue();
	}

	public void estimateDifferencesInReadDistributions(){
		//find genes with multiple exons and max number exons
		int maxNumberExons = -1;
		ArrayList<UCSCGeneLine> al = new ArrayList<UCSCGeneLine>();
		for (int i=0; i< geneLinesWithReads.length; i++){
			int numEx = geneLinesWithReads[i].getExons().length;
			if (numEx > 1 && geneLinesWithReads[i].getExonCounts() != null) {
				al.add(geneLinesWithReads[i]);
				if (numEx > maxNumberExons) maxNumberExons = numEx;
			}
		}
		UCSCGeneLine[] genesWithExonsAndReads = new UCSCGeneLine[al.size()];
		al.toArray(genesWithExonsAndReads);

		//collect counts
		int[][] treatment = new int[genesWithExonsAndReads.length][maxNumberExons];
		int[][] control = new int[genesWithExonsAndReads.length][maxNumberExons];

		for (int i=0; i< genesWithExonsAndReads.length; i++){
			float[][] tc = genesWithExonsAndReads[i].getExonCounts();
			Arrays.fill(treatment[i], -1);
			Arrays.fill(control[i], -1);
			int[] t = Num.convertToInt(tc[0]);
			int[] c = Num.convertToInt(tc[1]);
			System.arraycopy(t, 0, treatment[i], 0, t.length);
			System.arraycopy(c, 0, control[i], 0, c.length);
		}

		//estimate chi-square pvalues using R for resolution of extreemly small p-values, radiculously slow
		double[] pVals = Num.chiSquareIndependenceTest(treatment, control, saveDirectory, fullPathToR, true);

		//bonferroni correction
		float bc = (float)Num.minus10log10(genesWithExonsAndReads.length);

		//add back
		for (int i=0; i< genesWithExonsAndReads.length; i++){
			float[][] tc = genesWithExonsAndReads[i].getExonCounts();
			//set corrected p-value 
			//scores = pVal,upPVal,downPVal,skew,chiSqrPVal, bps, tSumPlus,tSumMinus,cSumPlus,cSumMinus, realND,mockND,upEmpFDR
			float[] scores = genesWithExonsAndReads[i].getScores();
			scores[4] = (float) pVals[i] + bc;
			if (scores[4] < 0) scores[4] = 0;
		}
	}

	public void scoreSpliceJunctions(){
		//collect #t #c expect
		float[][] toR = new float[totalNumberJunctionsScored][3];
		//for each gene
		int index = 0;
		for (int i=0; i< geneLinesWithReads.length; i++){
			SpliceJunction[] sjs = geneLinesWithReads[i].getSpliceJunctions();
			if (sjs == null) continue;
			//get scores   pVal, qValFDR, eFDR, pValSkew, chiSqrPVal, bps, Log2((sumT+1)/(sumC+1)), tSumPlus, tSumMinus, cSumPlus, cSumMinus
			float[] scores = geneLinesWithReads[i].getScores();
			//calculate expect, add one to control for zero T
			double numT = scores[7] + scores[8] +1;
			double numC = scores[9] + scores[10];
			float e = (float)(numT/(numT+numC));
			double scTC = numT/numC;
			double scCT = numC/numT;
			//for each junction
			for (int j=0; j< sjs.length; j++){
				toR[index++] = new float[] {sjs[j].getNumberTreatmentObservations(), sjs[j].getNumberControlObservations(), e};
				//calculate and set log2Ratio
				float log2Ratio = calculateLog2Ratio(sjs[j].getNumberTreatmentObservations(), sjs[j].getNumberControlObservations(), scTC, scCT);
				sjs[j].setLog2Ratio(log2Ratio);
			}
		}
		//calculate binomial pvalues in R for up and down
		double[][] pvalues = Num.binomialPValues(toR, saveDirectory, fullPathToR);
		//assign pvals
		index = 0;
		double bonCorr = Num.minus10log10(totalNumberJunctionsScored);
		for (int i=0; i< geneLinesWithReads.length; i++){
			SpliceJunction[] sjs = geneLinesWithReads[i].getSpliceJunctions();
			if (sjs == null) continue;
			//for each junction
			boolean nonZero = false;
			for (int j=0; j< sjs.length; j++){
				double pVal = pvalues[index][0];
				if (pVal< pvalues[index][1]) pVal = pvalues[index][1];
				index++;
				//multiple test correct 
				pVal = pVal + bonCorr;
				if (pVal < 0) continue;
				sjs[j].setBinomialPValue((float) pVal);
				nonZero = true;
			}
			if (nonZero == false) geneLinesWithReads[i].setSpliceJunctions(null);
			else Arrays.sort(sjs);
		}
	}

	public void printGeneModels(){
		try {
			File results = new File (saveDirectory, "drssResults.xls");
			PrintWriter out = new PrintWriter (new FileWriter( results));

			String misc = "GenomeVersion="+genomeVersion+ ", TotalTreatObs="+(int)numberTreatmentObservations;
			String url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";
			if (allGeneLines[0].getDisplayName() != null) out.print("#DisplayName\tName\tChr\tStrand\tStart\tStop\t");
			else out.print("#Name\tChr\tStrand\tStart\tStop\t");
			if (controlPointDirs != null){
				misc = misc+ ", TotalCtrlObs="+(int)numberControlObservations;
				out.print("pVal\tqValFDR\teFDR\tpValSkew\tpValDiffDist\tTotalRegionBPs\tLog2((sumT+1)/(sumC+1))\ttSumPlus\ttSumMinus\ttRPKM\tcSumPlus\tcSumMinus\tcRPKM\t");
				if (totalNumberJunctionsScored !=0) out.println("SpliceJunctions->\tName\tBonCorrBinPVal\tLog2((sumT+1)/(sumC+1))\ttSum\tcSum\t"+misc);
				else out.println(misc);
				for (int i=0; i< allGeneLines.length; i++){
					String name;
					if (allGeneLines[i].getDisplayName() !=null) name = allGeneLines[i].getDisplayName();
					else name = allGeneLines[i].getName();
					//url
					int start = allGeneLines[i].getTxStart() - 10000;
					if (start < 0) start = 0;
					int end = allGeneLines[i].getTxEnd() + 10000;
					out.print(url+allGeneLines[i].getChrom()+"&start="+start+"&end="+end+"\",\""+name+"\")\t");
					//print second text?
					if (allGeneLines[i].getDisplayName() !=null) out.print(allGeneLines[i].getName()+"\t");
					//scores 
					float[] s = allGeneLines[i].getScores();
					if (s.length == 1) out.print(allGeneLines[i].coordinates());
					else out.print(allGeneLines[i].coordinates()+"\t"+Misc.floatArrayToString(s, "\t"));
					//splices
					SpliceJunction[] sjs = allGeneLines[i].getSpliceJunctions();
					if (sjs!=null){
						out.print("\t");
						for (int x=0; x< sjs.length; x++){
							if (sjs[x].getBinomialPValue() < minSplicePVal) continue;
							out.print("\t");
							out.print(sjs[x]);
						}
					}
					out.println();
				}
			}
			else{
				out.println("Sum\tSum+\tSum-\tTotalRegionBPs\tRPKM\t"+misc);
				for (int i=0; i< allGeneLines.length; i++){
					String name;
					if (allGeneLines[i].getDisplayName() !=null) name = allGeneLines[i].getDisplayName();
					else name = allGeneLines[i].getName();
					//url
					int start = allGeneLines[i].getTxStart() - 10000;
					if (start < 0) start = 0;
					int end = allGeneLines[i].getTxEnd() + 10000;
					out.print(url+allGeneLines[i].getChrom()+"&start="+start+"&end="+end+"\",\""+name+"\")\t");
					//print second text?
					if (allGeneLines[i].getDisplayName() !=null) out.print(allGeneLines[i].getName()+"\t");
					float[] s = allGeneLines[i].getScores();
					out.println(allGeneLines[i].coordinates()+"\t"+s[0]+"\t"+ s[1]+"\t"+ s[2]+"\t"+ s[3]+"\t"+ s[4]);
				}
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}


	public void loadGeneModels(){
		//load gene models from refFlat for refSeq UCSC gene table
		UCSCGeneModelTableReader reader = null;
		if (refSeqFile != null){
			reader = new UCSCGeneModelTableReader(refSeqFile, 0);
			reader.splitByChromosome();
			geneModels = reader.getChromSpecificGeneLines();
			allGeneLines = reader.getGeneLines();
		}
		//or from bed file
		else if (bedFile != null) {
			Bed[] bed = Bed.parseFile(bedFile, 0, 0);
			allGeneLines = new UCSCGeneLine[bed.length];
			boolean addName = bed[0].getName().trim().equals("");
			for (int i=0; i< bed.length; i++){
				if (addName) bed[i].setName((i+1)+"");
				allGeneLines[i] = new UCSCGeneLine(bed[i]);
			}
			reader = new UCSCGeneModelTableReader();
			reader.setGeneLines(allGeneLines);
			reader.splitByChromosome();
			geneModels = reader.getChromSpecificGeneLines();
		}
		if (geneModels == null | geneModels.size() == 0) Misc.printExit("\nProblem loading your USCS gene model table or bed file.\n");
		//check ordering
		if (reader.checkStartStopOrder() == false) Misc.printExit("\nOne of your regions's coordinates are reversed. Check that each start is less than the stop.\n");
		//check number of regions
		if (reader.getGeneLines().length < 100) {
			Misc.printExit("\nToo few regions to scan! DRSS needs at least 100 regions to properly estimate FDRs.\n");
		}
	}

	/**Fetches the names of all the chromosomes in the data.*/
	public String[] fetchAllChromosomes(){
		HashSet<String> c = new HashSet<String>();
		Iterator<String> it = treatmentPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = treatmentMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		if(controlPointDirs != null) {
			it = controlPlusPointData.keySet().iterator();
			while (it.hasNext()) c.add(it.next());
			it = controlMinusPointData.keySet().iterator();
			while (it.hasNext()) c.add(it.next());
		}
		return Misc.hashSetToStringArray(c);
	}

	/**Fetchs the data for a particular chromosome.*/
	public void fetchData(){
		//merge treatment
		treatmentChromPlus = null;
		if (treatmentPlusPointData.containsKey(chromosome)) treatmentChromPlus = PointData.combinePointData(treatmentPlusPointData.get(chromosome), true);
		treatmentChromMinus = null;
		if (treatmentMinusPointData.containsKey(chromosome)) treatmentChromMinus = PointData.combinePointData(treatmentMinusPointData.get(chromosome), true);
		//merge control
		if (controlPointDirs != null){
			controlChromPlus = null;
			if (controlPlusPointData.containsKey(chromosome)) controlChromPlus = PointData.combinePointData(controlPlusPointData.get(chromosome), true);
			controlChromMinus = null;
			if (controlMinusPointData.containsKey(chromosome)) controlChromMinus = PointData.combinePointData(controlMinusPointData.get(chromosome), true);
			//load splice data hashes?
			if (treatmentSpliceFiles != null && controlSpliceFiles != null){
				spliceTreatmentCounts = loadAndMergeSplices(splitSpliceTreatmentFiles);
				spliceControlCounts = loadAndMergeSplices(splitSpliceControlFiles);
				if (spliceTreatmentCounts == null || spliceControlCounts == null) {
					spliceTreatmentCounts = null;
					spliceControlCounts = null;
				}
			}
		}
		//assign genome
		if (genomeVersion == null){
			if (treatmentChromPlus != null) genomeVersion = treatmentChromPlus.getInfo().getVersionedGenome();
			else if (treatmentChromMinus != null) genomeVersion = treatmentChromMinus.getInfo().getVersionedGenome();
		}
	}

	/**Window scans a chromosome collecting read count data and calculating binomial p-values.*/
	public void scanChromosome(){
		//System.out.print(".");
		//fetch data
		fetchData();
		//fetch, shift, and merge all positions from the treatment and control
		shiftStripPointData();
		//scan
		if (controlPointDirs !=null) {
			calculateBinomialPValues();
			if (findReducedRegions == false) calculateBinomialPValsForEmpiricalFDRs();
		}
		else sumTreatmentReadsForEnds(50, 1);
	}


	/**Calculates a log2( (tSum+1)/(cSum+1) ) on linearly scaled tSum and cSum based on the total observations.*/
	public static float calculateLog2Ratio( double tSum, double cSum, double scalarTC, double scalarCT){
		double t;
		double c;
		if (tSum !=0 ) {
			t = tSum * scalarCT;
			c = cSum;
		}
		else {
			c = cSum * scalarTC;
			t = tSum;
		}
		double ratio = (t+1)/(c+1);
		return (float)Num.log2(ratio);
	}



	/**Calculates the number of genes that exceed pVal thresholds for control vs control data.*/
	private void processControlsForEmpFDR (){
		//scores = pVal, upPVal,downPVal,skew, chiSqrPVal, bps, tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
		HashSet<Float> controlScoresHM = new HashSet<Float>();
		ArrayList<Float> allScoresAL = new ArrayList<Float>();
		for (int j=0; j< geneLinesWithReads.length; j++){
			float mockPVal = geneLinesWithReads[j].getScores()[11];
			if (mockPVal > 0) {
				controlScoresHM.add(new Float(mockPVal));
				allScoresAL.add(new Float(mockPVal));
			}
		}
		//Make and sort scores
		controlThresholds = Num.hashSetToFloat(controlScoresHM);
		Arrays.sort(controlThresholds);
		controlNumberGenes = new int[controlScoresHM.size()];
		float[] allScores = Num.arrayListOfFloatToArray(allScoresAL);
		Arrays.sort(allScores);

		//calculate number of genes that meet or exceed the thresholds
		for (int i=0; i<controlThresholds.length; i++){
			float test = controlThresholds[i];
			int index = Num.findClosestIndexToValue(allScores, test);
			int numGenes = allScores.length - index;
			if (numGenes <=0) break;
			controlNumberGenes[i] = numGenes;
		}
	}



	/**Converts binomial pvalues in SmoothingWindowInfo to qvalues using Storey's method.
	 * Adds random noise to windows with counts <21
	 **/
	private void convertPValuesToQValuesFuzzedUp(){
		//make up and down fuzzy pvalue generators
		FuzzyBinomialPValueGenerator upFuzz = new FuzzyBinomialPValueGenerator( saveDirectory, fullPathToR, expectedFractionUp);

		//collect -10Log10(pvalues)
		float[] up = new float[geneLinesWithReads.length+1];
		float maxPValUp = 0;
		for (int i=0; i< geneLinesWithReads.length; i++){
			//scores = pVal, upPVal,downPVal,skew,chiSqrPVal,bps,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
			float[] scores = geneLinesWithReads[i].getScores();
			//low number of obs? then fuzz pval
			int totalObs = Math.round(scores[6]+scores[7]+scores[8]+scores[9]);
			if (totalObs < 21){
				up[i] = (float)upFuzz.fetchFuzzyPVal((int)(scores[6]+scores[7]), (int)(scores[8]+scores[9]), true);
				//set max
				if (scores[1] > maxPValUp) maxPValUp = scores[1];
			}
			else{
				up[i] = scores[1];
			}
		}

		//set max as last
		up[geneLinesWithReads.length] = maxPValUp;

		//calc qvals		
		float[] upQVals = Num.qValueFDR(saveDirectory, up, fullPathToR, true, true);

		//transform total
		float transTotal = (float)Num.minus10log10(geneLinesWithReads.length);

		//assign to windows and multiple test correct the skew pvalue
		ArrayList<UCSCGeneLine> toSetQVal = new ArrayList<UCSCGeneLine>();
		for (int i=0; i< geneLinesWithReads.length; i++){
			//scores = pVal, upPVal,downPVal,skew,chiSqrPVal,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
			float[] scores = geneLinesWithReads[i].getScores();
			int totalObs = Math.round(scores[5]+scores[6]+scores[7]+scores[8]);
			if (totalObs < 21) {
				toSetQVal.add(geneLinesWithReads[i]);
				upFuzz.addAdjPValQVal(up[i],upQVals[i]);
			}
			else {
				scores[1] = upQVals[i];
				scores[2] = 0;
			} 
			//multiple test correct skew binomial pvalue
			scores[3] = scores[3] + transTotal;
			if (scores[3] < 0) scores[3] = 0;
		}

		//set max in FuzzyBinomialPValGenerators
		upFuzz.addAdjPValQVal(up[geneLinesWithReads.length],upQVals[geneLinesWithReads.length]);

		//make qvalue calculators
		upFuzz.generateSortedPAndQValArrays();
		//convert toSeqQval
		UCSCGeneLine[] toFixQVal = new UCSCGeneLine[toSetQVal.size()];
		toSetQVal.toArray(toFixQVal);
		//load qVals
		upFuzz.loadQValues(1, toFixQVal);
	}

	/**Converts binomial pvalues to qvalues using Storey's method.
	 * Adds random noise to windows with counts <21
	 **/
	private void convertPValuesToQValuesFuzzedUpDown(){
		//make up and down fuzzy pvalue generators
		FuzzyBinomialPValueGenerator upFuzz = new FuzzyBinomialPValueGenerator( saveDirectory, fullPathToR, expectedFractionUp);
		FuzzyBinomialPValueGenerator downFuzz = new FuzzyBinomialPValueGenerator( saveDirectory, fullPathToR, expectedFractionDown);

		//collect -10Log10(pvalues)
		float[] up = new float[geneLinesWithReads.length +1];
		float[] down = new float[geneLinesWithReads.length +1];
		float maxPValUp = 0;
		float maxPValDown = 0;
		for (int i=0; i< geneLinesWithReads.length; i++){
			//scores = pVal, upPVal,downPVal,skew,chiSqrPVal,bps,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
			float[] scores = geneLinesWithReads[i].getScores();
			//low number of obs? then fuzz pval
			int totalObs = Math.round(scores[6]+scores[7]+scores[8]+scores[9]);
			if (totalObs < 21){
				up[i] = (float)upFuzz.fetchFuzzyPVal((int)(scores[6]+scores[7]), (int)(scores[8]+scores[9]), true);
				down[i] = (float)downFuzz.fetchFuzzyPVal((int)(scores[8]+scores[9]), (int)(scores[6]+scores[7]), true);
				//set max
				if (scores[1] > maxPValUp) maxPValUp = scores[1];
				if (scores[2] > maxPValDown) maxPValDown = scores[2];
			}
			else{
				up[i] = scores[1];
				down[i] = scores[2];
			}
		}

		//set max as last
		up[geneLinesWithReads.length] = maxPValUp;
		down[geneLinesWithReads.length] = maxPValDown;

		//calc qvals
		float[] upQVals = Num.qValueFDR(saveDirectory, up, fullPathToR, true, true);
		float[] downQVals = Num.qValueFDR(saveDirectory, down, fullPathToR, true, true);

		//transform total
		float transTotal = (float)Num.minus10log10(geneLinesWithReads.length);

		//assign to windows and multiple test correct the skew pvalue
		ArrayList<UCSCGeneLine> toSetQVal = new ArrayList<UCSCGeneLine>();
		for (int i=0; i< geneLinesWithReads.length; i++){
			//scores = pVal, upPVal,downPVal,skew,chiSqrPVal,bps,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
			float[] scores = geneLinesWithReads[i].getScores();
			int totalObs = Math.round(scores[6]+scores[7]+scores[8]+scores[9]);
			if (totalObs < 21) {
				toSetQVal.add(geneLinesWithReads[i]);
				upFuzz.addAdjPValQVal(up[i],upQVals[i]);
				downFuzz.addAdjPValQVal(down[i],downQVals[i]);
			}
			else {
				scores[1] = upQVals[i];
				scores[2] = downQVals[i];
			} 
			//multiple test correct skew binomial pvalue
			scores[3] = scores[3] + transTotal;
			if (scores[3] < 0) scores[3] = 0;
		}

		//set max in FuzzyBinomialPValGenerators
		upFuzz.addAdjPValQVal(up[geneLinesWithReads.length],upQVals[geneLinesWithReads.length]);
		downFuzz.addAdjPValQVal(down[geneLinesWithReads.length],downQVals[geneLinesWithReads.length]);
		//make qvalue calculators
		upFuzz.generateSortedPAndQValArrays();
		downFuzz.generateSortedPAndQValArrays();
		//convert toSeqQval
		UCSCGeneLine[] toFixQVal = new UCSCGeneLine[toSetQVal.size()];
		toSetQVal.toArray(toFixQVal);
		//load qVals
		upFuzz.loadQValues(1, toFixQVal);
		downFuzz.loadQValues(2, toFixQVal);
	}


	/**Calculates binomial p-values for emp FDR estimation.
	 * Only call this after you are done with the PointData for a particular chromosome!
	 * This method subsamples and replaces the data!*/
	private void calculateBinomialPValsForEmpiricalFDRs(){
		//subsample PointData so T: Ca: Cb, on a chromosome basis not genome basis, theoretical problem here?! 
		PointData[] pdTCC = subSamplePointDataForEmpFDR();

		//make arrays to calculate binomial p-vals with R
		ArrayList<UCSCGeneLine> forBinom = new ArrayList<UCSCGeneLine>();

		//get UCSCGeneLine[]
		UCSCGeneLine[] genes = geneModels.get(chromosome);

		//for each gene 
		for (int i=0; i< genes.length; i++){
			float[] scores = genes[i].getScores();
			//calculate totals
			ExonIntron[] exons = genes[i].getExons();
			//fetch scores
			float tSum = 0; 
			float c1Sum = 0;
			float c2Sum = 0;
			for (int j=0; j< exons.length; j++){
				int start = exons[j].getStart();
				int stop = exons[j].getEnd();
				tSum += pdTCC[0].sumScoreBP(start, stop); 
				c1Sum += pdTCC[1].sumScoreBP(start, stop); 
				c2Sum += pdTCC[2].sumScoreBP(start, stop);  
			}
			//any scores?
			if ((tSum+c1Sum+c2Sum) == 0){
				scores[9] = 0;
				scores[10] = 0;
				scores[11] = 0;
			}
			else {
				//calc pval from cache?
				//yes
				if (tSum < numberCachedBinPVals && c1Sum < numberCachedBinPVals && c2Sum< numberCachedBinPVals){
					//calc pvals
					int t = Math.round(tSum);
					int c1 = Math.round(c1Sum);
					int c2 = Math.round(c2Sum);
					//calc realPVal
					float realPVal = pValsSkew[t][c2];
					//calc mockPVal
					float mockPVal = pValsSkew[c1][c2];
					//trim pVals sig fig to one, xxx.x
					realPVal = multBy10AndRound(realPVal);
					mockPVal = multBy10AndRound(mockPVal);
					//scores = pVal, upPVal,downPVal,skew,chiSqrPVal,bps,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR	
					scores[10] = realPVal;
					scores[11] = mockPVal;
					scores[12] = 0;	
				}
				//no so save for calculating with a call to R
				else {
					scores[10] = tSum;
					scores[11] = c1Sum;
					scores[12] = c2Sum;
					forBinom.add(genes[i]);
				}
			}
		}
		if (forBinom.size() != 0) calculatePValsFromRForEmpFDR(forBinom);
	}

	/**Calcs p-vals from R for EmpFDR that aren't in the look up tables.*/
	private void calculatePValsFromRForEmpFDR (ArrayList<UCSCGeneLine> al){
		//for each region calculate pvals
		int num = al.size();
		int[][] obs = new int[num*2][2];
		int index = 0;
		for (int i=0; i< num; i++){
			float[] scores = al.get(i).getScores();
			//scores = pVal,upPVal,downPVal,skew,chiSqrPVal,bps,tSumPlus,tSumMinus,cSumPlus,cSumMinus,tSum,c1Sum,c2Sum,upEmpFDR
			int t = Math.round(scores[10]);
			int c1 = Math.round(scores[11]);
			int c2 = Math.round(scores[12]);
			//t vs c2
			obs[index++] = new int[]{t,c2};
			//c1 vs c2
			obs[index++] = new int[]{c1,c2};

		}
		//fetch pvalues
		double[] pVals = Num.binomialPValues(0.5, saveDirectory, obs, fullPathToR, true, false,true);
		if (pVals == null || pVals.length != (2*num)) Misc.printErrAndExit("\nProblem fetching binomial pvalues from R for EmpFDR. Check that R really is at "+fullPathToR);
		//assign
		index = 0;
		for (int i=0; i< num; i++){
			float realPVal = (float) pVals[index++];
			float mockPVal = (float) pVals[index++];
			//trim pVals sig fig to one, xxx.x
			realPVal = multBy10AndRound(realPVal);
			mockPVal = multBy10AndRound(mockPVal);
			//scores = pVal, upPVal,downPVal,skew,chiSqrPVal,bps,tSumPlus,tSumMinus,cSumPlus,cSumMinus,realPVal,mockPVal,upEmpFDR
			float[] scores = al.get(i).getScores();	
			int t = (int)scores[10];
			int c1 = (int)scores[11];
			int c2 = (int)scores[12];
			scores[10] = realPVal;
			scores[11] = mockPVal;
			scores[12] = 0;
		}
	}

	/**Subsamples point data for empirical FDR calculations. t=c1=c2 based on whole genome, not individual chromosome.*/
	private PointData[] subSamplePointDataForEmpFDR(){
		//count observations
		int halfCtrlObs = (int)Math.round(numberControlObservations/2);
		//merge PointData
		PointData treatmentPD = PointData.combinePointData(new PointData[]{treatmentChromPlus,treatmentChromMinus}, true);
		PointData controlPD = PointData.combinePointData(new PointData[]{controlChromPlus,controlChromMinus}, true);
		//make c 2x t or t 1/2 c

		//match number of treatments to half the number of controls
		//subsample treatment?
		if (halfCtrlObs < numberTreatmentObservations) {
			double totalNumberToMatch = halfCtrlObs;
			double totalNumberInChrom = treatmentChromPlus.getInfo().getNumberObservations() + treatmentChromMinus.getInfo().getNumberObservations();
			int numToFetchForChrom = (int)Math.round(totalNumberToMatch * totalNumberInChrom/numberTreatmentObservations);
			if (numToFetchForChrom < treatmentPD.getInfo().getNumberObservations())treatmentPD = PointData.fetchRandomObservations(treatmentPD, numToFetchForChrom);
		}
		//subsample controls?
		else if (halfCtrlObs > numberTreatmentObservations) {
			double totalNumberToMatch = numberTreatmentObservations *2;
			double totalNumberInChrom = controlChromPlus.getInfo().getNumberObservations() + controlChromMinus.getInfo().getNumberObservations();
			int numToFetchForChrom = (int)Math.round(totalNumberToMatch * totalNumberInChrom/numberControlObservations);
			if (numToFetchForChrom < controlPD.getInfo().getNumberObservations()) controlPD = PointData.fetchRandomObservations(controlPD, numToFetchForChrom);
		}
		//split controls
		Point[] controlPts = Point.makePoints(controlPD.getPositions(), controlPD.getScores());
		Point[][] splitPts = Point.split(controlPts);
		//set
		PointData c1 = Point.extractPositionScores(splitPts[0]);
		PointData c2 = Point.extractPositionScores(splitPts[1]);
		return new PointData[]{treatmentPD,c1,c2};
	}

	/**For cases when just interested in treatment data.*/


	private void sumTreatmentReads(){
		//get UCSCGeneLine[]
		UCSCGeneLine[] genes = geneModels.get(chromosome);

		//for each gene 
		for (int i=0; i< genes.length; i++){
			//calculate totals
			ExonIntron[] exons = genes[i].getExons();
			//fetch scores
			float tSumPlus = 0; 
			float tSumMinus = 0;
			float totalBPs = 0;
			for (int j=0; j< exons.length; j++){
				int start = exons[j].getStart();
				int stop = exons[j].getEnd();
				totalBPs += stop - start;
				if (treatmentChromPlus != null) tSumPlus += treatmentChromPlus.sumScoreBP(start, stop); 
				if (treatmentChromMinus != null)tSumMinus += treatmentChromMinus.sumScoreBP(start, stop); 
			}
			float tSum = tSumPlus+ tSumMinus;
			float rpkm = calculateRPKM(millionMappedTreatmentReads, totalBPs, tSum);
			//scores = tSum, tSumPlus,tSumMinus, total region bps, rpkm
			float[] scores = new float[]{tSum,tSumPlus,tSumMinus, totalBPs, rpkm};
			//make window
			genes[i].setScores(scores);
		}
	}
	private void sumTreatmentReadsForEnds(int numberTerminalBPs, float minimumNumberReadsEachSide){
		int minimumLength = 2 * numberTerminalBPs;
		//get UCSCGeneLine[]
		UCSCGeneLine[] genes = geneModels.get(chromosome);

		//for each gene 
		for (int i=0; i< genes.length; i++){
			//calculate totals
			ExonIntron[] exons = genes[i].getExons();
			//see if minimum length requirement met
			float totalBPs = 0;
			for (int j=0; j< exons.length; j++){
				int start = exons[j].getStart();
				int stop = exons[j].getEnd();
				int length = stop-start;
				totalBPs += length; 
			}
			if (totalBPs < minimumLength) continue;
			
			//fetch scores
			float leftSum = 0;
			float rightSum = 0;

			//make container for bed lines
			StringBuilder sb = new StringBuilder();
			
			//first for left side
			int bpCovered = 0;		
			for (int j=0; j< exons.length; j++){
				int start = exons[j].getStart();
				int stop = exons[j].getEnd();
				int length = stop-start;
				//add all?
				if ((bpCovered + length) <= numberTerminalBPs){
					if (treatmentChromPlus != null) leftSum += treatmentChromPlus.sumScoreBP(start, stop); 
					if (treatmentChromMinus != null) leftSum += treatmentChromMinus.sumScoreBP(start, stop); 
					bpCovered += length;
					sb.append(genes[i].getChrom()+"\t"+start+"\t"+stop+"\n");
				}
				//add partial
				else {
					int diff = (bpCovered + length) - numberTerminalBPs;
					stop -= diff;
					if (treatmentChromPlus != null) leftSum += treatmentChromPlus.sumScoreBP(start, stop); 
					if (treatmentChromMinus != null) leftSum += treatmentChromMinus.sumScoreBP(start, stop);
					length = stop-start;
					bpCovered+=length;
					sb.append(genes[i].getChrom()+"\t"+start+"\t"+stop+"\n");
					break;
				}
			}	
			
			//continue?
			if (leftSum < minimumNumberReadsEachSide) continue;
			
			//then for right side
			bpCovered = 0;		
			for (int j=exons.length-1; j>= 0; j--){
				int start = exons[j].getStart();
				int stop = exons[j].getEnd();
				int length = stop-start;
				//add all?
				if ((bpCovered + length) <= numberTerminalBPs){
					if (treatmentChromPlus != null) rightSum += treatmentChromPlus.sumScoreBP(start, stop); 
					if (treatmentChromMinus != null) rightSum += treatmentChromMinus.sumScoreBP(start, stop); 
					sb.append(genes[i].getChrom()+"\t"+start+"\t"+stop+"\n");
					bpCovered += length;
				}
				//add partial
				else {
					int bpNeeded = numberTerminalBPs-bpCovered;
					start = stop - bpNeeded;
					if (treatmentChromPlus != null) rightSum += treatmentChromPlus.sumScoreBP(start, stop); 
					if (treatmentChromMinus != null) rightSum += treatmentChromMinus.sumScoreBP(start, stop);
					sb.append(genes[i].getChrom()+"\t"+start+"\t"+stop+"\n");
					length = stop-start;
					bpCovered+=length;
					break;
				}
			}		
			
			//print? 5' and 3' scores
			if (rightSum < minimumNumberReadsEachSide) continue;
			
			String strand = genes[i].getStrand();
			if (strand.equals("+")){
				System.out.println(genes[i].getDisplayName()+"\t"+genes[i].getName()+"\t"+genes[i].coordinates()+"\t"+leftSum+"\t"+rightSum+"\t"+(Num.log2((1+rightSum)/(1+leftSum))));
			}
			else if (strand.equals("-")) System.out.println(genes[i].getDisplayName()+"\t"+genes[i].getName()+"\t"+genes[i].coordinates()+"\t"+rightSum+"\t"+leftSum+"\t"+(Num.log2((1+leftSum)/(1+rightSum))));
			else System.out.println(genes[i].getName()+"\tUnstranded!");
			//System.out.println(sb);
		}
	}

	/**Main window scanner using binomial p-value as the score based on read counts.
	 * Be sure to replace the PointData scores with 1 if you want to use the real bin p-val.*/
	private void calculateBinomialPValues(){
		//make arrays to calculate binomial p-vals with R
		ArrayList<UCSCGeneLine> forBinom = new ArrayList<UCSCGeneLine>();

		//get UCSCGeneLine[]
		UCSCGeneLine[] genes = geneModels.get(chromosome);

		//for each gene 
		for (int i=0; i< genes.length; i++){
			//calculate totals
			ExonIntron[] exons = genes[i].getExons();
			//fetch scores
			float tSumPlus = 0; 
			float tSumMinus = 0;
			float cSumPlus = 0;
			float cSumMinus = 0;
			float[] tExonCounts = new float[exons.length];
			float[] cExonCounts = new float[exons.length];
			float totalExonicBP = 0;
			for (int j=0; j< exons.length; j++){
				int start = exons[j].getStart();
				int stop = exons[j].getEnd();
				totalExonicBP+= (stop - start);
				//treatment
				if (treatmentChromPlus != null) {
					float t = treatmentChromPlus.sumScoreBP(start, stop);
					tExonCounts[j] = (long)t;
					tSumPlus += t; 
				}
				if (treatmentChromMinus != null){
					float t = treatmentChromMinus.sumScoreBP(start, stop);
					tExonCounts[j] += (long)t;
					tSumMinus += t;
				}
				//control
				if (controlChromPlus != null){
					float c = controlChromPlus.sumScoreBP(start, stop);
					cExonCounts[j] = (long) c;
					cSumPlus += c; 
				}
				if (controlChromMinus != null){
					float c = controlChromMinus.sumScoreBP(start, stop);
					cExonCounts[j] += (long)c;
					cSumMinus += c; 
				}
			}

			//scores =                 pVal,upPVal,downPVal,skew,chiSqrPVal,bps,tSumPlus,tSumMinus, cSumPlus,cSumMinus,realND,mockND,upEmpFDR
			float[] scores = new float[]{0,   0,     0,       0, 0, totalExonicBP, tSumPlus,tSumMinus,cSumPlus, cSumMinus,  0,     0,      0};
			float tSum = tSumPlus+ tSumMinus;
			float cSum = cSumPlus+ cSumMinus;

			//save scores
			genes[i].setScores(scores);
			if ((tSum+cSum) == 0) continue;
			geneLinesWithReadsAL.add(genes[i]);

			//save exon counts?
			if ((tSum+cSum) >= 50) genes[i].setExonCounts(new float[][]{tExonCounts, cExonCounts});

			//calc diff binomial from cache? or directly?
			//yes
			if (tSum < numberCachedBinPVals && cSum < numberCachedBinPVals){
				//calc pvals
				int t = Math.round(tSum);
				int c = Math.round(cSum);
				float up = pValsDataUp[t][c];
				float down = pValsDataDown[c][t];
				//only test for skew in one direction, that minus is way less than plus
				float skew = pValsSkew[Math.round(tSumPlus)][Math.round(tSumMinus)];
				//set scores
				//scores = pVal,upPVal,downPVal,skew,chiSqrPVal,bps,tSumPlus,tSumMinus,cSum
				if (up > down) scores[0] = up;
				else scores[0] = -1* down;
				scores[1] = up;
				scores[2] = down;
				scores[3] = skew;
				genes[i].setScores(scores);
			}
			//no
			else forBinom.add(genes[i]);

			//score splices?
			if (spliceTreatmentCounts != null && (tSum + cSum) >= minNumObs2ScoreSpliceJunctions && exons.length > 1){
				//generate splice names
				String[] spliceNames = makeSpliceJunctions(exons);
				//for each junction, collect scores
				ArrayList<SpliceJunction> spAL = new ArrayList<SpliceJunction>();
				for (int x=0; x< spliceNames.length; x++){
					int numT =0;
					int numC =0;
					if (spliceTreatmentCounts.containsKey(spliceNames[x])) numT = spliceTreatmentCounts.get(spliceNames[x]);
					if (spliceControlCounts.containsKey(spliceNames[x])) numC = spliceControlCounts.get(spliceNames[x]);
					if (numT !=0 || numC !=0){
						spAL.add(new SpliceJunction(spliceNames[x], numT, numC));
					}
				}
				int numSpliceJunctions = spAL.size();
				if (numSpliceJunctions !=0){
					totalNumberJunctionsScored += numSpliceJunctions;
					SpliceJunction[] sjs = new SpliceJunction[numSpliceJunctions];
					spAL.toArray(sjs);
					genes[i].setSpliceJunctions(sjs);
				} 
			}
		}
		//calc binom from R
		if (forBinom.size()!=0) calculateBinomialsFromR(forBinom);
	}

	public static String[] makeSpliceJunctions (ExonIntron[] exons){
		int num = exons.length;
		String[] junctionNames = new String[num*(num-1)/2];
		//for each pairing
		int index = 0;
		for (int j=0; j< num; j++){
			//define the stop
			int end = exons[j].getEnd();
			//for all subsequent starts
			for (int k=j+1; k< num; k++){
				int start = exons[k].getStart();
				junctionNames[index++] = end+"_"+start;
			}
		}
		return junctionNames;
	}

	/**Takes UCSCGeneLines that had observations that exceeded cached pvalues and calculates
	 * directly from R.*/
	private void calculateBinomialsFromR(ArrayList<UCSCGeneLine> sm){
		int num = sm.size();
		//for each win, want to calculate up, down, skew
		int[][] skewObs = new int[num][2];
		int[][] upObs = new int[num][2];
		int[][] downObs = new int[num][2];
		for (int i=0; i< num; i++){
			//scores = pVal,upPVal,downPVal,skew,chiSqrPVal,bps,tSumPlus,tSumMinus,cSumPlus,cSumMinus, realND,mockND,upEmpFDR
			float[] scores = sm.get(i).getScores();
			int tSumPlus = Math.round(scores[6]);
			int tSumMinus = Math.round(scores[7]);
			int tSum = tSumPlus + tSumMinus;
			int cSum = Math.round(scores[8]+scores[9]);
			//set obs
			skewObs[i] = new int[]{tSumPlus, tSumMinus};
			upObs[i] = new int[]{tSum,cSum};
			downObs[i] = new int[]{cSum,tSum};
		}
		//fetch pvalues
		double[] skewPvals = Num.binomialPValues(0.5, saveDirectory, skewObs, fullPathToR, true, false, true);
		double[] upPvals = Num.binomialPValues(expectedFractionUp, saveDirectory, upObs, fullPathToR, true, false, true);
		double[] downPvals = Num.binomialPValues(expectedFractionDown, saveDirectory, downObs, fullPathToR, true, false,true);
		if (skewPvals == null || upPvals == null || downPvals == null || skewPvals.length != num || upPvals.length != num || downPvals.length != num) {
			Misc.printErrAndExit("\nProblem fetching binomial pvalues from R. Check that R really is at "+fullPathToR);
		}
		//assign scores
		for (int i=0; i< num; i++){
			UCSCGeneLine s = sm.get(i);
			//scores = pVal,upPVal,downPVal,skew,chiSqrPVal,bps,tSumPlus,tSumMinus,cSum,realND,mockND,upEmpFDR
			float[] scores = s.getScores();
			scores[1] = (float)upPvals[i];
			scores[2] = (float)downPvals[i];
			scores[3] = (float)skewPvals[i];
			//assign pVal score
			if (scores[1]> scores[2]) scores[0] = scores[1];
			else scores[0] = -1 * scores[2];
		}
	}

	/**Collects and calculates a bunch of stats re the PointData.*/
	private void calculateReadCountStatistics(){
		//fetch treatment PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (treatmentPointDirs);
		treatmentPlusPointData = PointData.convertArrayList2Array(combo[0]);
		treatmentMinusPointData = PointData.convertArrayList2Array(combo[1]);
		numberTreatmentObservations = PointData.totalObservationsMultiPointData(treatmentPlusPointData);
		numberTreatmentObservations += PointData.totalObservationsMultiPointData(treatmentMinusPointData);
		millionMappedTreatmentReads = numberTreatmentObservations/1000000;
		System.out.println("\t"+(int)numberTreatmentObservations+" Treatment Observations");

		//likewise control data
		if (controlPointDirs != null){
			combo = PointData.fetchStrandedPointDataNoMerge (controlPointDirs);
			controlPlusPointData = PointData.convertArrayList2Array(combo[0]);
			controlMinusPointData = PointData.convertArrayList2Array(combo[1]);
			numberControlObservations = PointData.totalObservationsMultiPointData(controlPlusPointData);
			numberControlObservations += PointData.totalObservationsMultiPointData(controlMinusPointData);
			System.out.println("\t"+(int)numberControlObservations+" Control Observations");
			//some stats for binomial pvalue calc and scaling the ratios
			expectedFractionUp = numberTreatmentObservations/(numberTreatmentObservations+numberControlObservations);
			expectedFractionDown = 1- expectedFractionUp;		
			scalarTC = numberTreatmentObservations/ numberControlObservations;
			scalarCT = numberControlObservations/numberTreatmentObservations;
			millionMappedControlReads = numberControlObservations/1000000;

			//make look up tables		
			pValsDataUp = Num.convertToFloat(Num.binomialPValMatrix(numberCachedBinPVals+1, expectedFractionUp, saveDirectory, fullPathToR, true));
			pValsDataDown = Num.convertToFloat(Num.binomialPValMatrix(numberCachedBinPVals+1, expectedFractionDown, saveDirectory, fullPathToR, true));
			pValsSkew = Num.convertToFloat(Num.binomialPValMatrix(numberCachedBinPVals+1, 0.5, saveDirectory, fullPathToR, true));		
		}
	}

	public HashMap<String, Integer> loadAndMergeSplices (HashMap<String, File> data){
		//any data
		File f = data.get(chromosome);
		if (f == null) return null;
		HashMap<String, Integer> spliceNameCounts = new HashMap<String, Integer>();
		try {
			String line;
			BufferedReader in = new BufferedReader ( new FileReader (f));
			while ((line = in.readLine()) != null){
				int index = line.indexOf(":");
				String name = line.substring(index+1);
				int count = Integer.parseInt(line.substring(0,index));
				if (spliceNameCounts.containsKey(name)) count += spliceNameCounts.get(name).intValue();
				spliceNameCounts.put(name, new Integer (count));
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem loading parsed splice data from "+f);
		}
		return spliceNameCounts;
	}

	/**Assumes bed files are sorted by chromosome.*/
	public HashMap<String, File> splitByChrom (File[] bedFiles){
		Pattern tab = Pattern.compile("\\t");
		HashMap<String, File> chromFiles = new HashMap<String, File>();
		try{
			HashMap<String, PrintWriter> chromPWs = new HashMap<String, PrintWriter>();
			PrintWriter out = null;
			String currChrom = "";			
			for (int i=0; i< bedFiles.length; i++){
				BufferedReader in = IO.fetchBufferedReader(bedFiles[i]);
				System.out.println("\t"+bedFiles[i].getName());
				String[] tokens;
				String line;
				while ((line=in.readLine()) != null){
					tokens = tab.split(line);
					if (tokens.length != 6) continue;
					//check chrom
					if (tokens[0].equals(currChrom) == false){
						currChrom = tokens[0];
						//look for pw
						if (chromPWs.containsKey(currChrom)) out = chromPWs.get(currChrom);
						else {
							File f = new File (saveDirectory, "tempFile_"+currChrom+"_"+Passwords.createRandowWord(7));
							out = new PrintWriter ( new FileWriter( f ));
							chromFiles.put(currChrom, f);
							chromPWs.put(currChrom, out);
						}
					}
					out.println(tokens[3]);
				}
			}
			//close PWs
			Iterator<String> it = chromPWs.keySet().iterator();
			while (it.hasNext()) chromPWs.get(it.next()).close();

		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing bed files.\n");
		}
		return chromFiles;
	}

	/**Shifts the positions halfPeakShift (+ for sense, - for antisense) sets the positions into the data.
	 * May replace all scores with 1 if stripScores == true.*/
	private void shiftStripPointData(){
		//fetch data from treatments
		if (treatmentChromPlus !=null) {
			int[] p = treatmentChromPlus.getPositions();
			addShift(p,halfPeakShift);
			if (useWeightedReads == false) treatmentChromPlus.stripScores();
		}
		if (treatmentChromMinus!=null) {
			int[] p = treatmentChromMinus.getPositions();
			addShift(p, -1*halfPeakShift);
			if (useWeightedReads == false) treatmentChromMinus.stripScores();
		}
		//fetch data from controls
		if (controlPointDirs != null) {
			if (controlChromPlus !=null) {
				int[] p = controlChromPlus.getPositions();
				addShift(p, halfPeakShift);
				if (useWeightedReads == false) controlChromPlus.stripScores();
			}
			if (controlChromMinus!=null){
				int[] p = controlChromMinus.getPositions();
				addShift(p, -1*halfPeakShift);
				if (useWeightedReads == false) controlChromMinus.stripScores();
			}
		}
	}

	/**Adds the toAdd to each int.*/
	public static void addShift(int[] positions, int toAdd){
		for (int i=0; i< positions.length; i++){
			positions[i] += toAdd;
			if (positions[i]<0) positions[i] = 0;
		}
	}

	/**Given a score threshold, returns the number of genes that would be generated using the
	 * precomputed values.  This is a way of minimizing the memory requirements.*/
	private int calculateNumberControlERs(float threshold){
		int index = Arrays.binarySearch(controlThresholds, threshold);
		if (index < 0){
			int conv = -index -1;
			if (conv >= controlThresholds.length) return 0;
			else return controlNumberGenes[conv];
		}
		else return controlNumberGenes[index];
	}


	/**Calculates empirical FDRs.*/
	private void estimateEmpiricalFDRs(UCSCGeneLine[] regionsWithReads, int scoreIndex, int fdrIndex){
		//sort by score
		Arrays.sort(regionsWithReads, new UCSCGeneLineComparatorScore(scoreIndex));

		//fetch scores
		float[] sortedScores = new float[regionsWithReads.length];
		for (int i=0; i< regionsWithReads.length; i++){
			float[] scores = regionsWithReads[i].getScores();
			sortedScores[i] = scores[scoreIndex];	
		}
		//for each score calc fdr = # control ERs/ #real ERs
		float testScore = 0;
		double numControls = 0;
		float maxFDR = 0;

		for (int i=0; i< regionsWithReads.length; i++){
			//calc num control windows
			testScore = regionsWithReads[i].getScores()[scoreIndex];
			numControls = calculateNumberControlERs(testScore);
			//any control windows?
			//yes, controls, calc FDR and assign to maxFDR
			if (numControls !=0) {				
				//calc num real
				double numReal = regionsWithReads.length - Num.findClosestIndexToValue(sortedScores, testScore);
				//calc fdr
				float fdr = 0;
				if (numReal !=0) fdr = new Double(-10 * Num.log10(numControls/numReal)).floatValue();
				if (fdr > maxFDR) maxFDR = fdr;
				//advance until score changes and assign max
				for (; i< regionsWithReads.length; i++){
					//different so break
					if (regionsWithReads[i].getScores()[scoreIndex] != testScore) break;
					//not different so set maxFDR
					else {						
						float[] s = regionsWithReads[i].getScores();
						s[fdrIndex] = maxFDR;
						//System.out.println("\t"+i+" Setting MaxFDR-> "+maxFDR+" score-> "+s[scoreIndex]);						
					}
				}
				//back off one so calc numcontrol on diff score
				if (i < regionsWithReads.length) i--;
			}
			//no  more controls so assign maxFDR to remainder
			else {
				for (; i< regionsWithReads.length; i++){
					float[] s = regionsWithReads[i].getScores();
					s[fdrIndex] = maxFDR;	
				}
			}
		}		
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new DefinedEndScanSeqs(args);
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
					case 't': treatmentPointDirs = IO.extractFiles(args[++i]); break;
					case 'c': controlPointDirs = IO.extractFiles(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'p': peakShift = Integer.parseInt(args[++i]); break;
					case 'f': findReducedRegions = false; break;
					case 'w': useWeightedReads = true; break;
					case 'u': refSeqFile = new File(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 'd': treatmentSpliceFiles = IO.extractFiles(args[++i]); break;
					case 'e': controlSpliceFiles = IO.extractFiles(args[++i]); break;
					case 'm': minNumObs2ScoreSpliceJunctions = Integer.parseInt(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//look for point directories
		if (treatmentPointDirs == null || treatmentPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your treatment PointData directories(s)!\n");
		//only one directory look deeper
		if (treatmentPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(treatmentPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) treatmentPointDirs = otherDirs;
		}
		//control data
		if (controlPointDirs != null){
			if (controlPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your control PointData directories(s)!\n");
			//only one directory look deeper
			if (controlPointDirs.length == 1){
				File[] otherDirs = IO.extractOnlyDirectories(controlPointDirs[0]);
				if (otherDirs != null && otherDirs.length > 0) controlPointDirs = otherDirs;
			}
		}
		//set half peak shift and windowSize
		halfPeakShift = (int)Math.round( ((double)peakShift)/2 );

		//need R
		if (fullPathToR == null || fullPathToR.exists() == false) Misc.printExit("\nPlease enter the full path to R, ie '-r /usr/bin/R'\n");

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		if (saveDirectory.exists() == false) saveDirectory.mkdir();

		//look for bed file
		if (refSeqFile == null && bedFile == null){
			Misc.printExit("\nPlease enter a regions file to use in scoring regions.\n");
		}
		//info
		if (useWeightedReads) System.out.println("Using read score probabilities...");
		if (findReducedRegions) System.out.println("Scanning for enriched and reduced regions, no EmpFDR estimation...");
		else System.out.println("Scanning for enriched regions...");
		System.out.println("\t"+peakShift+ "\tPeak shift");

		//check splice junctions
		if (treatmentSpliceFiles != null || controlSpliceFiles != null){
			if (treatmentSpliceFiles == null || controlSpliceFiles == null) Misc.printErrAndExit("\nSomething is wrong with one or both of your treatment or control splice junction files. Aborting\n");
		}
	}	

	public static float multBy10AndRound(float d){
		float comp = d*10.0f;
		comp = Math.round(comp);
		return comp;
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                           Defined Region Scan Seqs: Feb 2010                     **\n" +
				"**************************************************************************************\n" +
				"DRSS takes chromosome specific PointData xxx.bar.zip files and extracts scores under\n" +
				"each region to calculate several statistics including a binomial p-value, Storey\n" +
				"q-value FDR, an empirical FDR, a p-value for strand skew, and a chi-square test of\n" +
				"independence between the exon read count distributions between treatment and control\n" +
				"data (a test for alternative splicing). Several measures of read counts are provided\n" +
				"including counts for each strand, a normalized log2 ratio, and RPKMs (# reads per kb\n" +
				"of interrogated region per total million mapped reads). If a gene table is provided,\n" +
				"scores under each exon are summed to give a whole gene summary. It is also recommended\n" +
				"to run a gene table of introns (see the ExportIntronicRegions app) to look for\n" +
				"intronic retention and novel transfrags/ exons.  If one provides splice junction bed\n" +
				"files for treatment and control RNA-Seq data, see the NovoalignParser, splice\n" +
				"junctions will be scored for differential expression. This is an additional\n" +
				"calculation unrelated to the chi-square independance test. Lastly, if control\n" +
				"data is not provided, simple region sums are calculated.\n\n"+

				"Options:\n"+
				"-s Save directory, full path.\n"+
				"-t Treatment PointData directories, full path, comma delimited. These should\n" +
				"       contain unshifted stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories.\n" +
				"-c Control PointData directories, ditto. \n" +
				"-p Peak shift, average distance between + and - strand peaks for chIP-Seq data, see\n" +
				"       PeakShiftFinder. For RNA-Seq set to the smallest expected fragment size. Will\n" +
				"       be used to shift the PointData 3' by 1/2 the peak shift.\n"+
				"-r Full path to R loaded with Storey's q-value library, defaults to '/usr/bin/R'\n" +
				"       file, see http://genomics.princeton.edu/storeylab/qvalue/\n"+
				"-u UCSC RefFlat or RefSeq gene table file, full path. See,\n"+
				"       http://genome.ucsc.edu/cgi-bin/hgTables, (name1 name2(optional) chrom strand\n" +
				"       txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds)\n"+
				"-b (Or) a bed file (chr, start, stop,...), full path, See,\n" +
				"       http://genome.ucsc.edu/FAQ/FAQformat#format1\n"+

				"\nAdvanced Options:\n"+
				"-f Scan for just enriched regions, defaults to look for both. Only use with chIP-Seq\n" +
				"       datasets where the control is input. This turns on the empFDR estimation.\n"+
				"-d Treatment splice junction bed file(s) from the NovoalignParser, comma delimited,\n" +
				"       full path.\n"+
				"-e Control splice junction bed file(s), comma delimited, full path.\n"+
				"-m Minimum number of reads in associated gene before scoring splice junctions.\n"+
				"       Used in estimating the expected proportion of T and scaling the log2Ratio. \n" +
				"       Defaults to 100.\n"+
				"-w Use read score probabilities (assumes scores are > 0 and <= 1), defaults to\n" +
				"       assigning 1 to each read score. Experimental.\n"+
				"\n"+

				"Example: java -Xmx4G -jar pathTo/USeq/Apps/DefinedEndScanSeqs -t\n" +
				"      /Data/PolIIRep1/,/Data/PolIIRep2/ -c /Data/Input1/,Data/Input2/ -s\n" +
				"      /Data/PolIIResults -p 100 -b /Data/selectRegions.bed -f \n\n" +

		"**************************************************************************************\n");

	}
}
