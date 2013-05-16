package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import edu.utah.seq.useq.data.Region;
import trans.tpmap.*;
import org.apache.commons.math3.stat.inference.ChiSquareTest;


/**For identifing differentially methylated genomic regions from two condition bisulfite sequencing data
 * @author Nix
 * */
public class DefinedRegionBisSeq {

	//user defined fields
	private File[] tConPointDirs;
	private File[] tNonConPointDirs;
	private File[] cConPointDirs;
	private File[] cNonConPointDirs;
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private int minimumReadCoverage = 5;
	private String regionsFileName;

	private LinkedHashMap<String,SmoothingWindow[]> chromosomeRegion;

	//internal fields
	private HashMap<String,PointData[]> tConPlusPointData;
	private HashMap<String,PointData[]> tConMinusPointData;
	private HashMap<String,PointData[]> tNonConPlusPointData;
	private HashMap<String,PointData[]> tNonConMinusPointData;
	private HashMap<String,PointData[]> cConPlusPointData;
	private HashMap<String,PointData[]> cConMinusPointData;
	private HashMap<String,PointData[]> cNonConPlusPointData;
	private HashMap<String,PointData[]> cNonConMinusPointData;
	private int log2RatioIndex = 0;
	private int pValIndex = 1;
	private int numberRegions = 0;

	private double cacheNumber = 1000;
	private FisherExact fischerExact = new FisherExact((int)cacheNumber);
	private ChiSquareTest chiSquare = new ChiSquareTest();
	private String genomeVersion;

	//by chromosome
	private String chromosome;
	private PointData tConMergedChromPlus = null;
	private PointData tConMergedChromMinus = null;	
	private PointData tNonConMergedChromPlus = null;
	private PointData tNonConMergedChromMinus = null;
	private PointData cConMergedChromPlus = null;
	private PointData cConMergedChromMinus = null;	
	private PointData cNonConMergedChromPlus = null;
	private PointData cNonConMergedChromMinus = null;


	//constructors
	/**Stand alone.*/
	public DefinedRegionBisSeq(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//load hashes
		loadDataHashes();

		//for each chromosome of regions
		System.out.println("\nScanning regions by chromosome for differential methylation... ");

		for (String chr : chromosomeRegion.keySet()){
			chromosome = chr;
			System.out.print(chromosome + " ");

			//fetch data
			fetchData();

			//check
			if (tConMergedChromPlus == null || tConMergedChromMinus == null || tNonConMergedChromPlus == null || tNonConMergedChromMinus == null || 
					cConMergedChromPlus == null || cConMergedChromMinus == null || cNonConMergedChromPlus == null || cNonConMergedChromMinus == null) {
				System.out.print(" - Couldn't find all 8 datasets, skipping! ");
				continue;
			}
			//scan chromosome
			scanChromosome();
		}

		System.out.println();System.out.println();

		//B & H correct pvalues and make EnrichedRegions
		System.out.println("Applying a Benjamini and Hochberg FDR correction to the p-values...\n");
		correctPValuesBH();

		//print spreadsheet
		printSpreadSheet();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("Done! "+Math.round(diffTime)+" min\n");
	}



	/**Writes out an excel compatible tab delimited spreadsheet with hyperlinks for IGB.*/
	public void printSpreadSheet(){
		try{
			File file = new File(saveDirectory, regionsFileName+"_DRBSMinRC"+minimumReadCoverage+".xls");
			PrintWriter out = new PrintWriter (new FileWriter (file));
			//print header line
			out.println("#"+genomeVersion+"_IGBHyperLinks\tChr\tStart\tStop\tPseMedianLog2Ratio\tFDR\t#Obs\t#TCon\t#TNonCon\t#CCon\t#CNonCon\tPseMedianLog2TFractions\tPseMedianLog2CFractions");
			String url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";
			String tab = "\t";
			//for each region
			int index = 0;
			for (String chr: chromosomeRegion.keySet()){
				for (SmoothingWindow sw : chromosomeRegion.get(chr)){
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
					//scores
					float[] scores = sw.getScores();
					//	log2Ratio
					out.print(scores[log2RatioIndex]); out.print(tab);
					//	FDR
					out.print(scores[pValIndex]);
					//	counts
					for (int j=2; j< 7; j++){
						out.print(tab);
						out.print((int)scores[j]);
					}
					out.print(tab);
					out.print(scores[7]);
					out.print(tab);
					out.print(scores[8]);
					out.println();
				}
			}
			out.close();
		}catch (Exception e){
			System.out.println("\nError: problem printing spreadsheet report");
			e.printStackTrace();
		}
	}











	/**Fetchs the data for a particular chromosome.*/
	public void fetchData(){
		ArrayList<PointData> al = null;
		PointData[] pd;
		//merge converted
		tConMergedChromPlus = null;
		if (tConPlusPointData.containsKey(chromosome)) {
			pd = tConPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			tConMergedChromPlus = PointData.mergePointData(al, false, true);
			//set version
			if (genomeVersion == null) genomeVersion = pd[0].getInfo().getVersionedGenome();
		}
		tConMergedChromMinus = null;
		if (tConMinusPointData.containsKey(chromosome)) {
			pd = tConMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			tConMergedChromMinus = PointData.mergePointData(al, false, true);
		}
		//merge nonConverted
		tNonConMergedChromPlus = null;
		if (tNonConPlusPointData.containsKey(chromosome)) {
			pd = tNonConPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			tNonConMergedChromPlus = PointData.mergePointData(al, false, true);
		}
		tNonConMergedChromMinus = null;
		if (tNonConMinusPointData.containsKey(chromosome)) {
			pd = tNonConMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			tNonConMergedChromMinus = PointData.mergePointData(al, false, true);
		}

		//merge converted
		cConMergedChromPlus = null;
		if (cConPlusPointData.containsKey(chromosome)) {
			pd = cConPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			cConMergedChromPlus = PointData.mergePointData(al, false, true);
		}
		cConMergedChromMinus = null;
		if (cConMinusPointData.containsKey(chromosome)) {
			pd = cConMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			cConMergedChromMinus = PointData.mergePointData(al, false, true);
		}
		//merge nonConverted
		cNonConMergedChromPlus = null;
		if (cNonConPlusPointData.containsKey(chromosome)) {
			pd = cNonConPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			cNonConMergedChromPlus = PointData.mergePointData(al, false, true);
		}
		cNonConMergedChromMinus = null;
		if (cNonConMinusPointData.containsKey(chromosome)) {
			pd = cNonConMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			cNonConMergedChromMinus = PointData.mergePointData(al, false, true);
		}

		pd = null;
		al = null;
	}


	public void scanChromosome(){
		//merge strands
		ArrayList<PointData> al = new ArrayList<PointData>();
		al.add(tConMergedChromMinus);
		al.add(tConMergedChromPlus);
		PointData mergedTreatmentCon = PointData.mergePointData(al, false, true);
		al.clear();
		al.add(tNonConMergedChromMinus);
		al.add(tNonConMergedChromPlus);
		PointData mergedTreatmentNonCon = PointData.mergePointData(al, false, true);

		al.clear();
		al.add(cConMergedChromMinus);
		al.add(cConMergedChromPlus);
		PointData mergedControlCon = PointData.mergePointData(al, false, true);
		al.clear();
		al.add(cNonConMergedChromMinus);
		al.add(cNonConMergedChromPlus);
		PointData mergedControlNonCon = PointData.mergePointData(al, false, true);

		//fetch paired base observations meeting minimum read coverage in both T and C
		MethylatedBaseObservation[] mbo = fetchCommonBasesWithMinimumObservations (mergedTreatmentNonCon, mergedTreatmentCon, mergedControlNonCon, mergedControlCon, minimumReadCoverage);
		
		if (mbo == null) return;

		//scan
		scoreRegionsForDifferentialMethylation(mbo);

	}


	public void scoreRegionsForDifferentialMethylation(MethylatedBaseObservation[] mbo){

		SmoothingWindow[] regions = chromosomeRegion.get(chromosome);
		ArrayList<SmoothingWindow> forR = new ArrayList<SmoothingWindow>();
		numberRegions+= regions.length;

		int[] positions = MethylatedBaseObservation.fetchPositions(mbo);

		//for each region
		for (int i=0; i< regions.length; i++){

			int[] indexes = Num.findIndexes(regions[i].getStart(), regions[i].getStop(), positions);

			//for each index, calculate fraction methylated
			float[] fractions = new float[indexes[1]- indexes[0]];
			float pseLog2 = 0.0f;
			float pseTFrac =  0.0f;
			float pseCFrac =  0.0f;
			
			if (fractions.length !=0){
				float[] tFractions = new float [fractions.length];
				float[] cFractions = new float [fractions.length];
				int f = 0;
				
				for (int x=indexes[0]; x< indexes[1]; x++) {
					tFractions[f] = Num.log2(mbo[x].getFractionMethylatedT());
					cFractions[f] = Num.log2(mbo[x].getFractionMethylatedC());
					fractions[f] = Num.log2(mbo[x].getFractionMethylatedT()/mbo[x].getFractionMethylatedC());
					f++;
				}
				pseLog2 = (float) Num.pseudoMedian(fractions);
				pseTFrac = (float) Num.pseudoMedian(tFractions);
				pseCFrac = (float) Num.pseudoMedian(cFractions);
				
			}

			//build array for chiSquareTest
			long tNumCon = 0;
			long tNumNonCon = 0;
			long cNumCon = 0;
			long cNumNonCon = 0;
			for (int x=indexes[0]; x< indexes[1]; x++) {
				//nonConT, conT, nonConC, conC
				long[] c = mbo[x].getCounts();
				tNumCon+= c[1];
				tNumNonCon+= c[0];
				cNumCon+= c[3];
				cNumNonCon+= c[2];
			}

			//scores are log2, pval, num bp scored, tNumCon, tNumNonCon, cNumCon, cNumNonCon, pseTFrac, pseCFrac
			float[] scores = new float[]{pseLog2, 0, fractions.length, (float)tNumCon, (float)tNumNonCon, (float)cNumCon, (float)cNumNonCon, pseTFrac, pseCFrac};
			regions[i].setScores(scores);


			//watch out for no observations in one or the other datasets
			if ((tNumCon == 0 && tNumNonCon == 0) ||  (cNumCon == 0 && cNumNonCon == 0)) continue;
			
			long[][] counts = new long[][]{{tNumCon, tNumNonCon},{cNumCon,cNumNonCon}};

			//calculate p-value for differences in unstranded methylation
			boolean sendToR = false;
			double totalObservations = tNumCon + tNumNonCon + cNumCon + cNumNonCon;
			
			//use fischer's?
			if (totalObservations < cacheNumber) {
				double pNoLog = fischerExact.getTwoTailedP((int)tNumCon, (int)tNumNonCon, (int)cNumCon, (int)cNumNonCon);
				scores[pValIndex] = Num.minus10log10Float(pNoLog);
			}

			//use chi-square?
			else {	
				//too many for java chi-square?
				double chiSquareStat = Num.chiSquareTestStatistic(counts);
				if (chiSquareStat > 1400) sendToR = true;

				//use apache chi-square
				else{
					double pNoLog = -1;
					try {
						pNoLog = chiSquare.chiSquareTest(counts);
					} catch (Exception e) {}
					if (pValCheck(pNoLog) == false) sendToR = true;
					else scores[pValIndex] = Num.minus10log10Float(pNoLog);
				}
			}
			//assign
			if (sendToR) forR.add(regions[i]); 
		}

		//any difficult counts?
		if (forR.size() !=0) calculateChiSquareInRWithLookupBPs(forR);

		//check pvalues
		for (SmoothingWindow sm : regions){
			float[] s = sm.getScores();
			//need to watch out for slight neg -10log10(pval)s, these are insignificant
			if (pValCheck(s[pValIndex]) == false) s[pValIndex] = 0;	
		}
	}

	public void correctPValuesBH(){
		//load Point[] with all of the pvalues
		Point[] point = new Point[numberRegions];
		int index = 0;
		for (SmoothingWindow[] sw : chromosomeRegion.values()){
			for (int j=0; j<sw.length; j++){
				point[index] = new Point(index, sw[j].getScores()[pValIndex]);
				index++;
			}
		}
		//sort by score smallest -10Log10(pval) to largest
		Arrays.sort(point, new ComparatorPointAscendingScore());
		//correct
		Point.benjaminiHochbergCorrect(point, 0);
		//sort back to original position
		Arrays.sort(point, new ComparatorPointPosition());
		//assign FDRs to pVals
		index = 0;
		for (SmoothingWindow[] sw : chromosomeRegion.values()){
			for (int j=0; j<sw.length; j++){
				sw[j].getScores()[pValIndex] = point[index].getScore();
				index++;
			}
		}
	}

	/**Collects bases that have a minimum interrogation in both the T and C.*/
	public static MethylatedBaseObservation[] fetchCommonBasesWithMinimumObservations(PointData nonConT, PointData conT, PointData nonConC, PointData conC, int minimumReadCoverage){

		//fetch arrays
		int[] positionsNonConT = nonConT.getPositions();
		float[] readsNonConT = nonConT.getScores();
		int[] positionsConT = conT.getPositions();
		float[] readsConT = conT.getScores();

		int[] positionsNonConC = nonConC.getPositions();
		float[] readsNonConC = nonConC.getScores();
		int[] positionsConC = conC.getPositions();
		float[] readsConC = conC.getScores();

		//make containers for graph data
		ArrayList<MethylatedBaseObservation> ibAL = new ArrayList<MethylatedBaseObservation>();

		//collect all positions
		int[] allPositions = Num.returnUniques(new int[][]{positionsNonConT, positionsConT, positionsNonConC, positionsConC});

		//for each position 
		int indexNonConT =0;
		int indexConT =0;
		int indexNonConC =0;
		int indexConC =0;
		for (int i=0; i< allPositions.length; i++){
			int testPos = allPositions[i];
			//values for each
			float numNonConT =0;
			float numConT =0;
			float numNonConC =0;
			float numConC =0;

			//treatment
			//present in nonConT?
			for (int j=indexNonConT; j < positionsNonConT.length; j++){
				//match!
				if (testPos == positionsNonConT[j]){
					numNonConT = readsNonConT[j];
					indexNonConT++;
					break;
				}
				//less than
				if (testPos < positionsNonConT[j]) break;
				//greater than so keep advancing
				indexNonConT = j;
			}
			//present in conT?
			for (int j=indexConT; j < positionsConT.length; j++){
				//match!
				if (testPos == positionsConT[j]){
					numConT = readsConT[j];
					indexConT++;
					break;
				}
				//less than
				if (testPos < positionsConT[j]) break;
				//greater than so keep advancing
				indexConT = j;
			}
			//enough t obs?
			int numTObs = (int)(numConT+numNonConT);
			if (numTObs < minimumReadCoverage) continue;

			//control
			//present in nonConC?
			for (int j=indexNonConC; j < positionsNonConC.length; j++){
				//match!
				if (testPos == positionsNonConC[j]){
					numNonConC = readsNonConC[j];
					indexNonConC++;
					break;
				}
				//less than
				if (testPos < positionsNonConC[j]) break;
				//greater than so keep advancing
				indexNonConC = j;
			}
			//present in conC?
			for (int j=indexConC; j < positionsConC.length; j++){
				//match!
				if (testPos == positionsConC[j]){
					numConC = readsConC[j];
					indexConC++;
					break;
				}
				//less than
				if (testPos < positionsConC[j]) break;
				//greater than so keep advancing
				indexConC = j;
			}
			//enough t obs?
			int numCObs = (int)(numConC+numNonConC);
			if (numCObs < minimumReadCoverage) continue;

			//save it
			ibAL.add(new MethylatedBaseObservation (testPos, numNonConT, numConT, numNonConC, numConC));

		}
		MethylatedBaseObservation[] ibs = null;
		if (ibAL.size() !=0){
			ibs = new MethylatedBaseObservation[ibAL.size()];
			ibAL.toArray(ibs);
		}
		return ibs;
	}


	public void calculateChiSquareInRWithLookupBPs(ArrayList<SmoothingWindow> win){
		int num = win.size();

		//calculate chiSquareStatistics to one decimal
		LinkedHashMap<Integer,float[]> map = new LinkedHashMap<Integer,float[]>();
		int[] windowChiSquare10 = new int[num];
		for (int i=0; i< num; i++){
			//scores are log2, pval, obs, tNumCon, tNumNonCon, cNumCon, cNumNonCon
			float[] scores = win.get(i).getScores();
			long[][] counts = new long[][]{{(long)scores[3],(long)scores[4]} , {(long)scores[5],(long)scores[6]}};
			int chi = (int)Math.round(Num.chiSquareTestStatistic(counts) * 10);
			windowChiSquare10[i] = chi;
			Integer chiInteger = chi;
			if (map.containsKey(chiInteger) == false) map.put(chiInteger, scores);
		}

		int reducedNum = map.size();

		//build t and c arrays from reduced map
		int[][] t = new int[reducedNum][2];
		int[][] c = new int[reducedNum][2];
		int index = 0;
		for (float[] scores : map.values()){
			t[index][0] = (int)scores[3];
			t[index][1] = (int)scores[4];
			c[index][0] = (int)scores[5];
			c[index][1] = (int)scores[6];
			index++;
		}
		//calc logged pvalues, note that Double.MAX_VALUE is set for inf
		double[] pvalues = Num.chiSquareIndependenceTest(t, c, saveDirectory, fullPathToR, true);

		//assign pvalues to Points
		Point[] chiPoints = new Point[reducedNum];
		index =0;
		for (Integer integer: map.keySet()){
			float pval = (float)pvalues[index];
			chiPoints[index] = new Point(integer, pval);
			index++;
		}

		//sort by position (really chiSquareStatistic x 10)
		Arrays.sort(chiPoints, new ComparatorPointPosition());

		//split into ints and floats
		int[] chi10 = Point.extractPositions(chiPoints);
		float[] pvals = Point.extractScores(chiPoints);

		//assign pvalues to windows using binary search
		for (int i=0; i< num; i++){
			//fetch pvalue index
			int match = Arrays.binarySearch(chi10, windowChiSquare10[i]);
			//assign
			SmoothingWindow window = win.get(i);
			//scores are log2, pval, num bp scored, tNumCon, tNumNonCon, cNumCon, cNumNonCon
			float[] scores = window.getScores();
			scores[pValIndex] = pvals[match];			
		}	
	}


	public static boolean pValCheck(double pval){
		if (Double.isInfinite(pval) || Double.isNaN(pval) || pval <=0) return false;
		return true;
	}

	/**Collects and calculates a bunch of stats re the PointData.*/
	private void loadDataHashes(){	
		//null values, weird zip error.
		tConPlusPointData = null;
		tConMinusPointData = null;
		tNonConPlusPointData = null;
		tNonConMinusPointData = null;
		cConPlusPointData = null;
		cConMinusPointData = null;
		cNonConPlusPointData = null;
		cNonConMinusPointData = null;

		//fetch treatment PointData		
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (tConPointDirs);
		tConPlusPointData = PointData.convertArrayList2Array(combo[0]);
		tConMinusPointData = PointData.convertArrayList2Array(combo[1]);
		combo = PointData.fetchStrandedPointDataNoMerge (tNonConPointDirs);
		tNonConPlusPointData = PointData.convertArrayList2Array(combo[0]);
		tNonConMinusPointData = PointData.convertArrayList2Array(combo[1]);

		//garbage collect, getting weird zip error.
		System.gc();System.gc();System.gc();System.gc();

		//fetch control PointData		
		combo = PointData.fetchStrandedPointDataNoMerge (cConPointDirs);
		cConPlusPointData = PointData.convertArrayList2Array(combo[0]);
		cConMinusPointData = PointData.convertArrayList2Array(combo[1]);
		combo = PointData.fetchStrandedPointDataNoMerge (cNonConPointDirs);
		cNonConPlusPointData = PointData.convertArrayList2Array(combo[0]);
		cNonConMinusPointData = PointData.convertArrayList2Array(combo[1]);
	}



	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new DefinedRegionBisSeq(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-zA-Z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		File bedFile = null;
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i];
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': tConPointDirs = IO.extractFiles(args[++i]); break;
					case 'n': tNonConPointDirs = IO.extractFiles(args[++i]); break;
					case 'C': cConPointDirs = IO.extractFiles(args[++i]); break;
					case 'N': cNonConPointDirs = IO.extractFiles(args[++i]); break;
					case 'b': bedFile = new File(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'd': minimumReadCoverage = Integer.parseInt(args[++i]); break;
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
		if (tConPointDirs == null || tConPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your treatment converted PointData directories(s)!\n");
		if (cConPointDirs == null || cConPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your control converted PointData directories(s)!\n");
		if (tNonConPointDirs == null || tNonConPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your treatment non-converted PointData directories(s)!\n");
		if (cNonConPointDirs == null || cNonConPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your control non-converted PointData directories(s)!\n");

		//only one directory look deeper
		if (tConPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(tConPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) tConPointDirs = otherDirs;
		}
		if (cConPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(cConPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) cConPointDirs = otherDirs;
		}
		if (tNonConPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(tNonConPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) tNonConPointDirs = otherDirs;
		}
		if (cNonConPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(cNonConPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) cNonConPointDirs = otherDirs;
		}

		//check for R 
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}

		

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		else if (saveDirectory.exists() == false) saveDirectory.mkdirs();

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
		regionsFileName = Misc.removeExtension(bedFile.getName());

	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                DefinedRegionBisSeq: May 2012                     **\n" +
				"**************************************************************************************\n" +
				"Takes two condition (treatment and control) PointData from converted and non-converted\n" +
				"C bisulfite sequencing data parsed using the NovoalignBisulfiteParser and scores user\n" +
				"defined regions for differential methylation using either a fisher or chi-square test. \n" +
				"A Benjamini & Hockberg correction is applied to convert the pvalues to FDRs. Data is\n" +
				"only collected on Cs that meet the minimum read coverage threshold in both datasets. \n" +
				"The fraction differential methylation statistic is calculated by taking the\n" +
				"pseudomedian of all of the log2 paired base level fraction methylations in a given\n" +
				"region. To examine particular mC contexts (e.g. mCG), first filter your PointData\n" +
				"using the ParsePointDataContexts app.\n\n" +

				"Options:\n"+
				"-b A bed file of regions to score (tab delimited: chr start stop ...)\n"+
				"-s Save directory, full path.\n"+
				"-c Treatment converted PointData directories, full path, comma delimited. These should\n" +
				"       contain stranded chromosome specific xxx_-/+_.bar.zip files fro the NBP app.\n" +
				"       One can also provide a single directory that contains multiple PointData\n" +
				"       directories.\n" +
				"-C Control converted PointData directories, ditto. \n"+
				"-n Treatment non-converted PointData directories, ditto. \n" +
				"-N Control non-coverted PointData directories, ditto. \n"+

				"\nDefault Options:\n"+
				"-d Minimum per base read coverage, defaults to 5.\n"+
				"-r Full path to R, defaults to '/usr/bin/R'\n" +
				"\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/DefinedRegionBisStat -c /Sperm/Converted\n" +
				"      -n /Sperm/NonConverted -C /Egg/Converted -N /Egg/NonConverted -s /Res/DRBS\n" +
				"      -b /Res/CpGIslands.bed \n\n" +

		"**************************************************************************************\n");

	}
}
