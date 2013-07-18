package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import trans.tpmap.*;
import org.apache.commons.math3.stat.inference.*;;


/**For identifing differentially methylated genomic regions from two condition bisulfite sequencing data
 * @author Nix
 * */
public class BisSeq {

	//user defined fields
	private File[] tConPointDirs;
	private File[] tNonConPointDirs;
	private File[] cConPointDirs;
	private File[] cNonConPointDirs;
	private File saveDirectory;
	private File fullPathToR = new File ("/usr/bin/R");
	private int windowSize = 250;
	private int minNumObsInWindow = 5;
	private int minimumReadCoverage = 5;
	private float fdrThreshold = 30;
	private float log2RatioThreshold = 1.585f;
	private boolean printGraphs = true;
	private boolean scrambleControlData = false;

	//internal fields
	private int maxGap = 0;
	private File unstrandedPValueDirectory;
	private File windowLog2RatioDirectory;
	private File baseLog2RatioDirectory;
	private File windowsDirectory;
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
	private int numberSkippedWindows = 0;
	private int numberPassingWindows = 0;
	private WindowMaker windowMaker; 
	private int[][] windows;
	private SmoothingWindow[] smoothingWindow;
	private double cacheNumber = 1000;
	private FisherExact fischerExact = new FisherExact((int)cacheNumber);
	private ChiSquareTest chiSquare = new ChiSquareTest();
	private ArrayList<File> swiFile = new ArrayList<File>();
	private EnrichedRegion[] enrichedRegions = null;
	private String genomeVersion;

	//by chromosome
	private String[] chromosomes;
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
	public BisSeq(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);

		//load hashes
		loadDataHashes();

		fetchAllChromosomes();

		//for each chromosome
		System.out.println("\nWindow scanning chromosomes for differential methylation... ");
		for (int i=0; i< chromosomes.length; i++){
			chromosome = chromosomes[i];
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
			windowScanChromosomeBPs();

			if (smoothingWindow == null || smoothingWindow.length ==0) continue;

			//print graphs?
			if (printGraphs) writeBarFileGraphsOnCurrentSM();
			
			//filter windows for log2ratio and pval thresholds
			filterWindows();

			//save windows?
			if (smoothingWindow == null || smoothingWindow.length ==0) continue;

			File obFile = new File(windowsDirectory,chromosome);
			IO.saveObject(obFile, smoothingWindow);
			swiFile.add(obFile);
			obFile.deleteOnExit();

		}

		System.out.println();System.out.println();

		//B & H correct pvalues and make EnrichedRegions
		System.out.println("Applying a Benjamini and Hochberg FDR correction to the p-values...\n");
		correctPValues();

		//any enriched regions?
		if (enrichedRegions == null) System.out.println("No Enriched/Reduced Regions found that pass the pseLog2Ratio ("+ log2RatioThreshold +") and FDR ("+ fdrThreshold +") thresholds.\n");
		else {
			//sort by abs(log2Ratio)
			Arrays.sort(enrichedRegions, new ComparatorEnrichedRegionAbsoluteScore(log2RatioIndex));
			printEnrichedRegions();
			System.out.println(enrichedRegions.length+"\tEnriched/Reduced regions passed thresholds.\n");
		}

		//cleanup
		IO.deleteDirectory(windowsDirectory);

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/60000;
		System.out.println("Done! "+Math.round(diffTime)+" min\n");
	}

	public void printEnrichedRegions(){
		printEnrichedRegionSpreadSheet();
		//split by log2Ratio into enriched and reduced
		EnrichedRegion[] ers = null;
		ArrayList<EnrichedRegion> ersAL = new ArrayList<EnrichedRegion>();
		//enriched
		for (int i=0; i< enrichedRegions.length; i++){
			float[] scores = enrichedRegions[i].getBestWindow().getScores();
			if (scores[log2RatioIndex] > 0) ersAL.add(enrichedRegions[i]);
		}
		if (ersAL.size() != 0){
			ers = new EnrichedRegion[ersAL.size()];
			ersAL.toArray(ers);
			printBedFile(ers, new File (saveDirectory, "enrichedMethylatedRegions.bed"));
			ersAL.clear();
		}
		//reduced
		for (int i=0; i< enrichedRegions.length; i++){
			float[] scores = enrichedRegions[i].getBestWindow().getScores();
			if (scores[log2RatioIndex] < 0) ersAL.add(enrichedRegions[i]);
		}
		if (ersAL.size() != 0){
			ers = new EnrichedRegion[ersAL.size()];
			ersAL.toArray(ers);
			printBedFile(ers, new File (saveDirectory, "reducedMethylatedRegions.bed"));
			ersAL.clear();
		}



	}


	/**Writes out an bed file, used by IGB for graph display.*/
	public  void printBedFile (EnrichedRegion[] egrs, File file){
		try{
			PrintWriter out = new PrintWriter (new FileWriter (file));
			//genome version
			out.println("# genome_version = "+genomeVersion);

			//score names
			out.println("# name = best window (pseLog2Ratio : FDR : #Obs _ tNumCon : tNumNonCon : cNumCon : cNumNonCon : pseLog2TFractions, pseLog2CFractions)");
			out.println("# score = scaled 1:1000 best window abs(pseLog2Ratio)");

			//transform log2Ratio scores to 1:1000
			float[] tScores = new float[egrs.length];
			//for each ER
			for (int i=0; i< egrs.length; i++){
				float[] scores = egrs[i].getBestWindow().getScores();
				tScores[i] = scores[log2RatioIndex];
			}
			//transform
			Num.scale1To1000(tScores);

			//for each ER
			for (int i=0; i< egrs.length; i++){
				float[] scores = egrs[i].getBestWindow().getScores();
				out.print(egrs[i].getChromosome());
				out.print("\t");
				out.print(egrs[i].getStart());
				out.print("\t");
				out.print(egrs[i].getStop());
				out.print("\t");
				out.print(formatBestWindowScores(scores));
				out.print("\t");
				out.print(tScores[i]);
				out.println("\t.");
			}
			out.close();
		}catch (Exception e){
			e.printStackTrace();
		}
	}

	public String formatBestWindowScores(float[] scores){
		//scores are log2, pval, num bp scored, tNumCon, tNumNonCon, cNumCon, cNumNonCon
		StringBuilder sb = new StringBuilder();
		sb.append(Num.formatNumber(scores[log2RatioIndex], 2));
		sb.append(":");
		sb.append(Num.formatNumber(scores[pValIndex], 2));
		sb.append(":");
		sb.append((int)scores[2]);
		sb.append("_");
		sb.append((int)scores[3]);
		sb.append(":");
		sb.append((int)scores[4]);
		sb.append(":");
		sb.append((int)scores[5]);
		sb.append(":");
		sb.append((int)scores[6]);
		sb.append(":");
		sb.append(Num.formatNumber(scores[7], 4));
		sb.append(":");
		sb.append(Num.formatNumber(scores[8], 4));
		return sb.toString();
	}


	/**Writes out an excel compatible tab delimited spreadsheet with hyperlinks for IGB.*/
	public void printEnrichedRegionSpreadSheet(){
		try{
			File file = new File(saveDirectory, "differentialMethylatedRegions.xls");
			PrintWriter out = new PrintWriter (new FileWriter (file));
			//print header line
			out.println("#"+genomeVersion+"_IGBHyperLinks\tChr\tER_Start\tER_Stop\t#WindowsInER\t"+
			"BW_Start\tBW_Stop\tBW_PseMedianLog2Ratio\tBW_FDR\tBW_#Obs\tBW_#TCon\tBW_#TNonCon\tBW_#CCon\tBW_#CNonCon\tBW_PseMedianLog2TFractions\tBW_PseMedianLog2CFractions");

			String url = "=HYPERLINK(\"http://localhost:7085/UnibrowControl?version="+genomeVersion+"&seqid=";
			String tab = "\t";
			//for each ER
			for (int i=0; i< enrichedRegions.length; i++){
				//url
				int winStart = enrichedRegions[i].getStart() - 10000;
				if (winStart < 0) winStart = 0;
				int winEnd = enrichedRegions[i].getStop() + 10000;
				out.print(url+enrichedRegions[i].getChromosome()+"&start="+winStart+"&end="+winEnd+"\",\""+(i+1)+"\")\t");
				//chrom
				out.print(enrichedRegions[i].getChromosome()); out.print(tab);
				//start
				out.print(enrichedRegions[i].getStart()); out.print(tab);
				//stop
				out.print(enrichedRegions[i].getStop()); out.print(tab);
				//#Windows
				out.print(enrichedRegions[i].getNumberOfWindows()); out.print(tab);
				//best window data
				SmoothingWindow bw = enrichedRegions[i].getBestWindow();
				//coordinates
				out.print(bw.getStart()); out.print(tab);
				out.print(bw.getStop()); out.print(tab);
				//scores are log2, pval, num bp scored, tNumCon, tNumNonCon, cNumCon, cNumNonCon, pseTFrac, pseCFrac
				float[] scores = bw.getScores();
				//	log2Ratio
				out.print(scores[0]); 
				out.print(tab);
				//	FDR
				out.print(scores[1]);
				//	# obs and counts as ints
				for (int j=2; j< scores.length-2; j++){
					out.print(tab);
					out.print((int)scores[j]);
				}
				//fractions
				out.print(tab);
				out.print(scores[7]);
				out.print(tab);
				out.print(scores[8]);

				out.println();
			}
			out.close();
		}catch (Exception e){
			System.out.println("\nError: problem printing spreadsheet report");
			e.printStackTrace();
		}
	}

	private void filterWindows(){
		ArrayList<SmoothingWindow> goodWin = new ArrayList<SmoothingWindow>();

		for (int i=0; i< smoothingWindow.length; i++){
			float[] scores = smoothingWindow[i].getScores();
			if ( (Math.abs(scores[log2RatioIndex]) < log2RatioThreshold) || (scores[pValIndex]< fdrThreshold) ) numberSkippedWindows++;
			else goodWin.add(smoothingWindow[i]);
		}
		smoothingWindow = new SmoothingWindow[goodWin.size()];
		goodWin.toArray(smoothingWindow);
		numberPassingWindows+= smoothingWindow.length;
	}

	public void correctPValues(){
		if (numberPassingWindows == 0) return;

		//load Point[] with all of the pvalues
		Point[] point = new Point[numberPassingWindows];
		//for each chromosome of data
		int index = 0;
		for (int i=0; i< swiFile.size(); i++){
			SmoothingWindow[] sw = (SmoothingWindow[])IO.fetchObject(swiFile.get(i));
			for (int j=0; j<sw.length; j++){
				point[index] = new Point(index, sw[j].getScores()[pValIndex]);
				index++;
			}
		}
		//sort by score smallest -10Log10(pval) to largest
		Arrays.sort(point, new ComparatorPointAscendingScore());
		//correct
		Point.benjaminiHochbergCorrect(point, numberSkippedWindows);
		//sort back to original position
		Arrays.sort(point, new ComparatorPointPosition());

		//make ER generator threshold on fdr (already thresholded on log2ratio)
		EnrichedRegionMaker erm = new EnrichedRegionMaker(maxGap, new int[]{pValIndex}, new float[]{fdrThreshold}, log2RatioIndex);
		ComparatorSmoothingWindowPosition comparator = new ComparatorSmoothingWindowPosition();
		boolean noSWs = true;

		//for each chromosome of data
		index = 0;
		for (int i=0; i< swiFile.size(); i++){
			ArrayList<SmoothingWindow> goodSWEnriched = new ArrayList<SmoothingWindow>();
			ArrayList<SmoothingWindow> goodSWReduced = new ArrayList<SmoothingWindow>();
			SmoothingWindow[] sw = (SmoothingWindow[])IO.fetchObject(swiFile.get(i));
			chromosome = swiFile.get(i).getName();

			//swap pvals for fdrs
			for (int j=0; j<sw.length; j++){				
				float[] scores = sw[j].getScores();
				scores[pValIndex] = point[index++].getScore();
				//save it? threshold on fdr, already thresholded on log2ratio
				if (scores[pValIndex] >= fdrThreshold) {
					if (scores[log2RatioIndex] > 0 ) goodSWEnriched.add(sw[j]);
					else goodSWReduced.add(sw[j]);
				}
				else sw[j] = null;
			}
			//convert als to [] and make ERS
			//	enriched regions
			if (goodSWEnriched.size() !=0){
				noSWs = false;
				sw = new SmoothingWindow[goodSWEnriched.size()];
				goodSWEnriched.toArray(sw);
				Arrays.sort(sw, comparator);
				erm.addEnrichedRegions(sw, chromosome);
			}
			//	reduced regions
			if (goodSWReduced.size() !=0){
				noSWs = false;
				sw = new SmoothingWindow[goodSWReduced.size()];
				goodSWReduced.toArray(sw);
				Arrays.sort(sw, comparator);
				erm.addEnrichedRegions(sw, chromosome);
			}

		}
		//make enriched regions? otherwise leave null
		if (noSWs == false) {
			enrichedRegions = new EnrichedRegion[erm.getEnrichedRegionsAL().size()];
			erm.getEnrichedRegionsAL().toArray(enrichedRegions);

			//set best window log2 and FDR scores in ER
			for (EnrichedRegion er : enrichedRegions){
				float[] scores = er.getBestWindow().getScores();
				er.setLog2Ratio(scores[pValIndex]);
				er.setBinomialPValue(scores[pValIndex]);
			}
		}

	}

	/**Fetches the names of all the chromosomes in the data excluding lambda and phiX if present.*/
	public void fetchAllChromosomes(){
		HashSet<String> c = new HashSet<String>();
		Iterator<String> it = tConPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = tConMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = tNonConPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = tNonConMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = cConPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = cConMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = cNonConPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = cNonConMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());

		//any lambda, phiX, or adapter found?
		it = c.iterator();
		Pattern lambda = Pattern.compile(".*lambda.*", Pattern.CASE_INSENSITIVE);
		Pattern phiX = Pattern.compile(".*phix.*", Pattern.CASE_INSENSITIVE);
		Pattern adapter = Pattern.compile(".*adapt.*", Pattern.CASE_INSENSITIVE);
		ArrayList<String> toRemove = new ArrayList<String>();
		while (it.hasNext()){
			String chr = it.next();
			if (lambda.matcher(chr).matches() || phiX.matcher(chr).matches() || adapter.matcher(chr).matches()) toRemove.add(chr); 
		}
		for (int i=0; i< toRemove.size(); i++) c.remove(toRemove.get(i));
		chromosomes=  Misc.hashSetToStringArray(c);
		Arrays.sort(chromosomes);
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


	public void windowScanChromosomeBPs(){
//System.out.println("Here");
		smoothingWindow = null;

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
		
		if (scrambleControlData){
			//swap control observations, not positions
			Random rng = new Random();       
		    for (int n =0; n< mbo.length; n++) {
		        int k = rng.nextInt(mbo.length); 
		        MethylatedBaseObservation curr = mbo[n];
		        float currCNon = curr.getNonConC();
		        float currCCon = curr.getConC();
		        MethylatedBaseObservation rand = mbo[k];
		        float randCNon = rand.getNonConC();
		        float randCCon = rand.getConC();
		        curr.setNonConC(randCNon);
		        curr.setConC(randCCon);
		        rand.setNonConC(currCNon);
		        rand.setConC(currCCon);
		    }
		}
		
		//fetch the positions 
		int[] positions = MethylatedBaseObservation.fetchPositions(mbo);
		
		//write out ratio of ratios?
		if (printGraphs) {
			float[] diffMeth = MethylatedBaseObservation.fetchLog2DifferentialFractionMethylation(mbo);
			saveBaseDiffMethPointData(positions, diffMeth);
		}

		//make windows using all of the positions
		makeWindows(positions);
		if (windows.length == 0) return;
		
		//scan
		scoreWindowsForDifferentialMethylation(mbo, positions);

	}


	public void scoreWindowsForDifferentialMethylation(MethylatedBaseObservation[] mbo, int[] positions){

		ArrayList<SmoothingWindow> pValueWindows = new ArrayList<SmoothingWindow>();
		ArrayList<SmoothingWindow> forRWindows = new ArrayList<SmoothingWindow>();

		//for each window
		for (int i=0; i< windows.length; i++){

			int[] indexes = Num.findIndexes(windows[i][0], windows[i][1], positions);
			
			//for each index, calculate fraction methylated
			float[] fractions = new float[indexes[1]- indexes[0]];
			float[] tFractions = new float [fractions.length];
			float[] cFractions = new float [fractions.length];
			int f = 0;
			for (int x=indexes[0]; x< indexes[1]; x++) {
				tFractions[f] = Num.log2(mbo[x].getFractionMethylatedT());
				cFractions[f] = Num.log2(mbo[x].getFractionMethylatedC());
				fractions[f] = Num.log2(mbo[x].getFractionMethylatedT()/mbo[x].getFractionMethylatedC());
				f++;
			}
			float pseLog2 = (float) Num.pseudoMedian(fractions);
			float pseTFrac = (float) Num.pseudoMedian(tFractions);
			float pseCFrac = (float) Num.pseudoMedian(cFractions);

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

			long[][] counts = new long[][]{{tNumCon, tNumNonCon},{cNumCon,cNumNonCon}};

			//scores are log2, pval, num bp scored, tNumCon, tNumNonCon, cNumCon, cNumNonCon, pseTFrac, pseCFrac
			float[] scores = new float[]{pseLog2, 0, fractions.length, (float)tNumCon, (float)tNumNonCon, (float)cNumCon, (float)cNumNonCon, pseTFrac, pseCFrac};

			//make smoothing window
			SmoothingWindow win = new SmoothingWindow (windows[i][0], windows[i][1], scores);

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
			if (sendToR) forRWindows.add(win); 
			else pValueWindows.add(win);
		}

		//make []
		int numPValWindows = pValueWindows.size();
		int numRWindows = forRWindows.size();
		smoothingWindow = new SmoothingWindow[numPValWindows + numRWindows];
		for (int i=0; i< numPValWindows; i++) smoothingWindow[i] = pValueWindows.get(i);


		//any difficult counts?
		if (numRWindows !=0) {				
			//use R to calculate
			calculateChiSquareInRWithLookupBPs(forRWindows);
			//add to []
			int index = 0;
			for (int i=numPValWindows; i< smoothingWindow.length; i++) {
				smoothingWindow[i] = forRWindows.get(index++);	
			}				
			//sort			
			Arrays.sort(smoothingWindow, new ComparatorSmoothingWindowPosition());
		}

		//check pvalues
		for (SmoothingWindow sm : smoothingWindow){
			float[] s = sm.getScores();
			//need to watch out for slight neg -10log10(pval)s, these are insignificant
			if (pValCheck(s[pValIndex]) == false) s[pValIndex] = 0;	
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

	/**Makes a common set of windows using merged positions.*/
	public void makeWindows(int[] positions){
		windows = windowMaker.makeWindows(positions);
		//assign bp positions
		for (int i=0; i< windows.length; i++){
			windows[i][0] = positions[windows[i][0]];
			windows[i][1] = positions[windows[i][1]]+1;	//last base isn't included
		}
	}

	/**Writes stair step window bar graph files*/
	public void writeBarFileGraphsOnCurrentSM(){

		//build info object
		Info info;
		if (tConMergedChromPlus != null) info= tConMergedChromPlus.getInfo();
		else info = tConMergedChromMinus.getInfo();
		HashMap<String,String> notes = new HashMap<String,String>();
		notes.put(BarParser.DESCRIPTION_TAG, "Pseudo median of per base log2 ratios (fraction methylated in T/ fraction methylated in C) in window");
		notes.put(BarParser.UNIT_TAG, "Ratio");
		info.setNotes(notes);
		info.setStrand(".");	

		//save log2 ratio data
		saveSmoothedHeatMapData (0, smoothingWindow, info, windowLog2RatioDirectory, "#FF0000", true); //red

		//save pvalues
		notes.put(BarParser.DESCRIPTION_TAG, "Uncorrected p-values of unstranded differential methylation");
		notes.put(BarParser.UNIT_TAG, "-10Log10(pval)");
		saveSmoothedHeatMapData (1, smoothingWindow, info, unstrandedPValueDirectory, "#00FF00", false); //green

	}

	/**Saves bar heatmap/ stairstep graph files*/
	public void saveSmoothedHeatMapData (int scoreIndex, SmoothingWindow[] sm, Info info, File dir, String color, boolean posNeg){
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = info.getNotes();		
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		//color red
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, color);
		//what's window size
		map.put(BarParser.WINDOW_SIZE, windowSize+"");
		//save in info
		info.setNotes(map);
		//get heatmap positions and values
		PointData pd;		
		if (posNeg){
			HeatMapMakerPosNeg hm = new HeatMapMakerPosNeg(scoreIndex, 0, 0);
			pd = hm.makeHeatMapPositionValues(sm);
		}
		else {
			HeatMapMaker hm = new HeatMapMaker(scoreIndex,0);
			pd = hm.makeHeatMapPositionValues(sm, false);
		}

		pd.setInfo(info);
		pd.writePointData(dir);
		//clean up
		pd.nullPositionScoreArrays();
	}

	/**Saves bar point graph files*/
	public void saveBaseDiffMethPointData (int[] positions, float[] scores){
		//build info object
		Info info;
		if (tConMergedChromPlus != null) info= tConMergedChromPlus.getInfo();
		else info = tConMergedChromMinus.getInfo();
		info.setStrand(".");
		
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = info.getNotes();
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, "#A020F0"); //purple
		map.put(BarParser.DESCRIPTION_TAG, "Base level log2 ratios (fraction methylated in T/ fraction methylated in C)");
		map.put(BarParser.UNIT_TAG, "Log2Ratio");
		//save in info
		info.setNotes(map);
		//write bar
		PointData pd = new PointData();
		pd.setInfo(info);
		pd.setPositions(positions);
		pd.setScores(scores);
		pd.writePointData(baseLog2RatioDirectory);
		//clean up
		pd.nullPositionScoreArrays();
		positions = null;
		scores = null;
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BisSeq(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-zA-Z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i];
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': fdrThreshold = Float.parseFloat(args[++i]); break;
					case 'l': log2RatioThreshold = Float.parseFloat(args[++i]); break;
					case 'c': tConPointDirs = IO.extractFiles(args[++i]); break;
					case 'n': tNonConPointDirs = IO.extractFiles(args[++i]); break;
					case 'C': cConPointDirs = IO.extractFiles(args[++i]); break;
					case 'N': cNonConPointDirs = IO.extractFiles(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'w': windowSize = Integer.parseInt(args[++i]); break;
					case 'm': minNumObsInWindow = Integer.parseInt(args[++i]); break;
					case 'd': minimumReadCoverage = Integer.parseInt(args[++i]); break;
					case 'g': printGraphs = false; break;
					case 'a': scrambleControlData = true; break;
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

		//set half peak shift and windowSize
		if (windowSize == 0 ) Misc.printExit("\nPlease enter a positive length for the window size.\n");

		//make window maker 
		windowMaker = new WindowMaker(windowSize,minNumObsInWindow);

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		else if (saveDirectory.exists() == false) saveDirectory.mkdirs();
		if (printGraphs) {
			unstrandedPValueDirectory = new File (saveDirectory, "WindowUnCorrPValues");
			unstrandedPValueDirectory.mkdir();
			windowLog2RatioDirectory = new File (saveDirectory, "WindowPseLog2Ratio");
			windowLog2RatioDirectory.mkdir();
			baseLog2RatioDirectory = new File (saveDirectory, "BasePseLog2Ratio");
			baseLog2RatioDirectory.mkdir();
		}
		windowsDirectory = new File (saveDirectory, "BinaryWindowData");
		windowsDirectory.mkdir();

		//params
		System.out.println("Params:");
		System.out.println(fdrThreshold+ "\tFDR threshold");
		System.out.println(log2RatioThreshold+ "\tLog2Ratio threshold");
		System.out.println(windowSize+ "\tWindow size");
		System.out.println(minNumObsInWindow+"\tMinimum # obs in window");
		System.out.println(minimumReadCoverage+"\tMinimum per base read coverage");
		System.out.println(printGraphs+ "\tPrint graphs");
		System.out.println(scrambleControlData+ "\tScramble control data");
	}	



	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                  BisSeq: July 2013                               **\n" +
				"**************************************************************************************\n" +
				"Takes two condition (treatment and control) PointData from converted and non-converted\n" +
				"C bisulfite sequencing data parsed using the NovoalignBisulfiteParser and scores\n" +
				"regions for differential methylation using either a fisher exact or chi-square test \n" +
				"for changes in methylation.  A Benjamini & Hockberg correction is applied to convert\n" +
				"the pvalues to FDRs. Data is only collected on bases that meet the minimum\n" +
				"read coverage threshold in both datasets.  The fraction differential methylation\n" +
				"statistic is calculated by taking the pseudomedian of all of the log2 paired base level\n" +
				"fraction methylations in a given window. Overlapping windows that meet both the\n" +
				"FDR and pseLog2Ratio thresholds are merged when generating enriched and reduced\n" +
				"regions. BisSeq generates several tracks for browsing and lists of differentially\n" +
				"methlated regions. To examine only mCG contexts, first filter your PointData using the\n" +
				"ParsePointDataContexts app. \n\n" +

				"Options:\n"+
				"-s Save directory, full path.\n"+
				"-c Treatment converted PointData directories, full path, comma delimited. These should\n" +
				"       contain stranded chromosome specific xxx_-/+_.bar.zip files fro the NBP app.\n" +
				"       One can also provide a single directory that contains multiple PointData\n" +
				"       directories.\n" +
				"-C Control converted PointData directories, ditto. \n"+
				"-n Treatment non-converted PointData directories, ditto. \n" +
				"-N Control non-coverted PointData directories, ditto. \n"+
				"-a Scramble control data.\n"+

				"\nDefault Options:\n"+
				"-d Minimum per base read coverage, defaults to 5.\n"+
				"-w Window size, defaults to 250.\n"+
				"-m Minimum number reads in window, defaults to 5. \n" +
				"-f FDR threshold, defaults to 30 (-10Log10(0.01)).\n"+
				"-l Log2Ratio threshold, defaults to 1.585 (3x).\n"+
				"-r Full path to R, defaults to '/usr/bin/R'\n" +
				"-g Don't print graph files.\n"+


				"\n"+

				"Example: java -Xmx10G -jar pathTo/USeq/Apps/BisStat -c /Sperm/Converted -n \n" +
				"      /Sperm/NonConverted -C /Egg/Converted -N /Egg/NonConverted -s /Res/BisSeq\n" +
				"      -w 500 -m 10 -l 2 -f 50 \n\n" +

		"**************************************************************************************\n");

	}
}
