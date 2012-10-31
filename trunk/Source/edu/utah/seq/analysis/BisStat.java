package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;

import util.bio.cluster.HierarchicalClustering;
import util.bio.parsers.MultiFastaParser;
import util.bio.seq.Seq;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.parsers.*;
import trans.tpmap.*;

/**For processing bisulfite sequencing data
 * @author Nix
 * */
public class BisStat {

	//user defined fields
	private File[] convertedPointDirs;
	private File[] nonConvertedPointDirs;
	private File saveDirectory;
	private File fractionNonConvertedDirectory;
	private File fractionNonConvertedDirectoryPlus;
	private File fractionNonConvertedDirectoryMinus;
	private File windowFractionNonConvertedWindow;
	private File windowDensity1stQuartileWindow;
	private File windowDensity2_3QuartileWindow;
	private File windowDensity4thQuartileWindow;
	private File windowFractionNonConvertedPoint;
	private File serializedWindowObjectsDirectory;
	private File fdrNonConvertedDirectory;
	private int windowSize = 1000;
	private float minimumReadCoverage = 8;
	private HashMap<String, File> chromosomeFastaFiles = new HashMap<String, File>();
	private boolean printGraphs = true;
	private File fullPathToR = new File ("/usr/bin/R");
	private double expectedNonConverted = 0.005;
	private float fdrThreshold = 20;
	private boolean useLambda = false;
	private int minimumCsInWindow = 5;
	private boolean performStrandedAnalysis = true;
	private float firstQuartileThreshold = 0.25f;
	private float fourthQuartileThreshold = 0.75f;

	//internal fields
	private String[] chromosomes;
	private String genomicSequence = null;
	private int genomicSequenceLengthMinus3 = 0;
	private HashMap<String,PointData[]> convertedPlusPointData;
	private HashMap<String,PointData[]> convertedMinusPointData;
	private HashMap<String,PointData[]> nonConvertedPlusPointData;
	private HashMap<String,PointData[]> nonConvertedMinusPointData;
	private HashMap<String, BaseContext> baseContexts = null;
	private Pattern CG_Pattern = Pattern.compile("..CG.");
	private String lambdaChromosome = null;
	private String phiXChromosome = null;
	private ArrayList<float[]> baseFractionMethylation = new ArrayList<float[]>();

	//pvalue caching
	private float numberCachedBinPVals = 500;
	private Float[][] cachedPValues;
	private Float zeroFloat = new Float(0.0f);
	private ArrayList<Point> pointAL = new ArrayList<Point>();
	private Point[] point;
	private float minimumPValueToCorrect;
	private long pValueOffset = 0;
	private int fdrIndex = 0;

	//window scanning
	private WindowMaker windowMaker; 
	private int[][] windows;
	private SmoothingWindow[] smoothingWindow;
	private String[] scoreNames;
	private String[] scoreDescriptions;
	private String[] scoreUnits;

	//by chromosome
	private String chromosome;
	private PointData convertedMergedChromPlus = null;
	private PointData convertedMergedChromMinus = null;	
	private PointData nonConvertedMergedChromPlus = null;
	private PointData nonConvertedMergedChromMinus = null;	
	private PointData convertedChrom = null;
	private PointData nonConvertedChrom = null;


	//constructors
	/**Stand alone.*/
	public BisStat(String[] args){
		long startTime = System.currentTimeMillis();

		//set fields
		processArgs(args);
		setBaseContexts();

		//scan each chromosome
		scan();

		//print contexts
		printBaseContexts();
		printGenomicContextStats();

		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}

	//methods

	public void printBaseContexts(){
		//convert to an array
		Iterator<BaseContext> i = baseContexts.values().iterator();
		BaseContext[] bc = new BaseContext[baseContexts.size()];
		int index = 0;
		while (i.hasNext()) bc[index++] = i.next();
		for (int x=0; x< bc.length; x++) bc[x].setFractionNonConverted();
		Arrays.sort(bc);
		System.out.println("\nStats per xxCxx context for the number of sequenced observations (Obs) and genomic loci (Gen)");
		System.out.println("CBaseContext\tNonConObs\tConObs\tNonCon/TotObs\tNonConGen\tConGen\tNonConGen/TotGen");
		Misc.printArray(bc);

		//print histograms
		System.out.println("\nDistribution of the fraction mC for each individual genomic contexts meeting a minimum FDR threshold of "+fdrThreshold+" and minimum read coverage of "+minimumReadCoverage);
		for (int x=0; x< bc.length; x++){
			if (bc[x].getHistogram().getTotalBinCounts() == 0) continue;
			System.out.println();
			System.out.println(bc[x].getSequence());
			bc[x].getHistogram().printScaledHistogram();
		}

	}

	public void scan(){
		//fetch counts
		System.out.println("\nCalculating read count stats...");
		calculateReadCountStatistics();

		//this also removes the phiX and lambda data from that to be scanned if resent
		fetchAllChromosomes();

		//scan lambda
		if (useLambda) scanLambda();

		//make pvalues after potentially modifying the expect
		cachedPValues = Num.convertToFloatObjectArray(Num.binomialPValMatrix((int)numberCachedBinPVals+1, expectedNonConverted, saveDirectory, fullPathToR, true));

		//scan for pvalues
		String oldChrom = "";
		System.out.print("\nScanning data for pvalues... ");
		for (int i=0; i< chromosomes.length; i++){		
			chromosome = chromosomes[i];	

			//fetch the chrom sequence?
			if (oldChrom.equals(chromosome) ==false){
				File seqFile = chromosomeFastaFiles.get(chromosome);
				if (seqFile == null) {
					System.err.println("\n\tWarning, could not find a fasta file for "+chromosome+", skipping!");
					continue;
				}
				MultiFastaParser mfp = new MultiFastaParser(seqFile);				
				genomicSequence = mfp.getSeqs()[0].toUpperCase();
				genomicSequenceLengthMinus3 = genomicSequence.length()-3;
				oldChrom = chromosome;
			}

			System.out.print(chromosome+" ");
			//fetch data
			if (fetchDataAndRemove() == false) continue;
			//calculate base level pvalues
			if (performStrandedAnalysis){
				baseScanChromosomeContextsForPValues(nonConvertedMergedChromPlus, convertedMergedChromPlus);
				baseScanChromosomeContextsForPValues(nonConvertedMergedChromMinus, convertedMergedChromMinus);
			}
			else {
				convertedChrom = PointData.mergePairedPointDataNoSumming(convertedMergedChromPlus, convertedMergedChromMinus);
				nonConvertedChrom = PointData.mergePairedPointDataNoSumming(nonConvertedMergedChromPlus, nonConvertedMergedChromMinus);
				baseScanChromosomeContextsMergingStrandsForPValues();
			}
		}
		System.out.println();

		//correct using B&H	
		System.out.println("\nApplying the Benjamini & Hochberg FDR correction...");
		correctPValues();

		//refetch data pointers, reset index
		calculateReadCountStatistics();
		fetchAllChromosomes();
		fdrIndex = 0;

		//scan each chromosome for significance
		System.out.println("\nRescanning data...\n");
		System.out.println("Stranded chromosome per base fraction methylation stats:");
		System.out.println("Chr\t#NonCon\t#Con\tMean\tMedian\tStdDev\tMin\tMax\t10th\t90th");


		for (int i=0; i< chromosomes.length; i++){			
			chromosome = chromosomes[i];
			//fetch data
			if (fetchDataAndRemove() == false) continue;
			//fetch the chrom sequence?
			if (oldChrom.equals(chromosome) ==false){
				File seqFile = chromosomeFastaFiles.get(chromosome);
				if (seqFile == null) continue;
				MultiFastaParser mfp = new MultiFastaParser(seqFile);				
				genomicSequence = mfp.getSeqs()[0].toUpperCase();
				genomicSequenceLengthMinus3 = genomicSequence.length()-3;
				oldChrom = chromosome;
			}

			//merge strands?
			if (printGraphs || performStrandedAnalysis == false){
				convertedChrom = PointData.mergePairedPointDataNoSumming(convertedMergedChromPlus, convertedMergedChromMinus);
				nonConvertedChrom = PointData.mergePairedPointDataNoSumming(nonConvertedMergedChromPlus, nonConvertedMergedChromMinus);
			}

			//non stranded
			if (performStrandedAnalysis == false){
				PointData[] pd = baseScanChromosomeContextsMergingStrands();
				if (pd == null){
					System.err.println("\n"+chromosome+"\tNo base contexts. Skipping.");
				}
				else {
					baseFractionMethylation.add(pd[0].getScores());
					//graphs?
					if (printGraphs){
						//window scan?
						windowScanChromosome(pd[0]);
						pd[0].writePointData(fractionNonConvertedDirectory);
						pd[1].writePointData(fdrNonConvertedDirectory);
					}
					//chrom stats
					System.out.println();
					float[] fracs = pd[0].getScores();
					Arrays.sort(fracs);
					String stats = "";
					if (fracs.length > 10) stats = Num.statFloatArray(fracs);
					System.out.println(chromosome +".\t"+ (int)nonConvertedChrom.getInfo().getScoreTotal() +"\t"+(int)convertedChrom.getInfo().getScoreTotal() +"\t"+stats );
					Histogram histogram = new Histogram(0, 1.1, 11);
					histogram.countAll(fracs);
					histogram.printScaledHistogram();
				}
			}

			//stranded
			else{
				//calculate base level stats			
				PointData[] plusPD = baseScanChromosomeContexts(nonConvertedMergedChromPlus, convertedMergedChromPlus, false);		
				PointData[] minusPD = baseScanChromosomeContexts(nonConvertedMergedChromMinus, convertedMergedChromMinus, true);

				if (plusPD == null || minusPD == null){
					System.err.println("\n"+chromosome+"\tNo base contexts. Skipping.");
				}

				else {
					//save fractions
					baseFractionMethylation.add(plusPD[0].getScores());
					baseFractionMethylation.add(minusPD[0].getScores());

					//print graphs?
					if (printGraphs){
						//write stranded fractions
						plusPD[0].writePointData(fractionNonConvertedDirectoryPlus);
						minusPD[0].writePointData(fractionNonConvertedDirectoryMinus);
						//merge fractions for window scanning
						PointData merged = PointData.mergePairedPointDataNoSumming(plusPD[0], minusPD[0]);
						windowScanChromosome(merged);
						//write merged fdrs
						merged = PointData.mergePairedPointDataNoSumming(plusPD[1], minusPD[1]);
						merged.writePointData(fdrNonConvertedDirectory);
					}

					//print chromosome stats
					double numNonPlus = nonConvertedMergedChromPlus.getInfo().getScoreTotal();
					double numConPlus = convertedMergedChromPlus.getInfo().getScoreTotal();
					double numNonMinus = nonConvertedMergedChromMinus.getInfo().getScoreTotal();
					double numConMinus = convertedMergedChromMinus.getInfo().getScoreTotal();

					//calc for plus
					System.out.println();
					float[] fracs = plusPD[0].getScores();
					//enought to stat?
					if (fracs.length<10) System.out.println(chromosome+"+\tskipping, too few observations");
					else {
						Arrays.sort(fracs);
						System.out.println(chromosome +"+\t"+ (int)numNonPlus +"\t"+(int)numConPlus +"\t"+Num.statFloatArray(fracs) );
						Histogram histogram = new Histogram(0, 1.1, 11);
						histogram.countAll(fracs);
						histogram.printScaledHistogram();
					}
					//for minus
					fracs = minusPD[0].getScores();
					if (fracs.length<10) System.out.println(chromosome+"-\tskipping, too few observations");
					else {
						Arrays.sort(fracs);
						System.out.println(chromosome +"-\t"+ (int)numNonMinus +"\t"+(int)numConMinus +"\t"+Num.statFloatArray(fracs) );
						Histogram histogram = new Histogram(0, 1.1, 11);
						histogram.countAll(fracs);
						histogram.printScaledHistogram();
					}
				}
			}
		}
		System.out.println();

	}

	public void correctPValues(){
		//make Points
		point = new Point[pointAL.size()];
		pointAL.toArray(point);
		pointAL = null;
		//sort by score smallest to largest
		Arrays.sort(point, new ComparatorPointAscendingScore());
		//correct
		Point.benjaminiHochbergCorrect(point, pValueOffset);
		//sort back to original position
		Arrays.sort(point, new ComparatorPointPosition());
	}


	public void scanLambda(){
		if (lambdaChromosome ==null) {
			System.out.println("Didn't find any lambda sequence data?! Skipping.");
			return;
		}
		//fetch data
		chromosome = lambdaChromosome;
		//this removes the lambda data
		fetchDataAndRemove();
		//calc fraction non-converted, the expected non-converted
		double numberNonCon = Num.sumArrayReturnDouble(nonConvertedMergedChromPlus.getScores());
		double numberCon = Num.sumArrayReturnDouble(convertedMergedChromPlus.getScores());
		numberNonCon += Num.sumArrayReturnDouble(nonConvertedMergedChromMinus.getScores());
		numberCon += Num.sumArrayReturnDouble(convertedMergedChromMinus.getScores());
		expectedNonConverted = numberNonCon/(numberNonCon+numberCon);
		String fraction = Num.formatNumber(expectedNonConverted, 5);
		System.out.println("\nUsing Lambda data to set the expected fraction non-converted Cs to "+fraction+" ("+(int)numberNonCon+"/("+(int)numberNonCon+"+"+(int)numberCon+"))");

		return;
	}


	/**Fetches the names of all the chromosomes in the data excluding lambda and phiX if present.*/
	public void fetchAllChromosomes(){
		HashSet<String> c = new HashSet<String>();
		Iterator<String> it = convertedPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = convertedMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = nonConvertedPlusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());
		it = nonConvertedMinusPointData.keySet().iterator();
		while (it.hasNext()) c.add(it.next());

		//any lambda or phiX found?
		it = c.iterator();
		Pattern lambda = Pattern.compile(".*lambda.*", Pattern.CASE_INSENSITIVE);
		Pattern phiX = Pattern.compile(".*phix.*", Pattern.CASE_INSENSITIVE);
		lambdaChromosome = null;
		phiXChromosome = null;
		while (it.hasNext()){
			String chr = it.next();
			if (lambda.matcher(chr).matches()) lambdaChromosome = chr;
			else if (phiX.matcher(chr).matches()) phiXChromosome = chr;
		}
		if (lambdaChromosome != null) c.remove(lambdaChromosome);
		if (phiXChromosome != null) c.remove(phiXChromosome);

		chromosomes=  Misc.hashSetToStringArray(c);
	}

	/**Fetchs the data for a particular chromosome. Returns true if all four datasets were found.*/
	public boolean fetchDataAndRemove(){
		ArrayList<PointData> al = null;
		PointData[] pd;
		//merge converted
		convertedMergedChromPlus = null;
		if (convertedPlusPointData.containsKey(chromosome)) {
			pd = convertedPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			convertedMergedChromPlus = PointData.mergePointData(al, false, true);
		}
		convertedMergedChromMinus = null;
		if (convertedMinusPointData.containsKey(chromosome)) {
			pd = convertedMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			convertedMergedChromMinus = PointData.mergePointData(al, false, true);
		}
		//merge nonConverted
		nonConvertedMergedChromPlus = null;
		if (nonConvertedPlusPointData.containsKey(chromosome)) {
			pd = nonConvertedPlusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			nonConvertedMergedChromPlus = PointData.mergePointData(al, false, true);
		}
		nonConvertedMergedChromMinus = null;
		if (nonConvertedMinusPointData.containsKey(chromosome)) {
			pd = nonConvertedMinusPointData.remove(chromosome);
			al = PointData.convertArray2ArrayList(pd);
			nonConvertedMergedChromMinus = PointData.mergePointData(al, false, true);
		}
		pd = null;
		al = null;
		if ( convertedMergedChromPlus == null || nonConvertedMergedChromPlus == null || convertedMergedChromMinus == null || nonConvertedMergedChromMinus == null) return false;
		return true;
	}

	/**Window scans a chromosome.*/
	public void windowScanChromosome(PointData pd){
		//make windows 
		makeWindows(pd.getPositions());
		//any windows?
		if (windows.length == 0){
			System.out.println("\n\tSkipping "+chromosome+". No windows found with minimum reads of "+minimumCsInWindow+" within a window size of "+windowSize);
			return;
		}
		//scan
		scoreWindows(pd);
		//save array
		IO.saveObject( new File(serializedWindowObjectsDirectory, chromosome), smoothingWindow);
		//write 		
		writeBarFileGraphsOnCurrentSM();
	}

	/*Unstranded.*/
	private void scoreWindows(PointData pd){
		smoothingWindow = new SmoothingWindow[windows.length];

		//for each window
		for (int i=0; i< windows.length; i++){

			//fetch points 
			Point[] pts = pd.fetchPoints(windows[i][0], windows[i][1]);

			//calculate mean
			float[] f = Point.extractScores(pts);
			Arrays.sort(f);
			float fraction = (float)Num.mean(f);

			//count quartiles
			float first = 0.0f;
			float secondThird = 0.0f;
			float fourth = 0.0f;
			for (float test: f){
				if (test<= firstQuartileThreshold) first++;
				else if (test>= fourthQuartileThreshold) fourth++;
				else secondThird++;
			}
			float total = f.length;

			//make window
			smoothingWindow[i] = new SmoothingWindow (windows[i][0], windows[i][1], new float[]{fraction, first/total, secondThird/total, fourth/total});
		}
	}

	public static int[] fetchUniquePositions (Point[][] pts){
		ArrayList<int[]> intsAL = new ArrayList<int[]>();
		for (int i=0; i< pts.length; i++){
			if (pts[i] != null) intsAL.add(Point.extractPositions(pts[i]));
		}
		if (intsAL.size()==0) return null;
		return Num.returnUniques(intsAL);
	}


	/**Scores a chromosome for non-converted to total at base level.
	 * @return PointData[2] fraction and FDRs for a given strand*/
	private PointData[] baseScanChromosomeContexts(PointData nonCon, PointData con, boolean negativeStrand){
		if (nonCon == null || con == null) return null;
		//fetch arrays
		int[] positionsNonCon = nonCon.getPositions();
		float[] readsNonCon = nonCon.getScores();
		int[] positionsCon = con.getPositions();
		float[] readsCon = con.getScores();
		//make containers for graph data
		ArrayList<Float> fractionNonConAL = new ArrayList<Float>();
		ArrayList<Float> fdrAL = new ArrayList<Float>();
		ArrayList<Integer> positionsAL = new ArrayList<Integer>();

		//collect all positions
		int[] allPositions = Num.returnUniques(new int[][]{positionsNonCon, positionsCon});

		//for each position 
		int indexNonCon =0;
		int indexCon =0;
		for (int i=0; i< allPositions.length; i++){
			int testPos = allPositions[i];
			float numNonCon =0;
			float numCon =0;
			//present in nonCon?
			for (int j=indexNonCon; j < positionsNonCon.length; j++){
				//match!
				if (testPos == positionsNonCon[j]){
					numNonCon = readsNonCon[j];
					indexNonCon++;
					break;
				}
				//less than
				if (testPos < positionsNonCon[j]) break;
				//greater than so keep advancing
				indexNonCon = j;
			}
			//present in con?
			for (int j=indexCon; j < positionsCon.length; j++){
				//match!
				if (testPos == positionsCon[j]){
					numCon = readsCon[j];
					indexCon++;
					break;
				}
				//less than
				if (testPos < positionsCon[j]) break;
				//greater than so keep advancing
				indexCon = j;
			}


			float totalObservations = numCon+numNonCon;

			//calc fraction non-converted?
			float fnc = (numNonCon+1)/(totalObservations+2);
			float pValue = 0;
			if (numNonCon !=0) {
				//calculate pValue
				int nc = (int)numNonCon;
				int c = (int)numCon;
				//too big? reduce, this is a hack but very few should hit this ceiling
				if (totalObservations > numberCachedBinPVals){
					float multiplier = numberCachedBinPVals/totalObservations;
					nc = (int)(numNonCon * multiplier);
					c = (int)(numCon * multiplier);
				}
				pValue = cachedPValues[nc][c];	
			}

			//fetch fdr
			float fdr = 0;
			if (pValue > minimumPValueToCorrect) fdr = point[fdrIndex++].getScore();

			//watch out for out of bounds sequence due to partial matches to sequence termini
			if (testPos < 2 || testPos > genomicSequenceLengthMinus3) continue;
			//fetch genomic sequence
			String genSeq = genomicSequence.substring(testPos-2,testPos+3);
			if (negativeStrand) genSeq = Seq.reverseComplementDNA(genSeq);

			//fetch BaseContext and increment counters
			BaseContext bc = baseContexts.get(genSeq);
			if (bc != null) {
				bc.incrementNumberConvertedReads((long)numCon);
				bc.incrementNumberNonConvertedReads((long)numNonCon);
				//check fdr (it's in -10Log10 space) and min coverage
				if (fdr >= fdrThreshold) {
					bc.incrementNumberNonConvertedGenomicContexts();
					if (totalObservations >= minimumReadCoverage) {
						bc.addFractionNonConvertedToHistogram(fnc);
					}
				}
				//save for graphing?
				if (printGraphs && totalObservations >= minimumReadCoverage){
					positionsAL.add(new Integer(testPos));
					fractionNonConAL.add(new Float(fnc));
					fdrAL.add(new Float(fdr));
				}
				//TODO: fix! How!?
				//else bc.incrementNumberConvertedGenomicContexts();
			}

		}
		//clean up
		String genomeVersion = con.getInfo().getVersionedGenome();
		int readLength = con.getInfo().getReadLength();
		con = null;
		nonCon = null;
		positionsCon = null;
		positionsNonCon = null;
		readsCon = null;
		readsNonCon = null;

		PointData[] pds = null;
		if (printGraphs){
			//any positions?
			int[] positions = Num.arrayListOfIntegerToInts(positionsAL);
			if (positions == null || positions.length ==0) return null;


			//make PointData for fraction
			HashMap<String,String> map = new HashMap<String,String>();
			map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
			map.put(BarParser.GRAPH_TYPE_COLOR_TAG, "#E8E8E8");
			map.put(BarParser.DESCRIPTION_TAG, "Ratio (non-converted+1)/ (total reads over given C + 2)");
			PointData fractionPD = new PointData();
			Info info = new Info("Fraction non-converted Cs per base", genomeVersion, chromosome, ".", readLength, map);
			if (negativeStrand) info.setStrand("-");
			else info.setStrand("+");
			fractionPD.setInfo(info);
			fractionPD.setPositions(positions);
			float[] fractions = Num.arrayListOfFloatToArray(fractionNonConAL);
			fractionPD.setScores(fractions);

			//make PointData for fdr
			HashMap<String,String> mapFDR = new HashMap<String,String>();
			mapFDR.put(BarParser.GRAPH_TYPE_COLOR_TAG, "#C00000");
			mapFDR.put(BarParser.DESCRIPTION_TAG, "FDRs for methylated C");
			mapFDR.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
			PointData fdrPD = new PointData();
			Info infoFDR = new Info("Binomial p-values converted to FDRs using B & H", genomeVersion, chromosome, ".", readLength, mapFDR);
			if (negativeStrand) infoFDR.setStrand("-");
			else infoFDR.setStrand("+");
			fdrPD.setInfo(infoFDR);
			fdrPD.setPositions(positions);
			float[] fractionsFDR = Num.arrayListOfFloatToArray(fdrAL);
			fdrPD.setScores(fractionsFDR);

			pds = new PointData[]{fractionPD, fdrPD};
		}
		return pds;
	}

	/**Scores a chromosome for non-converted to total at base level. Combines stranded data from CG contexts.
	 * @return PointData[2] fraction and FDRs for a given strand*/
	private PointData[] baseScanChromosomeContextsMergingStrands(){
		//fetch arrays
		int[] positionsNonCon = nonConvertedChrom.getPositions();
		float[] readsNonCon = nonConvertedChrom.getScores();
		int[] positionsCon = convertedChrom.getPositions();
		float[] readsCon = convertedChrom.getScores();
		//make containers for graph data
		ArrayList<Float> fractionNonConAL = new ArrayList<Float>();
		ArrayList<Float> fdrAL = new ArrayList<Float>();
		ArrayList<Integer> positionsAL = new ArrayList<Integer>();

		//collect all positions
		int[] allPositions = Num.returnUniques(new int[][]{positionsNonCon, positionsCon});

		//for each position 
		int indexNonCon =0;
		int indexCon =0;
		for (int i=0; i< allPositions.length; i++){
			int testPos = allPositions[i];
			float numNonCon =0;
			float numCon =0;
			//present in nonCon?
			for (int j=indexNonCon; j < positionsNonCon.length; j++){
				//match!
				if (testPos == positionsNonCon[j]){
					numNonCon = readsNonCon[j];
					indexNonCon++;
					break;
				}
				//less than
				if (testPos < positionsNonCon[j]) break;
				//greater than so keep advancing
				indexNonCon = j;
			}
			//present in con?
			for (int j=indexCon; j < positionsCon.length; j++){
				//match!
				if (testPos == positionsCon[j]){
					numCon = readsCon[j];
					indexCon++;
					break;
				}
				//less than
				if (testPos < positionsCon[j]) break;
				//greater than so keep advancing
				indexCon = j;
			}

			//watch out for out of bounds sequence due to partial matches to sequence termini
			if (testPos < 2 || testPos > genomicSequenceLengthMinus3) continue;
			//fetch genomic sequence
			String genSeq = genomicSequence.substring(testPos-2,testPos+3);

			//is this a CG context?
			if (genomicSequence.subSequence(testPos, testPos+2).equals("CG")){
				testPos++;
				//look to see if any data on neg strand
				boolean incrementI = false;
				if (indexNonCon < positionsNonCon.length){
					int nextNonConPos = positionsNonCon[indexNonCon];
					if (testPos == nextNonConPos){
						numNonCon += readsNonCon[indexNonCon];
						indexNonCon++;
						incrementI = true;
					}
				}
				if (indexCon < positionsCon.length){
					int nextConPos = positionsCon[indexCon];
					if (testPos == nextConPos){
						numCon += readsCon[indexCon];
						indexCon++;
						incrementI = true;
					}
				}
				if (incrementI) i++;
			}

			//invert sequence?
			if (genSeq.charAt(2) == 'G') {
				genSeq = Seq.reverseComplementDNA(genSeq);
			}



			float totalObservations = numCon+numNonCon;

			//calc fraction non-converted?
			float fnc = (numNonCon+1)/(totalObservations+2);

			float pValue = 0;
			if (numNonCon !=0) {
				//calculate pValue
				int nc = (int)numNonCon;
				int c = (int)numCon;
				//too big? reduce, this is a hack but very few should hit this ceiling
				if (totalObservations > numberCachedBinPVals){
					float multiplier = numberCachedBinPVals/totalObservations;
					nc = (int)(numNonCon * multiplier);
					c = (int)(numCon * multiplier);
				}
				pValue = cachedPValues[nc][c];	
			}

			//fetch fdr
			float fdr = 0;
			if (pValue > minimumPValueToCorrect) fdr = point[fdrIndex++].getScore();

			//fetch BaseContext and increment counters
			BaseContext bc = baseContexts.get(genSeq);
			if (bc != null) {
				bc.incrementNumberConvertedReads((long)numCon);
				bc.incrementNumberNonConvertedReads((long)numNonCon);
				//check fdr (it's in -10Log10 space) and min coverage
				if (fdr >= fdrThreshold) {
					bc.incrementNumberNonConvertedGenomicContexts();
					if (totalObservations >= minimumReadCoverage) {
						bc.addFractionNonConvertedToHistogram(fnc);
					}
				}
				//save for graphing?
				if (printGraphs && totalObservations >= minimumReadCoverage){
					positionsAL.add(new Integer(testPos));
					fractionNonConAL.add(new Float(fnc));
					fdrAL.add(new Float(fdr));
				}

			}

		}
		//clean up
		String genomeVersion = convertedChrom.getInfo().getVersionedGenome();
		int readLength = convertedChrom.getInfo().getReadLength();
		positionsCon = null;
		positionsNonCon = null;
		readsCon = null;
		readsNonCon = null;

		PointData[] pds = null;
		if (printGraphs){
			//any positions?
			int[] positions = Num.arrayListOfIntegerToInts(positionsAL);
			if (positions == null || positions.length ==0) return null;

			//make PointData for fraction
			HashMap<String,String> map = new HashMap<String,String>();
			map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
			map.put(BarParser.GRAPH_TYPE_COLOR_TAG, "#E8E8E8");
			map.put(BarParser.DESCRIPTION_TAG, "Ratio (non-converted+1)/ (total reads over given C + 2), stranded counts for CG contexts are merged");
			PointData fractionPD = new PointData();
			Info info = new Info("Fraction non-converted Cs per base", genomeVersion, chromosome, ".", readLength, map);
			info.setStrand(".");
			fractionPD.setInfo(info);
			fractionPD.setPositions(positions);
			float[] fractions = Num.arrayListOfFloatToArray(fractionNonConAL);
			fractionPD.setScores(fractions);

			//make PointData for fdr
			HashMap<String,String> mapFDR = new HashMap<String,String>();
			mapFDR.put(BarParser.GRAPH_TYPE_COLOR_TAG, "#C00000");
			mapFDR.put(BarParser.DESCRIPTION_TAG, "FDRs for methylated C");
			mapFDR.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
			PointData fdrPD = new PointData();
			Info infoFDR = new Info("Binomial p-values converted to FDRs using B & H", genomeVersion, chromosome, ".", readLength, mapFDR);
			infoFDR.setStrand(".");
			fdrPD.setInfo(infoFDR);
			fdrPD.setPositions(positions);
			float[] fractionsFDR = Num.arrayListOfFloatToArray(fdrAL);
			fdrPD.setScores(fractionsFDR);

			pds = new PointData[]{fractionPD, fdrPD};
		}
		return pds;
	}


	private void baseScanChromosomeContextsMergingStrandsForPValues(){
		//fetch arrays
		int[] positionsNonCon = nonConvertedChrom.getPositions();
		float[] readsNonCon = nonConvertedChrom.getScores();
		int[] positionsCon = convertedChrom.getPositions();
		float[] readsCon = convertedChrom.getScores();

		//collect all positions
		int[] allPositions = Num.returnUniques(new int[][]{positionsNonCon, positionsCon});

		//for each position 
		int indexNonCon =0;
		int indexCon =0;
		for (int i=0; i< allPositions.length; i++){

			int testPos = allPositions[i];
			float numNonCon =0;
			float numCon =0;
			//present in nonCon?
			for (int j=indexNonCon; j < positionsNonCon.length; j++){
				//match!
				if (testPos == positionsNonCon[j]){
					numNonCon = readsNonCon[j];
					indexNonCon++;
					break;
				}
				//less than
				if (testPos < positionsNonCon[j]) break;
				//greater than so keep advancing
				indexNonCon = j;
			}
			//present in con?
			for (int j=indexCon; j < positionsCon.length; j++){
				//match!
				if (testPos == positionsCon[j]){
					numCon = readsCon[j];
					indexCon++;
					break;
				}
				//less than
				if (testPos < positionsCon[j]) break;
				//greater than so keep advancing
				indexCon = j;
			}

			//watch out for out of bounds sequence due to partial matches to sequence termini
			if (testPos < 2 || testPos > genomicSequenceLengthMinus3) continue;
			//is this a CG context?
			if (genomicSequence.subSequence(testPos, testPos+2).equals("CG")){
				//shift position one base downstream
				testPos++;
				//look to see if any data on neg strand
				boolean incrementI = false;
				if (indexNonCon < positionsNonCon.length){
					int nextNonConPos = positionsNonCon[indexNonCon];
					if (testPos == nextNonConPos){
						numNonCon += readsNonCon[indexNonCon];
						indexNonCon++;
						incrementI = true;
					}
				}
				if (indexCon < positionsCon.length){
					int nextConPos = positionsCon[indexCon];
					if (testPos == nextConPos){
						numCon += readsCon[indexCon];
						indexCon++;
						incrementI = true;
					}
				}
				if (incrementI) i++;
			}

			float totalObservations = numCon+numNonCon;
			float pValue = 0;
			if (numNonCon !=0) {
				//calculate pValue
				int nc = (int)numNonCon;
				int c = (int)numCon;
				//too big? reduce, this is a hack but very few should hit this ceiling
				if (totalObservations > numberCachedBinPVals){
					float multiplier = numberCachedBinPVals/totalObservations;
					nc = (int)(numNonCon * multiplier);
					c = (int)(numCon * multiplier);
				}
				pValue = cachedPValues[nc][c];	
			}

			//save it or increment the offset
			if (pValue > minimumPValueToCorrect) pointAL.add(new Point(fdrIndex++, pValue));
			else pValueOffset++;
		}
	}


	/**Scores a chromosome for non-converted to total at base level.*/
	private void baseScanChromosomeContextsForPValues(PointData nonCon, PointData con){
		//fetch arrays
		int[] positionsNonCon = nonCon.getPositions();
		float[] readsNonCon = nonCon.getScores();
		int[] positionsCon = con.getPositions();
		float[] readsCon = con.getScores();

		//collect all positions
		int[] allPositions = Num.returnUniques(new int[][]{positionsNonCon, positionsCon});

		//for each position 
		int indexNonCon =0;
		int indexCon =0;
		for (int i=0; i< allPositions.length; i++){
			int testPos = allPositions[i];
			float numNonCon =0;
			float numCon =0;
			//present in nonCon?
			for (int j=indexNonCon; j < positionsNonCon.length; j++){
				//match!
				if (testPos == positionsNonCon[j]){
					numNonCon = readsNonCon[j];
					indexNonCon++;
					break;
				}
				//less than
				if (testPos < positionsNonCon[j]) break;
				//greater than so keep advancing
				indexNonCon = j;
			}
			//present in con?
			for (int j=indexCon; j < positionsCon.length; j++){
				//match!
				if (testPos == positionsCon[j]){
					numCon = readsCon[j];
					indexCon++;
					break;
				}
				//less than
				if (testPos < positionsCon[j]) break;
				//greater than so keep advancing
				indexCon = j;
			}

			//calculate pvalue (in -10Log10Space)
			float totalObservations = numCon+numNonCon;
			Float pValue = zeroFloat;
			if (numNonCon !=0) {
				//calculate pValue
				int nc = (int)numNonCon;
				int c = (int)numCon;
				//too big? reduce, this is a hack but very few should hit this ceiling
				if (totalObservations > numberCachedBinPVals){
					float multiplier = numberCachedBinPVals/totalObservations;
					nc = (int)(numNonCon * multiplier);
					c = (int)(numCon * multiplier);
				}
				pValue = cachedPValues[nc][c];	
			}
			//save it or increment the offset
			if (pValue > minimumPValueToCorrect) pointAL.add(new Point(fdrIndex++, pValue));
			else pValueOffset++;
		}


	}


	public float[] baseScanLambda(PointData nonCon, PointData con){
		if (nonCon == null || con == null) return null;
		//fetch arrays
		int[] positionsNonCon = nonCon.getPositions();
		float[] readsNonCon = nonCon.getScores();
		int[] positionsCon = con.getPositions();
		float[] readsCon = con.getScores();

		float[] pValues = new float[positionsNonCon.length];

		//for each non-converted position 
		int indexCon =0;
		for (int i=0; i< positionsNonCon.length; i++){
			int testPos = positionsNonCon[i];
			float numNonCon =readsNonCon[i];
			float numCon =0;

			//present in con?
			for (int j=indexCon; j < positionsCon.length; j++){
				//match!
				if (testPos == positionsCon[j]){
					numCon = readsCon[j];
					indexCon++;
					break;
				}
				//less than
				if (testPos < positionsCon[j]){
					break;
				}
				//greater than so keep advancing
				indexCon = j;
			}

			//calculate pValue
			int nc = (int)numNonCon;
			int c = (int)numCon;
			float totalObservations = numNonCon+ numCon;
			//too big? reduce, this is a hack but very few should hit this ceiling
			if (totalObservations > numberCachedBinPVals){
				float multiplier = numberCachedBinPVals/totalObservations;
				nc = (int)(numNonCon * multiplier);
				c = (int)(numCon * multiplier);
			}
			float pValue = cachedPValues[nc][c];	
			pValues[i] = pValue;

		}
		return pValues;

	}



	public void setBaseContexts(){
		baseContexts = new HashMap <String,BaseContext>();
		String[] bases = {"G","A","T","C"};
		for (int a = 0; a< 4; a++){
			for (int b=0; b<4; b++){
				for (int c=0; c<4; c++){
					for (int d=0; d<4; d++){
						String bc = bases[a]+bases[b]+"C"+bases[c]+bases[d];
						baseContexts.put(bc, new BaseContext(bc));
					}
				}
			}
		}
	}


	/**Collects and calculates a bunch of stats re the PointData.*/
	private void calculateReadCountStatistics(){
		//fetch converted PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (convertedPointDirs);
		convertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		convertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
		combo = PointData.fetchStrandedPointDataNoMerge (nonConvertedPointDirs);
		nonConvertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		nonConvertedMinusPointData = PointData.convertArrayList2Array(combo[1]);

		//calc totals
		double totalConvertedGenomicContexts = PointData.totalObservationsMultiPointData(convertedPlusPointData);
		totalConvertedGenomicContexts += PointData.totalObservationsMultiPointData(convertedMinusPointData);
		double totalNonConvertedGenomicContexts = PointData.totalObservationsMultiPointData(nonConvertedPlusPointData);
		totalNonConvertedGenomicContexts += PointData.totalObservationsMultiPointData(nonConvertedMinusPointData);
	}

	private int[] fetchMergePositions(){
		ArrayList<int[]> posAL = new ArrayList<int[]>();
		posAL.add(nonConvertedChrom.getPositions());
		posAL.add(convertedChrom.getPositions());
		//merge
		int[][] toMerge = new int[posAL.size()][];
		for (int i=0; i< posAL.size(); i++) toMerge[i] = posAL.get(i);
		return Num.returnUniques(toMerge);
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

	/**Writes stair step and point window bar graph files*/
	public void writeBarFileGraphsOnCurrentSM(){
		//any windows?
		if (smoothingWindow == null || smoothingWindow.length == 0) return;
		//build info object
		Info info;
		if (convertedMergedChromPlus != null) info= convertedMergedChromPlus.getInfo();
		else info = convertedMergedChromMinus.getInfo();
		info = info.copy();
		HashMap<String,String> notes = new HashMap<String,String>();
		notes.put(BarParser.WINDOW_SIZE, windowSize+"");
		notes.put(BarParser.DESCRIPTION_TAG, Misc.stringArrayToString(scoreNames, ","));
		notes.put(BarParser.UNIT_TAG, Misc.stringArrayToString(scoreUnits, ","));
		info.setNotes(notes);
		info.setStrand(".");	

		//save data as heatmap and as points
		////fraction methylation
		saveSmoothedHeatMapData (0, smoothingWindow, info, windowFractionNonConvertedWindow, "#FF0000", false); //red
		////density of 1st quartile fraction methylation
		saveSmoothedHeatMapData (1, smoothingWindow, info, windowDensity1stQuartileWindow, "#0099FF", false); //blue
		////density of 2-3 quartile fraction methylation
		saveSmoothedHeatMapData (2, smoothingWindow, info, windowDensity2_3QuartileWindow, "#33FFFF", false); //
		////density of 4th quartile fraction methylation
		saveSmoothedHeatMapData (3, smoothingWindow, info, windowDensity4thQuartileWindow, "#CC3399", false); //



	}

	/**Saves bar heatmap/ stairstep graph files*/
	public void saveSmoothedHeatMapData (int scoreIndex, SmoothingWindow[] sm, Info info, File dir, String color, boolean posNeg){
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = info.getNotes();		
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_STAIRSTEP);
		//color red
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, color);
		//what's the source
		String fileNames = Misc.stringArrayToString(IO.fetchFileNames(convertedPointDirs),",");
		if (nonConvertedPointDirs !=null) fileNames = fileNames+" vs "+Misc.stringArrayToString(IO.fetchFileNames(nonConvertedPointDirs),",");
		map.put(BarParser.SOURCE_TAG, fileNames);
		//what's window size
		map.put(BarParser.WINDOW_SIZE, windowSize+"");
		//what's the unit on the scores
		map.put(BarParser.UNIT_TAG, scoreUnits[scoreIndex]);
		//description
		map.put(BarParser.DESCRIPTION_TAG, scoreDescriptions[scoreIndex]);
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
	public void saveSmoothedPointData (int scoreIndex, SmoothingWindow[] sm, Info info, File dir, String color){
		//add info to hashmap for writing to bar file
		HashMap<String,String> map = info.getNotes();
		//what graph type should be used to display it?
		map.put(BarParser.GRAPH_TYPE_TAG, BarParser.GRAPH_TYPE_BAR);
		//color
		map.put(BarParser.GRAPH_TYPE_COLOR_TAG, color);
		//what's the source
		String fileNames = Misc.stringArrayToString(IO.fetchFileNames(convertedPointDirs),",");
		if (nonConvertedPointDirs !=null) fileNames = fileNames+" vs "+Misc.stringArrayToString(IO.fetchFileNames(nonConvertedPointDirs),",");
		map.put(BarParser.SOURCE_TAG, fileNames);
		//what's window size
		map.put(BarParser.WINDOW_SIZE, windowSize+"");
		//what's the unit on the scores
		map.put(BarParser.UNIT_TAG, scoreUnits[scoreIndex]);
		//description
		map.put(BarParser.DESCRIPTION_TAG, scoreDescriptions[scoreIndex]);
		//save in info
		info.setNotes(map);
		//make centered position and
		int[] positions = new int[sm.length];
		float[] scores = new float[sm.length];
		for (int i=0; i< sm.length; i++){
			positions[i] = Num.calculateMiddleIntergenicCoordinates(sm[i].getStart(), sm[i].getStop());
			scores[i] = sm[i].getScores()[scoreIndex];
		}
		//write bar
		PointData pd = new PointData();
		pd.setInfo(info);
		pd.setPositions(positions);
		pd.setScores(scores);
		pd.writePointData(dir);
		//clean up
		pd.nullPositionScoreArrays();
		positions = null;
		scores = null;
		sm = null;
	}

	/**Sets score names/ descriptions/ units base on whether nonConverted data is present.*/
	public void setScoreStrings(){
		scoreNames = new String[]{
				"FractionNonConvertedCs",
				"Density1stQuartileFracNonCon",
				"Density2nd3rdQuartilesFracNonCon",
				"Density4thQuartileFracNonCon"

		};
		scoreDescriptions = new String[]{
				"Fraction of non coverted sequenced Cs in the window",
				"Density of base fractions <= 0.25 in the window",
				"Density of base fractions between 0.25 and 0.75 in the window",
				"Density of base fractions >= 0.75 in the window",

		};
		scoreUnits = new String[]{
				"fraction",	
				"density",
				"density",
				"density"
		};
	}

	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BisStat(args);
	}		

	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		File[] fastas = null;
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': fastas = IO.extractFiles(new File(args[++i])); break;
					case 'c': convertedPointDirs = IO.extractFiles(args[++i]); break;
					case 'n': nonConvertedPointDirs = IO.extractFiles(args[++i]); break;
					case 's': saveDirectory = new File(args[++i]); break;
					case 'w': windowSize = Integer.parseInt(args[++i]); break;
					case 'm': minimumCsInWindow = Integer.parseInt(args[++i]); break;
					case 'o': minimumReadCoverage = Float.parseFloat(args[++i]); break;
					case 'p': fdrThreshold = Float.parseFloat(args[++i]); break;
					case 'e': expectedNonConverted = Float.parseFloat(args[++i]); useLambda = false; break;
					case 'l': useLambda = true; break;
					case 'd': performStrandedAnalysis = false; break;
					case 'r': fullPathToR = new File(args[++i]); break;
					case 'h': printDocs(); System.exit(0);
					default: Misc.printExit("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		//look for fasta files
		if (fastas == null || fastas.length ==0) Misc.printErrAndExit("\nError: cannot find any fasta sequence files?\n");

		//look for point directories
		if (convertedPointDirs == null || convertedPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your converted PointData directories(s)!\n");
		//only one directory look deeper
		if (convertedPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(convertedPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) convertedPointDirs = otherDirs;
		}

		//nonConverted data
		if (nonConvertedPointDirs[0].isDirectory() == false) Misc.printExit("\nError: cannot find your nonConverted PointData directories(s)!\n");
		//only one directory look deeper
		if (nonConvertedPointDirs.length == 1){
			File[] otherDirs = IO.extractOnlyDirectories(nonConvertedPointDirs[0]);
			if (otherDirs != null && otherDirs.length > 0) nonConvertedPointDirs = otherDirs;
		}

		//set half peak shift and windowSize
		if (windowSize == 0 ) Misc.printExit("\nPlease enter a positive length for the window size.\n");

		//make window maker 
		windowMaker = new WindowMaker(windowSize,minimumCsInWindow);

		//look for and or create the save directory
		if (saveDirectory == null) Misc.printExit("\nError: enter a directory text to save results.\n");
		else if (saveDirectory.exists() == false) saveDirectory.mkdirs();


		if (printGraphs){
			//make dirs for base level fraction non-converted point data and their associated qvalues
			File baseDirectory = new File (saveDirectory, "Base");
			baseDirectory.mkdirs();
			fractionNonConvertedDirectory = new File (baseDirectory, "BaseFractionNonConverted");
			fractionNonConvertedDirectory.mkdirs();
			if (performStrandedAnalysis){
				fractionNonConvertedDirectoryPlus = new File (fractionNonConvertedDirectory, "Plus");
				fractionNonConvertedDirectoryPlus.mkdirs();
				fractionNonConvertedDirectoryMinus = new File (fractionNonConvertedDirectory, "Minus");
				fractionNonConvertedDirectoryMinus.mkdirs();
			}
			fdrNonConvertedDirectory = new File (baseDirectory, "BaseFDRsNonConverted");
			fdrNonConvertedDirectory.mkdirs();


			//make dir for window level fraction methylated point data
			File windowDirectory = new File (saveDirectory, "Window");
			windowDirectory.mkdirs();
			windowFractionNonConvertedWindow = new File(windowDirectory, "WindowFractionNonConverted");
			windowFractionNonConvertedWindow.mkdirs();
			//windowFractionNonConvertedPoint = new File(saveDirectory, "WindowFractionNonConverted_Point");
			//windowFractionNonConvertedPoint.mkdirs();

			//make dirs for density plots
			windowDensity1stQuartileWindow = new File(windowDirectory, "WindowDensity1stQuartileFracNonCon");
			windowDensity2_3QuartileWindow = new File(windowDirectory, "WindowDensity2nd3rdQuartilesFracNonCon");
			windowDensity4thQuartileWindow = new File(windowDirectory, "WindowDensity4thQuartileFracNonCon");
			windowDensity1stQuartileWindow.mkdirs();
			windowDensity2_3QuartileWindow.mkdirs();
			windowDensity4thQuartileWindow.mkdirs();

			serializedWindowObjectsDirectory = new File(saveDirectory, "SerializedWindowObjects");
			serializedWindowObjectsDirectory.mkdirs();
		}

		//load fasta files into hash
		chromosomeFastaFiles = new HashMap<String,File>();
		Pattern chrom = Pattern.compile("(.+)\\.fa.*");
		for (int i=0; i< fastas.length; i++){
			Matcher mat = chrom.matcher(fastas[i].getName());
			if (mat.matches()) chromosomeFastaFiles.put(mat.group(1), fastas[i]);
		}

		//set score items
		setScoreStrings();

		//check for R and required libraries
		if (fullPathToR == null || fullPathToR.canExecute()== false) {
			Misc.printExit("\nError: Cannot find or execute the R application -> "+fullPathToR+"\n");
		}

		//set minimum pval to minimum FDR
		minimumPValueToCorrect = fdrThreshold;

		System.out.println("\nParams:");
		System.out.println(minimumReadCoverage+"\tMinimum read coverage to call a base fraction methylation.");
		System.out.println(windowSize+"\tWindow size");
		System.out.println(minimumCsInWindow+"\tMinimum number Cs passing read coverage in window");
		System.out.println(firstQuartileThreshold + "\tFirst quartile window density threshold");
		System.out.println(fourthQuartileThreshold + "\tFourth quartile window density threshold");
		if (useLambda == false) System.out.println(expectedNonConverted+"\tExpected Non Converted");
		System.out.println(fdrThreshold+"\tFDR threshold");
		System.out.println(performStrandedAnalysis+"\tPerform standed analysis");
	}	

	/**Prints stats for various contexts in genomic space.*/
	public void printGenomicContextStats(){
		try {
			Pattern mCHG_Pattern = Pattern.compile("..C[ATC]G");
			Pattern mCHH_Pattern = Pattern.compile("..C[ATC][ATC]");
			//counts
			double mCG_Count = 0;
			double mCHG_Count = 0;
			double mCHH_Count = 0;
			double seqCs = 0;
			double seqCG = 0;
			double seqCHH = 0;
			double seqCHG = 0;
			double seqMethylCs = 0;
			double seqMethylCG = 0;
			double seqMethylCHH = 0;
			double seqMethylCHG = 0;

			//histograms
			Histogram mCGHistogram = null;
			Histogram mCHGHistogram = null;
			Histogram mCHHHistogram = null;
			Iterator<BaseContext> it = baseContexts.values().iterator();
			while (it.hasNext()){
				//increment counters
				BaseContext b = it.next();
				//read counts for c and mC
				seqCs += b.getNumberConvertedReads();
				seqMethylCs += b.getNumberNonConvertedReads();

				if (CG_Pattern.matcher(b.getSequence()).matches()){
					mCG_Count += b.getNumberNonConvertedGenomicContexts();
					if (mCGHistogram == null) mCGHistogram = b.getHistogram();
					else mCGHistogram.addCounts(b.getHistogram());
					//read counts
					seqCG += b.getNumberConvertedReads() ;
					seqMethylCG += b.getNumberNonConvertedReads();
				}
				else if (mCHG_Pattern.matcher(b.getSequence()).matches()){
					mCHG_Count += b.getNumberNonConvertedGenomicContexts();
					if (mCHGHistogram == null) mCHGHistogram = b.getHistogram();
					else mCHGHistogram.addCounts(b.getHistogram());
					//read counts
					seqCHG += b.getNumberConvertedReads();
					seqMethylCHG += b.getNumberNonConvertedReads();
				}
				else if (mCHH_Pattern.matcher(b.getSequence()).matches()){
					mCHH_Count += b.getNumberNonConvertedGenomicContexts();
					if (mCHHHistogram == null) mCHHHistogram = b.getHistogram();
					else mCHHHistogram.addCounts(b.getHistogram());
					//read counts
					seqCHH += b.getNumberConvertedReads();
					seqMethylCHH += b.getNumberNonConvertedReads();
				}
				else System.out.println("Warning, not counting "+b.getSequence());
			}

			//histograms
			System.out.println("\nDistribution of the fraction mC for composite genomic contexts meeting a minimum FDR threshold of "+fdrThreshold+" and minimum read coverage of "+minimumReadCoverage);

			if (mCGHistogram.getTotalBinCounts() !=0){
				System.out.println("mCG");
				mCGHistogram.printScaledHistogram();
			}
			if (mCHGHistogram.getTotalBinCounts() != 0){
				System.out.println("mCHG");
				mCHGHistogram.printScaledHistogram();
			}
			if (mCHHHistogram.getTotalBinCounts() !=0){
				System.out.println("mCHH");
				mCHHHistogram.printScaledHistogram();
			}


			double mC_Count = mCG_Count + mCHG_Count + mCHH_Count;

			System.out.println("\nStats based on aligned genomic contexts that meet a minimum FDR threshold of "+fdrThreshold+
			". WARNING: datasets must be subsampled to the same bp aligned for these stats to be cross dataset comparable.");

			printStatLine(mCG_Count, mC_Count, "mCG/mC");
			printStatLine(mCHG_Count, mC_Count, "mCHG/mC");
			printStatLine(mCHH_Count, mC_Count, "mCHH/mC");

			System.out.println("\nStats based on cumulative sums of all read sequences, no thresholds:");

			printStatLine(seqMethylCG, seqMethylCs, "mCG/mC");
			printStatLine(seqMethylCHG, seqMethylCs, "mCHG/mC");
			printStatLine(seqMethylCHH, seqMethylCs, "mCHH/mC");

			System.out.println();

			printStatLine(seqMethylCs, seqCs, "mC/C");
			printStatLine(seqMethylCG, seqCG, "mCG/CG");
			printStatLine(seqMethylCHG, seqCHG, "mCHG/CHG");
			printStatLine(seqMethylCHH, seqCHH, "mCHH/CHH");

			System.out.println();

			printStatLine(seqMethylCs, seqCs+seqMethylCs, "mC/(C+mC)");
			printStatLine(seqMethylCG, seqCG+seqMethylCG, "mCG/(CG+mCG)");
			printStatLine(seqMethylCHG, seqCHG+seqMethylCHG, "mCHG/(CHG+mCHG)");
			printStatLine(seqMethylCHH, seqCHH+seqMethylCHH, "mCHH/(CHH+mCHH)");

			float[] f = Num.collapseFloatArray(baseFractionMethylation);
			if (f != null && f.length > 2){
				System.out.println("\nStats based on the individual per base fraction methylations that pass the minimum read coverage, NO FDR threshold:");
				Arrays.sort(f);
				System.out.println("\tMean\tMedian\tStdDev\tMin\tMax\t10th\t90th");
				System.out.println("\t"+Num.statFloatArray(f));
				Histogram histogram = new Histogram(0, 1.1, 11);
				histogram.countAll(f);
				histogram.printScaledHistogram();
			}

		} catch (Exception e){
			e.printStackTrace();
		}

	}

	public static void printStatLine(double numerator, double denomenator, String name){
		double fraction = numerator/ denomenator;
		if (denomenator == 0) fraction = 0;
		System.out.println( "\t"+ Num.formatNumber(fraction, 3)+"\t("+(long)numerator+"/"+(long)denomenator+")\t"+name);
	}

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                  BisStat: August 2012                            **\n" +
				"**************************************************************************************\n" +
				"Takes PointData from converted and non-converted C bisulfite sequencing data parsed\n" +
				"using the NovoalignBisulfiteParser and generates several xxCxx context statistics and\n" +
				"graphs (bp and window level fraction converted Cs) for visualization in IGB.\n" +
				"BisStat estimates whether a given C is methylated using a binomial distribution where\n" +
				"the expect can be calculated using the fraction of non-converted Cs present in the\n" +
				"lambda data. Binomial p-values are converted to FDRs using the Benjamini & Hochberg\n" +
				"method. This app requires considerable RAM (10-64G).\n\n" +

				"Options:\n"+
				"-s Save directory, full path.\n"+
				"-c Converted PointData directories, full path, comma delimited. These should\n" +
				"       contain stranded chromosome specific xxx_-/+_.bar.zip files. One\n" +
				"       can also provide a single directory that contains multiple PointData\n" +
				"       directories.\n" +
				"-n Non-converted PointData directories, ditto. \n" +
				"-f Directory containing chrXXX.fasta(/.fa .zip/.gz OK) files for each chromosome.\n"+
				"\n"+
				"Default Options:\n"+
				"-p Minimimal FDR for non-converted C's to be counted as methylated, defaults to 20 a\n" +
				"       -10Log10(FDR = 0.01) conversion.\n"+
				"-e Expected fraction non-converted Cs due to partial bisulfite conversion and\n" +
				"       sequencing error, defaults to 0.005 .\n"+
				"-l Use the unmethylated lambda alignment data to set the expected fraction of\n" +
				"       non-converted Cs due to partial conversion and sequencing error. This is\n" +
				"       predicated on including a 'chrLambda' fasta sequence while aligning your data.\n"+
				"-o Minimum read coverage to count mC fractions, defaults to 8\n"+
				"-w Window size, defaults to 1000.\n"+
				"-m Minimum number Cs passing read coverage in window to score, defaults to 5. \n"+
				"-r Full path to R, defaults to '/usr/bin/R'\n" +
				"-d Merge stranded data, defaults to running a stranded analysis. Affects CG's.\n"+
				"-a First density quartile fraction methylation threshold, defaults to 0.25\n"+
				"-b Fourth density quartile fraction methylation threshold, defaults to 0.75\n"+


				"\n"+

				"Example: java -Xmx12G -jar pathTo/USeq/Apps/BisStat -c /Data/Sperm/Converted -n \n" +
				"      /Data/Sperm/NonConverted -s /Data/Sperm/BisSeq -w 5000 -m 10 -f\n" +
				"      /Genomes/Hg18/Fastas -o 10 \n\n" +

		"**************************************************************************************\n");

	}
}
