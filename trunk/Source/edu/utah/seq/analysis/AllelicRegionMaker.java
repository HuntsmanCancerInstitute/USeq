package edu.utah.seq.analysis;

import java.io.*;
import java.util.regex.*;
import java.util.*;
import util.gen.*;
import edu.utah.seq.data.*;
import edu.utah.seq.useq.data.RegionScoreText;

/**Indentifies regions with a given methylation profile.
 * @author Nix
 * */
public class AllelicRegionMaker {

	//user defined fields
	private File[] convertedPointDirs;
	private File[] nonConvertedPointDirs;
	private float minimumReadCoverage = 6;
	private float maximumReadCoverage = 10000;
	private int minimumCsInImprintBlock = 6;
	private int minimumLengthImprintBlock = 100;
	//private int numberBadCsInImprintBlock = 3;
	private float minimumFractionMethylation = 0.3f;
	private float maximumFractionMethylation = 0.7f;
	private int maxRunningBadBases =2;
	private float maximumFractionBadBases = 0.2f;

	//internal fields
	private String[] chromosomes;
	private HashMap<String,PointData[]> convertedPlusPointData;
	private HashMap<String,PointData[]> convertedMinusPointData;
	private HashMap<String,PointData[]> nonConvertedPlusPointData;
	private HashMap<String,PointData[]> nonConvertedMinusPointData;

	//by chromosome
	private String chromosome;
	private PointData convertedMergedChromPlus = null;
	private PointData convertedMergedChromMinus = null;	
	private PointData nonConvertedMergedChromPlus = null;
	private PointData nonConvertedMergedChromMinus = null;	
	private PointData convertedChrom = null;
	private PointData nonConvertedChrom = null;

	//constructor
	public AllelicRegionMaker(File[] convertedPointDirs, File[] nonConvertedPointDirs){
		this.convertedPointDirs = convertedPointDirs;
		this.nonConvertedPointDirs = nonConvertedPointDirs;
		//fetch counts
		fetchDataLinks();
		//this also removes the phiX and lambda data from that to be scanned if resent
		fetchAllChromosomes();
	}

	//methods

	public RegionScoreText[] scan(String chromosome){
			this.chromosome = chromosome;

			//fetch data
			fetchDataAndRemove();

			//are all four datasets present?
			if( convertedMergedChromPlus != null && nonConvertedMergedChromPlus != null && convertedMergedChromMinus != null &&  nonConvertedMergedChromMinus != null) {
				//merge strands
				convertedChrom = PointData.mergePairedPointDataNoSumming(convertedMergedChromPlus, convertedMergedChromMinus);
				nonConvertedChrom = PointData.mergePairedPointDataNoSumming(nonConvertedMergedChromPlus, nonConvertedMergedChromMinus);
				//create regions
				return scanChromForRegions2();

			}
			else return null;
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
		String lambdaChromosome = null;
		String phiXChromosome = null;
		while (it.hasNext()){
			String chr = it.next();
			if (lambda.matcher(chr).matches()) lambdaChromosome = chr;
			else if (phiX.matcher(chr).matches()) phiXChromosome = chr;
		}
		if (lambdaChromosome != null) c.remove(lambdaChromosome);
		if (phiXChromosome != null) c.remove(phiXChromosome);

		chromosomes=  Misc.hashSetToStringArray(c);
	}

	/**Fetchs the data for a particular chromosome.*/
	public void fetchDataAndRemove(){
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
	}

	public static int[] fetchUniquePositions (Point[][] pts){
		ArrayList<int[]> intsAL = new ArrayList<int[]>();
		for (int i=0; i< pts.length; i++){
			if (pts[i] != null) intsAL.add(Point.extractPositions(pts[i]));
		}
		if (intsAL.size()==0) return null;
		return Num.returnUniques(intsAL);
	}
	
	private class CBase {
		//fields
		float numberCon;
		float numberNonCon;
		int position;
		boolean goodBase;
		
		public CBase (float numberCon, float numberNonCon, int position, boolean goodBase){
			this.numberCon = numberCon;
			this.numberNonCon = numberNonCon;
			this.position = position;
			this.goodBase = goodBase;
		}
		
	}
	


	private class ImprintBlock {
		
		//fields
		int start;
		int end;
		int numberObs = 1;
		float numberNonCon;
		float numberCon;
		int numberBadCs = 0;

		//constructors
		public ImprintBlock (int start, float numberNonCon, float numberCon){
			this.start = start;
			this.end = start;
			this.numberNonCon = numberNonCon;
			this.numberCon = numberCon;
		}
		
		public ImprintBlock (ArrayList<CBase> bases){
			//remove trailing bad bases
			while (true){
				int lastIndex = bases.size()-1;
				if (bases.get(lastIndex).goodBase == false) bases.remove(lastIndex);
				else break;
			}
			
			for (CBase b: bases){
				numberObs++;
				numberNonCon += b.numberNonCon;
				numberCon += b.numberCon;
				if (b.goodBase == false) numberBadCs++;
			}
			start = bases.get(0).position;
			end = bases.get(bases.size()-1).position;
		}
		

		/**Returns null if thresholds aren't met*/
		public RegionScoreText fetchRegion(){
			//enough obs?
			if (numberObs < minimumCsInImprintBlock) return null;
			//long enough?
			if ((end - start) < minimumLengthImprintBlock) return null;
			//not too many bad bases?
			float fractionBadBases = (float)numberBadCs/(float)numberObs;
			if (fractionBadBases > maximumFractionBadBases) {
				//System.out.println("Failing fractionBadBases ");
				//print();
				return null;
			}
			//calc frac
			numberNonCon++;
			numberCon++;
			float fractionMeth = numberNonCon/(numberNonCon+numberCon);
			//chrom, start, stop, name, score, strand
			return new RegionScoreText(start, end+1, fractionMeth, numberObs+"_"+numberBadCs+"_"+(int)numberNonCon+"_"+(int)numberCon);
		}
		
		public void print(){
			System.out.println("\tstart " +start);
			System.out.println("\tend " +end);
			System.out.println("\tnumberObs " +numberObs);
			System.out.println("\tnumberBadCs " +numberBadCs);
			float fractBadBases = (float)numberBadCs/(float)numberObs;
			System.out.println("\tfractBadBases " +fractBadBases);
		}
	}
	

	/**Scores a chromosome for non-converted to total at base level.
	 * @return PointData[2] fraction and FDRs for a given strand*/
	private RegionScoreText[] scanChromForRegions2(){
		ArrayList<RegionScoreText> regionsAL = new ArrayList<RegionScoreText>();
		
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
		ArrayList<CBase> bases = new ArrayList<CBase>();
		int runningBadBases = 0;
		float numberBadBases = 0;
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
			float fnc = (numNonCon+1)/(totalObservations+2);
			
			//pass thresholds?
			boolean badBase = (totalObservations > maximumReadCoverage || totalObservations < minimumReadCoverage || fnc < minimumFractionMethylation || fnc > maximumFractionMethylation) ;		
			
			//good base, reset running bad
			if (badBase == false){
				bases.add(new CBase(numCon, numNonCon, testPos, true));
				runningBadBases = 0;
			}
			//bad base, only add if pre existing good base
			else if (bases.size() !=0) {
				runningBadBases++;
				numberBadBases++;
				//too many bad bases? 
				float fractionBadBases = numberBadBases/(float)bases.size();
				
				if (totalObservations < maximumReadCoverage && runningBadBases <= maxRunningBadBases && fractionBadBases <= maximumFractionBadBases){
					//ok, add bad base
					bases.add(new CBase(numCon, numNonCon, testPos, false));
				}
				//yes too many
				else {
					//attempt to fetch region
					int numObs = bases.size() - (runningBadBases -1);
					if (numObs >= minimumCsInImprintBlock){
						ImprintBlock im = new ImprintBlock (bases);
						RegionScoreText region = im.fetchRegion();
						if (region!= null) regionsAL.add(region);
					}
					runningBadBases = 0;
					numberBadBases = 0;
					bases.clear();
				}
			}
		}

		//close last block
		int numObs = bases.size() - (runningBadBases -1);
		if (numObs >= minimumCsInImprintBlock){
			ImprintBlock im = new ImprintBlock (bases);
			RegionScoreText region = im.fetchRegion();
			if (region!= null) regionsAL.add(region);
		}
		
		if (regionsAL.size() !=0) {
			RegionScoreText[] r = new RegionScoreText[regionsAL.size()];
			regionsAL.toArray(r);
			return r;
		}
		return null;
	}

	/**Scores a chromosome for non-converted to total at base level.
	 * @return PointData[2] fraction and FDRs for a given strand
	private RegionScoreText[] scanChromForRegions(){
		ArrayList<RegionScoreText> regionsAL = new ArrayList<RegionScoreText>();
		
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
		ImprintBlock im = null;
		int badBases = 0;
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
			float fnc = (numNonCon+1)/(totalObservations+2);
			

			//pass thresholds?
			boolean badBase = (totalObservations < minimumReadCoverage || fnc < minimumFractionMethylation || fnc > maximumFractionMethylation) ;
			if (badBase){
				//close it?
				if (im != null){
					if (im.numberBadCs < numberBadCsInImprintBlock) {
						im.numberBadCs++;
						badBases++;
					}
					else {
						im.numberBadCs -= badBases;
						RegionScoreText region = im.fetchRegion();
						if (region!= null) regionsAL.add(region);
						im = null;
						badBases = 0;
					}
				}
			}
			else {
				if (im == null) im = new ImprintBlock (testPos, numNonCon, numCon);
				else {
					im.numberObs++;
					im.numberCon += numCon;
					im.numberNonCon += numNonCon;
					im.end = testPos;
					badBases = 0;
				}
			}
		}

		//close last
		if (im != null){
			RegionScoreText region = im.fetchRegion();
			if (region!= null) regionsAL.add(region);
		}
		
		if (regionsAL.size() !=0) {
			RegionScoreText[] r = new RegionScoreText[regionsAL.size()];
			regionsAL.toArray(r);
			return r;
		}
		return null;
	}*/
	

	/**Collects and calculates a bunch of stats re the PointData.*/
	private void fetchDataLinks(){
		//fetch converted PointData and calculate total observations
		HashMap<String, ArrayList<PointData>>[] combo = PointData.fetchStrandedPointDataNoMerge (convertedPointDirs);
		convertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		convertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
		combo = PointData.fetchStrandedPointDataNoMerge (nonConvertedPointDirs);
		nonConvertedPlusPointData = PointData.convertArrayList2Array(combo[0]);
		nonConvertedMinusPointData = PointData.convertArrayList2Array(combo[1]);
	}

	public String[] getChromosomes() {
		return chromosomes;
	}
		
}
