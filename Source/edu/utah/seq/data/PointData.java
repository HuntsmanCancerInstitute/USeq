package edu.utah.seq.data;
import java.util.*;
import java.io.*;

import util.gen.*;
import edu.utah.seq.parsers.BarParser;


/**Container for holding genomic positions - scores data.
 * This data is split by versionedGenome, chromosome, and in some cases strand.
 * Methods should be created that return a HashMap<String-Chromosome_Strand,PointData> (ie chr1_+_ , PointData)
 * from various data sources.
 * Works in tandem with the BarParser to load data on demand.
 * The number of positions and scores are the same with direct index correspondence.
 * Positions can be repeated.
 * Final positions are in interbase coordinates.*/
public class PointData {

	//fields
	private Info info;
	private int[] positions;
	private	float[] scores;
	private BarParser barParser;

	//constructors
	public PointData (){
		info = new Info();
	}

	/**Reads in a barFile creating a PointData object.
	 * @param loadPositionScores - set to true to actually load the int[] and float[]s at the time of
	 * instantiation.  Otherwise they will load with the first call for either of the arrays.
	 * Call nullPositionScoreArrays() and System.gc() when done to save memory.*/
	public PointData (File barFile, boolean loadPositionScores){
		//make new BarParser
		barParser = new BarParser();
		//load header info but not necessarily the position scores
		barParser.readBarFile(barFile,loadPositionScores);
		info = new Info (barParser);
	}

	//methods
	/**Writes a zip compressed 'chromosome_strand_.bar.zip' file into the saveDirectory.
	 * ie. chr21_+_.bar.zip
	 * If strand is "." writes just chr21.bar.zip
	 * Will over write an existing file in the saveDirectory.*/
	public boolean writePointData(File saveDirectory) {
		//make text of file and check if it exists
		String fileName = info.getChromosome();
		if (info.getStrand().equals(".") == false) fileName = fileName+"_"+info.getStrand()+"_";
		fileName = fileName + ".bar";
		File barFile = new File (saveDirectory, fileName);
		if (barFile.exists()) System.out.println("\nWarning, over writing "+barFile);
		return writePointDataToFile(barFile);
	}

	/**Writes a zip compressed 'chromosome_strand_.bar.zip' file.
	 * Will over write an existing file in the saveDirectory.*/
	public boolean writePointDataToFile(File barFile) {
		//load data if needed
		this.getPositions();
		//make barParser?
		if (barParser == null) barParser = new BarParser();
		barParser.setZipCompress(true);
		//add notes to tagValues?
		HashMap<String,String> tagValues = barParser.getTagValues();
		if (info.getNotes() != null && info.getNotes().size() !=0) tagValues.putAll(info.getNotes());
		//set readLength in tag values, the strand info will be added by the BarParser
		if (info.getReadLength() !=0) tagValues.put(BarParser.READ_LENGTH_TAG, ""+info.getReadLength());
		//write bar	
		return barParser.writeBarFile(barFile, info.getChromosome(), info.getVersionedGenome(), 
				info.getStrand().charAt(0), positions, scores, tagValues);
	}

	/**Writes zip compressed 'chromosome_strand_.bar.zip' files into the saveDirectory.
	 * ie. chr21_+_.bar.zip
	 * If strand is "." writes just chr21.bar.zip
	 * Will over write existing files!*/
	public static boolean writePointData(HashMap<String,PointData> pd, File saveDirectory){
		Iterator<String> it = pd.keySet().iterator();
		while (it.hasNext()){
			boolean saved = pd.get(it.next()).writePointData(saveDirectory);
			if (saved == false) {
				return false;
			}
		}
		return true;
	}

	/*Will efficiently merge same chromosome PointData with minimal memory usage. Set convertScoresToHitCount to true to generate hit counts otherwise the scores are summed.
	 * If mixed strand is provided, the returned PointData's strand will be set to ".". Returns null if something bad happended.
	 * This is very slow!  Causing memory problems with the VM!
	public static PointData efficientMergeSummingScores(ArrayList<PointData> sameChromPointData, boolean convertScoresToHitCount){
		ArrayList<Point> ptsAL = null;
		boolean sameStrand = true;
		//assume that PointData are ordered as they should be 
		try{
			//load first file
			Info info = sameChromPointData.get(0).getInfo();
			BarParser ori = sameChromPointData.get(0).getBarParser();
			ori.setMergeIdenticalPositions(true);
			ori.setLoadPositionValues(true);
			ori.setReplaceWithHitCount(convertScoresToHitCount);
			ori.loadSimpleBarFile();			
			//convert to Point AL
			ptsAL = Point.makePointsAL(ori.getBasePositions(), ori.getValues());

			//for each BarParser
			for (int i=1; i< sameChromPointData.size(); i++) {
				//check strand
				if (info.getStrand().equals(sameChromPointData.get(i).getInfo().getStrand()) == false) sameStrand = false;
				BarParser bp = sameChromPointData.get(i).getBarParser();

				bp.setLoadPositionValues(false);
				bp.setCloseDataInputStream(false);
				bp.loadSimpleBarFile();

				Point testPoint;
				int index =0;
				while ((testPoint = bp.fetchNextPoint())!= null){
					int testPos = testPoint.getPosition();
					if (convertScoresToHitCount) testPoint.setScore(1);
					//find position in ArrayList of Points
					for (int x = index; x<ptsAL.size(); x++){
						Point alPoint = ptsAL.get(x);
						int alPos = alPoint.getPosition();						
						//same?
						if (alPos == testPos) {
							alPoint.incrementScore(testPoint.getScore());
							index = x;							
							break;
						}
						//test pos less than alPos
						else if (alPos > testPos){
							//insert testPoint before alPoint
							ptsAL.add(x, testPoint);
							index = x;
							//if (index != 0) index = x-1;							
							break;
						}
						else if (ptsAL.size()-1 == x){
							ptsAL.add(testPoint);
							Point currPoint;
							while ((currPoint = bp.fetchNextPoint()) != null) {
								if (convertScoresToHitCount) currPoint.setScore(1);
								if (testPoint.getPosition() == currPoint.getPosition()) testPoint.incrementScore(currPoint.getScore());
								else ptsAL.add(currPoint);
								testPoint = currPoint;
							}
							break;
						}
						//test pos must be greater than so do nothing
					}
				}
				//close data stream
				bp.getDis().close();
			}

			//build PointData
			Point[] finalP = new Point[ptsAL.size()];
			ptsAL.toArray(finalP);
			PointData pd = Point.extractPositionScores(finalP);
			info.setNumberObservations(finalP.length);
			info.setScoreTotal(Num.sumArray(pd.getScores()));
			info.getNotes().clear();
			if (sameStrand == false) info.setStrand(".");
			pd.setInfo(info);

			return pd;
		} catch (IOException e){
			e.printStackTrace();
			return null;
		}
	}*/


	/**Takes a boolean[] of masked/true positions and for every position, calculates the overlap of the read.
	 * If the overlap exceeds the max then it is tossed.  Leaves masked data in memory. If the read length is zero 
	 * then it just looks to see if it is in a masked base.*/
	public int filter(boolean[] maskedBases, double maximumFractionOverlap){
		//load data if not loaded
		getPositions();

		//vars
		ArrayList<Integer> goodPos = new ArrayList<Integer>(positions.length);
		ArrayList<Float> goodScr = new ArrayList<Float>(positions.length);
		double readLength = info.getReadLength();	
		int halfSize = (int)Math.round(readLength/2.0);
		int halfPlusOne = halfSize++;

		//zero read length? bisulfite data?
		if ((int)readLength == 0){
			//for each position, calculate fraction masked/true
			for (int i=0; i< positions.length; i++){
				int start = positions[i];
				//past the masking? or it isn't masked
				if (start >= maskedBases.length || maskedBases[start] == false) {
					goodPos.add(positions[i]);
					goodScr.add(scores[i]);
				}
			}
		}
		else {
			//for each position, calculate fraction masked/true
			for (int i=0; i< positions.length; i++){
				int start = positions[i]- halfSize;
				if (start < 0) start = 0;
				//past the masking?
				if (start > maskedBases.length) {
					goodPos.add(positions[i]);
					goodScr.add(scores[i]);
				}
				else {
					int stop = positions[i]+ halfPlusOne;
					if (stop > maskedBases.length) stop = maskedBases.length;
					double hits = 0;
					for (int j=start; j< stop; j++){
						if (maskedBases[j]) hits++;
					}
					double overlap = hits/readLength;
					if (overlap < maximumFractionOverlap) {
						goodPos.add(positions[i]);
						goodScr.add(scores[i]);
					}
				}
			}
		}
		//convert?
		int diff = positions.length - goodPos.size();
		if (diff == 0) return 0;
		positions = Num.arrayListOfIntegerToInts(goodPos);
		scores = Num.arrayListOfFloatToArray(goodScr);
		info.setNumberObservations(positions.length);
		info.setScoreTotal(Num.sumArrayReturnDouble(scores));
		return diff;
	}

	/**Takes a boolean[] of masked/true positions and for every position, calculates the overlap of the read.
	 * If the overlap >= the min then it is saved.  Leaves masked data in memory.*/
	public int fetch(boolean[] maskedBases, double minimumFractionOverlap){
		//load data if not loaded
		getPositions();

		//vars
		ArrayList<Integer> goodPos = new ArrayList<Integer>(positions.length);
		ArrayList<Float> goodScr = new ArrayList<Float>(positions.length);
		double readLength = info.getReadLength();
		int halfSize = (int)Math.round(readLength/2.0);
		int halfPlusOne = halfSize++;

		//zero read length? bisulfite data?
		if ((int) readLength == 0){
			//for each position, calculate fraction masked/true
			for (int i=0; i< positions.length; i++){
				//past the mask?
				if (positions[i]>= maskedBases.length) break;
				if (maskedBases[positions[i]] == true){
					goodPos.add(positions[i]);
					goodScr.add(scores[i]);
				}
			}
		}
		else {
			//for each position, calculate fraction masked/true
			for (int i=0; i< positions.length; i++){
				int start = positions[i]- halfSize;
				if (start < 0) start = 0;
				int stop = positions[i]+ halfPlusOne;
				if (stop > maskedBases.length) stop = maskedBases.length;
				double hits = 0;
				for (int j=start; j< stop; j++){
					if (maskedBases[j]) hits++;
				}
				double overlap = hits/readLength;
				if (overlap >= minimumFractionOverlap) {
					goodPos.add(positions[i]);
					goodScr.add(scores[i]);
				}
			}
		}
		//convert?
		int diff = positions.length - goodPos.size();
		if (diff == 0) return 0;
		positions = Num.arrayListOfIntegerToInts(goodPos);
		scores = Num.arrayListOfFloatToArray(goodScr);
		info.setNumberObservations(positions.length);
		info.setScoreTotal(Num.sumArrayReturnDouble(scores));
		return diff;
	}

	/**Gets the BarParser to (re)load the position and scores.*/
	public boolean loadPositionScores(){
		//load actual values
		if (barParser.readBarFile(barParser.getBarFile(),true) == false){
			System.out.println("ERROR: loading bar file position and values -> "+barParser.getBarFile());
			return false;
		}
		positions = barParser.getBasePositions();
		scores = barParser.getValues();
		return true;
	}

	/**Gets the BarParser to (re)load the position and scores.*/
	public boolean loadFlattenedPositionScores(boolean firstReplaceScoresWithOne){
		if (barParser == null) return false;
		//load actual values
		barParser.setMergeIdenticalPositions(true);
		barParser.setReplaceWithHitCount(firstReplaceScoresWithOne);
		if (barParser.readBarFile(barParser.getBarFile(),true) == false){
			System.out.println("ERROR: loading bar file position and values -> "+barParser.getBarFile());
			return false;
		}
		positions = barParser.getBasePositions();
		scores = barParser.getValues();
		info.setNumberObservations(positions.length);
		return true;
	}
	
	/**Takes an int[] of base positions, some duplicated, sorted and returns a minimal PointData object 
	 * containing unique positions where the associated score is the counts observed for each position.*/
	public static PointData flattenPositions(int[] positions){
		ArrayList<Point> points = new ArrayList<Point>();
		int count = 1;
		int oldPosition = positions[0];
		for (int i=1; i< positions.length; i++){
			//same position?
			if (positions[i] == oldPosition) count++;
			else {
				//add position
				points.add(new Point(oldPosition, (float)count));
				//reset
				count = 1;
				oldPosition = positions[i];
			}
		}
		//add last
		points.add(new Point(oldPosition, (float)count));
		
		return Point.extractPositionScores(points);
	}
	
	/**Expands a hit count PointData object returning an int[] of unique and duplicate 
	 * positions based on the score value at each position*/
	public static int[] expandCounts( PointData converted){
		//expand the converted
		int[] positions = converted.getPositions();
		float[] counts = converted.getScores();
		int[] exPos = new int[(int)converted.getInfo().getScoreTotal()];
		int index =0;
		//for each position
		for (int i=0; i< positions.length; i++){
			//for each count
			int count = (int)counts[i];
			for (int j=0; j< count; j++) exPos[index++] = positions[i];
		}
		return exPos;
	}

	/**Nulls the position and score arrays to save memory.
	 * These will reload from file if you call for them after nulling.*/
	public void nullPositionScoreArrays(){
		positions = null;
		scores = null;
	}

	/**Finds the min and maximum positions from all of the position arrays.
	 * Will load the data if not already loaded.
	 * @return int[2] {min, max}
	 * */
	public static int[] findMinMaxPosition( PointData[] pd){
		int min = 1000000000;
		int max = -1;
		for (int i=0; i< pd.length; i++){
			int[] pos = pd[i].getPositions();
			if (pos[0]< min) min = pos[0];
			if (pos[pos.length-1] > max) max = pos[pos.length-1];
		}
		return new int[] {min, max};
	}

	/**Given a start bp (included) and stop bp (not included), returns the sum of the associate scores.*/
	public float sumScoreBP (int startBp, int stopBp){
		//sum?
		int[] startStop = findIndexes(startBp, stopBp);
		if ((startStop[1] - startStop[0]) <=0) return 0;
		return sumScoreIndex(startStop[0], startStop[1]);
	}

	/**Given a start bp (included) and stop bp (not included), returns the sum of the associate positions, the hit count.*/
	public float sumPositionBP (int startBp, int stopBp){
		//sum?
		int[] startStop = findIndexes(startBp, stopBp);
		int hitCount = startStop[1] - startStop[0];
		if (hitCount <=0) return 0;
		return hitCount;
	}

	/**Given a start bp (included) and stop bp (not included), returns the number of reads and the sum of their associated scores.
	 * @return float[2]{numReads, sumScores}*/
	public float[] sumScoresPositionsBP (int startBp, int stopBp){
		int[] startStop = findIndexes(startBp, stopBp);
		float sumScores = 0;
		//num reads
		float numReads = startStop[1] - startStop[0];
		if (numReads !=0) sumScores = sumScoreIndex(startStop[0], startStop[1]);
		return new float[]{numReads, sumScores};
	}

	/**Returns an array of Point containing the slice defined by the start and stop(excluded).*/
	public Point[] fetchPoints (int startBp, int stopBp){
		int[] indexes = findIndexes (startBp, stopBp);	
		if (indexes == null || indexes[0] == indexes[1]) return null;
		int num = indexes[1] - indexes[0];
		if (num ==0) return null;
		int[] pos = new int[num];
		float[] vals = new float[num];
		int counter = 0;
		for (int i=indexes[0]; i< indexes[1]; i++){
			pos[counter] = positions[i];
			vals[counter++] = scores[i];
		}
		return Point.makePoints(pos, vals);
	}
	
	/**Returns an array of float containing the slice defined by the start and stop(excluded).*/
	public float[] fetchScores (int startBp, int stopBp){
		int[] indexes = findIndexes (startBp, stopBp);	
		if (indexes == null || indexes[0] == indexes[1]) return null;
		int num = indexes[1] - indexes[0];
		if (num ==0) return null;
		float[] vals = new float[num];
		int counter = 0;
		for (int i=indexes[0]; i< indexes[1]; i++){
			vals[counter++] = scores[i];
		}
		return vals;
	}
	
	/**Returns number of observations in the slice defined by the start and stop(excluded).*/
	public int countPoints (int startBp, int stopBp){
		int[] indexes = findIndexes (startBp, stopBp);	
		if (indexes == null) return 0;
		return indexes[1] - indexes[0];
	}

	/**Returns an array of Point containing the slices defined by the start and stop(excluded). Zeros the PointData scores for these.
	 * Returns null if no Points found.*/
	public Point[] fetchPointsZeroPointDataScores (int[][] startStopBp){
		ArrayList<Point> points = new ArrayList<Point>();
		//for each slice
		for (int i=0; i< startStopBp.length; i++){
			int startBp = startStopBp[i][0];
			int stopBp = startStopBp[i][1];
			int[] indexes = findIndexes (startBp, stopBp);	
			if (indexes == null || indexes[0] == indexes[1]) continue;
			int counter = 0;
			for (int j=indexes[0]; j< indexes[1]; j++){
				Point p = new Point (positions[j], scores[j]);
				points.add(p);
				scores[j] = 0;
			}
		}
		if (points.size() == 0) return null;
		Point[] pts = new Point[points.size()];
		points.toArray(pts);
		return pts;
	}

	/**Removes PointData with a score of zero returning the number removed.*/
	public int removeZeroScoringPoints(){
		ArrayList<Point> points = new ArrayList<Point>();
		int numberRemoved = 0;
		for (int i=0; i< positions.length; i++){
			if (scores[i] != 0) points.add(new Point(positions[i], scores[i]));
			else numberRemoved++;
		}
		Point[] p = new Point[points.size()];
		points.toArray(p);
		PointData pd = Point.extractPositionScores(p);
		positions = pd.getPositions();
		scores = pd.getScores();
		return numberRemoved;
	}

	/**Given a start bp (included) and stop bp (not included), returns start (included) and stop (not included) indexes.
	 * May return startIndex = endIndex, therefore nothing found.*/
	public int[] findIndexes(int startBp, int stopBp){
		//find start index, included
		int startIndex = Arrays.binarySearch(positions, startBp);
		if (startIndex < 0) {
			startIndex = (startIndex*-1) -1;
		}
		else {
			//find first instance of startBp
			while (true){
				int minOne = startIndex - 1;
				if (minOne < 0) break;
				if (positions[minOne] != startBp) break;
				startIndex = minOne;
			}
		}
		//find stop index, not included
		int stopBPMinOne = stopBp-1;
		int endIndex = Arrays.binarySearch(positions, stopBPMinOne);		
		if (endIndex < 0) {
			endIndex = (endIndex*-1) -1;

		}
		else {
			//find last instance of endBp
			while (true){
				int addOne = endIndex +1;
				if (addOne >= positions.length) break;
				if (positions[addOne] != stopBPMinOne) break;
				endIndex = addOne;
			}
			//add one to stop index, it's not included
			endIndex++;
		}
		return new int[]{startIndex, endIndex};
	}

	/**Given a start index (included) and stop index (NOT included), returns the sum of the associate scores.*/
	public float sumScoreIndex (int startIndex, int endIndex){
		float sum = 0;
		for (int i=startIndex; i< endIndex; i++) sum+= scores[i];
		return sum;
	}

	public Point[] windowSum(int startBp, int stopBp, int halfWindowSize){
		int[] startStopIndexes = findIndexes(startBp, stopBp);
		//check to see if anything was found
		if (startStopIndexes[0] == startStopIndexes[1]) return null;

		ArrayList<Point> points = new ArrayList<Point>();
		int currentPosition = -1;

		//for each position, sum down and up sorted positions array
		for (int i=startStopIndexes[0]; i< startStopIndexes[1]; i++){
			//check position
			if (currentPosition != positions[i]) currentPosition = positions[i];
			else continue;

			float sum = scores[i];

			//look down/ left?
			int left = i -1;
			while (left >=0){
				//check bp distance
				int distance = positions[i]- positions[left];
				if ( distance <= halfWindowSize) {
					sum += scores[left];
					left--;
				}
				else break;
			}
			//look up/right
			int right = i + 1;
			while (right < positions.length){
				int distance = positions[right] - positions[i];
				if (distance <= halfWindowSize){
					sum += scores[right];
					right++;
				}
				else break;
			}
			//make Point
			points.add(new Point(positions[i], sum));
		}
		Point[] pt = new Point[points.size()];
		points.toArray(pt);
		return pt;

	}



	/**Randomly divides data.
	 * @param splitIn should be 2 or 1/2, 3 for 1/3rds etc.*/
	public static HashMap<String,PointData>[] divide(HashMap<String,PointData> pd, int splitIn){
		//check numbers
		int totalObs = PointData.totalObservations(pd);
		int numInEachSplit = totalObs/splitIn;
		//combine  into an array
		InfoPoint[] ip = new InfoPoint[totalObs];
		int index =0;
		Iterator<String> it = pd.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			PointData pds = pd.get(chrom);
			int numObs = pds.getInfo().getNumberObservations();
			int[] pos = pds.getPositions();
			float[] scr = pds.getScores();
			Info info = pds.getInfo();
			for (int i=0; i< numObs; i++){
				ip[index++] = new InfoPoint(pos[i], scr[i], info);
			}
		}

		//randomize
		Misc.randomize(ip, 0);

		//for each splitIn
		HashMap<String,PointData>[] splits = new HashMap[splitIn];
		for (int x=0; x< splitIn; x++){
			int startInAll = x * numInEachSplit;
			//take numberWanted
			InfoPoint[] subIP = new InfoPoint[numInEachSplit];
			System.arraycopy(ip, startInAll, subIP, 0, numInEachSplit);
			//sort by chromosome
			ComparatorInfoPoint comp = new ComparatorInfoPoint();
			Arrays.sort(subIP, comp);

			//rebuild HashMap
			HashMap<String,PointData> split = new HashMap<String, PointData>();
			//make first AL
			ArrayList<Point> al = new ArrayList<Point>();
			Point p = new Point(subIP[0].getPosition(), subIP[0].getScore());
			al.add(p);
			Info info = subIP[0].getInfo();
			String chrom = info.getChromosome();
			for (int i=1; i< subIP.length; i++){
				//has chrom changed?
				String testChrom = subIP[i].getInfo().getChromosome();
				//no change just add new Point
				if (testChrom.equals(chrom)) al.add(new Point(subIP[i].getPosition(), subIP[i].getScore()));
				//yes change, add new entry to HashMap and reset
				else {
					//make PointData
					Point[] points = new Point[al.size()];
					al.toArray(points);
					Arrays.sort(points, new ComparatorPointPosition());
					PointData newPD = Point.extractPositionScores(points);
					info.setName("split");
					info.setNotes(null);
					newPD.setInfo(info);
					//add to hash
					split.put(chrom, newPD);				
					//reset
					info = subIP[i].getInfo();
					chrom = testChrom;
					al.clear();
					al.add(new Point(subIP[i].getPosition(), subIP[i].getScore()));
				}
			}
			//add last HashMap entry
			Point[] points = new Point[al.size()];
			al.toArray(points);
			Arrays.sort(points, new ComparatorPointPosition());
			PointData newPD = Point.extractPositionScores(points);
			info.setName("split");
			info.setNotes(null);
			newPD.setInfo(info);
			//add to hash
			split.put(chrom, newPD);

			//assign to array
			splits[x] = split;

			//run through and set number of points
			Iterator<String> itx = splits[x].keySet().iterator();
			while(itx.hasNext()){
				PointData px = splits[x].get(itx.next());
				int num= px.positions.length;
				px.getInfo().setNumberObservations(num);
			}
		}
		return splits;
	}

	/**Shifts bp positions numToShift 3' based on strand, if no strand info returns false.*/
	public boolean shiftPositions(int numToShift){
		String strand = info.getStrand();
		getPositions();
		if (strand.equals("+")) {
			for (int z=0; z< positions.length; z++) {
				positions[z] += numToShift;
			}
		}
		else if (strand.equals("-")) {
			for (int z=0; z< positions.length; z++) {
				positions[z] -= numToShift;
				if (positions[z] < 0) positions[z] = 0;
			}
		}
		else return false;
		return true;
	}

	/**Adds numToShift to each position, unstranded*/
	public void shiftPositionsUnstranded(int numToShift){
		for (int z=0; z< positions.length; z++) {			
			positions[z] += numToShift;
			if (positions[z] <0) positions[z] =0;
		}
	}

	/**Splits Stranded PointData in half.*/
	public static HashMap<String,PointData[]>[] splitStranded (HashMap<String,PointData[]> pd){
		//first split by strand
		HashMap<String,PointData>[] plusMinus = PointData.splitStrandedPointData(pd);
		//then split each strand in half
		HashMap<String,PointData>[] plus = PointData.divide(plusMinus[0], 2);
		HashMap<String,PointData>[] minus = PointData.divide(plusMinus[1], 2);
		//reassemble
		HashMap<String,PointData[]>[] split = new HashMap[2];
		split[0] = PointData.joinPointData(plus[0], minus[0]);
		split[1] = PointData.joinPointData(plus[1], minus[1]);
		return split;
	}



	/**Splits the data in half using random sampling.
	 * Will first call PointData.subSample() if not evenly divisable by 2.*/
	public static HashMap<String,PointData>[]  split(HashMap<String,PointData> pd){	
		//calculate the total number of observations
		int totalObs = totalObservations(pd);
		//need to drop one?
		if (totalObs % 2 != 0) {
			totalObs--;
			PointData.subSample(pd, totalObs);
		}
		int halfObs = totalObs/2;
		int numChroms = pd.size();

		//make array to track how to split, default is false
		boolean[][] chrToFlag = new boolean[numChroms][];

		//make array to track start and stops
		int[][] chrStartStops = new int[numChroms][];

		//for each chromosome, define start (included), stop (excluded)
		Iterator<PointData> it = pd.values().iterator();
		int total =0;
		int index =0;
		while(it.hasNext()){
			int start = total;
			PointData p = it.next();
			int numObs = p.getInfo().getNumberObservations();
			chrToFlag[index] = new boolean[numObs];
			int stop = numObs + total;
			int[] ss = {start, stop};
			chrStartStops[index++] = ss;
			total = stop;
		}

		//for half the number of observations
		Random random = new Random();
		int toPick = halfObs;
		while (toPick >0){
			//pick a random number from the total number observations
			int rndNum = random.nextInt(totalObs);
			//find the appropriate chromosome
			for (int j=0; j<numChroms; j++){
				int[] startStop = chrStartStops[j];
				if (rndNum >= startStop[0] && rndNum < startStop[1]) {
					//set a flag?
					index = rndNum - startStop[0];
					if (chrToFlag[j][index] == false){
						chrToFlag[j][index] = true;
						toPick--;
					}
					break;
				}
			}
		}
		//make hashes to hold results
		HashMap<String,PointData>[] splitData = new HashMap[2];
		splitData[0] = new HashMap();
		splitData[1] = new HashMap();

		//for each chromosome, split based on boolean array
		Iterator<String> it2 = pd.keySet().iterator();
		index = 0;
		while(it2.hasNext()){
			String chrom = it2.next();
			PointData p = pd.get(chrom);
			//split
			boolean[] flags = chrToFlag[index++];
			PointData[] split = p.split(flags);
			splitData[0].put(chrom, split[0]);
			splitData[1].put(chrom, split[1]);
		}
		return splitData;
	}

	/**Splits this.PointData's positions and scores based on the boolean[].
	 * @param b -must be same length as this.positions or returns null
	 * @return two new PointDatas with info of parent, watch it, these reference the parent's info.*/
	public PointData[] split(boolean[] b){
		if (b.length != positions.length) {
			System.out.println("Error: boolean length not equal to length of positions scores "+b.length+" "+info.getNumberObservations());
			return null;
		}
		//count number true
		int numTrue = Num.countNumberTrue(b);
		//make arrays
		Point[] p1 = new Point[numTrue];
		Point[] p2 = new Point[b.length - numTrue];
		//split
		int counter1 =0;
		int counter2 = 0;
		for (int i=0; i< positions.length; i++){
			if (b[i]) p1[counter1++] = new Point(positions[i], scores[i]);
			else p2[counter2++] = new Point(positions[i], scores[i]);
		}
		//make PointData
		PointData pd1 = Point.extractPositionScores(p1);
		PointData pd2 = Point.extractPositionScores(p2);
		//add info
		Info i1 = new Info("split1_"+info.getName(), info.getVersionedGenome(), info.getChromosome(), info.getStrand(), info.getReadLength(), info.getNotes());
		Info i2 = new Info("split2_"+info.getName(), info.getVersionedGenome(), info.getChromosome(), info.getStrand(), info.getReadLength(), info.getNotes());
		pd1.setInfo(i1);
		pd2.setInfo(i2);

		return new PointData[] {pd1, pd2};
	}

	/**Sub Samples stranded PointData[2]{+,-} to lowest number of observations replacing original PointData*/
	public static void subSamplePointData(PointData[] t, PointData[] c){
		int[] numObs= {
				t[0].getInfo().getNumberObservations(),
				t[1].getInfo().getNumberObservations(),
				c[0].getInfo().getNumberObservations(),
				c[1].getInfo().getNumberObservations(),
		};
		int lowest = Num.findSmallestInt(numObs);		
		if (numObs[0] != lowest) t[0] = PointData.fetchRandomObservations(t[0], lowest);
		if (numObs[1] != lowest) t[1] = PointData.fetchRandomObservations(t[1], lowest);
		if (numObs[2] != lowest) c[0] = PointData.fetchRandomObservations(c[0], lowest);
		if (numObs[3] != lowest) c[1] = PointData.fetchRandomObservations(c[1], lowest);
	}

	/**Fast. Be sure the numberToFetch is <= number of observations!*/
	public static PointData fetchRandomObservations (PointData pd, int numberToFetch){
		Point[] subSample = new Point[numberToFetch];
		Point[] points = Point.makePoints(pd.getPositions(), pd.getScores());
		Misc.randomize(points, 0);
		System.arraycopy(points, 0, subSample, 0, numberToFetch);
		Arrays.sort(subSample, new ComparatorPointPosition());
		PointData ss = Point.extractPositionScores(subSample);
		Info info = pd.getInfo();
		info.setNumberObservations(numberToFetch);
		ss.setInfo(info);
		return ss;
	}

	/**Lots of work to fetch random observation.
	 * Returns null if you ask for more observations that present.*/
	public static HashMap<String,PointData> fetchRandomObservations(HashMap<String,PointData> pd, int numWanted){
		//check numbers
		int totalObs = PointData.totalObservations(pd);
		if (numWanted > totalObs) return null;

		//combine  into an array
		InfoPoint[] ip = new InfoPoint[totalObs];
		int index =0;
		Iterator<String> it = pd.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			PointData pds = pd.get(chrom);
			int numObs = pds.getInfo().getNumberObservations();
			int[] pos = pds.getPositions();
			float[] scr = pds.getScores();
			Info info = pds.getInfo();
			for (int i=0; i< numObs; i++){
				ip[index++] = new InfoPoint(pos[i], scr[i], info);
			}
		}

		//randomize
		Misc.randomize(ip, 0);

		//take numberWanted
		InfoPoint[] subIP = new InfoPoint[numWanted];
		System.arraycopy(ip, 0, subIP, 0, numWanted);

		//sort by chromosome
		ComparatorInfoPoint comp = new ComparatorInfoPoint();
		Arrays.sort(subIP, comp);

		//rebuild HashMap
		HashMap<String,PointData> split = new HashMap<String, PointData>();
		//make first AL
		ArrayList<Point> al = new ArrayList<Point>();
		Point p = new Point(subIP[0].getPosition(), subIP[0].getScore());
		al.add(p);
		Info info = subIP[0].getInfo();
		String chrom = info.getChromosome();
		for (int i=1; i< subIP.length; i++){
			//has chrom changed?
			String testChrom = subIP[i].getInfo().getChromosome();
			//no change just add new Point
			if (testChrom.equals(chrom)) al.add(new Point(subIP[i].getPosition(), subIP[i].getScore()));
			//yes change, add new entry to HashMap and reset
			else {
				//make PointData
				Point[] points = new Point[al.size()];
				al.toArray(points);
				Arrays.sort(points, new ComparatorPointPosition());
				PointData newPD = Point.extractPositionScores(points);
				info.setName("sub sampled");
				info.setStrand(".");
				info.setNotes(null);
				newPD.setInfo(info);
				//add to hash
				split.put(chrom, newPD);				
				//reset
				info = subIP[i].getInfo();
				chrom = testChrom;
				al.clear();
				al.add(new Point(subIP[i].getPosition(), subIP[i].getScore()));
			}
		}
		//add last HashMap entry
		Point[] points = new Point[al.size()];
		al.toArray(points);
		Arrays.sort(points, new ComparatorPointPosition());
		PointData newPD = Point.extractPositionScores(points);
		info.setName("sub sampled");
		info.setStrand(".");
		info.setNotes(null);
		newPD.setInfo(info);
		//add to hash
		split.put(chrom, newPD);		
		return split;
	}

	/**Splits a stranded PointData hashmap by strand returning plus and minus hashmaps.*/
	public static HashMap<String,PointData>[] splitStrandedPointData (HashMap<String,PointData[]> strandedPD){
		HashMap<String,PointData>[] plusMinus = new HashMap[2];
		plusMinus[0] = new HashMap<String,PointData>();
		plusMinus[1] = new HashMap<String,PointData>();
		Iterator<String> it = strandedPD.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			PointData[] pd = strandedPD.get(chrom);			
			if (pd[0] != null) plusMinus[0].put(chrom, pd[0]);
			if (pd[1] != null) plusMinus[1].put(chrom, pd[1]);
		}
		return plusMinus;
	}

	/**Joins point data, assumes same chromosome set.*/
	public static HashMap<String,PointData[]> joinPointData(HashMap<String,PointData> plus, HashMap<String,PointData> minus){
		HashMap<String,PointData[]> joined = new HashMap<String,PointData[]>();
		Iterator<String> it = plus.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			PointData[] pd = new PointData[2];
			pd[0] = plus.get(chrom);
			pd[1] = minus.get(chrom);
			joined.put(chrom, pd);
		}
		return joined;
	}


	/**Reduces the size of the stranded point data.*/
	public static boolean subSampleStranded(HashMap<String,PointData[]> pd, int numPos){
		//split by strand	
		HashMap<String,PointData>[] plusMinus = splitStrandedPointData(pd);
		//split numToDrop by half
		int numToDrop = (int)Math.round((double)(numPos)/2);
		if(subSample(plusMinus[0], numToDrop)){
			if (subSample(plusMinus[1],numToDrop)) return true;
			else return false;
		}
		else return false;
	}

	/**Reduces the number of observations to the specified number by random selection. This does not work with 
	 * drastic sub samplings, ie > 1/2 the data.  Becomes inefficient, use fetchRandomObservations.*/
	public static boolean subSample(HashMap<String,PointData> pd, int numPos){		
		//calculate the total number of observations
		int totalObs = totalObservations(pd);	
		int toDrop = totalObs - numPos;
		if (toDrop <= 0) return false;
		int numChroms = pd.size();

		//make array to track how many to drop from each chromosome
		int[] chrToDrop = new int[numChroms];
		//make array to track start and stops
		int[][] chrStartStops = new int[numChroms][];

		//for each chromosome, define start (included), stop (excluded)
		Iterator<PointData> it = pd.values().iterator();
		int total =0;
		int index =0;
		while(it.hasNext()){
			int start = total;
			PointData p = it.next();
			int stop = p.getInfo().getNumberObservations() + total;
			int[] ss = {start, stop};
			chrStartStops[index++] = ss;
			total = stop;
		}

		//for each toDrop
		Random random = new Random();
		for (int i=0; i< toDrop; i++){
			//pick a random number
			int rndNum = random.nextInt(totalObs);
			//find the appropriate chromosome
			for (int j=0; j<numChroms; j++){
				int[] startStop = chrStartStops[j];
				if (rndNum >= startStop[0] && rndNum < startStop[1]) {
					chrToDrop[j]++;
					break;
				}
			}
		}

		//for each chromosome, fetch point data and drop values
		it = pd.values().iterator();

		index = 0;
		while(it.hasNext()){
			PointData p = it.next();
			//any to drop?
			if (chrToDrop[index]>0){
				p.drop(chrToDrop[index]);
			}
			index++;
		}

		return true;
	}

	/**Will draw the numPos wanted from the array of PointData sampling randomly within but equally throughout.*/
	public static HashMap<String, PointData> subSamplePercent(HashMap<String, PointData>[] pd, int numPos){
		//find total reads
		double total = numPos;
		double totalReads = PointData.totalObservations(pd);
		double fraction = total/totalReads;
		//for each HashMap, subsample
		int subSampleTotal = 0;
		for (int i=0; i<pd.length; i++){
			double numObs = PointData.totalObservations(pd[i]);
			int numPosNeeded = (int)Math.round(fraction * numObs);
			subSampleTotal += numPosNeeded;
			PointData.subSample(pd[i], numPosNeeded);
		}
		//make merge
		HashMap<String, PointData> merged = PointData.mergePointData(pd, false, true);
		//subsample
		if ((int)total != subSampleTotal) PointData.subSample(merged, (int)total);
		return merged;
	}

	/**Randomly removes the numToDrop from the position and score arrays.
	 * Use to shrink the size of the data.*/
	public void drop(int numToDrop){
		int numObs = info.getNumberObservations();
		//load data if not loaded
		getPositions();
		//use a boolean array to track indexes to drop, default is false
		boolean[] drop = new boolean[numObs];
		Random random = new Random();
		int counter = numToDrop;
		while (counter > 0){
			//get a random index
			int index = random.nextInt(numObs);
			//check if it has been changed, if not flag to drop
			if (drop[index] == false) {
				drop[index] = true;
				counter--;
			}
		}
		//make new arrays
		int newSize = numObs- numToDrop;
		int[] tPositions = new int[newSize];
		float[] tScores = new float[newSize];
		int index =0;
		for (int i=0; i< numObs; i++){
			if (drop[i] == false) {
				tPositions[index] = positions[i];
				tScores[index] = scores[i];
				index++;
			}
		}
		//reset objects
		positions = tPositions;
		scores = tScores;
		info.setNumberObservations(positions.length);
		info.setScoreTotal(Num.sumArrayReturnDouble(scores));
	}


	/**Randomly picks the numToPick from the position and score arrays and replaces.*/
	public boolean pick(int numToPick){
		if (numToPick > info.getNumberObservations()) {
			System.err.println("\nError: attempting to pick more observations that present in the PointData\n");
			return false;
		}
		//load data if not loaded
		getPositions();
		//make PointData array
		Point[] pda = Point.makePoints(positions, scores);
		//randomize
		Misc.randomize(pda, 0);
		//make new Point[]
		Point[] picks = new Point[numToPick];
		//copy
		System.arraycopy(pda, 0, picks, 0, numToPick);
		//sort
		Arrays.sort(picks, new ComparatorPointPosition());
		//assign 
		PointData p = Point.extractPositionScores(picks);
		positions = p.getPositions();
		scores = p.getScores();
		info.setNumberObservations(positions.length);
		info.setScoreTotal(Num.sumArrayReturnDouble(scores));
		return true;
	}

	/**Counts the total number of observations.
	 * This does not load the data.*/
	public static int totalObservationsMultiPointData(HashMap<String,PointData[]> pd){
		int total = 0;
		Iterator<PointData[]> it = pd.values().iterator();
		while(it.hasNext()){
			PointData[] p = it.next();
			for (int i=0; i< p.length; i++){
				//watch out for no stranded data
				if (p[i] == null) continue;
				total += p[i].getInfo().getNumberObservations();
			}
		}
		return total;
	}

	/**Sums the total scores. This does not load the data. Returns -1 if appropriate tag/value not found.*/
	public static double totalScoreMultiPointData(HashMap<String,PointData[]> pd){
		double total = 0;
		Iterator<PointData[]> it = pd.values().iterator();
		while(it.hasNext()){
			PointData[] p = it.next();
			for (int i=0; i< p.length; i++){
				//watch out for no stranded data
				if (p[i] == null) continue;
				double sum = p[i].getInfo().getScoreTotal();
				if (sum == 0) return -1;
				total += sum;
			}
		}
		return total;
	}

	/**Replaces all scores with 1.*/
	public void stripScores(){
		float[] scores = getScores();
		Arrays.fill(scores, 1);
	}

	/**Counts the number of observations by chromosome. This does not load the data.*/
	public static HashMap<String,Integer> totalObservationsByChromosome(HashMap<String,PointData[]> pd){
		HashMap<String,Integer> chromCount = new HashMap<String,Integer>();
		Iterator<String> it = pd.keySet().iterator();
		while(it.hasNext()){
			String chrom = it.next();
			PointData[] p = pd.get(chrom);
			int total = 0;
			for (int i=0; i< p.length; i++){
				//watch out for no stranded data
				if (p[i] == null) continue;
				total += p[i].getInfo().getNumberObservations();
			}
			chromCount.put(chrom, new Integer(total));
		}
		return chromCount;
	}

	/**Converts a HashMap<String,ArrayList<PointData>> to a HashMap<String,PointData[]>, note
	 * the PointData[] can contain any number of PointDatas.*/
	public static HashMap<String,PointData[]> convertArrayList2Array (HashMap<String,ArrayList<PointData>> pd){
		HashMap<String,PointData[]> array = new HashMap<String,PointData[]>();
		Iterator<String> it = pd.keySet().iterator();
		while (it.hasNext()){
			String chr = it.next();
			ArrayList<PointData> al = pd.get(chr);
			PointData[] p = new PointData[al.size()];
			al.toArray(p);
			array.put(chr, p);
		}
		return array;
	}

	public static ArrayList<PointData> convertArray2ArrayList(PointData[] pd){
		ArrayList<PointData> al = new ArrayList<PointData>(pd.length);
		for (int i=0; i< pd.length; i++) al.add(pd[i]);
		return al;
	}

	/**Counts the total number of observations.
	 * This does not load the data.*/
	public static int totalObservations(HashMap<String,PointData> pd){
		int total = 0;
		Iterator<PointData> it = pd.values().iterator();
		while(it.hasNext()){
			PointData p = it.next();
			total += p.getInfo().getNumberObservations();
		}
		return total;
	}

	/**Counts the total number of observations.
	 * This does not load the data.*/
	public static int totalObservations(HashMap<String,PointData>[] pd){
		int total = 0;
		for (int i=0; i< pd.length; i++) total += PointData.totalObservations(pd[i]);
		return total;
	}

	/**Combines without summing the scores, thus duplicate positions are likely.*/
	public static PointData combinePointData(PointData[] pd, boolean nullOriginalData){
		ArrayList<PointData> pdAL = new ArrayList<PointData>(pd.length);
		for (int i=0; i< pd.length; i++) {
			if (pd[i] != null) pdAL.add(pd[i]);
		}
		return PointData.combinePointData(pdAL, nullOriginalData);
	}

	/**Combines without summing the scores, thus duplicate positions are likely.*/
	public static PointData combinePointData(ArrayList<PointData> pdAL, boolean nullOriginalPosVal){
		//load data
		int num = pdAL.size();
		PointData[] pd = new PointData[num];
		pdAL.toArray(pd);

		//find total length and if mixed strands
		int totalPositions = 0;
		boolean allSameStrand = true;
		String strand = pd[0].getInfo().getStrand();		
		for (int i=0; i<num; i++) {
			totalPositions += pd[i].getInfo().getNumberObservations();
			//check strand?
			if (allSameStrand){
				String testStrand = pd[i].getInfo().getStrand();				
				if (strand.equals(testStrand) == false){
					allSameStrand = false;
				}
			}
		}

		//load Point[]
		Point[] points = new Point[totalPositions];	
		int counter = 0;
		for (int i=0; i<num; i++) {
			int[] pos = pd[i].getPositions();
			float[] scr = pd[i].getScores();
			//make Points
			for (int j=0; j< pos.length; j++) points[counter++] = new Point (pos[j],scr[j]);
			if (nullOriginalPosVal) pd[i].nullPositionScoreArrays();
		}

		//sort
		Arrays.sort(points, new ComparatorPointPosition());

		//split
		int[] positions = new int[totalPositions];
		float[] scores = new float[totalPositions];
		for (int i=0; i< totalPositions; i++){
			positions[i] = points[i].position;
			scores[i] = points[i].score;
		}
		points = null;

		//Info info = new Info()
		PointData pdSum = new PointData();
		pdSum.setInfo(pd[0].getInfo());
		pdSum.setPositions(positions);
		pdSum.setScores(scores);
		if (allSameStrand == false) pdSum.getInfo().setStrand(".");
		pdSum.getInfo().setName("merged");
		return pdSum;
	}

	//public static HashMap<String,PointData[]> combinePointData(ArrayList<HashMap<String,PointData[]>> al, boolean nullOriginalPosVal){
	/**Combines point data without summing scores. Thus duplicates likely.  Maintains strandedness*/
	public static HashMap<String,PointData[]> combinePointDataToHashMap(ArrayList<HashMap<String,PointData[]>> al, boolean nullOriginalPosVal){
		//split by strand
		HashMap<String,ArrayList<PointData>> plus = new HashMap<String,ArrayList<PointData>>();
		HashMap<String,ArrayList<PointData>> minus = new HashMap<String,ArrayList<PointData>>();
		int num = al.size();
		for (int i=0; i< num; i++){
			HashMap<String,PointData[]> pd = al.get(i);
			Iterator<String> it = pd.keySet().iterator();
			while (it.hasNext()){
				String chrom = it.next();
				PointData[] p = pd.get(chrom);
				//plus
				ArrayList<PointData> alPlus;
				if (plus.containsKey(chrom)) alPlus = plus.get(chrom);
				else {
					alPlus = new ArrayList<PointData>();
					plus.put(chrom, alPlus);
				}
				alPlus.add(p[0]);
				//plus
				ArrayList<PointData> alMinus;
				if (minus.containsKey(chrom)) alMinus = minus.get(chrom);
				else {
					alMinus = new ArrayList<PointData>();
					minus.put(chrom, alMinus);
				}
				alMinus.add(p[1]);
			}
		}
		//combine
		HashMap<String,PointData> plusCombine = combinePointData(plus,nullOriginalPosVal);
		HashMap<String,PointData> minusCombine = combinePointData(minus,nullOriginalPosVal);
		//combine strands
		HashMap<String,PointData[]> combo = PointData.joinPointData(plusCombine, minusCombine);
		return combo;
	}

	/**Combines without summing the scores, thus duplicate positions are likely.*/
	public static HashMap<String,PointData> combinePointData(HashMap<String,ArrayList<PointData>> hash, boolean nullOriginalPosVal){
		HashMap<String,PointData> combine = new HashMap<String,PointData>();
		Iterator<String> it = hash.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			ArrayList<PointData> pdAL = hash.get(chrom);
			PointData pd = PointData.combinePointData(pdAL, nullOriginalPosVal);
			combine.put(chrom, pd);
		}
		return combine;
	}

	/**Merges a HashMap<Chrom,PointData> by either summing same position scores or simply concatinating then sorting the arrays.*/
	public static HashMap<String,PointData> mergePointData(HashMap<String,PointData>[] hashes, boolean sumScores, boolean nullOriginalPosVal){
		HashMap<String,PointData> results = new HashMap<String, PointData>();

		//make set of all chroms
		HashSet<String> chroms = new HashSet<String>();
		for (int i=0; i< hashes.length; i++) chroms.addAll(hashes[i].keySet());

		//for each chromosome
		Iterator<String> it = chroms.iterator();
		while (it.hasNext()){
			String chrom = it.next();
			ArrayList<PointData> pds = new ArrayList<PointData>();
			//for each hash fetch chrom
			for (int i=0; i< hashes.length; i++){
				if (hashes[i].containsKey(chrom)) pds.add(hashes[i].get(chrom));
			}
			//merge?
			PointData merged;
			if (pds.size() == 1) merged = pds.get(0);
			else if (sumScores) merged = PointData.mergePointDataDepreciated(pds, nullOriginalPosVal);
			else merged = PointData.combinePointData(pds, nullOriginalPosVal);
			//add to results
			results.put(chrom, merged);
		}
		return results;
	}


	/**Merges PointData, nulls originals positions and scores, sets strand to '.' if mixed strand.
	 * Sums scores of identical positions. */
	public static PointData mergePointDataDepreciated(ArrayList<PointData> pdAL, boolean nullOriginalPosVal){
		//find total points and if all same strand
		int num = pdAL.size();
		boolean allSameStrand = true;
		PointData pd = pdAL.get(0);
		String strand = pd.getInfo().getStrand();
		int totalPoints = pd.getInfo().getNumberObservations();
		for (int i=1; i< num ; i++){
			pd = pdAL.get(i);
			totalPoints += pd.getInfo().getNumberObservations();
			//check strand?
			if (allSameStrand){
				if (strand.equals(pd.getInfo().getStrand()) == false){
					allSameStrand = false;
				}
			}
		}
		//make combine Point[]
		Point[] points = new Point[totalPoints];
		int index = 0;
		for (int i=0; i< num; i++){
			pd = pdAL.get(i);
			Point[] pts = Point.makePoints(pd);
			System.arraycopy(pts, 0, points, index, pts.length);
			index+= pts.length;
		}
		//sort
		Arrays.sort(points, new ComparatorPointPosition());
		//merge
		points = Point.sumIdenticalPositionScores(points);
		PointData pdMerged = Point.extractPositionScores(points);

		//Info info = new Info()
		PointData pdSum = Point.extractPositionScores(points);
		pdSum.setInfo(pd.getInfo());
		if (allSameStrand == false) pdSum.getInfo().setStrand(".");
		pdSum.getInfo().setName("merged");
		pdSum.getInfo().setNumberObservations(points.length);

		//null position and scores
		if (nullOriginalPosVal) for (int i=0; i<num; i++) pdAL.get(i).nullPositionScoreArrays();
		return pdSum;
	}

	/**This is a more memory efficient method for merging PointData, nulls originals positions and scores, sets strand to '.' if mixed strand.
	 * Sums scores of identical positions. This reloads the point data from file so anything done to it will not be used.*/
	public static PointData mergePointData (ArrayList<PointData> pdAL, boolean replaceScoresWithHitCount, boolean nullOriginals){
		PointData merged = pdAL.get(0);
		merged.loadFlattenedPositionScores(replaceScoresWithHitCount);
		for (int i=1; i< pdAL.size(); i++){
			PointData pd = pdAL.get(i);
			pd.loadFlattenedPositionScores(replaceScoresWithHitCount);
			merged = mergePairedPointData(merged, pd, nullOriginals);
		}
		//no need to sort, this done by mergePairedPointData
		return merged;
	}

	/**Merges two PointData arrays, nulls originals positions and scores, sets strand to '.' if mixed strand.
	 * Sums scores of identical positions.*/
	public static PointData mergePairedPointData(PointData one, PointData two, boolean nullOriginals){
		Point[] onePts = Point.makePoints(one);
		Point[] twoPts = Point.makePoints(two);
		//make combine Point[]
		Point[] points = new Point[onePts.length+ twoPts.length];
		for (int i=0; i< onePts.length; i++) points[i] = onePts[i];
		int index = onePts.length;
		for (int i=0; i< twoPts.length; i++) points[index++] = twoPts[i];
		//sort
		Arrays.sort(points, new ComparatorPointPosition());

		//merge
		points = Point.sumIdenticalPositionScores(points);
		PointData pdMerged = Point.extractPositionScores(points);

		//check strands
		String strand = one.getInfo().getStrand();
		if (strand.equals(two.getInfo().getStrand()) == false) strand = ".";

		//make info object
		Info info = new Info();
		info.setNumberObservations(points.length);
		info.setScoreTotal(Num.sumArray(pdMerged.getScores()));
		info.setStrand(strand);
		info.setName("Merged");
		info.setVersionedGenome(one.getInfo().getVersionedGenome());
		info.setChromosome(one.getInfo().getChromosome());
		info.setReadLength(one.getInfo().getReadLength());
		pdMerged.setInfo(info);

		//null objects
		if (nullOriginals){
			one.nullPositionScoreArrays();
			two.nullPositionScoreArrays();
		}
		onePts = null;
		twoPts = null;
		points = null;

		return pdMerged;
	}

	/**Merges two PointData arrays, sets strand to '.' if mixed strand. No sum of identical positions or nulling of originals.*/
	public static PointData mergePairedPointDataNoSumming(PointData one, PointData two){
		Point[] onePts = Point.makePoints(one);
		Point[] twoPts = Point.makePoints(two);
		//make combine Point[]
		Point[] points = new Point[onePts.length+ twoPts.length];
		for (int i=0; i< onePts.length; i++) points[i] = onePts[i];
		int index = onePts.length;
		for (int i=0; i< twoPts.length; i++) points[index++] = twoPts[i];
		//sort
		Arrays.sort(points, new ComparatorPointPosition());

		//merge
		PointData pdMerged = Point.extractPositionScores(points);

		//check strands
		String strand = one.getInfo().getStrand();
		if (strand.equals(two.getInfo().getStrand()) == false) strand = ".";

		//make info object
		Info info = new Info();
		info.setNumberObservations(points.length);
		info.setScoreTotal(Num.sumArray(pdMerged.getScores()));
		info.setStrand(strand);
		info.setName("Merged");
		info.setVersionedGenome(one.getInfo().getVersionedGenome());
		info.setChromosome(one.getInfo().getChromosome());
		info.setReadLength(one.getInfo().getReadLength());
		pdMerged.setInfo(info);

		onePts = null;
		twoPts = null;
		points = null;

		return pdMerged;
	}

	/**Loads a directory of xxx.bar PointData files into a HashMap where the
	 * key = */
	public static HashMap<String, ArrayList<PointData>> fetchPointData (File directory){
		//fetch files
		File[] files = IO.extractFiles(directory,".bar.zip");
		if (files == null || files.length == 0) return null;
		HashMap<String, ArrayList<PointData>> pointData = new HashMap<String, ArrayList<PointData>>();
		//for each file 
		for (int j=0; j< files.length; j++){
			//make a PointData object
			PointData pd = new PointData (files[j], false); 
			//add to hash
			String chr;
			if (pd.getInfo().getStrand().equals(".")) chr = pd.getInfo().getChromosome();
			else chr = pd.getInfo().getChromosome() +"_"+pd.getInfo().getStrand(); 
			if (pointData.containsKey(chr)) pointData.get(chr).add(pd);
			else {
				ArrayList<PointData> al = new ArrayList<PointData>();
				al.add(pd);
				pointData.put(chr, al);
			}
		}
		return pointData;
	}

	/**Multi directory version of fetchPointData(), set mergeDirectories to null to keep everything loaded in memory.*/
	public static HashMap<String, PointData>[] fetchPointData (File[] directories, File[] mergeDirectories, boolean sumScores){
		HashMap<String, PointData>[] hashes = new HashMap[directories.length];
		if (mergeDirectories == null) for (int i=0; i< directories.length; i++) {
			hashes[i] = fetchPointData(directories[i], null, sumScores);
		}
		else for (int i=0; i< directories.length; i++) hashes[i] = fetchPointData(directories[i], mergeDirectories[i], sumScores);
		return hashes;
	}

	/**Multi directory version of fetchStrandedPointData(). Fetches, but doesn't load or merge, PointData
	 * Returns null if no PointData files found.*/
	public static HashMap<String, PointData[]>[] fetchStrandedPointData (File[] directories){
		HashMap<String, PointData[]>[] hashes = new HashMap[directories.length];
		for (int i=0; i< directories.length; i++) {
			hashes[i] = fetchStrandedPointData(directories[i]);
			if (hashes[i] == null) return null;
		}
		return hashes;
	}

	/**Fetches PointData, returning two chromosome split hashmaps<String, ArrayList<PointData>>. 
	 * HashMap[0] is from the plus strand, 
	 * HashMap[1] minus strand. 
	 * Either may be null if single stranded.
	 * Will add unstranded data to plus strand.
	 * No loading of data.
	 * Returns null if no PointData files found.*/
	public static HashMap<String, ArrayList<PointData>>[] fetchStrandedPointDataNoMerge (File[] directories){
		//load hash map with stranded data
		HashMap<String,ArrayList<PointData>> plusPDAL = new HashMap<String,ArrayList<PointData>>();
		HashMap<String,ArrayList<PointData>> minusPDAL = new HashMap<String,ArrayList<PointData>>();
		for (int i = 0; i< directories.length; i++){
			//fetch files			
			File[] files = IO.extractFiles(directories[i], ".bar");
			if (files == null || files.length == 0) {
				files = IO.extractFiles(directories[i], ".bar.zip");
				if (files == null || files.length == 0) {
					System.err.println("\t\tError: no Point Data files found in -> "+directories[i]+"\n");
					return null;
				}	
			}
			//for each file 
			for (int j=0; j< files.length; j++){
				//make a PointData object
				PointData pd = new PointData (files[j], false);
				//get chrom and strand and hash
				String chr = pd.getInfo().getChromosome();
				String strand = pd.getInfo().getStrand();
				HashMap<String,ArrayList<PointData>> pdHash;
				if (strand.equals("-")) pdHash = minusPDAL;
				//add plus and no stranded data to plus
				else pdHash = plusPDAL;
				//add to hash
				if (pdHash.containsKey(chr)) pdHash.get(chr).add(pd);
				else {
					ArrayList<PointData> al = new ArrayList<PointData>();
					al.add(pd);
					pdHash.put(chr, al);
				}
			}
		}
		return new HashMap[] {plusPDAL, minusPDAL};
	}

	/**Fetches PointData, returning a chromosome split hashmap<String, PointData[2]>. 
	 * PointData[0] is from the plus strand, 
	 * PointData[1] minus strand. Scores are not summed.
	 * One or the other may be null if single stranded.
	 * Returns null if no PointData files found.
	 * Wil load data if multiple of the same strand found.*/
	public static HashMap<String, PointData[]> fetchStrandedCombinePointData(File[] directories){
		//fetch hash map with stranded data
		HashMap<String, ArrayList<PointData>>[] combo = fetchStrandedPointDataNoMerge (directories);
		HashMap<String,ArrayList<PointData>> plusPDAL = combo[0];
		HashMap<String,ArrayList<PointData>> minusPDAL = combo[1];

		//combine PointData ArrayLists convert to PointData[]s
		HashMap<String,PointData[]> strandedPD = new HashMap<String,PointData[]>();
		HashSet<String> combine = new HashSet<String>();
		combine.addAll(plusPDAL.keySet());
		combine.addAll(minusPDAL.keySet());
		Iterator<String> it = combine.iterator();
		while (it.hasNext()){
			String chrom = it.next();
			//fetch ALs
			ArrayList<PointData> alPlus = plusPDAL.get(chrom);
			ArrayList<PointData> alMinus = minusPDAL.get(chrom);
			//merge
			PointData plus = null;
			if (alPlus != null) plus = PointData.combinePointData(alPlus, true);
			PointData minus = null;
			if (alMinus != null) minus = PointData.combinePointData(alMinus, true);
			//add to strandedPD
			strandedPD.put(chrom, new PointData[]{plus, minus});

		}
		return strandedPD;
	}



	/**Fetches, but doesn't load or merge, PointData, returning a chromosome split hashmap.
	 * Returns null if no PointData files found.*/
	public static HashMap<String, PointData[]> fetchStrandedPointData(File directory){
		//load hash map with stranded data
		HashMap<String,ArrayList<PointData>> strandedPointData = new HashMap<String,ArrayList<PointData>>();
		//fetch files		
		File[] files = IO.extractFiles(directory, ".bar");
		if (files == null || files.length == 0) {
			files = IO.extractFiles(directory, ".bar.zip");
			if (files == null || files.length == 0) {
				System.out.println("\t\tError: no Point Data files found in -> "+directory+"\n");
				return null;
			}	
		}
		//for each file 
		for (int j=0; j< files.length; j++){
			//make a PointData object
			PointData pd = new PointData (files[j], false);
			//add to hash
			String chr = pd.getInfo().getChromosome();
			if (strandedPointData.containsKey(chr)) strandedPointData.get(chr).add(pd);
			else {
				ArrayList<PointData> al = new ArrayList<PointData>();
				al.add(pd);
				strandedPointData.put(chr, al);
			}
		}
		if (strandedPointData.size() == 0) return null;
		//convert to PointData[]s
		HashMap<String,PointData[]> realPD = new HashMap<String,PointData[]>();
		Iterator<String> it = strandedPointData.keySet().iterator();
		while (it.hasNext()){
			String chrom = it.next();
			ArrayList<PointData> al = strandedPointData.get(chrom);
			PointData[] pd = new PointData[al.size()];
			al.toArray(pd);
			realPD.put(chrom, pd);
		}
		return realPD;
	}

	/**Fetches PointData from a directory. If same chromosome found, ie from a stranded set, it will merge the data.
	 * If a mergeDirectory is provided, the merged data will be written to the directory otherwise it is kept in memory.
	 * @param mergeDirectory set to null to keep merged data in memory, otherwise it will be written to disk.*/
	public static HashMap<String, PointData> fetchPointData (File directory, File mergeDirectory, boolean sumScores){

		//load hash map with stranded data
		HashMap<String,ArrayList<PointData>> strandedPointData = new HashMap<String, ArrayList<PointData>>();

		//fetch files
		File[] files = IO.extractFiles(directory, ".bar");
		if (files == null || files.length == 0) {
			files = IO.extractFiles(directory, ".bar.zip");
			if (files == null || files.length == 0) {
				System.out.println("\t\tError: no Point Data files found in -> "+directory+"\n");
				return null;
			}	
		}
		//for each file 
		for (int j=0; j< files.length; j++){
			//make a PointData object
			PointData pd = new PointData (files[j], false);
			//add to hash
			String chr = pd.getInfo().getChromosome();
			if (strandedPointData.containsKey(chr)) strandedPointData.get(chr).add(pd);
			else {
				ArrayList<PointData> al = new ArrayList<PointData>();
				al.add(pd);
				strandedPointData.put(chr, al);
			}
		}
		//merge strands?, should be two for each chromosome
		HashMap<String,PointData> pointData = new HashMap<String, PointData>();
		Iterator<String> it = strandedPointData.keySet().iterator();
		while (it.hasNext()){
			String chr = it.next();
			ArrayList<PointData> pdToMerge = strandedPointData.get(chr);
			PointData merged;
			int numDataSets = pdToMerge.size();			
			if (numDataSets == 1) merged = pdToMerge.get(0);
			else if (sumScores) merged = PointData.mergePointDataDepreciated(pdToMerge, true);
			else merged = PointData.combinePointData(pdToMerge, true);
			//write to directory?
			if (mergeDirectory != null){
				merged.getScores();
				merged.writePointData(mergeDirectory);
				merged.nullPositionScoreArrays();
			}
			pointData.put(chr, merged);
		}
		return pointData;
	}

	//getters setters
	public Info getInfo() {
		return info;
	}
	public void setInfo(Info info) {
		this.info = info;
	}
	/**Will load positions and scores from file if null.*/
	public int[] getPositions() {
		if (positions == null) loadPositionScores();
		return positions;
	}
	public void setPositions(int[] position) {
		this.positions = position;
		info.setNumberObservations(positions.length);
	}
	/**Will load scores and positions from file if null.*/
	public float[] getScores() {
		if (scores == null) loadPositionScores();
		return scores;
	}
	public void setScores(float[] score) {
		this.scores = score;
		info.setScoreTotal(Num.sumArrayReturnDouble(scores));
	}
	public BarParser getBarParser() {
		return barParser;
	}

	/*
	public static void main(String[] args){
		PointData pd = new PointData (new File ("/Users/nix/HCI/SigSequencers/HighResProfHistMethInHumGenomeCell2007/PolII_Point/chr4_+_.bar.zip"), false);
		System.out.println(pd.getInfo());
	}*/

}
