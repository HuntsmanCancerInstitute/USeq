package edu.utah.seq.analysis;

import edu.utah.seq.data.Point;
import edu.utah.seq.data.PointData;

/**Container for bis seq data*/
public class BisSeqRegion {

		//fields
		Point[] nonConPlus;
		Point[] conPlus;
		Point[] nonConMinus;
		Point[] conMinus;
		boolean dataPresent = false;
			
		public BisSeqRegion(int start, int stop, PointData nonConvertedMergedChromPlus, PointData convertedMergedChromPlus, PointData nonConvertedMergedChromMinus, PointData convertedMergedChromMinus){
			//fetch arrays, might be null
			nonConPlus = nonConvertedMergedChromPlus.fetchPoints(start, stop);
			conPlus = convertedMergedChromPlus.fetchPoints(start, stop);
			nonConMinus = nonConvertedMergedChromMinus.fetchPoints(start, stop);
			conMinus = convertedMergedChromMinus.fetchPoints(start, stop);	
			
			if (nonConPlus !=null || conPlus !=null || nonConMinus!=null || conMinus!=null) dataPresent =true;
			else return;		
			
			//reset positions based on start
			if (nonConPlus!=null) Point.subtractConstantFromPositions(nonConPlus, start);
			if (conPlus!=null) Point.subtractConstantFromPositions(conPlus, start);
			if (nonConMinus!=null) Point.subtractConstantFromPositions(nonConMinus, start);
			if (conMinus!=null) Point.subtractConstantFromPositions(conMinus, start);
		}
		
		public void printMe(){
			float numNonCon = 0;
			float numCon =0;
			if (nonConPlus!=null) for (Point p: nonConPlus) numNonCon += p.getScore();
			if (conPlus!=null) for (Point p: conPlus) numCon += p.getScore();
			if (nonConMinus!=null) for (Point p: nonConMinus) numNonCon += p.getScore();
			if (conMinus!=null) for (Point p: conMinus) numCon += p.getScore();
			System.out.println("PrintMe "+numNonCon+"\t"+numCon);
		}

		public void merge(BisSeqRegion pts) {	
			if (dataPresent == false || pts.dataPresent == false) return;
			
			//this checks for nulls
			nonConPlus = Point.mergePoints(nonConPlus, pts.nonConPlus);
			conPlus = Point.mergePoints(conPlus, pts.conPlus);
			nonConMinus = Point.mergePoints(nonConMinus, pts.nonConMinus);
			conMinus = Point.mergePoints(conMinus, pts.conMinus);
			dataPresent = true;
		}

		public void scale(double scalar) {
			if (dataPresent == false) return;
			if (nonConPlus!=null) Point.scalePositions(nonConPlus, scalar);
			if (conPlus!=null) Point.scalePositions(conPlus, scalar);
			if (nonConMinus!=null)Point.scalePositions(nonConMinus, scalar);
			if (conMinus!=null) Point.scalePositions(conMinus, scalar);
		}
		
		public void sumIdentialPositions(){
			if (dataPresent == false) return;
			if (nonConPlus!=null) Point.sumIdenticalPositionScores(nonConPlus);
			if (conPlus!=null) Point.sumIdenticalPositionScores(conPlus);
			if (nonConMinus!=null) Point.sumIdenticalPositionScores(nonConMinus);
			if (conMinus!=null) Point.sumIdenticalPositionScores(conMinus);
		}
		
		public void invert(int lastBase){
			if (dataPresent == false) return;
			if (nonConPlus!=null) nonConPlus = Point.invert(nonConPlus, lastBase);
			if (conPlus!=null) conPlus = Point.invert(conPlus, lastBase);
			if (nonConMinus!=null) nonConMinus = Point.invert(nonConMinus, lastBase);
			if (conMinus!=null) conMinus = Point.invert(conMinus, lastBase);
		}
		
		public static String printCounts(float[][] counts){
			StringBuffer sb = new StringBuffer();
			for (int i=0; i<counts[0].length; i++){
				sb.append((i+1)+"\t"+counts[0][i]+"\t"+counts[1][i]+"\t"+counts[2][i]+"\t"+counts[3][i]+"\n");
			}
			return sb.toString();
		}
		
		/**Returns float[4][maxSize]: nonConPlus[0][] conPlus[1][] nonConMinus[2][] conMinus[3][]*/
		public float[][] fetchCounts(int maxSize){

			//Make PD 
			PointData pdNonConPlus = null;
			PointData pdConPlus = null;
			PointData pdNonConMinus = null;
			PointData pdConMinus = null;
			if (nonConPlus!=null) pdNonConPlus = Point.extractPositionScores(nonConPlus);
			if (conPlus!=null) pdConPlus = Point.extractPositionScores(conPlus);
			if (nonConMinus!=null) pdNonConMinus = Point.extractPositionScores(nonConMinus);
			if (conMinus!=null) pdConMinus = Point.extractPositionScores(conMinus);
			
			int lastBase = maxSize+1;
			float[] nonConPlusCounts = new float[lastBase];
			float[] conPlusCounts = new float[lastBase];
			float[] nonConMinusCounts = new float[lastBase];
			float[] conMinusCounts = new float[lastBase];
			
			//for each base
			for (int i=0; i< lastBase; i++){
				int start = i;
				int stop = i+1;
				Point[] pts = null;
				//counts
				float ncp = 0;
				float cp = 0;
				float ncm = 0;
				float cm = 0;
				
				//nonConPlus
				if (pdNonConPlus!=null) pts = pdNonConPlus.fetchPoints(start, stop);
				else pts = null;
				if (pts!=null) ncp = pts[0].getScore();
				
				//nonConPlus
				if (pdConPlus!=null) pts = pdConPlus.fetchPoints(start, stop);
				else pts = null;
				if (pts!=null) cp = pts[0].getScore();
				
				//nonConPlus
				if (pdNonConMinus!=null) pts = pdNonConMinus.fetchPoints(start, stop);
				else pts = null;
				if (pts!=null) ncm = pts[0].getScore();
				
				//nonConPlus
				if (pdConMinus!=null) pts = pdConMinus.fetchPoints(start, stop);
				else pts = null;
				if (pts!=null) cm = pts[0].getScore();
				
				nonConPlusCounts[i] = ncp;
				conPlusCounts[i] = cp;
				nonConMinusCounts[i] = ncm;
				conMinusCounts[i] = cm;
			}
			
			return new float[][]{nonConPlusCounts, conPlusCounts, nonConMinusCounts, conMinusCounts};
		}
		
		public int fetchLastBase(){
			int max = lastBase(nonConPlus);
			if (lastBase(conPlus) > max) max = lastBase(conPlus);
			if (lastBase(nonConMinus) > max) max = lastBase(nonConMinus);
			if (lastBase(conMinus) > max) max = lastBase(conMinus);
			return max;
		}
		
		public int lastBase(Point[] pt){
			if (pt == null) return -1;
			return pt[pt.length-1].getPosition();
		}
	}