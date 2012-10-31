package edu.utah.seq.data;

import java.util.*;
import java.io.*;
import util.gen.*;

public class Point implements Serializable{
	//fields
	int position;
	float score;
	
	//constructor
	public Point (int position, float score){
		this.position = position;
		this.score = score;
	}
	
	/**@return a near empty PointData object loaded with just positions and scores.*/
	public static PointData extractPositionScores(Point[] p){
		int[] positions = new int[p.length];
		float[] scores = new float[p.length];
		for (int i=0; i< p.length; i++){
			positions[i] = p[i].position;
			scores[i] = p[i].score;
		}
		PointData pd = new PointData();
		pd.setPositions(positions);
		pd.setScores(scores);
		return pd;
	}
	
	/**@return a near empty PointData object loaded with just positions and scores.*/
	public static PointData extractPositionScores(ArrayList<Point> al){
		Point[] p = new Point[al.size()];
		al.toArray(p);
		return extractPositionScores(p);
	}
	

	public static Point[] makePoints(int[] positions, float[] scores){
		Point[] p = new Point[positions.length];
		for (int i=0; i< p.length; i++) p[i] = new Point(positions[i], scores[i]);
		return p;
	}
	
	public static ArrayList<Point> makePointsAL(int[] positions, float[] scores){
		ArrayList<Point> p = new ArrayList<Point>(positions.length);
		for (int i=0; i< positions.length; i++) p.add(new Point(positions[i], scores[i]));
		return p;
	}
	
	public static Point[] makePoints(PointData pd){
		return makePoints(pd.getPositions(), pd.getScores());
	}
	/**Makes a Point[] using the index position in the al as the position.*/
	public static Point[] makePoints(ArrayList<Float> al){
		int num = al.size();
		Point[] p = new Point[num];
		for (int i=0; i< num; i++){
			p[i] = new Point(i,al.get(i).floatValue());
		}
		return p;
	}
	
	/**Applies the Benjamini & Hochberg FDR correction to the array. Use this method to correct threholded pvalues by setting the offset.
	 * @param sortedMin10Log10PVals -10Log10(pval) transformed and sorted in ascending order
	 * @param offset Number of pvalues with poor significance not included in the array
	 * @return nada, modifies the original array*/
	public static void benjaminiHochbergCorrect(Point[] sortedMin10Log10PVals, long offset){
		double num = sortedMin10Log10PVals.length + offset;
		double offsetDouble = (double)offset;
		double prior = 1;
		float maxFDR = 0;
		for (int i=1; i< sortedMin10Log10PVals.length; i++){
			double val = Num.antiNeg10log10(sortedMin10Log10PVals[i].score);
			val = val * num / (num-i-offsetDouble);
			if(val < prior) prior = val; 
			else val = prior;
			sortedMin10Log10PVals[i].score = Num.minus10log10Float(val);
			if (sortedMin10Log10PVals[i].score > maxFDR && Float.isInfinite(sortedMin10Log10PVals[i].score) == false) maxFDR = sortedMin10Log10PVals[i].score;
		}
		//look for Infinity and replace with max
		for (int i=1; i< sortedMin10Log10PVals.length; i++){
			if (Float.isInfinite(sortedMin10Log10PVals[i].score)) sortedMin10Log10PVals[i].score = maxFDR;
		}
		//check first value, if larger than second, replace with second
		if (sortedMin10Log10PVals[0].score > sortedMin10Log10PVals[1].score) sortedMin10Log10PVals[0].score = sortedMin10Log10PVals[1].score;
		
	}
	
	public static int[] extractPositions (Point[] pts){
		int[] p = new int[pts.length];
		for (int i=0; i< pts.length; i++) p[i] = pts[i].position;
		return p;
	}

	public static float[] extractScores (Point[] pts){
		float[] p = new float[pts.length];
		for (int i=0; i< pts.length; i++) p[i] = pts[i].score;
		return p;
	}
	
	/**Sums the scores of Points with the same position.*/
	public static Point[] sumIdenticalPositionScores(Point[] sortedPoints){
		ArrayList<Point> savedPoints = new ArrayList<Point>();
		Point p = sortedPoints[0];
		int pPos = p.getPosition();
		for (int i=1; i< sortedPoints.length; i++){
			int testPos = sortedPoints[i].getPosition();
			if (testPos == pPos) p.setScore(p.getScore()+sortedPoints[i].getScore());
			else {
				savedPoints.add(p);
				p = sortedPoints[i];
				pPos = p.getPosition();
			}
		}
		//add last
		savedPoints.add(p);
		//convert and assign
		Point[] pts = new Point[savedPoints.size()];
		savedPoints.toArray(pts);
		return pts;
	}
	
	/**Subtracts the constant from the position of each Point. Sets those < 0 to 0.*/
	public static void subtractConstantFromPositions(Point[] pts, int constant){
		for (Point pt : pts){
			int pos = pt.getPosition()-constant;
			if (pos < 0) pos = 0;
			pt.setPosition(pos);
		}
	}
	
	/**Multiplies each position by the scalar and rounds up.*/
	public static void scalePositions(Point[] pts, double scalar){
		for (Point pt : pts){
			double pos = pt.getPosition();
			int newPos = (int)Math.round(pos*scalar);
			pt.setPosition(newPos);
		}
	}
	
	/**Merges and sorts by position.*/
	public static Point[] mergePoints (Point[] p1, Point[] p2){
		//check for nulls
		if (p2 == null || p1 == null) {
			if (p1 != null) return p1;
			if (p2 != null) return p2;
			return null;
		}
		Point[] merged = new Point[p1.length+p2.length];
		System.arraycopy(p1, 0, merged, 0, p1.length);
		System.arraycopy(p2, 0, merged, p1.length, p2.length);
		Arrays.sort(merged, new ComparatorPointPosition());
		return sumIdenticalPositionScores(merged);
	}
	
	/**Splits the Point[] in half.*/
	public static Point[][] split(Point[] p){
		Misc.randomize(p, 0);
		int half = p.length/2;
		Point[] one = new Point[half];
		Point[] two = new Point[half];
		System.arraycopy(p, 0, one, 0, half);
		System.arraycopy(p, half, two, 0, half);
		ComparatorPointPosition comp = new ComparatorPointPosition();
		Arrays.sort(one, comp);
		Arrays.sort(two,comp);
		return new Point[][]{one,two};
	}
	public void incrementScore(float add){
		score += add;
	}
	public String toString(){
		return position+"\t"+score;
	}

	public int getPosition() {
		return position;
	}

	public void setPosition(int position) {
		this.position = position;
	}

	public float getScore() {
		return score;
	}

	public void setScore(float score) {
		this.score = score;
	}

	/**Subtracts the each position from the lastBase and inverts the array.*/
	public static Point[] invert(Point[] pts, int lastBase) {		
		int index = pts.length -1;	
		Point[] inverse = new Point[pts.length];
		for (int i=0; i< pts.length; i++) {
			pts[i].setPosition(lastBase - pts[i].getPosition());
			inverse[index] = pts[i];
			index--;
		}
		return inverse;
	}
}
