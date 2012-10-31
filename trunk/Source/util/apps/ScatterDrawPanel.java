package util.apps;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.util.HashSet;
import java.util.Iterator;
import javax.swing.JPanel;
import util.gen.*;


/**
 * Helper class for {@link ScatterPlot}, the guts of drawing the scatter plot.
 *
 */
public class ScatterDrawPanel extends JPanel{
	//fields
	private Line2D[] lines;
	private int numLines;
	private Rectangle2D[] recs;
	private int numPts;
	private boolean skipZeros = false;
	private boolean logValues = false;  //this doesn't graph properly
	
	public ScatterDrawPanel(String[] args){
		//args
		File first = new File (args[0]);
		File second = new File (args[1]);
		if (args.length > 2){
			skipZeros = true;
			System.out.print("\nIgnoring ");
		}
		
		//fetch float arrays
		float[] expA = null;
		float[] expB = null;
		try {
			expA = (float[])IO.fetchObject(first);
			expB = (float[])IO.fetchObject(second);	
		}catch(Exception e){
			try {
				expA = Num.convertToFloat((int[])IO.fetchObject(first));
				expB = Num.convertToFloat((int[])IO.fetchObject(second));
			}catch(Exception f){
				System.out.println("\nSomething is wrong with your files? Correct and restart.\n");
				f.printStackTrace();
			System.exit(1);
			}
		}
		//Strip zeros?
		if (skipZeros){
			int size = expA.length;
			float[][] cc = new float[2][];
			cc[0] = expA;
			cc[1] = expB;
			float[][] zeroed = Num.removeZeroValues(cc);
			expA = zeroed[0];
			expB = zeroed[1];
			size = size - expA.length;
			System.out.println(size+ " zero values!");
		}

		
		//make grid lines
		lines = new Line2D[3];
		Point2D ptTL = new Point2D.Double(5,5);
		Point2D ptTR = new Point2D.Double(690,5);
		Point2D ptBL = new Point2D.Double(5,690);
		Point2D ptBR = new Point2D.Double(690,690);
		lines[0] = new Line2D.Double(ptTL, ptBL);
		lines[1] = new Line2D.Double(ptBL, ptBR);
		lines[2] = new Line2D.Double(ptBL, ptTR);
		numLines = lines.length;
		//print stats
		System.out.println((float)Num.standardDeviation(expA)+ " Standard Deviation "+first.getName());
		System.out.println((float)Num.standardDeviation(expB)+ " Standard Deviation "+second.getName());
		System.out.println((float)PearsonCorrelation.correlationCoefficient(expA, expB)+ " Pearson Correlation Coefficient.");
		
		//log2 em?
		if (logValues){
			int num = expA.length;
			double log2 = Math.log(2);
			for (int i=0; i< num; i++){
				//return Math.log(number)/Math.log(2);
				expA[i] = (float)(Math.log(expA[i])/log2);
				expB[i] = (float)(Math.log(expB[i])/log2);
			}
		}
		
		//use a hash to clear dups x:y, round to ints
		int numValues = expA.length;
		if (numValues!=expB.length){
			System.out.println("\nThe files do not contain the same number of values!? Aborting...\n");
			System.exit(1);
		}
		System.out.println("Processing "+numValues*2+" values");
		HashSet hash = new HashSet(numValues/100);
		for (int i=0; i<numValues; i++) {
			hash.add(Math.round(expA[i])+":"+ Math.round(expB[i]));
		}
		
		//convert hash to float[]
		expA = new float[hash.size()];
		expB = new float[hash.size()];
		Iterator it = hash.iterator();
		String pair;
		String[] items;
		int counter =0;
		while (it.hasNext()){
			pair = (String)it.next();
			items = pair.split(":");
			expA[counter] = Float.parseFloat(items[0]);
			expB[counter] = Float.parseFloat(items[1]);
			counter++;
		}
		hash = null;
		
		//get max range and scalar
		float max = Num.findHighestFloat(expA);
		float maxB = Num.findHighestFloat(expB);
		
		if (maxB>max) max = maxB;
		double scalar = (double)max/685.0;
		
		System.out.println("Max value: "+max);

		
		//make pts
		numPts = expA.length;

		recs = new Rectangle2D[numPts];
		Point2D pt;
		double x;
		double y;
		for (int i=0; i<numPts; i++){
			x = ((double)expA[i]/scalar) + 3.5;
			y = 688.5- ((double)expB[i]/scalar);
			recs[i] = new Rectangle2D.Double(x,y,3,3);
		}
		
	}
	
	public void paintComponent(Graphics g){
		super.paintComponent(g);
		Graphics2D g2 = (Graphics2D)g;
		
		g2.setPaint(Color.BLACK);
		//draw rectangles
		for (int i=0; i<numPts; i++) g2.fill(recs[i]);
		
		g2.setPaint(Color.BLUE);
		//draw lines
		for (int i=0; i<numLines; i++) g2.draw(lines[i]);
	}
	public static double[] logEm(int[] ints){
		int num = ints.length;
		double[] doubles = new double[num];
		for (int i=0; i<num; i++){
			doubles[i] = Math.log(ints[i]);
		}
		return doubles;
	}
}