package gata.geneGlyphs;
import gata.main.*;

import java.util.*;
import java.awt.geom.*;

/**
 * @author Nix
Makes a Scale bar or ruler with intervals of between 10 and 100 pixels in length at 1, 10, 100 bp or kb or mb sizes
based on pixPerNt.
 */
public class ScaleBar {
	//fields
	private ArrayList scaleBarLines = new ArrayList();
	private ArrayList scaleRulerLines = new ArrayList();
	private String scaleBarLabel;
	private double scaleBarWidth;
	private AnnoSpecParams annoParams;
	private GATAParams gataParams;
	private ArrayList labels; //master list of all lables
	

	public ScaleBar(AnnoSpecParams annoParams, GATAParams gataParams) {
		this.annoParams = annoParams;
		this.gataParams = gataParams;
		labels = annoParams.getLabels();
		annoParams.setScaleBar(this);
		
	}
	/**Makes a scale ruler with intervals between 10 and 100 pixels in length at 1, 10, 100 bp or kb or mb sizes,
		change the scaleRulerStartStopX and scaleRulerY coordinates to relocate the bar*/
	public void makeScaleRuler(){
		scaleRulerLines.clear();
		double pixPerNuc =gataParams.getPixPerNt();
		double[] scaleRulerStartStop= new double[] {1,(double)gataParams.getZoomedWidth()};
		
		double scaleRulerY = annoParams.getScaleRulerY();
	
		String unit;
		double scalar = 1;
		for (double x = 0.01; x < 10000000; x *= 10) {
			double pix = pixPerNuc * x;
			if (pix > 10 && pix <= 100) {
				if (x < 1000)
					unit = "bp";
				else if (x >= 1000 && x <= 100000) {
					unit = "kb";
					scalar = 1000;
				} else {
					unit = "mb";
					scalar = 1000000;
				}
				x = x / scalar;
				//main line
				scaleRulerLines.add(
				new Line2D.Double(
					scaleRulerStartStop[0],
					scaleRulerY,
					scaleRulerStartStop[1],
					scaleRulerY));
				double counter = 0;
				//double halfPix = pix/2;
				double yMinus = scaleRulerY-1;
				double yPlus = scaleRulerY+1;
				while (counter <= scaleRulerStartStop[1]) {
				//caps
					double xCorr = scaleRulerStartStop[0] + counter;
					scaleRulerLines.add(
					new Line2D.Double(xCorr, yMinus, xCorr, yPlus));
					counter += pix;
				}
				//set global variables
				scaleBarLabel = (int) x + " " + unit;
				scaleBarWidth = pix;
				//make label
				double xMidPtPix = pix/2 +10;
				labels.add(scaleBarLabel);
				labels.add(new double[]{xMidPtPix,scaleRulerY+10}); //fudge factor to shove down
				break;
			}
		}
	}
	//Not using in final GATA
	/**Makes a scale bar between 10 and 100 pixels in length at 1, 10, 100 bp or kb or mb sizes,
	change the scaleBarX and scaleBarY coordinates to relocate the bar*/
	public void makeScaleBar(){
		scaleBarLines.clear();
		double pixPerNuc = gataParams.getPixPerNt();
		double scaleBarX = annoParams.getScaleBarX();
		double scaleBarY = annoParams.getScaleBarY();
		String unit;
		double scalar = 1;
		for (double x = 0.01; x < 10000000; x *= 10) {
			double pix = pixPerNuc * x;
			if (pix > 10 && pix <= 100) {
				if (x < 1000)
					unit = "bp";
				else if (x >= 1000 && x <= 100000) {
					unit = "kb";
					scalar = 1000;
				} else {
					unit = "mb";
					scalar = 1000000;
				}
				x = x / scalar;
				scaleBarLines.add(
					new Line2D.Double(
						scaleBarX,
						scaleBarY,
						scaleBarX + pix,
						scaleBarY));
				scaleBarLines.add(
					new Line2D.Double(
						scaleBarX,
						scaleBarY - 2,
						scaleBarX,
						scaleBarY + 2));
				scaleBarLines.add(
					new Line2D.Double(
						scaleBarX + pix,
						scaleBarY - 2,
						scaleBarX + pix,
						scaleBarY + 2));
				//set global variables
				scaleBarLabel = (int) x + " " + unit;
				scaleBarWidth = pix;
				//make label
				double xMidPtPix = (pix/2) + scaleBarX;
				labels.add(scaleBarLabel);
				labels.add(new double[]{xMidPtPix,scaleBarY+10}); //fudge factor to shove down
				break;
			}
		}
	}
	public String getScaleBarLabel() {
		return scaleBarLabel;
	}
	public double getScaleBarWidth() {
		return scaleBarWidth;
	}
	public Line2D.Double[] getScaleBarLines() {
		return GATAUtil.ALToLine2DArray(scaleBarLines);
	}
	public Line2D.Double[] getScaleRulerLines(){
		return GATAUtil.ALToLine2DArray(scaleRulerLines);
	}
	public Line2D.Double[] getScaleLines(){
		ArrayList x = (ArrayList)scaleRulerLines.clone();
		x.addAll(scaleBarLines);
		return GATAUtil.ALToLine2DArray(x);
	}
	/**Used in hiding scale ruler bar*/
	public void clearScaleLineArrays(){
		scaleBarLines.clear();
		scaleRulerLines.clear();
	}
}
