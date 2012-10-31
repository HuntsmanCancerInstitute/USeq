package trans.graphics;

import java.awt.*;
import java.awt.geom.*;
import java.awt.font.*;

import util.gen.*;



/**
 * Graph for {@link IntervalPlotter}.
 *
 */
public class Graph {
	private Line2D[] graphLines;
	private Line2D[] dataLines;
	private double h;			//height of graph
	private double hOffSet = 1;	//offSet from bottom of graph, only used if no negative values
	private double w;			//width of graph
	private double wOffSet = 2;	//offSet from left side of graph
	private double x;			//X-coor top left
	private double y;			//y-coor top left
	private double hPlusY;
	private int textSpacer = 1;	//spacer between labels and vertical graph bar
	private double[] scaledVals; 
	private double[] scaledPositions;
	private double[] values;				 //={2,5,3,8,10,22,21,24,19,18,5, 2, 3 , 2 , 6 ,1};
	private double[] positions; 			//={5,6,7,8,9, 12,14,15,16,17,25,26,30, 31, 32,33};
	private int numDataLines;
	private int numGraphLines;
	private Color graphLineColor = Color.WHITE;
	private Color dataLineColor;
	private double[] minMaxY;
	private double[] minMaxX;
	private double rangeY;
	private double rangeX;
	private double scalarY = 0;
	private double zeroY;
	private Rectangle2D labelBoundsTop;
	private Rectangle2D labelBoundsBottom;
	private String topLabel;
	private String bottomLabel;
	private float xCoorTopLabel;
	private float yCoorTopLabel;
	private float xCoorBottomLabel;
	private float yCoorBottomLabel;
	private double pixelsPerBase;
	
	public Graph (double[] values, double[] positions, double width, double height, Color dataLineColor){
		this.values = values;
		this.positions = positions;
		this.dataLineColor = dataLineColor;
		numDataLines = values.length;
		
		w = width;
		h = height;
		//initial settings can reset y by calling setMinMaxY with a combine data set to have their y the same scale
		minMaxY = Num.findMinMaxDoubleValues(values);
		setMinMaxY(minMaxY);
		minMaxX = Num.findMinMaxDoubleValues(positions);
		rangeX = minMaxX[1]-minMaxX[0];
	}
	
	/**Use to set a common y scale, pass it a minMax for the entire dataset.*/
	public void setMinMaxY(double[] minMaxY){
		this.minMaxY = minMaxY;
		rangeY = minMaxY[1]-minMaxY[0];
		//any negative values?
		if (minMaxY[0]<0) hOffSet = 0;
		//set labels
		if (minMaxY[1]<100)topLabel = Num.formatNumber(minMaxY[1],1);
		else topLabel = (int)Math.round(minMaxY[1]) +"";
		if (minMaxY[0]<100) bottomLabel = Num.formatNumber(minMaxY[0],1);
		else bottomLabel = (int)Math.round(minMaxY[0]) +"";
	}
	
	/**Returns bottom y coordinate of graph.*/
	public double makeGlyphs(double topLeftXCoor, double topLeftYCoor, double widthScalar, double heigthScalar){
		w *= widthScalar;
		h *= heigthScalar;
		x = topLeftXCoor;
		y = topLeftYCoor;
		scaleValuesPositions();
		makeGridLines();
		makeDataLines();
		labelBoundsTop = null;
		hPlusY = h+y;
		return hPlusY;
	}
	
	public void draw(Graphics2D g2, FontRenderContext fontContext, Font font){
		//draw graphlines
		g2.setPaint(graphLineColor);
		for (int i=0; i<numGraphLines; i++) g2.draw(graphLines[i]);
		
		//draw labels
		if (labelBoundsTop == null){
			labelBoundsTop = font.getStringBounds(topLabel, fontContext);
			labelBoundsBottom = font.getStringBounds(bottomLabel, fontContext);
			xCoorTopLabel = (float)(x-labelBoundsTop.getWidth()-textSpacer);
			yCoorTopLabel = (float)(y-labelBoundsTop.getY());
			xCoorBottomLabel = (float)(x-labelBoundsBottom.getWidth()-textSpacer);
			yCoorBottomLabel = (float)(h+y);
		}
		g2.drawString(topLabel,xCoorTopLabel, yCoorTopLabel);
		g2.drawString(bottomLabel,xCoorBottomLabel, yCoorBottomLabel);
		
		
		//draw dataLines
		g2.setPaint(dataLineColor);
		for (int i=0; i<numDataLines; i++) g2.draw(dataLines[i]);
	}
	
	public void makeGridLines(){
		zeroY = y+h;
		//any negative values?
		if (minMaxY[0]<0) zeroY += scalarY*minMaxY[0];
		//make graph grind lines
		graphLines = new Line2D[2];
		//left side
		graphLines[0] = new Line2D.Double(x,y, x,y+h);
		//zero line 
		graphLines[1] = new Line2D.Double(x,zeroY, x+w,zeroY);
		numGraphLines = graphLines.length;
	}
	public void scaleValuesPositions(){
		//get scalars and set base pix conversions
		if (rangeY!=0) scalarY = (h-hOffSet)/rangeY;
		pixelsPerBase = (w-wOffSet)/rangeX;
		
		//take min and subtract from all values then multiply by h or w/range to convert to 0 to h or w
			//if neg values present in values then want to not set base as 0
				double minY;
				if (minMaxY[0]<0) minY=0;
				else minY = minMaxY[0];
		int numVals = values.length;
		scaledVals = new double[numVals];
		scaledPositions = new double[numVals];
		
		for (int i=0; i<numDataLines; i++){
			scaledVals[i] = ((double)values[i]-minY)*scalarY;
			scaledPositions[i] = ((double)positions[i]-minMaxX[0])*pixelsPerBase;
		}
	}
	
	public void makeDataLines(){
		//make data lines
		dataLines = new Line2D[numDataLines];
		double modX = wOffSet+ x;
		double modY = y + h - hOffSet;
		//adjust for negative values? must shove either side of zero
		double finalX;
		if (minMaxY[0]<0) {
			modY += scalarY*minMaxY[0];	
			for (int i=0; i<numDataLines; i++){
				finalX = scaledPositions[i]+modX;
				if (scaledVals[i]>0) dataLines[i] =new Line2D.Double(finalX,modY-1,finalX,modY-scaledVals[i]);
				else if (scaledVals[i]<0) dataLines[i] =new Line2D.Double(finalX,modY+1,finalX,modY-scaledVals[i]);
				else dataLines[i] =new Line2D.Double(finalX,modY,finalX,modY-scaledVals[i]);
			}
		}
		//no just positive values
		else {
			for (int i=0; i<numDataLines; i++){
				finalX = scaledPositions[i]+modX;
				dataLines[i] =new Line2D.Double(finalX,modY,finalX,modY-scaledVals[i]);
			}
		}
	}

	public double getWOffSet() {
		return wOffSet;
	}
	public double getX() {
		return x;
	}
	public double getPixelsPerBase() {
		return pixelsPerBase;
	}
	public double[] getPositions() {
		return positions;
	}
}
