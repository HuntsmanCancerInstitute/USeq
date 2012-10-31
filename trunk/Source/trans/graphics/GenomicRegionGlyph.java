package trans.graphics;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.font.FontRenderContext;
import java.awt.geom.*;

public class GenomicRegionGlyph {
	
	//fields
	private int bpLength;
	private Rectangle2D rectangle; 
	private Point2D point;
	private String label;
	private Rectangle2D labelBounds;
	private float xCoorLabel;
	private float yCoorLabel;
	private static double labelGlyphSpacer = 5;
	private static double yLabelFudge = 4;
	private boolean rightSide;
	private Color rectangleColor;
	
	//constructor
	public GenomicRegionGlyph(String label, int bpLength){
		this.label = label;
		this.bpLength = bpLength;
	}
	
	//methods
	public void makeRectangleAndPoint(boolean rightSide, double topXCoor, double topYCoor, double widthScalar, double height,  Color rectangleColor){
		//make rectangle
		double pixelWidth = bpLength*widthScalar;
		double topLeftXCoor;
		if (rightSide) topLeftXCoor = topXCoor;
		else topLeftXCoor = topXCoor - pixelWidth;
		rectangle = new Rectangle2D.Double(topLeftXCoor, topYCoor, pixelWidth, height);
		//make point connector point
		double y = height/2 + topYCoor;
		point = new Point2D.Double(topXCoor, y); 
		//set boolean for drawing label
		this.rightSide = rightSide;
		this.rectangleColor = rectangleColor;
	}
	
	public void draw(Graphics2D g2, FontRenderContext fontContext, Font font){
		//draw label
		if (labelBounds == null){
			labelBounds = font.getStringBounds(label, fontContext);
			//y
			yCoorLabel = new Double(rectangle.getCenterY() + yLabelFudge).floatValue();
			//x
			double widthLabel = labelBounds.getWidth();
			if (rightSide) xCoorLabel = new Double(rectangle.getX() + rectangle.getWidth() + labelGlyphSpacer).floatValue();
			else xCoorLabel = new Double(rectangle.getX()- widthLabel - labelGlyphSpacer).floatValue();	
		}
		g2.setColor(Color.BLACK);
		g2.drawString(label,xCoorLabel, yCoorLabel);
		//draw rectangle
		g2.setColor(rectangleColor);
		g2.draw(rectangle);
	}

	public Point2D getPoint() {
		return point;
	}

	public Rectangle2D getRectangle() {
		return rectangle;
	}

	public String getLabel() {
		return label;
	}
	
	
	
	
	
}
