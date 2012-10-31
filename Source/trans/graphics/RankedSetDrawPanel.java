package trans.graphics;
import gata.main.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.font.FontRenderContext;
import java.awt.geom.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.*;
import javax.imageio.ImageIO;
import javax.swing.JPanel;
import util.gen.*;

/**
 * Panel with the meat for {@link RankedSetAnalysis}
 *
 */
public class RankedSetDrawPanel extends JPanel{
	
	//fields
	private Font font = new Font("Dialog",Font.PLAIN,10);
	private FontRenderContext context;
	private TextFrame textFrame;
	//dimensions
	private double maxSideWidth = 200;
	private double alignLeftX = 220;
	private double widthMiddle = 200;
	private double alignRightX = alignLeftX+ widthMiddle;
	private double spacer = 5;
	private double glyphHeight = 10;
	private double baseScalar;
	private double panelWidth = (2*alignLeftX) + widthMiddle;
	private double panelHeight;
	//regions
	private GenomicRegion[] one;
	private GenomicRegion[] two;
	//graphics
	Line2D[] connectorLines;
	Color rectangleIntersectingColor = Color.BLACK;
	Color rectangleNonIntersectingColor = Color.RED;
	
	public RankedSetDrawPanel(RankedSetAnalysis main){
		//set TextFrame reference
		textFrame = main.getTextFrame();
		
		//get regions
		one = main.getRegionsOne();
		two = main.getRegionsTwo();
		
		//find longest region and set scalar
		int longest = GenomicRegion.findBiggestGenomicRegion(one);
		int longestTwo = GenomicRegion.findBiggestGenomicRegion(two);
		if (longestTwo > longest) longest = longestTwo;
		baseScalar = maxSideWidth/longest;
		
		//make GenomicRegionGlyphs
		makeGlyphs();
		
		//make Rectangles in GRGs
		makeRectangles();
		
		//make connector lines between GRGs
		makeConnectorLines();
		
		//add listeners
		addMouseListener(new MouseHandler());
	}
	
	public void makeConnectorLines(){
		ArrayList lines = new ArrayList(one.length);
		//for each left side GenomicRegion
		for (int i=0; i<one.length; i++){
			Point2D onePoint = one[i].getGlyph().getPoint();
			//fetch ArrayList of intersecting GenomicRegions
			ArrayList overlappingRegions = one[i].getIntersectingRegions();
			for (int j=0; j< overlappingRegions.size(); j++){
				//make line using Points
				GenomicRegion two = (GenomicRegion)overlappingRegions.get(j);
				Point2D twoPoint = two.getGlyph().getPoint();
				lines.add(new Line2D.Double(onePoint, twoPoint));
			}
		}
		//convert ArrayList to Line2D[]
		connectorLines = new Line2D[lines.size()];
		lines.toArray(connectorLines);
	}
	
	public void makeRectangles(){
		Color color;
		//left side
		double runningY = spacer;
		for (int i=0; i<one.length; i++){
			GenomicRegionGlyph g = one[i].getGlyph();
			if (one[i].getIntersectingRegions().size() !=0 ) color = rectangleIntersectingColor;
			else color = rectangleNonIntersectingColor;
			g.makeRectangleAndPoint(false, alignLeftX, runningY, baseScalar, glyphHeight,color);
			runningY = runningY + glyphHeight + spacer;
		}
		//right side
		double runningYTwo = spacer;
		for (int i=0; i<two.length; i++){
			GenomicRegionGlyph g = two[i].getGlyph();
			if (two[i].getIntersectingRegions().size() !=0 ) color = rectangleIntersectingColor;
			else color = rectangleNonIntersectingColor;
			g.makeRectangleAndPoint(true, alignRightX, runningYTwo, baseScalar, glyphHeight, color);
			runningYTwo = runningYTwo + glyphHeight + spacer;
		}
		//set panel height
		if (runningY >= runningYTwo ) panelHeight = runningY+spacer;
		else panelHeight = runningYTwo = runningYTwo + spacer;
	}
	
	/**Makes and places a GenomicRegionGlyph in each region using the rank number as the label.*/
	public void makeGlyphs(){
		for (int i=0; i<one.length; i++){
			GenomicRegionGlyph g = new GenomicRegionGlyph((i+1)+"", one[i].getStop()-one[i].getStart());
			one[i].setGlyph(g);
		}
		for (int i=0; i<two.length; i++){
			GenomicRegionGlyph g = new GenomicRegionGlyph((i+1)+"", two[i].getStop()-two[i].getStart());
			two[i].setGlyph(g);
		}
	}
	
	public void paintComponent(Graphics g){
		super.paintComponent(g);
		Graphics2D g2 = (Graphics2D)g;
		g2.setFont(font);
		context = g2.getFontRenderContext();
		g2.setColor(Color.BLACK);
		//draw rectangles
		for (int i=0; i<one.length; i++) {
			one[i].getGlyph().draw(g2,context, font);
		}
		for (int i=0; i<two.length; i++) {
			two[i].getGlyph().draw(g2,context, font);
		}
		//draw lines
		g2.setColor(Color.BLUE);
		for (int i=0; i<connectorLines.length; i++){
			g2.draw(connectorLines[i]);
		}
		
	}

	public double getPanelHeight() {
		return panelHeight;
	}

	public double getPanelWidth() {
		return panelWidth;
	}
	
	private class MouseHandler extends MouseAdapter {
		public void mousePressed(MouseEvent event){
			//where was it pressed
			Point2D mousePressedPoint = event.getPoint();
			//find which rectangle
			String side = "L";
			GenomicRegion clicked = findRegion(one, mousePressedPoint);
			if (clicked == null) {
				clicked = findRegion(two, mousePressedPoint);
				side = "R";
			}
			if (clicked !=null){
				//write out info
				textFrame.append(side+formattedInfo(clicked));
			}
		}	
	}
	
	public String formattedInfo(GenomicRegion g){
		StringBuffer sb = new StringBuffer();
		sb.append(g.getGlyph().getLabel()); //rank
		sb.append(" ");
		sb.append(g.getScore());
		sb.append(" ");
		sb.append(g.getChromosome());
		sb.append(":");
		sb.append(g.getStart());
		sb.append("-");
		sb.append(g.getStop());
		sb.append(" ");
		if (g.getNotes()!=null) sb.append(g.getNotes());
		sb.append("\n");
		return sb.toString();
	}

	/**Scans an array of GenomicRegions Rectangles to see if the Point is present, ie they clicked it.*/
	public GenomicRegion findRegion(GenomicRegion[] regions, Point2D point){
		for (int i=0; i< regions.length; i++){
			if (regions[i].getGlyph().getRectangle().contains(point)) return regions[i];
		}
		return null;
	}
	
	public  void saveBufferedImage(double scale, File file, boolean antiAlias){
		BufferedImage bufferedImage=null;
		Graphics2D g2=null;
		double width = (panelWidth+10)*scale; 
		double height = (panelHeight+10)*scale;
		
		bufferedImage = new BufferedImage ((int)Math.round(width), (int)Math.round(height), BufferedImage.TYPE_INT_ARGB);
		g2 = bufferedImage.createGraphics();
		g2.addRenderingHints(GATAUtil.fetchRenderingHints());
		//kill antialiasing to make sharp lines and text for computer display, not good for printed figures
		if (antiAlias==false){
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
			g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF);
		}
		//set scale
		if (scale >1 ) g2.scale(scale,scale);
		//draw background
		g2.setColor(Color.WHITE);
		g2.fill(new Rectangle2D.Double(0,0,width+50, height+50));
		//set font
		g2.setFont(font);
		context = g2.getFontRenderContext();
		g2.setColor(Color.BLACK);
		//draw rectangles
		for (int i=0; i<one.length; i++) {
			one[i].getGlyph().draw(g2,context, font);
		}
		for (int i=0; i<two.length; i++) {
			two[i].getGlyph().draw(g2,context, font);
		}
		//draw lines
		g2.setColor(Color.BLUE);
		for (int i=0; i<connectorLines.length; i++){
			g2.draw(connectorLines[i]);
		}
		
		try{
			ImageIO.write(bufferedImage, "PNG", file);
		} catch (IOException e){
			e.printStackTrace();
		}
	}
}