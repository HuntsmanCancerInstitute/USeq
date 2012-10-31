package trans.graphics;

import java.awt.*;
import java.awt.geom.*;

/**
 * Box for {@link IntervalPlotter}.
 *
 */
public class Box {
	//fields
	private Rectangle2D box;
	private Color color;
	private double baseStart;
	private double baseStop;
	private double baseFirstPosition;
	private double height = 10;
	
	public Box(double baseStart, double baseStop, double baseFirstPosition, Color color){
		this.baseStart = baseStart;
		this.baseStop = baseStop;
		this.baseFirstPosition = baseFirstPosition;
		this.color = color;
	}
	
	public double makeGlyphs(Graph graph, double y){
		double pixelsPerBase = graph.getPixelsPerBase();
		double xBase = graph.getX() + graph.getWOffSet();
		double xStart = (baseStart-baseFirstPosition)*pixelsPerBase + xBase;
		double xStop = (baseStop-baseFirstPosition)*pixelsPerBase + xBase;
		box = new Rectangle2D.Double(xStart,y, 1+xStop-xStart, height);
		return y+height;
	}
	
	public void draw(Graphics2D g2){
		g2.setPaint(color);
		g2.fill(box);
	}

}
