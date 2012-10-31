package gata.geneGlyphs;

import gata.main.*;
import gata.plotter.*;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.awt.event.*;
import java.util.*;
import java.awt.image.*;

/**
 * @author nix
 * Swing Panel holding the gene annotation shapes/ glyphs.
 */
public class GlyphPanel extends JPanel {
	//fields
	private GATAParams gataParams;
	private AnnoSpecParams annoParams;
	private Graphics2D g2;
	private ScaleBar scaleBar;
	private GenericGlyph[][] genericGlyphs;
	private GeneRepGlyph[] geneGlyphs;
	private Console console;
	private JViewport viewPort;

	//private ScaleBar scaleBar;
	private BasicStroke thinStroke;
	private BasicStroke normalStroke;
	private BasicStroke fatStroke;
	private BasicStroke dashStroke;
	private BasicStroke scaleBarStroke = new BasicStroke(1,BasicStroke.CAP_BUTT,BasicStroke.JOIN_BEVEL);
	private BasicStroke[] trackStrokes;
	private Font font;

	//fields
	private Color DNAColor;
	private Color RNAColor;
	private Color proteinColor;
	private Color labelColor;
	private Color scaleLineColor;
	private Line2D.Double[] normalDNALines;
	private Line2D.Double[] thinDNALines;
	private Line2D.Double[] fatRNALines;
	private Line2D.Double[] thinRNALines;
	private Line2D.Double[] dottedProteinLines;
	private Line2D.Double[] thinProteinLines;
	private Line2D.Double[] scaleLines;
	private ArrayList labels;
	private ArrayList[] trackLabels;
	private boolean anyGenericGlyphs = false;
	private boolean[] trackVis;
	private boolean[] trackLabelVis;
	private int numberTracks;

	public GlyphPanel(AnnoSpecParams annoParams, GATAParams gataParams) {
		this.annoParams = annoParams;
		this.gataParams = gataParams;
		geneGlyphs = annoParams.getGeneGlyphs();
		scaleBar = annoParams.getScaleBar();
		anyGenericGlyphs = annoParams.getAnyGenericGlyphs();
		setTrackParams();
		
		//to listen for mouse clicks
		addMouseListener(new MouseHandler());
	}
	public void setConsole(Console console){
		this.console = console;
	}
	
	public void setTrackParams(){
		trackVis = annoParams.getTrackVis();
		trackLabelVis = annoParams.getTrackLabelVis();
		if (trackVis !=null) numberTracks = trackVis.length;
	}
	
	/**Fetches glyphParams to speed up drawing, call this to initialize variables forGlyphPanel*/
	public void setDrawParams() {
		thinStroke = gataParams.getThinStroke();
		normalStroke = gataParams.getNormalStroke();
		fatStroke = gataParams.getFatStroke();
		dashStroke = gataParams.getDashStroke();
		font = gataParams.getFont();
		DNAColor = gataParams.getDNAColor();
		RNAColor = gataParams.getRNAColor();
		proteinColor = gataParams.getProteinColor();
		labelColor = gataParams.getLabelColor();
		scaleLineColor = annoParams.getScaleLineColor();
		if (anyGenericGlyphs){
			trackStrokes = annoParams.getTrackStrokes();
		}
	}
	/**Call this to set or reset the array of Line2D[]'s representing the glyphs.
	 * Also resets labels and scale bar.  These will change with zooming/ resizing.*/
	public void setLine2DArrays() {
		normalDNALines = GATAUtil.ALToLine2DArray(annoParams.getNormalDNALines());
		thinDNALines = GATAUtil.ALToLine2DArray(annoParams.getThinDNALines());
		fatRNALines = GATAUtil.ALToLine2DArray(annoParams.getFatRNALines());
		thinRNALines = GATAUtil.ALToLine2DArray(annoParams.getThinRNALines());
		dottedProteinLines =
		GATAUtil.ALToLine2DArray(annoParams.getDottedProteinLines());
		thinProteinLines = GATAUtil.ALToLine2DArray(annoParams.getThinProteinLines());
		scaleLines = scaleBar.getScaleLines();
		labels = annoParams.getLabels();
		if (anyGenericGlyphs) {
			genericGlyphs = annoParams.getTrackedGenericGlyphs();
			trackLabels = annoParams.getTrackLabels();
		} 
	}

	public void drawGlyphPanel(Graphics2D g2){
		//draw labels
		//Couldn't extract these calculations from paintComponent method since one
		//needs the Graphics2D to call context and it can't be instantiated!  Bastards!
		g2.setFont(font);
		g2.setColor(labelColor);
		GATAUtil.drawLabels(labels, font, g2);
		
		//track labels
		if (anyGenericGlyphs){
			for (int i = 0; i < numberTracks; i++) {
				if (trackLabelVis[i]==false) continue;
				GATAUtil.drawLabels(trackLabels[i], font, g2);
			}
		}
		//draw scaleBar/Ruler lines
		g2.setStroke(scaleBarStroke);
		g2.setColor(scaleLineColor);
		GATAUtil.drawLines(scaleLines, g2);

		//draw GenericGlyphs
		if (anyGenericGlyphs) {
			for (int i = 0; i < numberTracks; i++) {
				if (trackVis[i]==false) continue; 
				g2.setStroke(trackStrokes[i]);
				int len2 = genericGlyphs[i].length;
				for (int j = 0; j < len2; j++) {
					g2.setColor(genericGlyphs[i][j].getColor());
					g2.draw(genericGlyphs[i][j].getPath());
				}
			}
		}

		//DNA
		g2.setColor(DNAColor);
		g2.setStroke(normalStroke);
		GATAUtil.drawLines(normalDNALines, g2);
		g2.setStroke(thinStroke);
		GATAUtil.drawLines(thinDNALines, g2);
		//Protein
		g2.setColor(proteinColor);
		GATAUtil.drawLines(thinProteinLines, g2);
		g2.setStroke(dashStroke);
		GATAUtil.drawLines(dottedProteinLines, g2);
		//RNA
		g2.setStroke(thinStroke);
		g2.setColor(RNAColor);
		GATAUtil.drawLines(thinRNALines, g2);
		g2.setStroke(fatStroke);
		GATAUtil.drawLines(fatRNALines, g2);
		
	}

	/**Watch order of drawing, optimize to avoid changing stroke*/
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		g2 = (Graphics2D) g;
		drawGlyphPanel(g2);
	}

	public  BufferedImage makeBufferedImage(double scale, boolean saveVisible){
		BufferedImage bufferedImage=null;
		Graphics2D g2buff=null;
		double width = getWidth()*scale;
		double height = (annoParams.getBiggestGlyphPanelY()-annoParams.getSmallestGlyphPanelY())*scale;
		
		//create image matching what is shown
		if (saveVisible){
			int partialWidth = (int)Math.round(((double)viewPort.getWidth())*scale);
			int partialHeight = (int)Math.round(((double)viewPort.getHeight())*scale);
			bufferedImage = new BufferedImage(partialWidth, partialHeight, BufferedImage.TYPE_INT_ARGB);
			g2buff = bufferedImage.createGraphics();
			Point xy = viewPort.getViewPosition();
			int xPt = (int)Math.round(xy.getX()*scale*-1);
			int yPt = (int)Math.round(xy.getY()*scale*-1);
			g2buff.translate(xPt,yPt);  	
		}
		else {//full sized image
			bufferedImage = new BufferedImage ((int)Math.round(width), (int)Math.round(height), BufferedImage.TYPE_INT_ARGB);
			g2buff = bufferedImage.createGraphics();
			g2buff.translate(0,-200*scale); //to kill 200 pix buffer, trailing is also lost
		}
		g2buff.addRenderingHints(GATAUtil.fetchRenderingHints());
		//kill antialiasing to make sharp scale ruler
		g2buff.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
		g2buff.scale(scale,scale);
		g2buff.setColor(gataParams.getBackgroundColor());
		g2buff.fill(new Rectangle2D.Double(0,0,getWidth(), getHeight()+400));
		drawGlyphPanel(g2buff);
		g2buff.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		return bufferedImage;
	}


	private class MouseHandler extends MouseAdapter {
		public void mousePressed(MouseEvent event) {
			Point2D pt = event.getPoint();							
			//scan GeneRepGlyphs
			for (int i = geneGlyphs.length - 1; i >= 0; i--) {
				if ((geneGlyphs[i].getBoundingBox()).contains(pt)) {
			console.printToTextArea(geneGlyphs[i].toString());
					return;
				}
			}
			//scan GenericGlyphs
			if (anyGenericGlyphs) {
				for (int i = genericGlyphs.length - 1; i >= 0; i--) {
					if (genericGlyphs[i]==null)continue;
					for (int j = genericGlyphs[i].length - 1; j >= 0; j--) {
						if ((genericGlyphs[i][j].getBoundingBox())
							.contains(pt)) {
			console.printToTextArea(genericGlyphs[i][j].toString());
							return;
						}
					}
				}
			}
			console.printToTextArea("----------------------------------\n");
		}
	}
	public JViewport getViewPort() {
		return viewPort;
	}
	public void setViewPort(JViewport viewport) {
		viewPort = viewport;
	}

}
