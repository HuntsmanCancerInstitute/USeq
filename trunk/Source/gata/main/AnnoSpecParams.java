package gata.main;

import gata.geneGlyphs.*;

import java.util.*;
import java.awt.*;
import java.awt.geom.*;
import java.io.*;

/**
 * @author Nix
Contains fields specific for each annotation panel, either for the Reference Sequence or Comparative Sequence
 */
public class AnnoSpecParams {
	
//	FIELDS, many are initialize using the PlotterPreferences object
	  //location and scale
	  private double referenceY=200; //pixel y-coor where annotation panel should start DNA glyph
	  private int annoPanelHeight=400; //height of annotation panel
	private double scaleRulerY=210; //where to initially draw scale Ruler, this is reset
	private double scaleBarX=10;
	private double scaleBarY=200;
	
	//params to figure size of drawn panels
	private int smallestGlyphPanelY;
	private int biggestGlyphPanelY;

	//arrayLists for annotation panel
	private ArrayList labels = new ArrayList(); //arraylist containing a text, and an int[2] midptX,ycoord
	private ArrayList[] trackLabels; //ditto
	private ArrayList thinDNALines = new ArrayList();
	private ArrayList normalDNALines = new ArrayList();
	private ArrayList thinRNALines = new ArrayList();
	private ArrayList fatRNALines = new ArrayList();
	private ArrayList thinProteinLines = new ArrayList();
	private ArrayList dottedProteinLines = new ArrayList();

	private GenericGlyph[][] trackedGenericGlyphs;
	private GeneRepGlyph[] geneGlyphs;
	private Annotation annotation;
	private Graphics2D g2;
	private ScaleBar scaleBar;
	private Color scaleLineColor = Color.CYAN;
	private boolean scaleBarVis = true;
	private boolean ScaleRulerVis = true;

	private GlyphPanel glyphPanel;
	private File gffFile; //full path file text for gff file
	
	private double ntAtXZero;

	private boolean anyGenericGlyphs;
	private boolean[] trackVis;
	private boolean[] trackLabelVis;
	private BasicStroke trackStroke[];
	private int trackNumber; //number of generic tracks
	private float[] trackThickness;

	//	adder methods
	public void addLabel(ArrayList list) {
		labels.addAll(list);
	}
	public void addNormalDNALines(ArrayList al) {
		normalDNALines.addAll(al);
	}
	public void addProteinLines(
		ArrayList thinLinesProt,
		ArrayList dottedLinesProt) {
		thinProteinLines.addAll(thinLinesProt);
		dottedProteinLines.addAll(dottedLinesProt);
	}
	public void addRNALines(Line2D.Double baseLine, ArrayList fatLines) {
		thinRNALines.add(baseLine);
		fatRNALines.addAll(fatLines);
	}
	public void addThinDNALines(ArrayList al) {
		thinDNALines.addAll(al);
	}
	//misc methods
	/**Wipes all the Line2D.Double ArrayLists, call before reDrawing
	 * the Glyphs when params have been changed.*/
	public void clearLines() {
		thinDNALines.clear();
		normalDNALines.clear();
		thinRNALines.clear();
		fatRNALines.clear();
		thinProteinLines.clear();
		dottedProteinLines.clear();
		labels.clear();
	}
	public int getAnnoPanelHeight() {
		return annoPanelHeight;
	}
	public Annotation getAnnotation() {
		return annotation;
	}
	public boolean getAnyGenericGlyphs() {
		return anyGenericGlyphs;
	}	
	public ArrayList getDottedProteinLines() {
		return dottedProteinLines;
	}	
	public ArrayList getFatRNALines() {
		return fatRNALines;
	}
	public GeneRepGlyph[] getGeneGlyphs() {
		return geneGlyphs;
	}
	public GlyphPanel getGlyphPanel() {
		return glyphPanel;
	}
	public Graphics2D getGraphics2DN() {
		return g2;
	}	
	public ArrayList getLabels() {
		return labels;
	}
	public ArrayList getNormalDNALines() {
		return normalDNALines;
	}
	
	public double getReferenceY() {
		return referenceY;
	}
	public ArrayList getThinDNALines() {
		return thinDNALines;
	}
	public ArrayList getThinProteinLines() {
		return thinProteinLines;
	}
	public ArrayList getThinRNALines() {
		return thinRNALines;
	}
	public boolean[] getTrackLabelVis() {
		return trackLabelVis;
	}
	public boolean[] getTrackVis() {
		return trackVis;
	}
	public GenericGlyph[][] getTrackedGenericGlyphs() {
		return trackedGenericGlyphs;
	}
	public ArrayList[] getTrackLabels() {
		return trackLabels;
	}
	public int getTrackNumber() {
		return trackNumber;
	}		
	public void setAnnotation(Annotation annotation) {
		this.annotation = annotation;
	}
	public void setAnyGenericGlyphs(boolean x) {
		anyGenericGlyphs = x;
	}	
	public void setGeneGlyphs(GeneRepGlyph[] glyphs) {
		geneGlyphs = glyphs;
	}
	public void setGlyphPanel(GlyphPanel panel) {
		glyphPanel = panel;
	}
	//setter methods
	public void setGraphics2DN(Graphics2D g2) {
		this.g2 = g2;
	}
	public void setTrackLabelVis(boolean[] bs) {
		trackLabelVis = bs;
	}
	public void setTrackVis(boolean[] bs) {
		trackVis = bs;
	}
	public void setTrackedGenericGlyphs(GenericGlyph[][] glyphs) {
		trackedGenericGlyphs = glyphs;
	}
	public void setTrackLabels(ArrayList[] lists) {
		trackLabels = lists;
	}
	public void setTrackNumber(int i) {
		trackNumber = i;
	}
	public void setAnnoPanelHeight(int i) {
		annoPanelHeight = i;
	}	
	public int getSmallestGlyphPanelY() {
		return smallestGlyphPanelY;
	}
	public void setSmallestGlyphPanelY(int i) {
		smallestGlyphPanelY = i;
	}
	public int getBiggestGlyphPanelY() {
		return biggestGlyphPanelY;
	}
	public void setBiggestGlyphPanelY(int i) {
		biggestGlyphPanelY = i;
	}
	public void setReferenceY(double d) {
		referenceY = d;
	}
	public double getScaleRulerY() {
		return scaleRulerY;
	}
	public void setScaleRulerY(double d) {
		scaleRulerY = d;
	}
	public ScaleBar getScaleBar() {
		return scaleBar;
	}
	public void setScaleBar(ScaleBar bar) {
		scaleBar = bar;
	}
	public double getNtAtXZero() {
		return ntAtXZero;
	}
	public void setNtAtXZero(double d) {
		ntAtXZero = d;
	}
	public File getGffFile() {
		return gffFile;
	}
	public void setGffFile(File file) {
		gffFile = file;
	}
	public BasicStroke[] getTrackStrokes() {
		return trackStroke;
	}
	public void setTrackThickness(float[] fs) {
		trackThickness = fs;
		int len = fs.length;
		trackStroke = new BasicStroke[len];
		for (int i = 0; i < len; i++) {
			trackStroke[i] =
				new BasicStroke(
					trackThickness[i],
					BasicStroke.CAP_BUTT,
					BasicStroke.JOIN_MITER,
					15.0F);
		}
	}
	public float[] getTrackThickness() {
		return trackThickness;
	}
	public Color getScaleLineColor() {
		return scaleLineColor;
	}
	public void setScaleLineColor(Color color) {
		scaleLineColor = color;
	}
	public boolean isScaleBarVis() {
		return scaleBarVis;
	}
	public boolean isScaleRulerVis() {
		return ScaleRulerVis;
	}
	public void setScaleBarVis(boolean b) {
		scaleBarVis = b;
	}
	public void setScaleRulerVis(boolean b) {
		ScaleRulerVis = b;
	}
	public double getScaleBarX() {
		return scaleBarX;
	}
	public double getScaleBarY() {
		return scaleBarY;
	}
	public void setScaleBarX(double d) {
		scaleBarX = d;
	}
	public void setScaleBarY(double d) {
		scaleBarY = d;
	}

}
