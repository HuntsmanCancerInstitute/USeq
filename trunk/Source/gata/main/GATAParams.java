package gata.main;

import gata.aligner.*;
import gata.menu.*;
import gata.plotter.*;
import gata.plotter.Console;

import java.awt.*;
import java.io.*;

/**
 * @author Nix
Master params file for the GATAPlotter program.
 */
public class GATAParams {
	
	//drawing paramemters 
	private int width; //width of annotation and alignment panel
	private int zoomedWidth;  //modified panel width needed for zooming
	private int height; //height of alignment panel
	
	//alignment panel spacing of dna shapes as a % of panel width or height
	private double A; //%gap btw top of panel and top DNA bar
	private double B; //%gap btw left of panel and furthest left DNA bar
	private double C; //%width largest DNA bar glyphs
	private double D; //%thickness of DNA bars
	private double E; //% gap btw top of panel and bottom DNA bar
	private double F; //% gap btw multiple alignments
	
	//params to figure size of drawn panels
	private int smallestAlignPanelY;
	private int biggestAlignPanelY;

	//spacers for annotation panel
	private double pixBtwProtRNA;//pixels between protein and RNA glyphs (translation and transcription/exons)
	private double pixBtwTransGrps; //pixels between TransGrps
	private double pixBtwTransGrpsGene;//pixels between last/first TransGroup and gene
	private double pixBtwTransGrpLabel;//pixels between TransGroup and Label
	private double pixBtwTracks;
	private double pixBtwLabelTrack;
	private double pixBtwDNATracks;

	//thickness for annotation panel
	private double pixHeightArrows;//pixel height of arrow and stop line cap
	private float pixThinLineThickness;
	private float pixFatLineThickness; //thickness of fat lines 
	private float pixNormalLineThickness;
	private double heightLabels=10;//guestimate on the height of rendered text
	private float defaultTrackThickness = 5f;
	
	//colors for annotation panel
	private Color backgroundColor;
	private Color labelColor;
	private Color DNAColor; 
	private Color RNAColor; 
	private Color proteinColor; 
		
	//file names
	private File objectsFile; //contains the GATAligner arraylist objects
	
	//object references
	private Console console;
	private AlignPanel aAlignPanel;
	private AlignParams AP;
	private ToolsFrame TF;
	private GATAPlotter gata;
	private GATAFrame gataFrame;
	private PlotterPreferences plotterPreferences;
	private AnnoSpecParams refAnnoSpecParams;
	private AnnoSpecParams compAnnoSpecParams;

	//booleans to determine what is drawn in annotation panel
	private boolean DNAVis = true;
	private boolean RNAVis = true;
	private boolean proteinVis = true;
	private boolean labelsVis = true;
	private boolean genericLabelsVis = true;
	private boolean gffRefPresent = false;
	private boolean gffCompPresent = false;
	
	//misc
	private Font font = new Font("Dialog", 1, 10); //for labels
	private String trackScoreConverter = "linear";  //used to recalculate score for color shading (linear, log,lnx)
	private double pixPerNt;	//set by width of panel and size of DNA sequences
	private double ntPerPixel;
	
	//Strokes for drawing lines
	private float[] dash = { 2.0f };

  // boolean methods	  
  public boolean areGenericLabelsVis() {
	  return genericLabelsVis;
  }
  public boolean areLabelsVis() {
	  return labelsVis;
  }
  public boolean isDNAVis() {
	  return DNAVis;
  }
  public boolean isGffRefPresent() {
	  return gffRefPresent;
  }
  public boolean isGffCompPresent() {
		return gffCompPresent;
  }
  public boolean isProteinVis() {
	  return proteinVis;
  }
  public boolean isRNAVis() {
	  return RNAVis;
  }
  

  //getter methods	  
  public double getA() {
	  return A;
  }
  public AlignPanel getAlignPanel() {
	  return aAlignPanel;
  }
  public AlignParams getAlignParams() {
	  return AP;
  }

  public double getB() {
	  return B;
  }
  public Color getBackgroundColor() {
	  return backgroundColor;
  }
  public double getC() {
	  return C;
  }
  public Console getConsole() {
	  return console;
  }
  public double getD() {
	  return D;
  }
  public BasicStroke getDashStroke() {
	  return new BasicStroke(
	pixThinLineThickness,
	BasicStroke.CAP_BUTT,
	BasicStroke.JOIN_MITER,
	10.0f,
	dash,
	0.0f);
  }
  public float getDefaultTrackThickness() {
	  return defaultTrackThickness;
  }
  public Color getDNAColor() {
	  return DNAColor;
  }

  public double getE() {
	  return E;
  }
  public double getF() {
	  return F;
}
  public double getFatLineThicknes() {
	  return (double) pixFatLineThickness;
  }

  public BasicStroke getFatStroke() {
	  return new BasicStroke(
	pixFatLineThickness,
	BasicStroke.CAP_BUTT,
	BasicStroke.JOIN_BEVEL);
  }
  public Font getFont() {
	  return font;
  }
  public int getHeight() {
	  return height;
  }
  public double getHeightLabels() {
	  return heightLabels;
  }
  public Color getLabelColor() {
	  return labelColor;
  }

  public BasicStroke getNormalStroke() {
	  return new BasicStroke(
	pixNormalLineThickness,
	BasicStroke.CAP_BUTT,
	BasicStroke.JOIN_BEVEL);
  }

  public double getNtPerPixel() {
	  return ntPerPixel;
  }
  public double getPixBtwDNATracks() {
	  return pixBtwDNATracks;
  }
  public double getPixBtwLabelTrack() {
	  return pixBtwLabelTrack;
  }
  public double getPixBtwProtRNA() {
	  return pixBtwProtRNA;
  }
  public double getPixBtwTracks() {
	  return pixBtwTracks;
  }
  public double getPixBtwTransGrpLabel() {
	  return pixBtwTransGrpLabel;
  }
  public double getPixBtwTransGrps() {
	  return pixBtwTransGrps;
  }
  public double getPixBtwTransGrpsGene() {
	  return pixBtwTransGrpsGene;
  }
  public float getPixFatLineThickness() {
	  return pixFatLineThickness;
  }
  public double getPixHeightArrows() {
	  return pixHeightArrows;
  }
  public double getPixNormalLineThickness() {
	  return (double) pixNormalLineThickness;
  }
  public double getPixPerNt() {
	  return pixPerNt;
  }
  public Color getProteinColor() {
	  return proteinColor;
  }

  public Color getRNAColor() {
	  return RNAColor;
  }

  public BasicStroke getThinStroke() {
	  return new BasicStroke(
	pixThinLineThickness,
	BasicStroke.CAP_BUTT,
	BasicStroke.JOIN_BEVEL);
  }
  public ToolsFrame getToolsFrameRef() {
	  return TF;
  }

  public int getWidth() {
	  return width;
  }
  public int getZoomedWidth() {
	  return zoomedWidth;
  }

  //setter methods
  public void setAlignPanel(AlignPanel a) {
	  aAlignPanel = a;
  }
  public void setAlignParams(AlignParams ap) {
	  AP = ap;
  }

  public void setConsole(Console console) {
	  this.console = console;
  }

  public void setHeightLabels(double d) {
	  heightLabels = d;
  }
  public void setNtPerPixel(double d) {
	  ntPerPixel = d;
  }
  public void setPixBtwProtRNA(double d) {
	  pixBtwProtRNA = d;
  }
  public void setPixBtwTransGrps(double d) {
	  pixBtwTransGrps = d;
  }
  public void setPixBtwTransGrpsGene(double d) {
	  pixBtwTransGrpsGene = d;
  }
  public void setPixFatLineThickness(float f) {
	  pixFatLineThickness = f;
  }
  public void setPixHeightArrows(double d) {
	  pixHeightArrows = d;
  }
  public void setPixNormalLineThickness(float f) {
	  pixNormalLineThickness = f;
  }
  public void setPixPerNt(double d) {
	  pixPerNt = d;
  }
  public void setToolsFrameRef(ToolsFrame aTF) {
	  TF = aTF;
  }

  public void setZoomedWidth(int i) {
	  zoomedWidth = i;
  }
	public File getObjectsFile() {
		return objectsFile;
	}
	public void setObjectsFile(File file) {
		objectsFile = file;
	}
	public GATAPlotter getGataPlotter() {
		return gata;
	}
	public void setGataPlotter(GATAPlotter gata) {
		this.gata = gata;
	}
	public void setGffRefPresent(boolean b) {
		gffRefPresent = b;
	}
	public void setGffCompPresent(boolean b) {
		gffCompPresent = b;
	}
	
	public AlignPanel getAAlignPanel() {
		return aAlignPanel;
	}
	public boolean isLabelsVis() {
		return labelsVis;
	}
	public float getPixThinLineThickness() {
		return pixThinLineThickness;
	}
	public void setA(double d) {
		A = d;
	}
	public void setAAlignPanel(AlignPanel panel) {
		aAlignPanel = panel;
	}

	public void setB(double d) {
		B = d;
	}
	public void setBackgroundColor(Color color) {
		backgroundColor = color;
	}
	public void setC(double d) {
		C = d;
	}
	public void setD(double d) {
		D = d;
	}
	public void setDNAColor(Color color) {
		DNAColor = color;
	}
	public void setDNAVis(boolean b) {
		DNAVis = b;
	}
	public void setE(double d) {
		E = d;
	}
	public void setF(double d) {
		F = d;
	}
	public void setHeight(int i) {
		height = i;
	}
	public void setLabelColor(Color color) {
		labelColor = color;
	}
	public void setLabelsVis(boolean b) {
		labelsVis = b;
	}
	public void setPixBtwDNATracks(double d) {
		pixBtwDNATracks = d;
	}
	public void setPixBtwLabelTrack(double d) {
		pixBtwLabelTrack = d;
	}
	public void setPixBtwTracks(double d) {
		pixBtwTracks = d;
	}
	public void setPixBtwTransGrpLabel(double d) {
		pixBtwTransGrpLabel = d;
	}
	public void setPixThinLineThickness(float f) {
		pixThinLineThickness = f;
	}
	public void setProteinColor(Color color) {
		proteinColor = color;
	}
	public void setProteinVis(boolean b) {
		proteinVis = b;
	}
	public void setRNAColor(Color color) {
		RNAColor = color;
	}
	public void setRNAVis(boolean b) {
		RNAVis = b;
	}
	public void setWidth(int i) {
		width = i;
	}
	public GATAFrame getGataFrame() {
		return gataFrame;
	}
	public void setGataFrame(GATAFrame frame) {
		gataFrame = frame;
	}

	public int getBiggestAlignPanelY() {
		return biggestAlignPanelY;
	}
	public int getSmallestAlignPanelY() {
		return smallestAlignPanelY;
	}
	public void setBiggestAlignPanelY(int i) {
		biggestAlignPanelY = i;
	}
	public void setSmallestAlignPanelY(int i) {
		smallestAlignPanelY = i;
	}
	public void setDefaultTrackThickness(float f) {
		defaultTrackThickness = f;
	}
	public boolean isGenericLabelsVis() {
		return genericLabelsVis;
	}
	public void setGenericLabelsVis(boolean b) {
		genericLabelsVis = b;
	}
	public PlotterPreferences getPlotterPreferences() {
		return plotterPreferences;
	}
	public void setPlotterPreferences(PlotterPreferences preferences) {
		plotterPreferences = preferences;
	}
	public AnnoSpecParams getCompAnnoSpecParams() {
		return compAnnoSpecParams;
	}
	public AnnoSpecParams getRefAnnoSpecParams() {
		return refAnnoSpecParams;
	}
	public void setCompAnnoSpecParams(AnnoSpecParams params) {
		compAnnoSpecParams = params;
	}
	public void setRefAnnoSpecParams(AnnoSpecParams params) {
		refAnnoSpecParams = params;
	}

}
