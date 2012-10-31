package gata.main;
import java.awt.*;
import java.io.*;
/**
 * @author nix
 Used to set GATAPlotter preferences upon booting, saved upon close.
 */
public class PlotterPreferences implements Serializable {
	//FIELDS
	//location and scale
	//drawing paramemter defaults, can be over ridden by args for panels
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
	
	//spacers for annotation panel
	private double pixBtwProtRNA;//pixels between protein and RNA glyphs (translation and transcription/exons)
	private double pixBtwTransGrps; //pixels between TransGrps
	private double pixBtwTransGrpsGene;//pixels between last/first TransGroup and gene
	private double pixBtwTransGrpLabel;//pixels between TransGroup and Label
	private double pixBtwTracks;
	private double pixBtwLabelTrack;
	private double pixBtwDNATracks;

	//thickness for annotation panel
	private double pixHeightArrows ;//pixel height of arrow and stop line cap
	private float pixThinLineThickness;
	private float pixFatLineThickness ; //thickness of fat lines 
	private float pixNormalLineThickness ;
	
	
	//colors for annotation panel
	private Color backgroundColor ;
	private Color labelColor ;
	private Color scaleLineColor;
	private Color DNAColor; 
	private Color RNAColor; 
	private Color proteinColor; 
	
	//METHODS
	public void fetchAndAssignGATAParams(GATAParams p){
		width= p.getWidth();
		zoomedWidth = p.getZoomedWidth();
		height= p.getHeight();
		A= p.getA();
		B= p.getB();
		C= p.getC();
		D= p.getD();
		E= p.getE();
		F= p.getF();
		pixBtwProtRNA= p.getPixBtwProtRNA();
		pixBtwTransGrps= p.getPixBtwTransGrps();
		pixBtwTransGrpsGene= p.getPixBtwTransGrpsGene();
		pixBtwTransGrpLabel= p.getPixBtwTransGrpLabel();
		pixBtwTracks= p.getPixBtwTracks();
		pixBtwLabelTrack= p.getPixBtwLabelTrack();
		pixBtwDNATracks= p.getPixBtwDNATracks();
		pixHeightArrows= p.getPixHeightArrows();
		pixThinLineThickness= p.getPixThinLineThickness();
		pixFatLineThickness= p.getPixFatLineThickness();
		pixNormalLineThickness= (float)p.getPixNormalLineThickness();
		backgroundColor= p.getBackgroundColor();
		labelColor= p.getLabelColor();
		DNAColor= p.getDNAColor();
		RNAColor= p.getRNAColor();
		proteinColor= p.getProteinColor();
	}
	
	public void setAnnoScaleDefaults(AnnoSpecParams asp){
		asp.setReferenceY(200);
		asp.setAnnoPanelHeight(400);
		asp.setScaleRulerY(210);	
		asp.setScaleLineColor(Color.CYAN);
		asp.setScaleBarVis(true);
		asp.setScaleRulerVis(true);
	}	
	
	/**Contains all the default preferences for GATAParams, alter these to change the default view*/
	public void setPlotterDefaults(GATAParams p){
		if (p.isGffRefPresent()) setAnnoScaleDefaults(p.getRefAnnoSpecParams());
		if (p.isGffCompPresent()) setAnnoScaleDefaults(p.getCompAnnoSpecParams());
		p.setWidth(1000);
		p.setZoomedWidth(1000);
		p.setHeight(200);
		p.setA(0.1);
		p.setB(0.025);
		p.setC(0.95);
		p.setD(0.1);
		p.setE(0.8);
		p.setF(0.01);
		p.setPixBtwProtRNA(2);
		p.setPixBtwTransGrps(4);
		p.setPixBtwTransGrpsGene(5);
		p.setPixBtwTransGrpLabel(0);
		p.setPixBtwTracks(4);
		p.setPixBtwLabelTrack(0);
		p.setPixBtwDNATracks(17);
		p.setPixHeightArrows(4);
		p.setPixThinLineThickness(1);
		p.setPixFatLineThickness(4);
		p.setPixNormalLineThickness(4);
		p.setBackgroundColor(Color.BLACK);
		p.setLabelColor(Color.WHITE);
		p.setDNAColor(Color.BLUE);
		p.setRNAColor(Color.RED);
		p.setProteinColor(Color.YELLOW);
		p.setDNAVis(true);
		p.setRNAVis(true);
		p.setProteinVis(true);
		p.setLabelsVis(true);
		p.setGenericLabelsVis(true);
	}
	/**Contains all the default preferences for GATAParams*/
	public void setCurrentPlotterPreferences(GATAParams p){
		if (p.isGffRefPresent()) setAnnoScaleDefaults(p.getRefAnnoSpecParams());
		if (p.isGffCompPresent()) setAnnoScaleDefaults(p.getCompAnnoSpecParams());
		p.setWidth(width);
		p.setZoomedWidth(zoomedWidth);
		p.setHeight(height);
		p.setA(A);
		p.setB(B);
		p.setC(C);
		p.setD(D);
		p.setE(E);
		p.setF(F);
		p.setPixBtwProtRNA(pixBtwProtRNA);
		p.setPixBtwTransGrps(pixBtwTransGrps);
		p.setPixBtwTransGrpsGene(pixBtwTransGrpsGene);
		p.setPixBtwTransGrpLabel(pixBtwTransGrpLabel);
		p.setPixBtwTracks(pixBtwTracks);
		p.setPixBtwLabelTrack(pixBtwLabelTrack);
		p.setPixBtwDNATracks(pixBtwDNATracks);
		p.setPixHeightArrows(pixHeightArrows);
		p.setPixThinLineThickness(pixThinLineThickness);
		p.setPixFatLineThickness(pixFatLineThickness);
		p.setPixNormalLineThickness(pixNormalLineThickness);
		p.setBackgroundColor(backgroundColor);
		p.setLabelColor(labelColor);
		p.setDNAColor(DNAColor);
		p.setRNAColor(RNAColor);
		p.setProteinColor(proteinColor);
	}	
}
