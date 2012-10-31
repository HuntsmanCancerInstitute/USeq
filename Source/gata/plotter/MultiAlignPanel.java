package gata.plotter;

import gata.aligner.*;
import gata.geneGlyphs.*;
import gata.main.*;

import java.util.*;
import java.awt.*;
import java.awt.geom.*;
import javax.swing.*;
import java.awt.event.*;
import java.io.*;

/**
 * @author Nix
 * Panel holding the DNA box-line-box alignment
 * */
public class MultiAlignPanel extends JPanel {
	/*
	//global user parameters to apply to all alignments
	//location params for bounding boxes
	private double ntPerPixel;
	
	private double minScore;
	private double maxScore;
	private double maxScoreColor;
	private double minScoreColor;
	private Alignment[][] alignments;
	private Alignment[][] visableAlignments;
	private AlignParams[] alignParams;
	private GATAParams params;
	
	private Rectangle2D[] sequenceBoxes;

	private Rectangle2D lineRec;
	private Rectangle2D selectionRec = new Rectangle2D.Double();
	
	//misc stuff
	private ArrayList shapes = new ArrayList(3);
	private int orientation = 0;
	private double score = 0;
	private Console console;
	private  ToolsFrame tools;
	private JViewport viewPort; //for zooming
	private int widthBuffImage;  //for saving buffered image
	private	int heightBuffImage;  //for saving buffered image
	
	//params for either ref or comp annotation panels
	private AnnoSpecParams refParams;
	private boolean gffRefPresent;
	private AnnoSpecParams compParams;    
	private boolean gffCompPresent;
	
	private Point mousePressedPoint;
	private boolean selectRecVis = false;
	private int lenVisAligns; //length visable alignments arrayList
	
	public MultiAlignPanel(GATAParams params){
		//get plotterparam object and set variables
		this.params=params;
		setPresenceOfGFFs();
		
		//get saved objects
		//ArrayList objs = GATAUtil.fetchArrayList(params.getObjectsFile());
		String[] files = {"/Users/nix/Desktop/delme/delmeYesExtract1>pseudoobsc.gata","/Users/nix/Desktop/delme/delmeYesExtract2>erecta_eve.gata"};
		ArrayList objs1 = GATAUtil.fetchArrayList(new File(files[0]));
		ArrayList objs2 = GATAUtil.fetchArrayList(new File(files[0]));
		alignments = new Alignment[2][];
		alignments[0] = (Alignment[])objs1.get(1);
		alignments[1] = (Alignment[])objs2.get(1);
		
		//get param object and set fields
		alignParams = new AlignParams[2];
		alignParams[0] = (AlignParams)objs1.get(0);
		alignParams[1] = (AlignParams)objs2.get(0);
		
		//check that alignParams are all the same
		boolean alignParamsOK = checkAlignParams(alignParams);
		System.out.println("AlignParams :"+alignParamsOK);
		
		//set min and max alignment scores
		maxScore = alignParams[0].getMAX_SCORE();
		minScore = fetchMinScore(alignParams);
		
		//make array of visable alignments
		makeVisableAlignments();
		
		//get refs
		console = params.getConsole();
		
		//set AP in plotterParams, should all be the same except the comp seq info
		params.setAlignParams(alignParams[0]);
		
		//make shapes using width and height from plotterParams
		makeAllShapes(params.getWidth(), params.getHeight());
		
		//add listeners to monitor box alignment selector and base position
		addMouseListener(new MouseHandler());
		addMouseMotionListener(new MouseMotionHandler());
		
	}
	
	public void setToolsRef(ToolsFrame aTF){tools = aTF;} //needed since tf is created last, and needed in mouse motion listener, don't want to call it there due to speed penalties
	
	public void setPresenceOfGFFs(){
		refParams = params.getRefAnnoSpecParams();
		gffRefPresent = params.isGffRefPresent();
		compParams = params.getCompAnnoSpecParams();
		gffCompPresent = params.isGffCompPresent();
		
	}
	
	public void zoom(int zScalar){
		double newAlignFrameWidth = (double)alignFrameWidth + (double)alignFrameWidth*(double)zScalar/100; //set new width
		//for centering viewPort
		int paramWidth = params.getWidth();
		int newX;
		
		//takes a percent (25%) and adds that to the width and redraws the objects
		if (zScalar == 0 || newAlignFrameWidth <= paramWidth ) {
			alignFrameWidth = paramWidth; 
			newX =0;
		}
		else {
			Point oldPoint = viewPort.getViewPosition();
			double oldViewX = oldPoint.getX();
			double halfWidthViewPort = 0.5*(double)viewPort.getWidth();
			double m = (oldViewX + halfWidthViewPort)/(double)alignFrameWidth;
			newX = (int)Math.round(m*newAlignFrameWidth - halfWidthViewPort);
			alignFrameWidth= (int)newAlignFrameWidth;
		} 
		
		makeAllShapes(alignFrameWidth, alignFrameHeight);//this changes plotter params too
		//set size of panel for scroll pane
		setPreferredSize(new Dimension(alignFrameWidth,alignFrameHeight+400)); //400px buffer
		//hide then show to kill flicker
		//redraw GlyphPanelRef
		if(gffRefPresent) {
			refParams.getGlyphPanel().setVisible(false);
			redrawGlyphPanel(refParams);
		} 
		if(gffCompPresent) {
			compParams.getGlyphPanel().setVisible(false);
			redrawGlyphPanel(compParams);
		} 
		setVisible(false);
		revalidate(); //simply repaint() ing doesn't work use revalidate();
		viewPort.setViewPosition(new Point(newX,(int)viewPort.getViewPosition().getY()));
		setVisible(true);
		if(gffRefPresent) refParams.getGlyphPanel().setVisible(true);
		if(gffCompPresent) compParams.getGlyphPanel().setVisible(true);
		
	}
	
	
	
	
	public void redrawGlyphPanel(AnnoSpecParams asp){
		//clearlines
		asp.clearLines();
		//set modified width
		ScaleBar sb = asp.getScaleBar();	
		//if (asp.isScaleBarVis()) sb.makeScaleBar();
		if (asp.isScaleRulerVis()) {
			params.setZoomedWidth(alignFrameWidth);
			sb.makeScaleRuler();
		}
		Annotation anno= asp.getAnnotation();	 
		anno.drawGlyphs();
		if (asp.getAnyGenericGlyphs()) {
			anno.drawGenericGlyphs();
		} 
		GlyphPanel gp = asp.getGlyphPanel();	
		gp.setLine2DArrays();
		gp.setPreferredSize(new Dimension(alignFrameWidth,asp.getAnnoPanelHeight()+400)); //for buffer
		gp.revalidate();    
	}
	
	/**Makes all shapes for faster drawing.*/
/*	public void makeAllShapes(double height, double width){
		//make ref seq fill panel
		double refX = width*params.getB(); //first base of reference sequence
		double refY = height*params.getA()+200; //buffer for making larger panel
		double refW = width*params.getC();
		double refCompH = height*params.getD();
				
		double compX;
		double compY = height*params.getE()+200; //buffer for making larger panel
		double compW;


		double scalar;
		double C= width*params.getC();
		double D= height*params.getD();
		double F = height*params.getF();
		
		ntPerPixel = (double)alignParams[0].getLENGTH_REFSEQ()/refW;
		params.setNtPerPixel(ntPerPixel);
		params.setPixPerNt(1/ntPerPixel);
		if (gffRefPresent)refParams.setNtAtXZero(alignParams[0].getSTART_INDEX_REFSEQ() - (refX*ntPerPixel));

		//run thru each alignment
		int numAlignments = alignParams.length;
		ArrayList seqBoxes = new ArrayList();
		double alignHeight;
		for (int i=0; i<numAlignments; i++){
			//create sequenceBoxes and alignment shapes
			seqBoxes.add(new Rectangle2D.Double(refX,refY,refW,refCompH));
			if (alignParams[i].isCompVisible()){
				//calculate compX compW
				compW = alignParams[i].getLENGTH_COMPSEQ()/ntPerPixel;
				//center comp
				double diff = width-compW;
				compX = diff/2;
				seqBoxes.add(new Rectangle2D.Double(compX,compY,compW,refCompH));
				alignHeight = (compY+D)-refY;
			}
			else alignHeight = D;
			//create shapes for each alignment
			createShapes(alignments[i], refX, refY, compX, compY, alignParams[0].getLENGTH_REFSEQ(), C, D);
			//advance Y coords for next alignment
			refY = refY + alignHeight + F;
			compY = compY +alignHeight + F;
		}
		
check again!
		
		//set Y's for calculating minimum panel size
		params.setSmallestAlignPanelY((int)(height*params.getA()+200-1));
		params.setBiggestAlignPanelY((int)(compY-F+1));        
		
		//create rectangle for mouse line selection
		//lineRec = new Rectangle2D.Double(0,refY+refH,C*2,compY-(refY+refH));
	}
	
	public static void createShapes(Alignment[] alignments,
			double refX, double refY, double compX, double compY, double ntLenBig,
			double C, double D){
		int len = alignments.length;
		//make shapes for each alignment object based on params
		for (int i =0; i<len; i++){
			alignments[i].makeShapes(refX, refY, compX, compY, ntLenBig, C, D);
		}
	}
	
	public void paintComponent(Graphics g){
		super.paintComponent(g);
		Graphics2D g2 = (Graphics2D)g;
		//Run through array of visible alignments based on score
		for (int i =0; i<lenVisAligns; i++){
			visableAlignments[i].drawAlignment(g2);
		}
		g2.setColor(Color.BLUE);
		g2.draw(ref);
		if (drawCompBox) g2.draw(comp);
		
		//draw selection rectangle
		if (selectRecVis){
			g2.setColor(Color.CYAN);
			g2.draw(selectionRec);
		}
	}
	
	/**Makes arrays of visable Alignment based on min and max scores, these are reset by sliders.*/
/*	public void makeVisableAlignments(){
		int len = alignments.length;
		visableAlignments = new Alignment[len][];
		ArrayList al = new ArrayList(alignments[0].length);
		for (int j=0; j<len; j++){
			al.clear();
			int numAligns = alignments[j].length;
			for (int i =0; i<numAligns; i++){
				score = alignments[j][i].getScore();
				if (score < minScore || score > maxScore) continue;
				al.add(alignments[j][i]);
			}
			int numVisAligns = al.size();
			visableAlignments[j] = new Alignment[numVisAligns];
			al.toArray(visableAlignments[j]);    	
		}
	}
	
	public void setMinScore(double min){
		minScore = min;
		makeVisableAlignments();
		repaint();
	}
	public void setMaxScore(double max){
		maxScore = max;
		makeVisableAlignments();
		repaint();
	}
	*/
	/**Checks to see if the same params were used used in making the alignments, note score cut off is not checked*/
	public static boolean checkAlignParams (AlignParams[] alignParams){
		int numAPs = alignParams.length;
		if (numAPs==1) return true;
		int lengthRefSeq = alignParams[0].getLENGTH_REFSEQ();
		String refName = alignParams[0].getNAME_REFSEQ();
		int mismatch = alignParams[0].getMISMATCH();
		int gapCreate = alignParams[0].getGAP_CREATE();
		int gapExit = alignParams[0].getGAP_EXT();
		int windowSize = alignParams[0].getWIN_SIZE();
		int match = alignParams[0].getMATCH();
		int startIndexRef = alignParams[0].getSTART_INDEX_REFSEQ();
		
		for (int i=1; i<numAPs; i++){
			if (lengthRefSeq != alignParams[i].getLENGTH_REFSEQ()) return false;
			if (refName.equals(alignParams[i].getNAME_REFSEQ())== false) return false;
			if (mismatch != alignParams[0].getMISMATCH()) return false;
			if (gapCreate != alignParams[0].getGAP_CREATE()) return false;
			if (gapExit != alignParams[0].getGAP_EXT()) return false;
			if (windowSize != alignParams[0].getWIN_SIZE()) return false;
			if (match != alignParams[0].getMATCH()) return false;
			if (startIndexRef != alignParams[0].getSTART_INDEX_REFSEQ()) return false;
		}
		return true;
	}
	
	/**Finds lowest minScore in an array of AlignParams.*/
	public static int fetchMinScore (AlignParams[] alignParams){
		int numAPs = alignParams.length;
		int minScore = alignParams[0].getMIN_SCORE();
		for (int i=1; i<numAPs; i++){
			if (minScore> alignParams[i].getMIN_SCORE()) minScore = alignParams[i].getMIN_SCORE();
		}
		return minScore;
	}
	/*
	private class MouseHandler extends MouseAdapter {
		private Point2D pt;
		private Rectangle2D rec1;
		private Rectangle2D rec2;
		private Line2D ln;
		private ConservedSeqs conservedSeqs;
		private double dist;
		private double pixels;
		private long bases;
		
		public void mouseClicked(MouseEvent event){
			int clicks = event.getClickCount();
			pt = event.getPoint();
			for (int i =0; i<lenVisAligns; i++){
				//get shapes array for a particular alignment
				shapes = visableAlignments[i].getShapes();
				rec1 = (Rectangle2D)shapes.get(0);
				rec2 = (Rectangle2D)shapes.get(1);
				//check to see if the selected a box
				if (rec1.contains(pt) || rec2.contains(pt)) {
					if (clicks == 2) console.printToTextArea(visableAlignments[i].getParentLocAlignString());
					else console.printToTextArea(visableAlignments[i].getAlignmentString());
				}
				//check to see if they selected a line
				//first check to see if it's within the line zone
				else if (lineRec.contains(pt)){
					ln = (Line2D)shapes.get(2);
					dist = ln.ptLineDist(pt);
					//is the mouse within 5 pixels of a line
					if (dist <5){
						if (clicks == 2) console.printToTextArea(visableAlignments[i].getParentLocAlignString());
						else console.printToTextArea(visableAlignments[i].getAlignmentString());
					}
				}
			}
			console.printToTextArea("----------------------------------\n");
		}
		public void mousePressed(MouseEvent event){
			mousePressedPoint = event.getPoint();
			selectRecVis = true;
			selectionRec.setFrame(0,0,0,0);
		}
		
		public void mouseReleased(MouseEvent event){
			//the following is used to return a selected dna sequence given a selection rectangle
			//	the kicker is that if alignment boxes are selected then the other orthologous seq is returned
			//check if they made a selection by dragging
			if (selectionRec.getWidth()>4){
				boolean refSelected = false;
				boolean compSelected = false;
				int leftBaseRef=0;
				int rightBaseRef=0;
				int leftBaseComp=0;
				int rightBaseComp=0;
				
				//find out if they selected top or bottom sequence box
				if (selectionRec.intersects(ref)) {
					refSelected = true;
					//get flanks ref
					pixels = selectionRec.getX() - refX;
					leftBaseRef = (int)Math.round((ntPerPixel * pixels));
					pixels = selectionRec.getX() + selectionRec.getWidth() - refX;
					rightBaseRef = (int)Math.round((ntPerPixel * pixels));
					if (leftBaseRef<0 || leftBaseRef>ntLenRefSeq) leftBaseRef=0;
					if (rightBaseRef<0 || rightBaseRef>ntLenRefSeq) rightBaseRef = (int)ntLenRefSeq;
				} 
				if (selectionRec.intersects(comp)){ 
					compSelected = true;
					//get flanks comp
					pixels = selectionRec.getX() - compX;
					leftBaseComp = (int)Math.round((ntPerPixel * pixels));
					pixels = selectionRec.getX() + selectionRec.getWidth() - compX;
					rightBaseComp = (int)Math.round((ntPerPixel * pixels));
					if (leftBaseComp<0 || leftBaseComp>ntLenCompSeq) leftBaseComp=0;
					if (rightBaseComp<0 || rightBaseComp>ntLenCompSeq) rightBaseComp = (int)ntLenCompSeq;	
				}		
				if (refSelected && compSelected){ // both selected fire ConservedSeqs with coords
					if (conservedSeqs == null) conservedSeqs = new ConservedSeqs(params);
					conservedSeqs.fetchNewRegions(leftBaseRef, rightBaseRef, leftBaseComp, rightBaseComp);		
				}
				else if (refSelected || compSelected){ //if one or the other
					int startRef=(int)ntLenRefSeq;
					int stopRef=0;
					int startComp=(int)ntLenCompSeq;
					int stopComp=0;
					boolean recFound = false;
					int recNum = 1; //for refSelected
					if (compSelected) recNum = 0; //for compSelected
					//find selected rectangles and expand start stops
					for (int i =0; i<lenVisAligns; i++){
						//get shapes array for a particular alignment
						shapes = visableAlignments[i].getShapes();
						rec1 = (Rectangle2D)shapes.get(recNum); 
						//check to see if they selected a box
						if (selectionRec.intersects(rec1)) {
							recFound = true;
							int[] coor = visableAlignments[i].getRelativeCoordinates();
							//expand sizes
							if (coor[0]<startRef) startRef=coor[0];
							if (coor[1]>stopRef) stopRef=coor[1];
							if (coor[2]<startComp) startComp = coor[2];
							if (coor[3]>stopComp) stopComp= coor[3];
						}
					}
					if(recFound && refSelected){
						//calculate distances from selection rectangle to first block
						int F = startRef-leftBaseRef; //could be negative
						int R = rightBaseRef-stopRef; //ditto
						//modify comp bps
						startComp -=F;
						stopComp +=R;
						if (startComp<0) startComp=0;
						if (stopComp>ntLenCompSeq) stopComp = (int)ntLenCompSeq;
						leftBaseComp = startComp;
						rightBaseComp = stopComp;
					}
					else if(recFound && compSelected){
						//calculate distances from selection rectangle to first block
						int F = startComp-leftBaseComp; //could be negative
						int R = rightBaseComp-stopComp; //ditto
						//modify comp bps
						startRef -=F;
						stopRef +=R;
						if (startRef<0) startRef=0;
						if (stopRef>ntLenRefSeq) stopRef = (int)ntLenRefSeq;
						leftBaseRef = startRef;
						rightBaseRef = stopRef;
					}				
					else if (refSelected && recFound==false)  {
						leftBaseComp = 0;
						rightBaseComp = 0;
					}
					else if (compSelected && recFound==false)  {
						leftBaseRef = 0;
						rightBaseRef = 0;
					}
					//fire conservedSeqs if an intersecting rectangle was found
					if (conservedSeqs == null) conservedSeqs = new ConservedSeqs(params);
					conservedSeqs.fetchNewRegions(leftBaseRef, rightBaseRef, leftBaseComp, rightBaseComp);
				}
			}
			selectRecVis = false;
			repaint();
		}
	}
	private class MouseMotionHandler implements MouseMotionListener{
		private Point2D pt;
		private double pixels;
		private long bases;    
		public void mouseMoved (MouseEvent event){
			//if in box get base number
			pt = event.getPoint();
			/*if (ref.contains(pt)){
				pixels = pt.getX() - refX;
				bases = Math.round((ntPerPixel * pixels) + refRealStart); 
				tools.setNtRefSeq((int)bases);    
			}
			else if (comp.contains(pt)){
				pixels = pt.getX() - compX;
				bases = Math.round((ntPerPixel * pixels) + compRealStart);
				tools.setNtCompSeq((int)bases);
			}
		}
		public void mouseDragged(MouseEvent mouseEvent) {
			selectionRec.setFrameFromDiagonal(mousePressedPoint, mouseEvent.getPoint());
			repaint();
		}   
	}
	public void setConsole(Console console) {
		this.console = console;
	}
	public void setViewPort(JViewport viewport) {
		viewPort = viewport;
	}
	public int getHeightBuffImage() {
		return heightBuffImage;
	}
	public int getWidthBuffImage() {
		return widthBuffImage;
	}
	public JViewport getViewPort() {
		return viewPort;
	}
	*/
}