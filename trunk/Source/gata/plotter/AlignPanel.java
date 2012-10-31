package gata.plotter;

import gata.aligner.*;
import gata.geneGlyphs.*;
import gata.main.*;

import java.util.*;
import java.awt.*;
import java.awt.geom.*;
import javax.swing.*;
import java.awt.event.*;
import java.awt.image.*;

/**
 * @author Nix
 * Panel holding the DNA box-line-box alignment
 * */
public class AlignPanel extends JPanel {
    
    //user parameters
    private double A;
    private double B;
    private double C;
    private double D;
    private double E;
    private double ntLenBig;
    private double ntLenSmall;
    private double ntLenRefSeq;
    private double ntLenCompSeq;
    private double ntPerPixel;
    private double refRealStart;   //nucleotide ref
    private double compRealStart;  //nucleotide ref
   
    private double windowSize;
    private double maxScore;
    private double minScore;
    private double maxScoreColor;
    private double minScoreColor;
    private Alignment[] alignments;
	private Alignment[] visableAlignments;
    private AlignParams AP;
    private GATAParams params;
    
    //make two sequence rectangles centering both in the frame
	private double refX; //first base of reference sequence
	private double refY;
	private double refW;
	private double refH;
	private double compX;
	private double compY;
	private double compW;
	private double compH;
    
	private Rectangle2D ref;
	private Rectangle2D comp;
	private Rectangle2D lineRec;
	private Rectangle2D selectionRec = new Rectangle2D.Double();
	
    //misc stuff
	private ArrayList shapes = new ArrayList(3);
	private int orientation = 0;
	private double score = 0;
	private int lenAlignments;
	private Console console;
	private int alignFrameWidth;  //dimensions of JPanel and Frame
	private int alignFrameHeight;
	private  ToolsFrame tools;
	private JViewport viewPort; //for zooming
	private int widthBuffImage;  //for saving buffered image
	private	int heightBuffImage;  //for saving buffered image
	private boolean drawCompBox = true; //draw blue outline box representing comparative sequence
	
	//params for either ref or comp annotation panels
	private AnnoSpecParams refParams;
    private boolean gffRefPresent;
	private AnnoSpecParams compParams;    
	private boolean gffCompPresent;
 
    private Point mousePressedPoint;
    private boolean selectRecVis = false;
    private int lenVisAligns; //length visable alignments arrayList
    
    public AlignPanel(GATAParams params){
		//get plotterparam object and set variables
		this.params=params;
		setPresenceOfGFFs();
		
    	//get saved objects
    	ArrayList objs = GATAUtil.fetchArrayList(params.getObjectsFile());

        alignments = (Alignment[])objs.get(1);
        lenAlignments = alignments.length;
        
        //get param object and set fields
        AP = (AlignParams)objs.get(0);
        windowSize = (double)AP.getWIN_SIZE();
        ntLenRefSeq = (double)AP.getLENGTH_REFSEQ();
        ntLenCompSeq = (double)AP.getLENGTH_COMPSEQ();
        refRealStart= (double)AP.getSTART_INDEX_REFSEQ();
        compRealStart = (double)AP.getSTART_INDEX_COMPSEQ();
        minScore = (double)AP.getMIN_SCORE();
        maxScore = (double)AP.getMATCH() * windowSize;
        
        //make array of visable alignments
		makeVisableAlignments();

        //get refs
        console = params.getConsole();
        
        //set AP in plotterParams
		params.setAlignParams(AP);
        
        alignFrameHeight = params.getHeight();
        alignFrameWidth = params.getWidth();
             
        //make shapes using width and height from plotterParams
        makeAllShapes(alignFrameWidth, alignFrameHeight);
        
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
    
    
    /**Generates a BufferedImage from the panel.  Can be scaled or cropped to visible.*/
    public  BufferedImage makeBufferedImage(double scale, boolean saveVisible){
		BufferedImage bufferedImage;
		Graphics2D g2buff;
				
    	// create image matching what is shown
    	if (saveVisible){
			widthBuffImage = (int)Math.round(((double)viewPort.getWidth())*scale);
			heightBuffImage = (int)Math.round(((double)viewPort.getHeight())*scale);
			bufferedImage = new BufferedImage(widthBuffImage, heightBuffImage, BufferedImage.TYPE_INT_ARGB);
			g2buff = bufferedImage.createGraphics();
			Point xy = viewPort.getViewPosition();
			int xPt = (int)Math.round(xy.getX()*scale*-1);
			int yPt = (int)Math.round(xy.getY()*scale*-1);
			g2buff.translate(xPt,yPt);    		
    	}
		// create full size image
    	else {
			widthBuffImage = (int)Math.round(alignFrameWidth*scale);
			heightBuffImage = (int)Math.round(alignFrameHeight*scale);
			bufferedImage = new BufferedImage(widthBuffImage, heightBuffImage, BufferedImage.TYPE_INT_ARGB);
			g2buff = bufferedImage.createGraphics();
			g2buff.translate(0,-200*scale); //to kill 200 pix buffer, trailing is also lost
    	}
		g2buff.addRenderingHints(GATAUtil.fetchRenderingHints());
		g2buff.scale(scale, scale);
		g2buff.setColor(Color.WHITE);
		g2buff.fill(new Rectangle2D.Double(0,0,getWidth(), getHeight()+400)); //to create solid background
    	
		//Run through array of visible alignments
		for (int i =0; i<lenVisAligns; i++){
			visableAlignments[i].drawAlignment(g2buff);
		}
		g2buff.setColor(Color.BLUE);
		g2buff.draw(ref);
		if (drawCompBox) g2buff.draw(comp);
		return bufferedImage;
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
    
    public void makeAllShapes(int width, int height){
        //method for zooming
        double scalar;
        A= height*params.getA()+200; //buffer for making larger panel
        B= width*params.getB();
        C= width*params.getC();
        D= height*params.getD();
        E= height*params.getE()+200; //buffer for making larger panel
        refY = A;
        refH = D;
        compY = E;
        compH = D;

        //make switches for whether ref or comp are biggest
        if (ntLenRefSeq >= ntLenCompSeq) {
            ntLenBig = ntLenRefSeq;
            ntLenSmall = ntLenCompSeq;
            scalar = ntLenSmall / ntLenBig;
            refX = B;
            compX = B + (0.5 * (C-(C * scalar)));
            refW = C;
            compW = C * scalar;
        }
        else {
            
            ntLenBig = ntLenCompSeq;
            ntLenSmall = ntLenRefSeq;
            scalar = ntLenSmall / ntLenBig;
            refX = B + (0.5 * (C-(C * scalar)));
            compX = B;
            refW = C * scalar;
            compW = C;
        }
        
        //ntPerPixel
        ntPerPixel = ntLenBig/C;
		params.setNtPerPixel(ntPerPixel);
		params.setPixPerNt(1/ntPerPixel);
        
        //set X zero position for annotation
        if (gffRefPresent)refParams.setNtAtXZero(refRealStart - (refX*ntPerPixel));
		if (gffCompPresent)compParams.setNtAtXZero(compRealStart - (compX*ntPerPixel));
		
        //create big rectangles
        ref = new Rectangle2D.Double(refX,refY,refW,refH);
        comp = new Rectangle2D.Double(compX,compY,compW,compH);
        //set Y's for calculating minimum panel size
        params.setSmallestAlignPanelY((int)(refY-1));
        params.setBiggestAlignPanelY((int)(compY+compH+1));        
        //create shapes for each alignment
        createShapes(alignments, refX, refY, compX, compY, ntLenBig, C, D);
        //create rectangle for mouse line selection
        lineRec = new Rectangle2D.Double(0,refY+refH,C*2,compY-(refY+refH));
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
    
    /**Resets an array of visable alignments*/
    public void makeVisableAlignments(){
    	int len = alignments.length;
    	ArrayList al = new ArrayList(len);
		for (int i =0; i<len; i++){
			score = alignments[i].getScore();
			if (score < minScore || score > maxScore) continue;
			al.add(alignments[i]);
		}
		lenVisAligns = al.size();
		visableAlignments = new Alignment[lenVisAligns];
		al.toArray(visableAlignments);    	
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
            if (ref.contains(pt)){
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
	public Alignment[] getVisableAlignments() {
		return visableAlignments;
	}
	public void setViewPort(JViewport viewport) {
		viewPort = viewport;
	}
	public void setAlignFrameWidth(int i) {
		alignFrameWidth = i;
	}
	public int getAlignFrameHeight() {
		return alignFrameHeight;
	}
	public int getAlignFrameWidth() {
		return alignFrameWidth;
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
	public void setDrawCompBox(boolean drawCompBox) {
		this.drawCompBox = drawCompBox;
	}
	public Alignment[] getAlignments() {
		return alignments;
	}
}