package gata.aligner;

import gata.main.*;

import java.util.*;
import java.io.*;
import java.awt.geom.*;
import java.awt.*;

/**
 * @author nix
 * Class to represent a windowed sub alignment, a portion of a BLAST local alignment.
 * A serialize array of these are saved and made availible for the GATAPlotter Program.*/
public class Alignment implements Serializable, Comparable{
    private String seq1;    //refseq!
    private int startSeq1;
    private int stopSeq1;
    private String seq2;  //compseq
    private int startSeq2;
    private int stopSeq2;
    private int score;
    private double bitScore;
    private String expect;
    private int orientation;
    private final AlignParams AP;
    private int[] coordinates;
    private String oriWord = "+/+";
    private ArrayList shapes = new ArrayList();
    private LocalAlignment LA;
	private double realStartRef;	//start index
	private double realStartComp;	//ditto
	private double relStartSeq1;
	private double relStartSeq2;
	private double relStopSeq1; 
	private double relStopSeq2;
	private int[] relativeCoordinates;
	private final double LAMBDA;
	private final double K;
	private final double EFF_M;
	private final double EFF_N;
	private Color color; //used in shading alignment block rectangle
	private Color lineColor; //might not be used, for drawing +/-
	private Rectangle2D refSeqRec;
	private Rectangle2D compSeqRec;
    private Line2D line;
    private boolean compVisible=true; //switch as to whether to draw comp box and connecting line
    
	//stroke parameters
	private static final BasicStroke LINE_BS = new BasicStroke(1.0f);
	private static final float dash[] = {5.0f};
	private static final BasicStroke LINE_DASH = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,
		  BasicStroke.JOIN_MITER,10.0f,dash,0.0f);

    public Alignment(String s1, int start1, int stop1, String s2, int start2, int stop2, 
    	int aScore, int ori, AlignParams ap, LocalAlignment la){
        AP = ap;
        orientation = ori;
        seq1 = s1;
        startSeq1 = start1;
        stopSeq1 = stop1;
        seq2 = s2;
        startSeq2 = start2;
        stopSeq2 = stop2;
        score = aScore;
        LA = la;
        if (ori==1) oriWord = "+/-";
        coordinates = new int[] {startSeq1, stopSeq1, startSeq2, stopSeq2};	
		realStartRef = AP.getSTART_INDEX_REFSEQ();
		realStartComp = AP.getSTART_INDEX_COMPSEQ();
		//get relative start stop
		relStartSeq1 =1+startSeq1-realStartRef;
		relStartSeq2 =1+startSeq2-realStartComp;
		relStopSeq1 =1+stopSeq1-realStartRef;
		relStopSeq2 =1+stopSeq2-realStartComp;
		relativeCoordinates = new int[] {(int)relStartSeq1,(int)relStopSeq1,(int)relStartSeq2,(int)relStopSeq2};
		if (relStartSeq2>relStopSeq2){
			relativeCoordinates[2]= (int)relStopSeq2;
			relativeCoordinates[3]= (int)relStartSeq2;
    	}
			  
		//calculate bitScore and E
		LAMBDA = AP.getLAMBDA();
		K = AP.getK();
		EFF_M = AP.getEFF_M();
		EFF_N = AP.getEFF_N();
		expect = GATAUtil.calculateEValue(LAMBDA,K,score,EFF_N,EFF_M);
		bitScore = GATAUtil.convertRawScoreToBit(LAMBDA,K,score);
		//set colors based on score and orientation
		double maxScoreColor = (double)AP.getMATCH() * AP.getWIN_SIZE();
		double minScoreColor = (double)AP.getMIN_SCORE();
		int convScore = 200 - (int)((score - minScoreColor) * 200/(maxScoreColor-minScoreColor));
		color = new Color(convScore,convScore,convScore);
		if(orientation==1) lineColor = new Color(255,convScore,convScore);
			
			  
    }
     
    //getter methods
    public int getScore(){return score;}
	public double getBitScore(){return bitScore;}
    public int getOrientation(){return orientation; }
    public int[] getCoordinates(){return coordinates;}
    public ArrayList getShapes(){return shapes; }
    public String getSeq1(){return seq1;}
    public String getSeq2(){return seq2;}
    //setter methods
    public void setSeqs(String s1, String s2){
        seq1=s1;
        seq2=s2;
    }
    
    /**Called every time the alignPanel is repaint() ed*/
    public void drawAlignment(Graphics2D g2){
		//check to see if it's inverse, if so then use dashed line and diff color
		if (orientation == 1) {
			g2.setColor(lineColor);
			g2.fill(refSeqRec);
			if (compVisible){
				g2.fill(compSeqRec);
				g2.setStroke(LINE_DASH);
				g2.draw(line);
				g2.setStroke(LINE_BS);
			}
		}
		else {
			g2.setColor(color);
			g2.fill(refSeqRec);
			if (compVisible){
				g2.fill(compSeqRec);
				g2.draw(line);
			}
		}
    }
    
    public void printAlignment(){
        System.out.println("Sub Alignment");
        System.out.println("        " + startSeq1 + GATAUtil.spaces(seq1,startSeq1) + stopSeq1);
        System.out.println("refseq: " + seq1);
        System.out.println("        " + GATAUtil.genDashes(seq1, seq2));
        System.out.println("cmpseq: " + seq2);
        System.out.println("        " + startSeq2 + GATAUtil.spaces(seq2,startSeq2) + stopSeq2);
        System.out.println(" Score: " + bitScore+"("+score+") bits,  Expect: "+expect+",  Ori: " + oriWord);
        System.out.println();
    }
    
    public String getAlignmentString(){
        return 
        "\nSub Alignment\n"+
        "        " + startSeq1 + GATAUtil.spaces(seq1,startSeq1) + stopSeq1 +"\n"+
        "refseq: " + seq1 +"\n"+
        "        " + GATAUtil.genDashes(seq1, seq2)+"\n"+
        "cmpseq: " + seq2+"\n"+
        "        " + startSeq2 + GATAUtil.spaces(seq2,startSeq2) + stopSeq2+"\n"+
        " Score: " + bitScore+"("+score+") bits,  Expect: "+expect+",  Ori: " + oriWord+"\n";
    }
    
     public String getParentLocAlignString(){
        return LA.getLocalAlignString(startSeq1, seq1.length());
    }
    public void makeShapes(double refX, double refY, double compX,
    double compY, double ntLenBig, double widthC, double heightD){
        //clear shapes array
        shapes.clear();
        //create comp rectangle
        int lenCompSeq = (int)(relStopSeq2 - relStartSeq2);
        
        double xCompSeq;
        if (orientation == 1) {
            //Width
            lenCompSeq = lenCompSeq * -1;
            //X
            xCompSeq = (widthC * (relStopSeq2 -1) / ntLenBig) + compX;
        }
        else {
            xCompSeq = (widthC * (relStartSeq2 -1)/ ntLenBig) + compX;
        }
        lenCompSeq +=1;
        
        //Width
        double widthCompSeq = widthC * lenCompSeq / ntLenBig;
        
        compSeqRec = new Rectangle2D.Double(xCompSeq,compY,widthCompSeq,heightD);
        shapes.add(compSeqRec);
        //create ref rectagle
        //X
        int lenRefSeq = (int)((relStopSeq1 - relStartSeq1) + 1);
        double xRefSeq = (widthC * (relStartSeq1 - 1)/ntLenBig) + refX;
        
        //Width
        double widthRefSeq = widthC * lenRefSeq/ ntLenBig;
        refSeqRec = new Rectangle2D.Double(xRefSeq, refY, widthRefSeq, heightD);
        shapes.add(refSeqRec);
        
        //create line
        //find midpoints
        double midPtComp = ((lenCompSeq/2) * widthC / ntLenBig) + xCompSeq;
        double midPtRef = ((lenRefSeq/2) * widthC / ntLenBig) + xRefSeq;
        line = new Line2D.Double(midPtRef, refY+ heightD, midPtComp, compY);
        shapes.add(line);
    } 
    
    public int compareTo(Object otherObject){
		Alignment other = (Alignment) otherObject;
        	//for sorting first on score, second on position with regards to refseq
        	int x;
        	if (score < other.score) return -1;
        	if (score > other.score) return 1;
        	//scores are equal then do the following
        	if (startSeq1 < other.startSeq1) return -1;
        	if (startSeq1 > other.startSeq1) return 1;
        	return 0;
    }
	public int[] getRelativeCoordinates() {
		return relativeCoordinates;
	}

	public void setCompVisible(boolean compVisible) {
		this.compVisible = compVisible;
	}
	public boolean isCompVisible() {
		return compVisible;
	}
}