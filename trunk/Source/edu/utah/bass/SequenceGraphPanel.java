package edu.utah.bass;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.BufferedImage;
import java.io.*;

import javax.imageio.ImageIO;
import javax.swing.JLabel;
import javax.swing.JPanel;

import util.gen.Num;
import util.gen.Swing;

import java.util.*;

public class SequenceGraphPanel extends JPanel{

	//fields
	private static final long serialVersionUID = 1L;
	private SequenceGraph sequenceGraph;
	private Line2D[] graphLines; 
	private Line2D[] matchLines;
	private Line2D[] legendLines;
	private Color matchColor = Color.BLUE;
	private int numBasesTop;
	private int totalNumberBases;
	private JLabel[] bases;
	private JLabel[] legends;
	//private Font baseLabelFont = new Font("Courier",1,36);
	private Font baseLabelFont = new Font("Lucida Sans Typewriter",Font.BOLD,28);
	private Font legendLabelFont = new Font("Dialog",1,18);
	private BasicStroke graphLineStroke = new BasicStroke (3.0f);
	private BasicStroke matchLineStroke = new BasicStroke (2.0f);
	private float legendLineStrokeValue = 3.0f;
	private BasicStroke legendLineStroke = new BasicStroke (legendLineStrokeValue);
	private boolean savePNG;
	private double panelWidth = 0;
	private double panelHeight = 0;
	private double runningX = 0;
	private double runningY = 0;
	private double masterStartX;
	private double masterStartY;
	private int maxLegendLabelHeight;
	private double[] scores;
	private double[] scaledScores;
	private double scoreScalar;
	private char[] bpsTop;
	private char[] bpsBottom;
	private double maxLabelHeight = 0;
	private double maxScore;

	//defaults
	private double topAndSideSpacer = 10;
	private double graphHeight = 150;
	private double graphTextSpacer = 5;
	private double textTextSpacer = 5;
	private double matchSpacer = 0;
	
	private double legendLineSpacer = 2;
	private double legendLineBaseSpacer = 3;
	private double legendLineWidth = 4;


	//constructors
	public SequenceGraphPanel (SequenceGraph sequenceGraph){
		//misc
		savePNG = sequenceGraph.isSavePNG();
		this.sequenceGraph = sequenceGraph;
		setLayout(null);
		setBackground(Color.white);
		
		
		//initialize Line and label arrays
		numBasesTop = sequenceGraph.getSequenceTop().length;
		if (sequenceGraph.getSequenceBottom() != null) totalNumberBases = numBasesTop * 2;
		else totalNumberBases = numBasesTop;		
		bases = new JLabel[totalNumberBases];
		bpsTop = sequenceGraph.getSequenceTop();
		bpsBottom = sequenceGraph.getSequenceBottom();
		
		//scalar
		scaleScores();

		//make labels and calculate maxLabelHeight
		makeBaseLabels();
		
		//make legend
		makeLegend();
		
		//set top base labels and make and set top graph lines
		runningY = masterStartY + graphHeight + graphTextSpacer;
		double graphStartY = masterStartY + graphHeight;
		ArrayList<Line2D.Double> linesAL = new ArrayList<Line2D.Double>();
		for (int i=0; i<numBasesTop; i++){
			//get it's dimensions
			Dimension dim = bases[i].getPreferredSize();
			double width = dim.getWidth();
			double height = dim.getHeight();
			//increment panelWidth
			panelWidth += width;
			//set label in panel
			bases[i].setBounds((int)Math.round(runningX), (int)Math.round(runningY), (int)Math.round(width),  (int)Math.round(height));
			runningX += width;
			//make line above label centered on label x,y start; x,y end?
			if (scores[i] > 0){
				double x= runningX - (width/2);
				double yEnd = graphStartY - scaledScores[i];
				linesAL.add(new Line2D.Double (x, graphStartY, x, yEnd));
			}
		}

		//set bottom base labels and make and set bottom graph lines?
		if (numBasesTop != totalNumberBases){
			ArrayList<Line2D.Double> matchLinesAL = new ArrayList<Line2D.Double>();
			runningX = masterStartX;
			runningY = masterStartY + graphHeight + graphTextSpacer + maxLabelHeight+ textTextSpacer;
			graphStartY = runningY + maxLabelHeight + graphTextSpacer;
			int index = 0;
			double matchLineLength = textTextSpacer - (matchSpacer *2);
			double yMatchBegin = runningY - matchSpacer;
			double yMatchStop = yMatchBegin - matchLineLength;
			
			for (int i=numBasesTop; i<totalNumberBases; i++){
				//get it's dimensions
				Dimension dim = bases[i].getPreferredSize();
				double width = dim.getWidth();
				double height = dim.getHeight();
				//set label in panel
				bases[i].setBounds((int)Math.round(runningX), (int)Math.round(runningY), (int)Math.round(width),  (int)Math.round(height));
				runningX += width;
				
				//make line below label centered on label x,y start; x,y end
				double x= runningX - (width/2);
				if (scores[i] > 0){
					double yEnd = graphStartY + scaledScores[i];
					linesAL.add( new Line2D.Double (x, graphStartY, x, yEnd));
				}
				
				//make line above label indicating a match
				if (basesMatch(bpsTop[index], bpsBottom[index])){
					matchLinesAL.add(new Line2D.Double (x, yMatchBegin, x, yMatchStop));
				}
				//increment index
				index++;
			}
			//convert matches
			matchLines = new Line2D.Double[matchLinesAL.size()];
			matchLinesAL.toArray(matchLines);
		}

		//convert graphLines
		graphLines = new Line2D.Double[linesAL.size()];
		linesAL.toArray(graphLines);

		//set panel size
		setSize((int)panelWidth, (int)panelHeight);

		//save png?
		if (savePNG){
			saveBufferedImage(1, sequenceGraph.getPngFile());
		}

	}


	//methods
	/**Saves a png for this panel, don't paint(g2), this requires a head and will throw a headless error.*/
	public void saveBufferedImage(double scale, File file){
		setBounds(0,0, new Double(panelWidth).intValue(), new Double(panelHeight).intValue());
		double width = (getPanelWidth())*scale; 
		double height = (getPanelHeight())*scale;
		
		BufferedImage bufferedImage = new BufferedImage ((int)Math.round(width), (int)Math.round(height), BufferedImage.TYPE_INT_RGB);
		Graphics2D g2 = bufferedImage.createGraphics();
		//set scale
		if (scale >1 ) g2.scale(scale,scale);
		//make and set background rectangle 
		Rectangle2D rec = new Rectangle2D.Double(0,0,width,height);
		g2.setColor(Color.WHITE);
		g2.fill(rec);
		
		//draw graph lines
		g2.setStroke(graphLineStroke);
		g2.setColor(Color.BLACK);
		for (int i=0; i<graphLines.length; i++) g2.draw(graphLines[i]);
		
		//draw legend lines
		g2.setStroke(legendLineStroke);
		g2.setColor(matchColor);
		for (int i=0; i< legendLines.length; i++) {
			g2.draw(legendLines[i]);
		}

		//draw match lines?
		if (matchLines!=null){
			g2.setStroke(matchLineStroke);
			g2.setColor(matchColor);
			for (int i=0; i<matchLines.length; i++) g2.draw(matchLines[i]);
		}
		
		//set rendering hints for text
		g2.addRenderingHints(fetchRenderingHints());
		
		//drawLabels, for some damn reason the JLabels don't draw, something to do with a null layout manger, must tune the offset to match the size of the letters! Hmm...
		for (int i=0; i< bases.length; i++){
			Swing.drawJLabel(g2, bases[i], 1.4f);
		}
		
		for (int i=0; i< legends.length; i++){
			Swing.drawJLabel(g2, legends[i], 1.5f);
		}
		
		g2.dispose();
		try{
			ImageIO.write(bufferedImage, "png", file);
		} catch (IOException e){
			e.printStackTrace();
		}
	}
	
	/** Map of hints for high quality/ slow rendering*/
	public static RenderingHints fetchRenderingHints () {
		RenderingHints hints = new RenderingHints(null);
		hints.put(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		hints.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		hints.put(RenderingHints.KEY_DITHERING, RenderingHints.VALUE_DITHER_DISABLE);
		hints.put(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		hints.put(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);
		hints.put(RenderingHints.KEY_ALPHA_INTERPOLATION, RenderingHints.VALUE_ALPHA_INTERPOLATION_QUALITY);
		hints.put(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY);
		return hints;
	}
	
	public void makeLegend(){
		double halfShift = 0;
		if (savePNG == false) {
			halfShift = legendLineStrokeValue/2.0f;
			legendLineSpacer += halfShift;
		}
		
		//how many to make?
		if (bpsBottom == null) {
			legends = new JLabel[4];
			legendLines = new Line2D.Double[4];
		}
		else {
			legends = new JLabel[7];
			legendLines = new Line2D.Double[8];
		}

		//first make labels
		String maxScoreLabel = Integer.toString((int)Math.round(maxScore));
		String halfMaxScoreLabel = Integer.toString((int)Math.round(maxScore/2));
		String[] s = new String[]{maxScoreLabel, halfMaxScoreLabel, "0", "0", halfMaxScoreLabel, maxScoreLabel};
		for (int x=0; x< legends.length-1; x++)legends[x]= makeLegendLabel(s[x]);
		legends[legends.length-1] = makeLegendLabel("Scaled to % editing at 20% A->I with "+sequenceGraph.getMatrixName()+ " preferences.");
		
		//find max size set placement
		Dimension dim = legends[0].getPreferredSize();
		int maxWidth = (int)dim.getWidth();
		maxLegendLabelHeight = (int)dim.getHeight();
		
		//place 1st label 100 and line
		double startX = topAndSideSpacer+ maxWidth + legendLineSpacer - halfShift;
		double x100 = topAndSideSpacer;
		legends[0].setBounds((int)x100, (int) topAndSideSpacer, maxWidth, maxLegendLabelHeight);		
		masterStartY = maxLegendLabelHeight/2 +topAndSideSpacer;
		masterStartX = startX + legendLineWidth+ legendLineBaseSpacer;
		double endX = startX+legendLineWidth + halfShift;
		legendLines[0] = new Line2D.Double(startX, masterStartY, endX, masterStartY);
		
		//place 2nd label 50 and line
		double y =  masterStartY+(graphHeight/2);
		legendLines[1] = new Line2D.Double(startX, y, endX, y);
		dim = legends[1].getPreferredSize();
		double y50 = y - (dim.getHeight()/2);
		double x50 = topAndSideSpacer + maxWidth - dim.getWidth();
		legends[1].setBounds((int)x50, (int)y50, (int)dim.getWidth(), (int)dim.getHeight());
		
		//place 3rd label 0 and line
		y= masterStartY+graphHeight;
		legendLines[2] = new Line2D.Double(startX, y, endX, y); 
		dim = legends[2].getPreferredSize();
		double y0 = y - (dim.getHeight()/2);
		double x0 = topAndSideSpacer + maxWidth - dim.getWidth();
		legends[2].setBounds((int)x0, (int)y0, (int)dim.getWidth(), (int)dim.getHeight());
		
		//add connector line
		legendLines[3] = new Line2D.Double(startX + halfShift, masterStartY, startX + halfShift, masterStartY+graphHeight); 
		
		//increment panel width
		panelWidth += (masterStartX + topAndSideSpacer);
		runningX = masterStartX;
		
		//make bottom legend
		if (bpsBottom == null) {
			y = y0 + dim.getHeight()+ maxLabelHeight+ textTextSpacer;
			dim = legends[3].getPreferredSize();
			legends[3].setBounds((int)topAndSideSpacer, (int)y, (int)dim.getWidth(), (int)dim.getHeight());
			panelHeight = y + dim.getHeight() + topAndSideSpacer;
			return;
		}
		
		//place 0 bottom label and line
		y = y + (graphTextSpacer*2) +(maxLabelHeight*2) + textTextSpacer;
		legendLines[4] = new Line2D.Double(startX, y, endX, y);
		dim = legends[3].getPreferredSize();
		y0 = y - (dim.getHeight()/2);
		legends[3].setBounds((int)x0, (int)y0, (int)dim.getWidth(), (int)dim.getHeight());
		
		//add connector
		legendLines[5] = new Line2D.Double(startX + halfShift, y, startX + halfShift, y+graphHeight);
		
		//place 50 bottom label and line
		y = y+ graphHeight/2;
		legendLines[6] = new Line2D.Double(startX, y, endX, y);
		dim = legends[4].getPreferredSize();
		y50 = y - (dim.getHeight()/2);
		legends[4].setBounds((int)x50, (int)y50, (int)dim.getWidth(), (int)dim.getHeight());
		
		//place 100 bottom label
		y = y+ graphHeight/2;
		legendLines[7] = new Line2D.Double(startX, y, endX, y);
		dim = legends[5].getPreferredSize();
		double y100 = y - (dim.getHeight()/2);
		legends[5].setBounds((int)x100, (int)y100, (int)dim.getWidth(), (int)dim.getHeight());
		
		//place legend
		y = y100 + dim.getHeight()+ textTextSpacer;
		dim = legends[6].getPreferredSize();
		legends[6].setBounds((int)topAndSideSpacer, (int)y, (int)dim.getWidth(), (int)dim.getHeight());
		panelHeight = y + dim.getHeight() + topAndSideSpacer;
	}
	
	public void makeBaseLabels(){
		for (int i=0; i<numBasesTop; i++){
			//make the label
			bases[i] = makeBaseLabel (bpsTop[i]);
			//shade to indicate unscored?
			if (scores[i] == -1){
				bases[i].setForeground(Color.LIGHT_GRAY);
			}
			//get it's dimensions
			Dimension dim = bases[i].getPreferredSize();
			double height = dim.getHeight();
			//look for max height
			if (maxLabelHeight < height) maxLabelHeight = height;
		}
		
		if (numBasesTop != totalNumberBases){
			int index = 0;
			for (int i=numBasesTop; i<totalNumberBases; i++){
				//make the label
				bases[i] = makeBaseLabel (bpsBottom[index++]);
				//shade to indicate unscored?
				if (scores[i] == -1){
					bases[i].setForeground(Color.LIGHT_GRAY);
				}
				//get it's dimensions
				Dimension dim = bases[i].getPreferredSize();
				double height = dim.getHeight();
				//look for max height
				if (maxLabelHeight < height) maxLabelHeight = height;
			}
		}
	}
	
	public boolean basesMatch(char a, char b){
		char x = Character.toLowerCase(a);
		char y = Character.toLowerCase(b);
		if (x == y) return false;
		if (x == 'g' && y == 'c') return true;
		if (x == 'c' && y == 'g') return true;
		if (x == 'a' && (y == 't' || y =='u')) return true;
		if ((x == 't' || x == 'u') && y == 'a') return true;
		return false;
	}
	
	public void scaleScores(){
		double[] scoresTop = sequenceGraph.getScoresTop();
		double[] scoresBottom = sequenceGraph.getScoresBottom();
		if (scoresBottom == null) scores = scoresTop;
		else {
			scores = new double[totalNumberBases];
			System.arraycopy(scoresTop, 0, scores, 0, numBasesTop);
			System.arraycopy(scoresBottom, 0, scores, numBasesTop, numBasesTop);
		}
		maxScore = Num.findMinMaxDoubleValues(scores)[1];
		scoreScalar = graphHeight / maxScore;
		scaledScores = new double[totalNumberBases];
		for (int i=0; i< totalNumberBases; i++){
			if (scores[i] == -1) scaledScores[i] = 0;
			else scaledScores[i] = scores[i] * scoreScalar;
		}
	}

	public JLabel makeBaseLabel (char text){
		JLabel lab = new JLabel(Character.toString(text));
		lab.setFont(baseLabelFont);
		add(lab);
		return lab;
	}
	
	public JLabel makeLegendLabel (String text){
		JLabel lab = new JLabel(text);
		lab.setFont(legendLabelFont);
		lab.setForeground(matchColor);
		add(lab);
		return lab;
	}
	
	public void paintComponent(Graphics g){
		super.paintComponent(g);
		Graphics2D g2 = (Graphics2D)g;
				
		//draw graph lines
		g2.setStroke(graphLineStroke);
		for (int i=0; i<graphLines.length; i++) g2.draw(graphLines[i]);
		
		//draw legend lines
		g2.setStroke(legendLineStroke);
		g2.setColor(matchColor);
		for (int i=0; i< legendLines.length; i++) {
			g2.draw(legendLines[i]);
		}

		//draw match lines?
		if (matchLines!=null){
			g2.setStroke(matchLineStroke);
			g2.setColor(matchColor);
			for (int i=0; i<matchLines.length; i++) g2.draw(matchLines[i]);
		}
	}


	public double getPanelWidth() {
		return panelWidth;
	}


	public double getPanelHeight() {
		return panelHeight;
	}


}
