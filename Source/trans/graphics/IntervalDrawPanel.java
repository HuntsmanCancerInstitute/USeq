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
import trans.main.*;
import util.bio.seq.*;
import util.gen.*;


/**
 * Panel with the meat for {@link IntervalPlotter}
 *
 */
public class IntervalDrawPanel extends JPanel{
	//fields
	private Font font = new Font("Dialog",Font.PLAIN,10);
	private FontRenderContext context;
	private double h = 50;		//height of graph
	private double w = 101;		//width of graph
	private double x = 35;		//X-coor top left
	private double y = 30;		//y-coor top left
	private double ySpacer = 10;	//spacer to put between graphs
	private int peakFlanks = 5;	//add bp to either side of bp peak;
	private Graph[] intGraphs;
	private Graph[] matrixMismatchGraphs;
	private Box[] boxes;
	private Box[] bindingPeaks = null;
	private int numIntGraphs;
	private int numMatrixMismatchGraphs;
	private int numBoxes;
	private int numBindingPeaks =0;
	private String title;
	private String chromosome;
	private double runningY;				//furthest south
	private double runningX;		//furthest east
	private TextFrame textFrame;
	private Interval interval;
	private int intervalNumber;
	
	//for selection rectangle
	private Rectangle2D selectionRec = new Rectangle2D.Double();
	private Point mousePressedPoint;
	private boolean selectRecVis = false;
	
	public IntervalDrawPanel(Interval interval,int intervalNumber, TextFrame textFrame, boolean showPSPM){
		this.textFrame = textFrame;
		this.interval = interval;
		this.intervalNumber = intervalNumber;
		//assuming that interval is fully loaded
		chromosome = interval.getChromosome();
		title = intervalNumber+" "+chromosome+": "+interval.getStart1stOligo()+"-"+ 
		(interval.getStartLastOligo()+interval.getSizeOfOligoMinusOne())+ ", "+(interval.getStartLastOligo()+interval.getSizeOfOligoMinusOne()-
				interval.getStart1stOligo())+"bp";
		ArrayList graphsAL = new ArrayList();
		
		//make treatment graphs
		Oligo[] oligos = interval.getOligos();
		int numOligos = oligos.length;
		
		float[][] treatments = new float[numOligos][];
		double[] positions = new double[numOligos];
		double[] matches = new double[numOligos];
		double[] aveTreatments = new double[numOligos]; 
		float[][] controls = new float[numOligos][];
		double[] aveControls = new double[numOligos]; 
		
		for (int i=0; i<numOligos; i++){
			treatments[i] = oligos[i].getTreatmentIntensities(interval.getNumberTreatmentIntensities());
			positions[i] = oligos[i].getStart();
			matches[i] = oligos[i].getMatches();
			aveTreatments[i] = Num.averageFloatArray(treatments[i]);
			controls[i] = oligos[i].getControlIntensities(interval.getNumberControlIntensities());
			aveControls[i] = Num.averageFloatArray(controls[i]);
		}
		
		//calculate minMax
		float[] minMaxTreatments = Num.findMinMaxFloatArrays(treatments);
		float[] minMaxControls = Num.findMinMaxFloatArrays(controls);
		double[] minMax = {minMaxTreatments[0], minMaxTreatments[1]};
		if (minMax[0]>minMaxControls[0]) minMax[0]=minMaxControls[0];
		if (minMax[1]<minMaxControls[1]) minMax[1]=minMaxControls[1];
		
		//calculate width, 
		//double width = positions[numOligos-1]-positions[0];
		//w = width/2;
		w = numOligos * 12;
		runningX = w + x;
		
		Graph graph;
		//make treatment average
		graph = new Graph(aveTreatments, positions, w, h, Color.RED);
		graph.setMinMaxY(minMax);
		graphsAL.add(graph);
		
		//make control average
		graph = new Graph(aveControls, positions, w, h, Color.GREEN);
		graph.setMinMaxY(minMax);
		graphsAL.add(graph);
		
		//make average difference
		double[] differences = new double[numOligos];
		for (int i=0; i<numOligos; i++){
			differences[i] = aveTreatments[i]-aveControls[i];
		}
		graphsAL.add(new Graph(differences, positions, w, h, Color.YELLOW));
		
		//make fold difference
		double[] ratios = new double[numOligos];
		for (int i=0; i<numOligos; i++){
			if (aveControls[i] !=0) ratios[i] = Math.log(aveTreatments[i]/aveControls[i])/Num.log2;
			else ratios[i] = 0 ;
		}
		graphsAL.add(new Graph(ratios, positions, w, h, Color.BLUE));
		
		//make smoothed scores	
		double[] smoothedScores = new double[numOligos];
		for (int k=0; k<numOligos; k++){
			smoothedScores[k] = oligos[k].getSmoothedScore();
		}
		graphsAL.add(new Graph(smoothedScores, positions, w, h, Color.LIGHT_GRAY));
		
		//make Binding Peaks, ignoring length of oligo for plotting purposes
		BindingPeak[] peaks = interval.getBindingPeaks();
		if (peaks != null){
			numBindingPeaks = peaks.length * 2;
			bindingPeaks = new Box[numBindingPeaks];
			int counter =0;
			for (int k=0; k<numBindingPeaks; k++){
				//(double baseStart, double baseStop, double baseFirstPosition, Color color)
				int bpLeft = oligos[peaks[counter].getLeftFlankingIndex()].getStart();
				int bpRight =  oligos[peaks[counter].getRightFlankingIndex()].getStart();
				//make region
				bindingPeaks[k] = new Box(bpLeft, bpRight, positions[0], Color.DARK_GRAY);
				//make peak +/- peakFlanks bp
				k++;
				int modPeakBp = peaks[counter].getPeakBP() - (int)Math.round(((double)interval.getSizeOfOligoMinusOne())/2);
				bpLeft = modPeakBp -peakFlanks;
				bpRight = modPeakBp +peakFlanks;
				bindingPeaks[k] = new Box(bpLeft, bpRight, positions[0], Color.RED);
				//
				counter++;
			}
		}
		
		//make best windows
		SubWindow sub = interval.getBestSubWindow();
		if (sub!= null) boxes = new Box[2];
		else {
			boxes = new Box[1];
		}
		trans.main.Window win = interval.getBestWindow();
		boxes[0] = new Box(win.getStart1stOligo(), win.getStartLastOligo(), positions[0], Color.ORANGE);
		
		//make best sub window, if it exists
		if (sub!=null){
			Oligo[] subOligos = sub.getOligos();
			boxes[1] = new Box(subOligos[0].getStart(), subOligos[subOligos.length-1].getStart(), positions[0], Color.PINK);
		}
		numBoxes = boxes.length;
		
		//make matrix hits and mismatch graphs
		if (showPSPM) {
			matrixMismatchGraphs = new Graph[2];
			int numBases = interval.getSizeOfOligoMinusOne()+ interval.getStartLastOligo()-interval.getStart1stOligo();
			double[] bases = new double[numBases];
			int base = interval.getStart1stOligo();
			for (int i= 0; i<numBases; i++) bases[i] = base++;
			matrixMismatchGraphs[0] = new Graph(interval.getBaseScores(), bases, w, h, Color.CYAN);
			matrixMismatchGraphs[1] = new Graph(matches, positions, w, h, Color.WHITE);
			numMatrixMismatchGraphs = 1;
		}
		else {
			matrixMismatchGraphs = new Graph[1];
			matrixMismatchGraphs[0] = new Graph(matches, positions, w, h, Color.WHITE);
		}
		numMatrixMismatchGraphs = matrixMismatchGraphs.length;		
		
		
		//make treatment graphs
		for (int i=0; i<interval.getNumberTreatmentIntensities(); i++){
			double[] intensities = new double[numOligos];
			for (int j=0; j<numOligos; j++){
				intensities[j]= treatments[j][i];
			}			
			graph = new Graph(intensities, positions, w, h, Color.RED);
			graph.setMinMaxY(minMax);
			graphsAL.add(graph);
		}
		//make control graphs
		for (int i=0; i<interval.getNumberControlIntensities(); i++){
			double[] intensities = new double[numOligos];
			for (int j=0; j<numOligos; j++){
				intensities[j]= controls[j][i];
			}
			graph = new Graph(intensities, positions, w, h, Color.GREEN);
			graph.setMinMaxY(minMax);
			graphsAL.add(graph);
		}
		
		//convert ArrayLists to arrays
		numIntGraphs = graphsAL.size();
		intGraphs = new Graph[numIntGraphs];
		graphsAL.toArray(intGraphs);
		
		//make initial lines  (double topLeftXCoor, double topLeftYCoor, double widthScalar, double heigthScalar)  returns bottom y Coordinate
		runningY = y;
		for (int i=0; i<5; i++){
			runningY = intGraphs[i].makeGlyphs(x, runningY,1,1) + ySpacer;
		}
		
		//make boxes to represent binding regions, print on one line
		double height = 0;
		for (int i=0; i<numBindingPeaks; i++){
			height = bindingPeaks[i].makeGlyphs(intGraphs[0], runningY);
		}
		if (height !=0 ) runningY = height + ySpacer;
		
		//make boxes to represent the best window and best sub window
		for (int i=0; i<numBoxes; i++){
			runningY = boxes[i].makeGlyphs(intGraphs[0], runningY) + ySpacer;
		}
		//make matrix mismatch graphs
		for (int i=0; i<numMatrixMismatchGraphs; i++){
			runningY = matrixMismatchGraphs[i].makeGlyphs(x, runningY,1,1) + ySpacer;
		}
		for (int i=5; i<numIntGraphs; i++){
			runningY = intGraphs[i].makeGlyphs(x, runningY,1,1) + ySpacer;
		}
		
		//add listeners to monitor box alignment selector and base position
		addMouseListener(new MouseHandler());
		addMouseMotionListener(new MouseMotionHandler());
	}
	
	/**Uses pixelsPerBase and X from first Graph to convert a given X Cooridinate to base pairs.*/
	public int convertPixelXCoorToBP(double x){
		double deltaPixels = x-intGraphs[0].getX();;
		double deltaBP = deltaPixels * (1/ intGraphs[0].getPixelsPerBase());
		return (int)Math.round(intGraphs[0].getPositions()[0] + deltaBP);
	}
	
	/**Finds the closes Oligo to a bp point and returns the indexed score.*/
	public static double fetchClosestOligoScore(Oligo[] oligos, int position){
		int num = oligos.length;
		int distance = 1000000;
		int index = -1; //will throw array out of bounds exception if not corrected
		int testDistance;
		for (int i=0; i<num; i++){
			testDistance = Math.abs(position-oligos[i].getStart());
			if (testDistance < distance) {
				distance = testDistance;
				index = i;
			}
		}
		return oligos[index].getSmoothedScore();
	}
	
	/**Scans oligos for a max score when given a bp range.*/
	public static double fetchMaxOligoScore(Oligo[] oligos, int startBp, int stopBp){
		int num = oligos.length;
		double maxScore = 0;
		for (int i=0; i<num; i++){
			//within region?
			int pos = oligos[i].getStart();
			if ((pos >=startBp) && (pos <= stopBp)){
				//check if better than max
				if (oligos[i].getSmoothedScore()>maxScore) maxScore = oligos[i].getSmoothedScore();
			}
		}
		return maxScore;
	}
	
	
	private class MouseHandler extends MouseAdapter {
		private StringBuffer line;
		private int start;
		private int stop;
		
		public void mousePressed(MouseEvent event){
			mousePressedPoint = event.getPoint();
			selectRecVis = true;
			selectionRec.setFrame(0,0,0,0);
		}
		public void mouseReleased(MouseEvent event){
			start = convertPixelXCoorToBP(mousePressedPoint.getX());
			stop = convertPixelXCoorToBP(event.getPoint().getX());
			if (stop<start){
				int holder = start;
				start = stop;
				stop = holder; 
			}
			//rank
			line = new StringBuffer(intervalNumber+" ");
			//best subWindow score
			line.append(" ");
			line.append(Num.formatNumberOneFraction(interval.getBestSubWindow().getMedianRatio()));
			//closest smoothed score to click point if just one or
			line.append(" ");
			if (stop == start) {
				double score = fetchClosestOligoScore(interval.getOligos(), start);
				line.append(Num.formatNumberOneFraction(score));
			}
			//max smoothed score withing click and release
			else {
				double score = fetchMaxOligoScore(interval.getOligos(), start, stop);
				line.append(Num.formatNumberOneFraction(score));
			}
			line.append(" ");
			line.append(chromosome);
			line.append(" ");
			line.append(start);
			if ((stop-start)>0){
				line.append(" ");
				line.append(stop);
				if (Misc.isNotEmpty(interval.getSequence())){
					line.append(" ");
					line.append(Seq.fetchSubSequence(start, stop, interval.getStart1stOligo(), interval.getSequence()));
				}
			}
			line.append ("\n");
			textFrame.append(line.toString());
			//selectRecVis = false;
			repaint();
		}
	}
	private class MouseMotionHandler implements MouseMotionListener{
		public void mouseMoved (MouseEvent event){
		}
		public void mouseDragged(MouseEvent mouseEvent) {
			selectionRec.setFrameFromDiagonal(mousePressedPoint, mouseEvent.getPoint());
			repaint();
		}   
	}
	
	public void paintComponent(Graphics g){
		super.paintComponent(g);
		Graphics2D g2 = (Graphics2D)g;
		g2.setFont(font);
		context = g2.getFontRenderContext();
		//draw title
		g2.setColor(Color.WHITE);
		g2.drawString(title, 5,10);
		//draw intensity graphs
		for (int i=0; i<numIntGraphs; i++) intGraphs[i].draw(g2,context, font);
		//draw binding peaks
		for (int i=0; i<numBindingPeaks; i++) bindingPeaks[i].draw(g2);
		//draw boxes
		for (int i=0; i<numBoxes; i++) boxes[i].draw(g2);
		//draw matrix mismatch graphs
		for (int i=0; i<numMatrixMismatchGraphs; i++) matrixMismatchGraphs[i].draw(g2,context, font);
		//draw selection rectangle
		if (selectRecVis){
			g2.setColor(Color.CYAN);
			g2.draw(selectionRec);
		}
	}
	
	public  void saveBufferedImage(double scale, File file, boolean antiAlias){
		BufferedImage bufferedImage=null;
		Graphics2D g2=null;
		double width = (runningX+10)*scale; 
		double height = (runningY+10)*scale;
		
		bufferedImage = new BufferedImage ((int)Math.round(width), (int)Math.round(height), BufferedImage.TYPE_INT_ARGB);
		g2 = bufferedImage.createGraphics();
		g2.addRenderingHints(GATAUtil.fetchRenderingHints());
		//kill antialiasing to make sharp lines and text for computer display, not good for printed figures
		if (antiAlias==false){
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
			g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF);
		}
		if (scale >1 ) g2.scale(scale,scale);
		g2.setColor(Color.BLACK);
		g2.fill(new Rectangle2D.Double(0,0,width+50, height+50));
		//draw glyphs
		g2.setFont(font);
		context = g2.getFontRenderContext();
		//draw title
		g2.setColor(Color.WHITE);
		g2.drawString(title, 5,10);
		//draw  intensity graphs
		for (int i=0; i<numIntGraphs; i++) intGraphs[i].draw(g2,context, font);
		//draw binding peaks
		for (int i=0; i<numBindingPeaks; i++) bindingPeaks[i].draw(g2);
		//draw boxes
		for (int i=0; i<numBoxes; i++) boxes[i].draw(g2);
		//draw matrix mismatch graphs
		for (int i=0; i<numMatrixMismatchGraphs; i++) matrixMismatchGraphs[i].draw(g2,context, font);
		
		try{
			ImageIO.write(bufferedImage, "PNG", file);
		} catch (IOException e){
			e.printStackTrace();
		}
	}
	
	public double getRunningX() {
		return runningX;
	}
	public double getRunningY() {
		return runningY;
	}
	public String getTitle() {
		return title;
	}
}