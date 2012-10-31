package util.bio.cluster;
import java.awt.*;
import java.awt.geom.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.*;
import javax.imageio.ImageIO;
import javax.swing.*;

import util.gen.*;



/**
 * Helper class for {@link HierarchicalClustering}, the guts of drawing the clustering plot.
 * For something so 'simple,' this was a bitch to write, hardest coding yet!
 *
 */
public class ClusterDrawPanel extends JPanel{
	//fields
	private Line2D[] lines;
	private int numLines;
	private JLabel[] labels;
	private int labelIndex = 0;
	private Cluster[] clusters;
	private HierarchicalClustering vhc;
	//defaults
	private double fileNameLineSpacer = 5;
	private double fileNameSpacer = 15;
	private double startingXCoordinate = 5;
	private double xCoordinateIncrement = 25;
	private Color failColor = Color.RED;
	private Font labelFont = new Font("Dialog",1,12);
	private double minimalCCPercent = 50;
	//set params
	private double panelWidth;
	private double panelHeight;
	private double runningX;
	private double runningY;
	
	
	public ClusterDrawPanel(HierarchicalClustering vhc){
		this.vhc = vhc;
		setLayout(null);
		setBackground(Color.WHITE);
		//get clusters
		ArrayList al = vhc.getClusters();
		clusters = new Cluster[al.size()];
		al.toArray(clusters);
		
		//initialize Line and label arrays
		numLines = (clusters.length-1) * 3;
		lines = new Line2D.Double[numLines];
		labels = new JLabel[(clusters.length *2) -1];
		
		//make file text labels
		double maxLabelWidth = 0;
		for (int i=0; i<clusters.length; i++){
			JLabel name = makeLabel (clusters[i].getName());
			clusters[i].setNameLabel(name);
			labels[labelIndex++] = name;
			double width = name.getPreferredSize().getWidth();
			if (width > maxLabelWidth) maxLabelWidth = width;
		}
		
		//calculate panel width and runningX
		runningX = ((clusters.length-1) * xCoordinateIncrement) + startingXCoordinate + fileNameLineSpacer;
		panelWidth = runningX + maxLabelWidth + fileNameLineSpacer;	
		
		//set file text labels, calculate center Y coordinate
		runningY = 10;
		for (int i=0; i<clusters.length; i++){
			JLabel name = clusters[i].getNameLabel();
			double lableHeight = setLabel (name, runningX, runningY);
			double halfHeightY = (lableHeight/ 2.0) + runningY;
			clusters[i].setNameLabelCenterYCoordinate(halfHeightY);			
			//set new Y coordinate for next label
			runningY = runningY + lableHeight +fileNameSpacer;
		}
		
		//assign panel height
		panelHeight = runningY;
		
		//build lines
		makeConnectorLines();

		//save png?
		if (vhc.savePNG()){
			saveBufferedImage(1, vhc.getPngResultFile());
		}
	}
	
	/**Builds tree and associated corr coeff JLabels.*/
	public void makeConnectorLines(){
		//draw clusters by scanning and spliting parents
		ArrayList clustersAL = new ArrayList();
		clustersAL.add(vhc.getMegaCluster());
		double x = startingXCoordinate;
		int index = 0;
		while (true){ 
			boolean noParents = true;
			for (int i=0; i<clustersAL.size(); i++){
				Cluster parentCluster = (Cluster)clustersAL.get(i);
				//parent found?
				if (parentCluster.getTotalClusters() != 1){
					noParents = false;
					//fetch top and bottom
					Cluster top = parentCluster.getTopCluster();
					Cluster bottom = parentCluster.getBottomCluster();
					
					//get number of files in top
					int numClustersInTop = (int)top.getTotalClusters();
					int numClustersInBottom = (int)bottom.getTotalClusters();
					
					//get position of top
					double upperBound;
					if (numClustersInTop >1){					
						upperBound = fetchClusterCenterY(top);
					}
					else {
						upperBound = top.getNameLabelCenterYCoordinate();
						//draw horizontal connector
						lines[index++] = new Line2D.Double(x, upperBound, runningX -fileNameLineSpacer , upperBound);
					}
					double lowerBound;

					if (numClustersInBottom >1){					
						lowerBound = fetchClusterCenterY(bottom);
					}
					else {
						lowerBound = bottom.getNameLabelCenterYCoordinate();
						//draw horizontal
						lines[index++] = new Line2D.Double(x, lowerBound, runningX -fileNameLineSpacer, lowerBound);
					}
					//draw vertical line
					Point2D topPoint = new Point2D.Double(x, upperBound);
					Point2D bottomPoint = new Point2D.Double(x, lowerBound);
					lines[index++] = new Line2D.Double(topPoint, bottomPoint);
					
					//set points in children for latter reference in drawing connector lines
					top.setParentPoint(topPoint);
					bottom.setParentPoint(bottomPoint);
					
					//draw a connector line to parent?
					Point2D parentPoint = parentCluster.getParentPoint();
					double midpoint = ( (lowerBound-upperBound)/2.0 ) + upperBound;
					if (parentPoint !=null){
						Point2D middlePoint = new Point2D.Double(x,midpoint);
						lines[index++] = new Line2D.Double(parentPoint, middlePoint);
					}
					
					//place label for correlation coefficient
					double num = parentCluster.getCorrelationCoefficient();
					//set R as a percent to shrink size of number
					num = 100 * num * num;
					JLabel cc = makeLabel (""+(int)Math.round(num));
					labels[labelIndex++] = cc;
					if (num < minimalCCPercent) cc.setForeground(failColor);
					Dimension dim = cc.getPreferredSize();
					double halfHeight = dim.getHeight()/2.0;
					cc.setBounds((int)Math.round(x+fileNameLineSpacer), (int)Math.round(midpoint-halfHeight), (int)Math.round(dim.getWidth()),  (int)Math.round(dim.getHeight())); 
				
					//replace parent with children
					clustersAL.remove(i);
					clustersAL.add(i, bottom);
					clustersAL.add(i,top);
					
					x+=xCoordinateIncrement;
				}
			}

			if (noParents){
				break;
			}
		}
	}
	
	/**Returns the 'middle' Y coordinate for the cluster.*/
	public double fetchClusterCenterY(Cluster c){
		//last in line?
		Cluster top = c.getTopCluster();
		if (top == null) return c.getNameLabelCenterYCoordinate();
		Cluster bottom = c.getBottomCluster();
		double upper = fetchTopCenterY(c);
		double lower = fetchBottomCenterY(c);
		//symetrical?
		if (top.getTotalClusters() == bottom.getTotalClusters()) return (lower - upper)/2.0 + upper;
		//asymetrical
		//find top
		double topY = fetchClusterCenterY(top);
		//find bottom
		double bottomY = fetchClusterCenterY(bottom);
		return (bottomY - topY)/2.0 + topY;
	}
	
	/**Recurse through clusters always going higher to find center*/
	public double fetchTopCenterY(Cluster cluster){
		if (cluster.getTotalClusters() == 1) return cluster.getNameLabelCenterYCoordinate();
		return fetchTopCenterY(cluster.getTopCluster());
	}
	
	/**Recurse through clusters always going lower to find center*/
	public double fetchBottomCenterY(Cluster cluster){
		if (cluster.getTotalClusters() == 1) return cluster.getNameLabelCenterYCoordinate();
		return fetchBottomCenterY(cluster.getBottomCluster());
	}
	
	public static double setLabel(JLabel lab, double x, double y) {
		Dimension dim = lab.getPreferredSize();
		lab.setBounds((int)Math.round(x), (int)Math.round(y), (int)Math.round(dim.getWidth()),  (int)Math.round(dim.getHeight()));
		return lab.getHeight();
	}
	
	public JLabel makeLabel (String text){
		JLabel lab = new JLabel(text);
		lab.setFont(labelFont);
		add(lab);
		return lab;
	}
	
	public void paintComponent(Graphics g){
		super.paintComponent(g);
		Graphics2D g2 = (Graphics2D)g;
		//draw lines
		for (int i=0; i<numLines; i++) g2.draw(lines[i]);
	}

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
		//draw lines
		g2.setColor(Color.BLACK);
		for (int i=0; i<numLines; i++) g2.draw(lines[i]);
		//set rendering hints for text
		g2.addRenderingHints(fetchRenderingHints());
		//drawLabels, for some damn reason the JLabels don't draw, something to do with a null layout manger
		for (int i=0; i< labels.length; i++){
			Swing.drawJLabel(g2, labels[i], 1.05f);
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
	
	public double getPanelHeight() {
		return panelHeight;
	}

	public double getPanelWidth() {
		return panelWidth;
	}

	public JLabel[] getLabels() {
		return labels;
	}
}