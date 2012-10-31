package trans.graphics;
import gata.main.GATAUtil;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.BufferedImage;
import java.io.*;
import javax.imageio.ImageIO;
import javax.swing.JPanel;
import java.util.*;

import util.gen.*;

/**Panel and calcs for {@link VirtualCel}.*/
public class VirtualCelPanel extends JPanel{
	//fields
	private VirtualCelFrame frame;
	private Color[] colors;
	private double colorScalar;
	private float[][] intensities;
	private int lengthX;
	private int lengthY;
	private BufferedImage bufferedImage;
	private Ellipse2D[] maskedEllipses = null;
	private double zoom = 1;
	private boolean writePNG = false;
	private File pngFile;
	private int maxIntensityValue = 20000;
	private float maskValue = 0;
	
	//for selection rectangle
	private Ellipse2D selectionEllipse = new Ellipse2D.Float();
	private Point mousePressedPoint;
	private boolean selectionVisible = false;
	
	public VirtualCelPanel (float[][] intensities, int maxIntensityValue, float maskValue, VirtualCelFrame frame){
		this.intensities = intensities;
		this.frame = frame;
		this.maxIntensityValue = maxIntensityValue;
		this.maskValue = maskValue;
		lengthX = intensities.length;
		lengthY = intensities[0].length;
		
		//make array of colors
		colors = makeColors();
		
		//calculate scalor to use in converting intensity to color index
		float[] minMax = Num.findMinMaxFloatArrays(intensities);
		if (minMax[1] > maxIntensityValue) minMax[1] = maxIntensityValue;
		float range = minMax[1]-minMax[0];
		colorScalar = (double)colors.length / (double)range;
		
		//make BufferedImage
		makeBufferedImage();
		
		//add listeners to monitor elipse selector and keyboard
		addMouseListener(new MouseHandler());
		addMouseMotionListener(new MouseMotionHandler());
		addKeyListener(new KeyHandler());
		setFocusable(true);
		
	}
	
	/**For printing pngs*/
	public VirtualCelPanel (float[][] intensities, int maxIntensityValue, File pngFile){
		this.intensities = intensities;
		writePNG = true;
		this.pngFile = pngFile;
		lengthX = intensities.length;
		lengthY = intensities[0].length;
		
		//make array of colors
		colors = makeColors();
		
		//calculate scalor to use in converting intensity to color index
		float[] minMax = Num.findMinMaxFloatArrays(intensities);
		if (minMax[1] > maxIntensityValue) minMax[1] = maxIntensityValue;
		float range = minMax[1]-minMax[0];
		colorScalar = (double)colors.length / (double)range;
		
		//make BufferedImage
		makeBufferedImage();
	}
	
	private class MouseHandler extends MouseAdapter {
		public void mousePressed(MouseEvent event){
			mousePressedPoint = event.getPoint();	
			selectionVisible = true;
		}
		public void mouseClicked(MouseEvent event){
			//double clicked?
			if (event.getClickCount()>1){
				//set new maskedEcllipse
				if (maskedEllipses == null){
					maskedEllipses = new Ellipse2D[1];
					maskedEllipses[0] = new Ellipse2D.Double(selectionEllipse.getX(), selectionEllipse.getY(), 
							selectionEllipse.getWidth(), selectionEllipse.getHeight());
				}
				else {
					Ellipse2D[] addOne = new Ellipse2D[maskedEllipses.length + 1];
					addOne[maskedEllipses.length] = new Ellipse2D.Double(selectionEllipse.getX(), selectionEllipse.getY(), 
							selectionEllipse.getWidth(), selectionEllipse.getHeight());
					for (int i=0; i< maskedEllipses.length; i++){
						addOne[i] = maskedEllipses[i];
					}
					maskedEllipses = addOne;
				}
				selectionVisible = false;
				repaint();
			}
		}
	}
	private class MouseMotionHandler implements MouseMotionListener{
		public void mouseMoved (MouseEvent event){
		}
		public void mouseDragged(MouseEvent mouseEvent) {
			//reset points based on scalar
			double xP = mousePressedPoint.getX() / zoom;
			double yP = mousePressedPoint.getY() / zoom;
			double xR = mouseEvent.getPoint().getX() / zoom;
			double yR = mouseEvent.getPoint().getY() / zoom;		
			selectionEllipse.setFrameFromDiagonal(xP,yP,xR,yR);
			repaint();
		}   
	}
	private class KeyHandler extends KeyAdapter {
		public void keyPressed(KeyEvent event){
			int keyCode = event.getKeyCode();
			if (event.isShiftDown()){
				//arrows
				if (keyCode == KeyEvent.VK_LEFT){modifySelection(-1,0);}
				else if (keyCode == KeyEvent.VK_RIGHT){modifySelection(1,0);}
				else if (keyCode == KeyEvent.VK_UP){modifySelection(0,-1);}
				else if (keyCode == KeyEvent.VK_DOWN){modifySelection(0,1);}
			}
			else {
				//arrows
				if (keyCode == KeyEvent.VK_LEFT){nudgeSelection(-1,0);}
				else if (keyCode == KeyEvent.VK_RIGHT){nudgeSelection(1,0);}
				else if (keyCode == KeyEvent.VK_UP){nudgeSelection(0,-1);}
				else if (keyCode == KeyEvent.VK_DOWN){nudgeSelection(0,1);}
				//zooming
				else if (keyCode == KeyEvent.VK_I){zoom(1);}
				else if (keyCode == KeyEvent.VK_O){zoom(-1);}
				//mask
				else if (keyCode == KeyEvent.VK_M){mask();}
				else if (keyCode == KeyEvent.VK_Q){close();}
			}
		}
		//free resources and return to CelMapper
		public void close(){
			frame.hide();
			setEllipses();
			frame.getCelMasker().writeCelFile();
			frame.dispose();
		}
		
		public void mask(){
			if (maskedEllipses != null){
				setEllipses();
				//redraw
				maskedEllipses = null;
				makeBufferedImage();
				repaint();
			}
		}
		public void setEllipses(){
			if (maskedEllipses != null){
				for (int i=0; i< maskedEllipses.length; i++){
					//get bounds to limit search
					Rectangle2D rec = maskedEllipses[i].getBounds2D();
					int startX = (int)Math.round(rec.getX())-1;
					int stopX = (int)Math.round(rec.getX()+ rec.getWidth())+ 1;
					int startY = (int)Math.round(rec.getY())-1;
					int stopY = (int)Math.round(rec.getY() + rec.getHeight())+ 1;
					//check that they don't overrun array
					if (startX < 0) startX = 0;
					if (stopX > lengthX) stopX = lengthX;
					if (startY < 0) startY = 0;
					if (stopY > lengthY) stopY = lengthY;
					//look for points and set to maskValue, masked value
					for (int x=startX; x< stopX; x++){
						for (int y=startY; y< stopY; y++){
							if (maskedEllipses[i].contains(x,y)) intensities[x][y] = maskValue;
						}
					}
				}
			}
		}
		
		public void zoom(int s){
			if (s == -1){ 
				if (zoom >1) zoom--;
				else zoom = zoom/2;
				
			}
			else if (s == 1){
				if (zoom <1) zoom = zoom *2;
				else zoom++;
			}
			int x = (int)Math.round((zoom * (double)lengthX));
			int y = (int)Math.round((zoom * (double)lengthY));
			setPreferredSize(new Dimension(x, y));
			revalidate();
			repaint();
			
		}
		
		public void modifySelection(int x, int y){
			if (selectionVisible){
				double w = selectionEllipse.getWidth();
				double h = selectionEllipse.getHeight();
				selectionEllipse.setFrame(selectionEllipse.getX(), selectionEllipse.getY(), w+x, h+y);
				repaint();
			}
		}
		public void nudgeSelection(int x, int y){
			if (selectionVisible){
				double xC = selectionEllipse.getX();
				double yC = selectionEllipse.getY();
				selectionEllipse.setFrame(xC+x, yC+y, selectionEllipse.getWidth(), selectionEllipse.getHeight());
				repaint();
			}
		}
	}
	
	
	public void makeBufferedImage(){
		System.out.print("Making Buffered Image: ");
		Graphics2D g2=null;
		bufferedImage = new BufferedImage (lengthX, lengthY, BufferedImage.TYPE_INT_ARGB);
		g2 = bufferedImage.createGraphics();
		g2.addRenderingHints(GATAUtil.fetchRenderingHints());
		//kill antialiasing to make sharp lines and text for computer display, not good for printed figures
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF);
		g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF);
		
		//draw features
		int counter = 0;
		Ellipse2D elipse;
		int colorIndex;
		for (int x=0; x< lengthX; x++){
			counter++;
			if (counter == 100){
				counter = 0;
				System.out.print((lengthX-x)+" ");
			}
			for (int y=0; y< lengthX; y++){
				//set params
				colorIndex = (int)Math.round(colorScalar * intensities[x][y]);
				//if zero index, black, don't draw
				if (colorIndex == 0) continue;
				if (colorIndex >= 768) colorIndex = 767;
				elipse = new Ellipse2D.Float(x,y,1,1);
				//draw it
				g2.setColor(colors[colorIndex]);
				g2.draw(elipse);
			}
		}
		if (writePNG){
			System.out.print("Writing ");
			try{
				ImageIO.write(bufferedImage, "PNG", pngFile);
			} catch (IOException e){
				e.printStackTrace();
			}
		}
		System.out.println("Done!");
	}
	
	/**Generates 768 diff colors along a black red yellow white gradient.*/
	public static Color[] makeColors(){
		Color[] colors = new Color[256 + 256 + 256];
		for (int i=0; i<256; i++){
			colors[i] = new Color(i,0,0);
		}
		int counter = 0;
		for (int i=256; i<512; i++){
			colors[i] = new Color(255,counter++,0);
		}
		counter = 0;
		for (int i=512; i<768; i++){
			colors[i] = new Color(255,255,counter++);
		}
		return colors; 
	}
	
	public void paintComponent(Graphics g){
		super.paintComponent(g);
		Graphics2D g2 = (Graphics2D)g;
		g2.scale(zoom, zoom);
		//draw image
		g2.drawImage(bufferedImage,null,0,0);
		//draw any masked Ellipses
		if (maskedEllipses != null){
			g2.setColor(Color.BLUE);
			for (int i=0; i< maskedEllipses.length; i++){
				g2.fill(maskedEllipses[i]);
			}
		}
		//draw selection rectangle
		if (selectionVisible){
			g2.setColor(Color.CYAN);
			g2.draw(selectionEllipse);			
		}
	}
}