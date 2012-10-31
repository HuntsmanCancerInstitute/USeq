package gata.menu;

import gata.aligner.*;
import gata.geneGlyphs.*;
import gata.main.*;
import gata.plotter.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.imageio.*;
import java.awt.image.*;
import java.awt.geom.*;

/**
 * @author nix
 * Panel and methods to capture and save a GATAPlot and a high resolution PNG file.
 */
public class SaveImagePanel extends JPanel implements ActionListener {
	//fields
	private int x = 10;
	private int y = 10;
	private JTextField saveFolderField;
	private JTextField nameField;
	private JTextField widthField;
	private JFrame frame;
	private JFileChooser chooser;
	private int maxWidth;
	private int maxHeight;
	private GATAParams params;
	private AlignPanel alignPanel;
	private boolean saveVisible =true;
	
	public SaveImagePanel(GATAParams params, JFrame frame) {
		//modify panel 
		setBackground(Color.WHITE);
		setLayout(null);
		chooser = new JFileChooser();
		chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		this.params = params;
		this.frame = frame;
		alignPanel = params.getAlignPanel();
		
		//get Aligner saved folder
		File objectsFile = params.getObjectsFile();
		String fileName = "";
		try{
			fileName = objectsFile.getCanonicalPath();
		}
		catch(IOException e){
			e.printStackTrace();
		}
		chooser.setCurrentDirectory(objectsFile);
		//save folder
		JLabel saveFolderLabel = GATAUtil.makeLabel(this,"Select a folder for saving the GATAPlot image.",14,0);
		saveFolderLabel.setForeground(Color.BLUE);
		int[] saveFolderLabelSize = GATAUtil.setLabel(saveFolderLabel, x, y);
		y += saveFolderLabelSize[1]+5;
		
		saveFolderField = GATAUtil.makeFieldSimple(this,objectsFile.getParent());
		int[] alignerFileSize = GATAUtil.setField(saveFolderField, x, y, 400);
		x+= alignerFileSize[0]+5;
		
		JButton saveFolderBrows = GATAUtil.makeButton("Browse", 12, new BrowseButtonAction(saveFolderField,frame,chooser),this);
		int[] saveFolderBrowsSize = GATAUtil.setButton(saveFolderBrows, x,y);
		maxWidth = saveFolderBrowsSize[0]+10+x;
		y+=saveFolderBrowsSize[1]+10;
		x= 10;
		
		//text
		JLabel nameLabel = GATAUtil.makeLabel(this,"Enter a text for the GATAPlot image.",14,0);
		nameLabel.setForeground(Color.BLUE);
		int[] nameLabelSize = GATAUtil.setLabel(nameLabel, x, y);
		x += nameLabelSize[0]+5;
		
		nameField = GATAUtil.makeFieldSimple(this, objectsFile.getName()+".png");
		int[] nameFieldSize = GATAUtil.setField(nameField, x, y, 200);
		y += nameFieldSize[1]+10;
		x=10;
		
		//Crop image question
		int I =5 ;//spacer
		JLabel cropLabel = GATAUtil.makeLabel(this,"Crop to visible?", 14, 0);
		cropLabel.setForeground(Color.BLUE);
		int[] cropLabelSize = GATAUtil.setLabel(cropLabel, x, y);
		x += cropLabelSize[0] + I;
		ButtonGroup cropGroup = new ButtonGroup();
		JRadioButton yes = new JRadioButton("Yes", false);
		yes.addActionListener(this);
		yes.setBackground(Color.WHITE);
		cropGroup.add(yes);
		yes.setSelected(true);
		JRadioButton no = new JRadioButton("No", true);
		no.addActionListener(this);
		no.setBackground(Color.WHITE);
		cropGroup.add(no);
		Dimension dimYes = yes.getPreferredSize();
		Dimension dimNo = no.getPreferredSize();
		yes.setBounds(
			x,
			y,
			(int) dimYes.getWidth(),
			(int) dimYes.getHeight());
		add(yes);
		x += dimYes.getWidth() + I;
		no.setBounds(
			x,
			y,
			(int) dimNo.getWidth(),
			(int) dimNo.getHeight());
		add(no);
		y += dimYes.getHeight() + 10;
		x=10;
		//pixel width
		JLabel widthLabel = GATAUtil.makeLabel(this,"Enter a size in pixels for the GATAPlot image.",14,0);
		widthLabel.setForeground(Color.BLUE);
		int[] widthLabelSize = GATAUtil.setLabel(widthLabel, x, y);
		x += widthLabelSize[0]+5;
		
		widthField = GATAUtil.makeFieldSimple(this, alignPanel.getViewPort().getWidth() +"");
		int[] widthFieldSize = GATAUtil.setField(widthField, x, y, 100);
		y += widthFieldSize[1]+10;
		x=10;

		
		//seperator
		JSeparator sep = new JSeparator();
		sep.setBounds(0,y,maxWidth,5);
		add(sep);
		y+=10;
		
		//Go
		JButton go = GATAUtil.makeButton(" Save GATAPlot ", 12, new GoButtonAction(),this);
		Dimension goDim = go.getPreferredSize();
		JButton cancel = GATAUtil.makeButton(" Cancel ", 12, new CancelButtonAction(),this);
		Dimension cancelDim = cancel.getPreferredSize();
		int wid = (int)(goDim.getWidth()+5+cancelDim.getWidth());
		x = maxWidth/2 - wid/2;
		int[] goSize = GATAUtil.setButton(go,x,y);
		x+= goSize[0]+5;
		GATAUtil.setButton(cancel,x,y);
		maxHeight = goSize[1]+y+35;

	}
	//used to monitor yes no buttons on cropping
	public void actionPerformed(ActionEvent e){
		String label = ((JRadioButton)e.getSource()).getText();
		//set pixel width in widthField
		if (label.equals("Yes")){
			widthField.setText(alignPanel.getViewPort().getWidth() +"");
			saveVisible = true;
		}
		else {
			widthField.setText(alignPanel.getWidth()+"");
			saveVisible = false;			
		} 
	}
	
	public int getMaxHeight() {
		return maxHeight;
	}
	public int getMaxWidth() {
		return maxWidth;
	}
	private class GoButtonAction implements ActionListener {
		
		public void actionPerformed(ActionEvent e) {
			//look for directory
			if (GATAUtil.checkFile(saveFolderField.getText())==false) return;
			
			//check that there is a text
			if (GATAUtil.isStringEmpty(nameField.getText())) {
				GATAUtil.throwWarning(frame, nameField,"Please enter a text for the image!");
				return;
			}
			//check that there is a size
			int size;
			try{
				size = Integer.parseInt(widthField.getText());
			}catch (NumberFormatException ex){
				GATAUtil.throwWarning(frame, widthField,"Please enter an integer width in pixels.");
				return;				
			}
			//calculate scalar
			double scalar;
			if(saveVisible) scalar = (double)size/(double)alignPanel.getViewPort().getWidth();
			else scalar = (double)size/(double)alignPanel.getWidth();
			
			File file = new File(saveFolderField.getText(), nameField.getText());
			//make png
			try{
				makePNG(file, scalar);
				frame.dispose();
			} catch (OutOfMemoryError o){
				GATAUtil.throwWarning(frame, widthField,
				"Sorry, your requested image is too big, reduce the size.\n  Alternatively, start GATAPlotter on the command line \n  with the following 'java -Xmx256m -jar GATAPlotter'");
			}	
		}
	}
	private class CancelButtonAction implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			frame.dispose();
		}
	}
	public void makePNG(File fileName, double scaleMultiplier){
		//booleans
		boolean gffRefPres = params.isGffRefPresent();
		boolean gffCompPres = params.isGffCompPresent();
		GlyphPanel refPanel =null;
		GlyphPanel compPanel =null;
		JViewport viewPort = alignPanel.getViewPort();
		if (gffRefPres) refPanel=params.getRefAnnoSpecParams().getGlyphPanel();
		if (gffCompPres) compPanel=params.getCompAnnoSpecParams().getGlyphPanel();
		
		//calculate size of combo
		int alignW;
		int alignH;
		int comboW;
		int comboH;
		int refH=0;
		int compH=0;
		if (saveVisible){ //crop panel
			alignW = (int)Math.round(((double)viewPort.getWidth())*scaleMultiplier);
			alignH = (int)Math.round(((double)viewPort.getHeight())*scaleMultiplier);
			comboW = alignW+2; //for 1pix blue border
			comboH = alignH+2;
			if (gffRefPres){ 
				refH = (int)Math.round(((double)refPanel.getViewPort().getHeight())*scaleMultiplier);
				comboH +=refH+1;
			}
			if (gffCompPres){
				compH = (int)Math.round(((double)compPanel.getViewPort().getHeight())*scaleMultiplier);
				comboH +=compH+1;
			}
		}
		else { //full size panel
			alignW = (int)Math.round(alignPanel.getAlignFrameWidth()*scaleMultiplier);
			alignH = (int)Math.round(alignPanel.getAlignFrameHeight()*scaleMultiplier);	
			comboW = alignW+2; //for 1pix blue border
			comboH = alignH+2;
			if (gffRefPres){ 
				refH = (int)Math.round(((double)refPanel.getViewPort().getHeight())*scaleMultiplier);
				comboH +=refH+1;
			}
			if (gffCompPres){
				compH = (int)Math.round(((double)compPanel.getViewPort().getHeight())*scaleMultiplier);
				comboH +=compH+1;
			}
		}

		//combine panels
		BufferedImage combo = new BufferedImage(comboW, comboH, BufferedImage.TYPE_INT_ARGB);
		Graphics2D g2 = combo.createGraphics();
		g2.setColor(Color.BLUE);
		g2.fill(new Rectangle2D.Double(0,0,comboW, comboH)); //to create solid background
		
		//all this null business is to try and keep the heap memory low and prevent out of memory errors
		BufferedImage buff;

		if (gffRefPres) {
			buff=refPanel.makeBufferedImage(scaleMultiplier,saveVisible);
			g2.drawImage(buff,1,1, (ImageObserver)new Canvas());
			buff=alignPanel.makeBufferedImage(scaleMultiplier,saveVisible); 
			g2.drawImage (buff, 1, refH+2, (ImageObserver)new Canvas());
			buff=null;
			if (gffCompPres) {
				buff = compPanel.makeBufferedImage(scaleMultiplier,saveVisible);
				g2.drawImage(buff,1,refH+alignH+3,(ImageObserver)new Canvas()); 
				buff=null;
			}
		}
		else if (gffCompPres){
			buff=alignPanel.makeBufferedImage(scaleMultiplier,saveVisible);
			g2.drawImage (buff, 1, 1, (ImageObserver)new Canvas());
			buff = compPanel.makeBufferedImage(scaleMultiplier,saveVisible);
			g2.drawImage(buff,1,alignH+2,(ImageObserver)new Canvas());
			buff=null;
		}
		else {
			buff=alignPanel.makeBufferedImage(scaleMultiplier,saveVisible);
			g2.drawImage (buff, 1, 1, (ImageObserver)new Canvas()); 
			buff=null;
		}
		try{
			ImageIO.write(combo, "PNG", fileName);
			combo=null;
		} catch (IOException e){
			e.printStackTrace();
		}
	}

}
