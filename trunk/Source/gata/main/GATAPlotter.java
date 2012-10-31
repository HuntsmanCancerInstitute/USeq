package gata.main;

import gata.aligner.*;
import gata.geneGlyphs.*;
import gata.menu.*;
import gata.plotter.*;
import gata.plotter.Console;

import java.awt.*;
import javax.swing.*;
import java.io.*;

/**
 * @author Nix
 * Launcher class for GATAPlotter.  Call main.
 * */

public class GATAPlotter {
	//fields
	private GATAParams gataParams;
	private GATAFrame gataFrame;
	private AnnoSpecParams refParams;
	private AnnoSpecParams compParams;
	private JScrollPane compScrollPane;
	private JScrollPane alignScrollPane;
	private JSplitPane innerSplitPane;
	//private ToolBar toolBar;
	
	public static void main(String[] args) {
		new GATAPlotter();
	}
	public GATAPlotter(){
		//make new GATAParams object
		gataParams = new GATAParams();
		gataParams.setGataPlotter(this);
		//fetch plotter preferences or make new one and fill GATAParams
		PlotterPreferences pp = (PlotterPreferences)GATAUtil.fetchObject(new File("GATAPlotterPreferences"));
		if (pp==null) {pp = new PlotterPreferences();
			pp.setPlotterDefaults(gataParams);
		}
		else pp.setCurrentPlotterPreferences(gataParams);
		gataParams.setPlotterPreferences(pp);
		//set objectsFile
		File objFile = getObjectsFile();
		gataParams.setObjectsFile(objFile);

		//make Frame
		gataFrame = new GATAFrame(gataParams);
		gataParams.setGataFrame(gataFrame);
		gataFrame.setResizable(true);
		gataFrame.setTitle("GATA Plotter");
		gataFrame.setSize(600, 100);
		gataFrame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		
		
		gataFrame.show();
		//fire request for new files
		new FileInputForm(gataParams);

	}
	public File getObjectsFile(){
		AlignerPreferences ap;
		try {
			ObjectInputStream in = new ObjectInputStream( new FileInputStream("GATAlignerPreferences"));
			ap = (AlignerPreferences) in.readObject();
			in.close();
			return new File(ap.getStorageLoc(),ap.getBaseName());
		} catch (Exception e) {
			System.out.println("Cannot find/open GATAlignerPreferences. Using defaults");
			return new File("");
		}
				
	}
	public void launchGATAPlotter(){			
		//alignment panel
		AlignPanel alignPanel = new AlignPanel(gataParams);
		gataParams.setAlignPanel(alignPanel);
		gataParams.getGataFrame().setAlignPanel(alignPanel);
		// save ref to panel in params for sliders
		alignPanel.setBackground(Color.WHITE);
		alignPanel.setLayout(new FlowLayout());
		alignPanel.setPreferredSize(
			new Dimension(gataParams.getWidth(), gataParams.getHeight()+400)); // extra 400 is a buffer for making bigger panel
		
		//make annotation
		boolean gffRefPres = gataParams.isGffRefPresent();
		boolean gffCompPres = gataParams.isGffCompPresent();
		
		GlyphPanel refAnnoPanel =null;
		GlyphPanel compAnnoPanel =null;
					
		if (gffRefPres|| gffCompPres) {
			if (gffRefPres){
				//make ref annotation
				refParams = gataParams.getRefAnnoSpecParams();				
				makeAnnotation(refParams);
				gataFrame.setRefAnnotation(refParams.getAnnotation());
				refAnnoPanel = refParams.getAnnotation().getGlyphPanel();
				refParams.setGlyphPanel(refAnnoPanel);
				gataFrame.setRefGlyphPanel(refAnnoPanel);
			}
			if (gffCompPres){
				//make comp annotation
				compParams = gataParams.getCompAnnoSpecParams();
				makeAnnotation(compParams);
				gataFrame.setCompAnnotation(compParams.getAnnotation());
				compAnnoPanel = compParams.getAnnotation().getGlyphPanel();
				compParams.setGlyphPanel(compAnnoPanel);
				gataFrame.setCompGlyphPanel(compAnnoPanel);
			}			

			//adjust annotation menues for possible generic tracks
			gataFrame.completeAnnotationMenu(gffRefPres, gffCompPres, refParams, compParams);

			//make scrollPanes for each panel
			JScrollPane refScrollPane=null;
			if(gffRefPres) {
				refScrollPane = new JScrollPane(refAnnoPanel);
				refScrollPane.setPreferredSize(new Dimension(gataParams.getWidth(), refParams.getAnnoPanelHeight()+20)); //for buffer
				JViewport viewPort = refScrollPane.getViewport();
				refAnnoPanel.setViewPort(viewPort);
				viewPort.setViewPosition(new Point (0,200)); //provides buffer				
				GATAUtil.setScrollPaneScrollBarIncrement(refScrollPane, 10);
			} 
			alignScrollPane = new JScrollPane(alignPanel);
				alignScrollPane.setPreferredSize(new Dimension(gataParams.getWidth(), gataParams.getHeight()+20)); //provides buffer			
				JViewport port = alignScrollPane.getViewport();
				alignPanel.setViewPort(port); //needed for center zooming
				port.setViewPosition(new Point ((int)port.getViewPosition().getX(),200)); //provides buffer
				GATAUtil.setScrollPaneScrollBarIncrement(alignScrollPane, 10);
			compScrollPane=null;
			if(gffCompPres){
				compScrollPane = new JScrollPane(compAnnoPanel); 
				compScrollPane.setPreferredSize(new Dimension(gataParams.getWidth(), compParams.getAnnoPanelHeight()+20)); //provides buffer
				JViewport viewPort = compScrollPane.getViewport();
				compAnnoPanel.setViewPort(viewPort);
				viewPort.setViewPosition(new Point (0,200)); //provides buffer				
				GATAUtil.setScrollPaneScrollBarIncrement(compScrollPane, 10);
			}
		
			//set scrollbars and splitPanes
			JSplitPane outerSplitPane =null;
			
			//all three present
			if (gffRefPres && gffCompPres){
				alignScrollPane.setHorizontalScrollBar(refScrollPane.getHorizontalScrollBar());
				compScrollPane.setHorizontalScrollBar(alignScrollPane.getHorizontalScrollBar());
				compScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
				innerSplitPane= new JSplitPane(JSplitPane.VERTICAL_SPLIT,alignScrollPane,compScrollPane);
				innerSplitPane.setDividerSize(1);
				outerSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, refScrollPane,innerSplitPane);
			}
			//only two present
			else {
				if (gffCompPres==false) {
					alignScrollPane.setHorizontalScrollBar(refScrollPane.getHorizontalScrollBar());
					alignScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
					outerSplitPane =new JSplitPane(JSplitPane.VERTICAL_SPLIT, refScrollPane, alignScrollPane);
				} 
				if (gffRefPres==false) {
					compScrollPane.setHorizontalScrollBar(alignScrollPane.getHorizontalScrollBar());
					compScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
					outerSplitPane =new JSplitPane(JSplitPane.VERTICAL_SPLIT, alignScrollPane, compScrollPane);
					innerSplitPane = outerSplitPane;  //to enable Divider management upon hiding/showing
				} 
			}
			
			//add scrollPanes to a splitPane
			outerSplitPane.setDividerSize(1);
	
			//make and load frame with splitPane	
			Container contentPane = gataFrame.getContentPane();
			contentPane.add(outerSplitPane, BorderLayout.CENTER);
		}
		//if no annotation then default to regular align frame, no splitframe
		else {
			Container contentPane = gataFrame.getContentPane();
			alignScrollPane = new JScrollPane(alignPanel);
			alignScrollPane.setPreferredSize(new Dimension(gataParams.getWidth(), gataParams.getHeight()+20)); //provides buffer
			GATAUtil.setScrollPaneScrollBarIncrement(alignScrollPane, 10);
			JViewport port = alignScrollPane.getViewport();
			alignPanel.setViewPort(port); //needed for center zooming
				port.setViewPosition(new Point ((int)port.getViewPosition().getX(),200)); //provides buffer
				alignScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		
			contentPane.add(alignScrollPane, BorderLayout.CENTER);
		}
		gataFrame.pack();
		Rectangle gfBounds = gataFrame.getBounds();

		//fire new toolframe
		ToolsFrame tf = new ToolsFrame((int)gfBounds.getX(), (int)(gfBounds.getY()+gfBounds.getHeight()), gataParams);
		tf.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
		alignPanel.setToolsRef(tf);
		Rectangle tfBounds = tf.getBounds();
		int x = (int)(tfBounds.getX()+tfBounds.getWidth());
		int y = (int)(tfBounds.getY());
		int w = (int)(gfBounds.getWidth()-tfBounds.getWidth());
		int h = (int)(tfBounds.getHeight());
		
		//make new text Console  
		Console console = new Console(x, y, w, h,  "Text Console");
		console.printToTextArea("Click on a box or line to retrieve the alignment.\n"+
			  "Double click to see the alignment in the context of the larger Blast alignment.\n\n");
		gataParams.setConsole(console);
		gataParams.getAlignPanel().setConsole(console);
		if (gffRefPres) refAnnoPanel.setConsole(console);
		if (gffCompPres) compAnnoPanel.setConsole(console);
		
		//set frame title
		AlignParams alignParams = gataParams.getAlignParams();
		String refSeqName = alignParams.getNAME_REFSEQ();
		String compSeqName = alignParams.getNAME_COMPSEQ();
		if (refSeqName.length()>30) refSeqName = refSeqName.substring(0,30)+"...";
		if (compSeqName.length()>30) compSeqName = compSeqName.substring(0,30)+"...";
		gataFrame.setTitle("GATAPlotter: "+refSeqName+" vs. "+compSeqName);
			  
		//show frames
		gataFrame.show();
		tf.show();
		console.show();
		
	}
	public void makeAnnotation(AnnoSpecParams annoParams){			
		//make annotation panel
		Annotation anno = new Annotation(annoParams, gataParams);
		//recalculate height and resize panel size
		annoParams.clearLines();
		int smallestY = annoParams.getSmallestGlyphPanelY();
		int biggestY= annoParams.getBiggestGlyphPanelY();
		annoParams.setAnnoPanelHeight(biggestY-smallestY);
		annoParams.setReferenceY(200+ annoParams.getReferenceY()-smallestY);	//100pix buffer
		annoParams.setScaleRulerY(200+ annoParams.getScaleRulerY()-smallestY);	//100pix buffer
		anno.redrawAnnotation();
		annoParams.setAnnotation(anno);
		
	}		
	public JScrollPane getCompScrollPane() {
		return compScrollPane;
	}
	public JScrollPane getAlignScrollPane() {
		return alignScrollPane;
	}
	public JSplitPane getInnerSplitPane() {
		return innerSplitPane;
	}
}