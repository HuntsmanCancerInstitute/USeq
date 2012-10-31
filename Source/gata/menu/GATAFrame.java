package gata.menu;

import gata.aligner.*;
import gata.geneGlyphs.*;
import gata.main.*;
import gata.plotter.*;
import gata.plotter.Console;

import java.awt.event.*;
import javax.swing.*;
import java.awt.*;

import java.io.*;

/**
 * @author Nix
 * Frame holding all of the menues for running the GATAPlotter App. Is also a container for the
 * aligner panel and any glyph/annotation panels.  Sort of ugly!
 * */
public class GATAFrame extends JFrame {
	//fields
	private GATAParams gataParams;
	private Annotation refAnnotation;
	private Annotation compAnnotation;
	private AnnoSpecParams refParams;
	private AnnoSpecParams compParams;
	private boolean gffRefPres = false;
	private boolean gffCompPres = false;
	private boolean refGenericGlyphs = false;
	private boolean compGenericGlyphs = false;
	private GlyphPanel refGlyphPanel;
	private GlyphPanel compGlyphPanel;	
	private JMenu mAnnotation;
	private ConservedSeqs conservedSeqs; //might be null
	private AlignPanel alignPanel;
	private PlotterPreferences plotterPreferences;
	private GATAFrame gataFrame;
	private AnnoSpecParams genParams = new AnnoSpecParams(); //used as a dummy to pass along in conditional statements, set to ref or compParams
	private TextWindow documentWindow;  //for html help docs
	private int innerDividerLoc = 0; //for keeping track of where the divider for the align/comp is
	
	/**needed to add later since annotation doesn't exist when ToolFrame is created*/
	public void setRefAnnotation(Annotation a){
		refAnnotation = a;
	}
	/**needed to add later since annotation doesn't exist when ToolFrame is created*/
	public void setCompAnnotation(Annotation a){
		compAnnotation = a;
	}	
	/**must be called after the initial creation when the user has possibly loaded generic annotation*/
	public void completeAnnotationMenu(boolean gffRefPres, boolean gffCompPres, 
		AnnoSpecParams refParams, AnnoSpecParams compParams){
		//set params for frame
		this.gffRefPres = gffRefPres;
		this.gffCompPres = gffCompPres;
		this.refParams = refParams;
		this.compParams = compParams;
		
		//build annotation menue
		addShowHideLabel(mAnnotation, " Gene Groups ");
		addShowHide(mAnnotation, " Protein ",true);
		addShowHide(mAnnotation, " RNA ",true);
		addShowHide(mAnnotation, " DNA ",true);
		
		if (gffRefPres) addScaleRuler(mAnnotation, " Scale Ruler Ref");
		if (gffCompPres) addScaleRuler(mAnnotation, " Scale Ruler Comp");
		if (gffRefPres || gffCompPres) addTrack(mAnnotation, " Tracks ",false);
		//add tracks
		if (gffRefPres){
			GenericGlyph[][] gg = refParams.getTrackedGenericGlyphs();
			if (gg != null){
				int len = gg.length;
				for (int i=0; i<len; i++){
					addTrack(mAnnotation, " RefTrk "+(i+1)+": "+gg[i][0].getGenericFeature().getFeatureType(), true);
				}
			}
		}
		if (gffCompPres){
			GenericGlyph[][] gg = compParams.getTrackedGenericGlyphs();
			if (gg != null){
				int len = gg.length;
				for (int i=0; i<len; i++){
					addTrack(mAnnotation, " CompTrk "+(i+1)+": "+gg[i][0].getGenericFeature().getFeatureType(),true);
				}
			}
		}			
	
		mAnnotation.add(new JSeparator());
		mAnnotation.add(new MenuAction(" Annotation ", " Pix Btw Labels & Trans Groups "));
		mAnnotation.add(new MenuAction(" Annotation ", " Pix Btw Trans Groups "));
		mAnnotation.add(new MenuAction(" Annotation ", " Pix Btw Trans Groups & DNA "));
		mAnnotation.add(new MenuAction(" Annotation ", " Pix Btw Protein & RNA "));
		//check and assign booleans for generic glyphs
		if (gffRefPres) refGenericGlyphs = refParams.getAnyGenericGlyphs();
		if (gffCompPres) compGenericGlyphs = compParams.getAnyGenericGlyphs();
		if (refGenericGlyphs|| compGenericGlyphs){
			mAnnotation.add(new MenuAction(" Annotation ", " Pix Btw DNA & Tracks "));
			mAnnotation.add(new MenuAction(" Annotation ", " Pix Btw Labels & Tracks "));
			mAnnotation.add(new MenuAction(" Annotation ", " Pix Btw Tracks "));
		}
		mAnnotation.add(new JSeparator());
		mAnnotation.add(new MenuAction(" Annotation ", " Thin Line Thickness "));
		mAnnotation.add(new MenuAction(" Annotation ", " Normal Line Thickness "));
		mAnnotation.add(new MenuAction(" Annotation ", " Fat Line Thickness "));
		mAnnotation.add(new MenuAction(" Annotation ", " Arrow Thickness "));
		mAnnotation.add(new JSeparator());
		mAnnotation.add(new MenuAction(" Annotation ", " Background Color "));
		mAnnotation.add(new MenuAction(" Annotation ", " Label Color "));
		if (refGenericGlyphs|| compGenericGlyphs) addScaleTrackColors(mAnnotation);
		
	}
	
	//constructor
	public GATAFrame(GATAParams gataParams) {
		this.gataParams = gataParams;
		plotterPreferences = gataParams.getPlotterPreferences();
		gataFrame = gataParams.getGataFrame();
		
		setSize(600, 100);
		JMenuBar bar = new JMenuBar();
		setJMenuBar(bar);
		//File
		JMenu mFile = new JMenu("    File  ");
		mFile.add(new MenuAction(" File ", " Open New Alignment "));
		mFile.add(new MenuAction(" File ", " Close Alignment "));
		mFile.add(new MenuAction(" File ", " Save GATAPlot Image "));
		mFile.add(new MenuAction(" File ", " Save GATAPlotter Settings "));
		mFile.add(new JSeparator());
		mFile.add(new MenuAction(" File ", " Quit GATAPlotter "));
		bar.add(mFile);

		//Alignment
		JMenu mAlignment = new JMenu("  Alignment  ");
		mAlignment.add(new MenuAction(" Alignment ", " Window Width "));
		mAlignment.add(new MenuAction(" Alignment ", " Window Height "));
		mAlignment.add(new MenuAction(" Alignment ", " Dist Btw Top & Ref Seq "));
		mAlignment.add(new MenuAction(" Alignment ", " Dist Btw Top & Comp Seq "));
		mAlignment.add(new MenuAction(" Alignment ", " Thickness of Seq Bars "));
		mAlignment.add(new MenuAction(" Alignment ", " Dist Btw Left & Seq Bar "));
		mAlignment.add(new MenuAction(" Alignment ", " Width of Largest Seq Bar "));
	    mAlignment.add(new MenuAction(" Alignment ", " Set Nucleotides Per Pixel "));
		
		mAlignment.add(new JSeparator());
		mAlignment.add(new MenuAction(" Alignment "," Fetch Conserved Sequences "));
		mAlignment.add(new MenuAction(" Alignment ", " GATAligner Parameters "));
		
		//for each aligment .. finish
		mAlignment.add(new JSeparator());
		addShowHide(mAlignment, " Comp Seq Vis ",false);
		bar.add(mAlignment);

		//Annotation
		mAnnotation = new JMenu("  Annotation  ");
		bar.add(mAnnotation);

		//Windows
		JMenu mWindows = new JMenu("  Windows  ");

		mWindows.add(new MenuAction(" Windows ", " Show All "));
		mWindows.add(new MenuAction(" Windows ", " Hide All "));
		mWindows.add(new MenuAction(" Windows ", " Redraw Using Defaults "));
		bar.add(mWindows);
		
		//Docs
		bar.add(Box.createHorizontalGlue());
		JMenu mDocs = new JMenu(" Documentation      ");
		mDocs.add(new MenuAction(" Documentation ", " About "));
		mDocs.add(new MenuAction(" Documentation ", " Plotter Help "));
		mDocs.add(new MenuAction(" Documentation ", " Aligner Help"));
		mDocs.add(new MenuAction(" Documentation ", " Examples "));
		bar.add(mDocs);
		
		//add window listener to close or open other windows
		addWindowListener(new WindowStateMonitor());
				
	}
	public JMenuItem addShowHide(JMenu menu, String label,boolean color) {
		String menuText = menu.getText().trim() +":"+label.trim();
		JRadioButtonMenuItem show = new JRadioButtonMenuItem(" Show");
		show.addActionListener(new RadioAction(menuText, "Show"));
		show.setSelected(true);
		JRadioButtonMenuItem hide = new JRadioButtonMenuItem(" Hide");
		hide.addActionListener(new RadioAction(menuText, "Hide"));
		ButtonGroup groupSH = new ButtonGroup();
		groupSH.add(show);
		groupSH.add(hide);
		JMenu showHide = new JMenu(label);
		showHide.add(show);
		showHide.add(hide);
		if(color) showHide.add(new MenuAction(menuText, "Color"));
		return menu.add(showHide);
	}
	public JMenuItem addScaleRuler(JMenu menu, String label) {
		String menuText = menu.getText().trim() +":"+label.trim();
		JRadioButtonMenuItem show = new JRadioButtonMenuItem(" Show");
		show.addActionListener(new RadioAction(menuText, "Show"));
		show.setSelected(true);
		JRadioButtonMenuItem hide = new JRadioButtonMenuItem(" Hide");
		hide.addActionListener(new RadioAction(menuText, "Hide"));
		ButtonGroup groupSH = new ButtonGroup();
		groupSH.add(show);
		groupSH.add(hide);
		JMenu showHide = new JMenu(label);
		showHide.add(show);
		showHide.add(hide);
		showHide.add(new MenuAction(menuText, "Color"));
		showHide.add(new MenuAction(menuText, "Location"));
		return menu.add(showHide);
	}
	
	public JMenuItem addScaleTrackColors(JMenu menu) {
		String label = " Scale Track Colors By Score ";
		String menuText = menu.getText().trim() +":"+label.trim();
		JRadioButtonMenuItem linear = new JRadioButtonMenuItem(" Linear");
		linear.addActionListener(new RadioAction(menuText, "Linear"));

		JRadioButtonMenuItem log10 = new JRadioButtonMenuItem(" Log(10)");
		log10.addActionListener(new RadioAction(menuText, "Log(10)"));
		
		JRadioButtonMenuItem log2 = new JRadioButtonMenuItem(" Log(2)");
		log2.addActionListener(new RadioAction(menuText, "Log(2)"));
		
		JRadioButtonMenuItem lnx = new JRadioButtonMenuItem(" Lnx");
		lnx.addActionListener(new RadioAction(menuText, "Lnx"));
				
		JRadioButtonMenuItem none = new JRadioButtonMenuItem(" None");
		none.addActionListener(new RadioAction(menuText, "None"));
		
		
		ButtonGroup groupSH = new ButtonGroup();
		groupSH.add(linear);
		groupSH.add(log10);
		groupSH.add(log2);
		groupSH.add(lnx);
		groupSH.add(none);
		
		JMenu showHide = new JMenu(label);
		showHide.add(linear);
		showHide.add(log10);
		showHide.add(log2);
		showHide.add(lnx);
		showHide.add(none);
		return menu.add(showHide);
	}	
	
	public JMenuItem addShowHideWithPath(JMenu menu, String label,String path) {
			String menuText = path+":"+menu.getText().trim() +":"+label.trim();
			JRadioButtonMenuItem show = new JRadioButtonMenuItem(" Show");
			show.addActionListener(new RadioAction(menuText, "Show"));
			show.setSelected(true);
			JRadioButtonMenuItem hide = new JRadioButtonMenuItem(" Hide");
			hide.addActionListener(new RadioAction(menuText, "Hide"));
			ButtonGroup groupSH = new ButtonGroup();
			groupSH.add(show);
			groupSH.add(hide);
			JMenu showHide = new JMenu(label);
			showHide.add(show);
			showHide.add(hide);
			return menu.add(showHide);
	}
	public JMenuItem addShowHideLabel(JMenu menu, String label) {
			String menuText = menu.getText().trim() +":"+label.trim();
			JRadioButtonMenuItem show = new JRadioButtonMenuItem(" Show");
			show.addActionListener(new RadioAction(menuText, "Show"));
			show.setSelected(true);
			JRadioButtonMenuItem hide = new JRadioButtonMenuItem(" Hide");
			hide.addActionListener(new RadioAction(menuText, "Hide"));
			ButtonGroup groupSH = new ButtonGroup();
			groupSH.add(show);
			groupSH.add(hide);
			JMenu track = new JMenu(label);
			track.add(show);
			track.add(hide);
			track.add(addShowHideWithPath(track,"Label",menu.getText().trim()));
			return menu.add(track);
		}	

	public JMenuItem addTrack(JMenu menu, String label, boolean color) {
		String annoLabel = "Annotation:"+label.trim();
		JRadioButtonMenuItem show = new JRadioButtonMenuItem(" Show");
		show.addActionListener(new RadioAction(annoLabel, "Show"));
		show.setSelected(true);
		JRadioButtonMenuItem hide = new JRadioButtonMenuItem(" Hide");
		hide.addActionListener(new RadioAction(annoLabel, "Hide"));
		ButtonGroup groupSH = new ButtonGroup();
		groupSH.add(show);
		groupSH.add(hide);
		JMenu track = new JMenu(label);
		track.add(show);
		track.add(hide);
		if (color) track.add(new MenuAction(menu.getText().trim()+":"+label.trim(), "Color"));
		track.add(new MenuAction(menu.getText().trim()+":"+label.trim(), "Thickness"));
		track.add(addShowHideWithPath(track,"Label",menu.getText().trim()));
		
		return menu.add(track);
	}	
	
	/**Used to listen for show and hide menu selections for various annotation glyphs*/
	private class RadioAction implements ActionListener {
		String path;
		String name;
		
		public RadioAction(String path, String name) {
			this.path = path.trim();
			this.name = name.trim();
		}
		public void actionPerformed(ActionEvent e) {
			//System.out.println("Radio_ path: "+path+ "\nname: "+ text +"\n\t(path is split on: into path[0]...)");
			String[] paths = path.split(":");
			if (paths[0].equals("Annotation")){
				if(paths[1].equals("Gene Groups")){
					if(paths.length==3){
						if (name.equals("Show")){
							if (gataParams.isLabelsVis()==false){
								gataParams.setLabelsVis(true);
								redrawAnnotationPanels();
							}
						}
						else { //by default hide
							if (gataParams.isLabelsVis()){
								gataParams.setLabelsVis(false);
								redrawAnnotationPanels();
							}
						}
					}
					else {
						if (name.equals("Show")){
							gataParams.setProteinVis(true);
							gataParams.setRNAVis(true);
							gataParams.setDNAVis(true);
							gataParams.setLabelsVis(true);
							redrawAnnotationPanels();						
						}	
						else { //hide
							gataParams.setProteinVis(false);
							gataParams.setRNAVis(false);
							gataParams.setDNAVis(false);
							gataParams.setLabelsVis(false);
							redrawAnnotationPanels();
						}
					}
				}
				else if(paths[1].equals("Protein")){
					if (name.equals("Show")){
						if (gataParams.isProteinVis()==false){
							gataParams.setProteinVis(true);
							redrawAnnotationPanels();
						}						
					}	
					else { //hide
						if (gataParams.isLabelsVis()){
							gataParams.setProteinVis(false);
							redrawAnnotationPanels();
						}
					
					}
				}
				else if(paths[1].equals("RNA")){
					if (name.equals("Show")){
						if (gataParams.isRNAVis()==false){
							gataParams.setRNAVis(true);
							redrawAnnotationPanels();
						}						
					}	
					else {
						if (gataParams.isRNAVis()){
							gataParams.setRNAVis(false);
							redrawAnnotationPanels();
						}
					
					}
				}
				else if(paths[1].equals("DNA")){
					if (name.equals("Show")){
						if (gataParams.isDNAVis()==false){
							gataParams.setDNAVis(true);
							redrawAnnotationPanels();
						}						
					}	
					else {
						if (gataParams.isDNAVis()){
							gataParams.setDNAVis(false);
							redrawAnnotationPanels();
						}
					
					}	
				}
				else if (paths[1].startsWith("Scale Ruler")){
					//set which params should be used
					if (paths[1].endsWith("Ref")) genParams = refParams;
					else genParams = compParams;

					if (name.equals("Show")){
						if (genParams.isScaleRulerVis()==false){
							genParams.setScaleRulerVis(true);
							genParams.setScaleBarVis(true);
							reallyRedrawAnnotationPanel(genParams);
						}						
					}	
					else {
						if (genParams.isScaleRulerVis()){
							genParams.setScaleRulerVis(false);
							genParams.setScaleBarVis(false);
							genParams.getScaleBar().clearScaleLineArrays();
							reallyRedrawAnnotationPanel(genParams);
						}
					
					}
				}
				else if(paths[1].equals("Scale Track Colors By Score")){
				//0 = linear, 1 = log10, 2 = lnx, 3= reset to solid, 4 = log2
					if (name.equals("Linear")) colorizeGenericGlyphs(0);
					else if (name.equals("Log(10)")) colorizeGenericGlyphs(1);
					else if (name.equals("Log(2)")) colorizeGenericGlyphs(4);
					else if (name.equals("Lnx")) colorizeGenericGlyphs(2);
					else colorizeGenericGlyphs(3);
					if (gffRefPres) refGlyphPanel.repaint();
					if (gffCompPres) compGlyphPanel.repaint();
				}

				else if(paths[1].equals("Tracks")){
					//show hide all tracks
					if (paths.length==2){
						if(name.equals("Show")){
							if (gffRefPres) GATAUtil.setBooleanArray(refParams.getTrackVis(), true);
							if (gffCompPres) GATAUtil.setBooleanArray(compParams.getTrackVis(), true);
							redrawAnnotationPanels();
						}
						else { //hide tracks 
							if (gffRefPres) GATAUtil.setBooleanArray(refParams.getTrackVis(), false);
							if (gffCompPres) GATAUtil.setBooleanArray(compParams.getTrackVis(), false);
							redrawAnnotationPanels();
						}
					}
					//show hide all track labels
					else {
						if(name.equals("Show")){
							if (gffRefPres) GATAUtil.setBooleanArray(refParams.getTrackLabelVis(), true);
							if (gffCompPres) GATAUtil.setBooleanArray(compParams.getTrackLabelVis(), true);
							redrawAnnotationPanels();
						}
						else { //hide track labels 
							if (gffRefPres) GATAUtil.setBooleanArray(refParams.getTrackLabelVis(), false);
							if (gffCompPres) GATAUtil.setBooleanArray(compParams.getTrackLabelVis(), false);
							redrawAnnotationPanels();
						}
					}
					//see menu for code to modify all track thickness
				}
				
				//modify individual tracks
				else {
					//get track number and ref/comp
					String[] items = paths[1].split(" "); //RefTrk or CompTrk is in items[0]
					int trackNum = Integer.parseInt(items[1]) - 1;
					
					//set which params should be used
					if (items[0].equals("RefTrk")) genParams = refParams;
					else genParams = compParams;
					
					//check if they want to modify the label vis or the track vis
					if (paths.length>3){
						boolean[] trackLabelVis = genParams.getTrackLabelVis();
						if (name.equals("Show")){
							if (trackLabelVis[trackNum]==false){
								trackLabelVis[trackNum]=true;
								reallyRedrawAnnotationPanel(genParams);
							}						
						}	
						else {
							if (trackLabelVis[trackNum]==true){
								trackLabelVis[trackNum]=false;
								reallyRedrawAnnotationPanel(genParams);
							}
						}
					}
					else{				
						//get whether they want to show or hide a track
						boolean[] trackVis = genParams.getTrackVis();
						if(name.equals("Show")){
							if (trackVis[trackNum]==false){						
								trackVis[trackNum]=true;
								reallyRedrawAnnotationPanel(genParams);
							}
						}
						else { //hide a particular track 
							if (trackVis[trackNum]==true){
								trackVis[trackNum]=false;
								reallyRedrawAnnotationPanel(genParams);
							}							
						}
					}				
				}
			}
			//by default Alignment menu was called
			else{
				//get alignments from AlignPanel, get compPane
				Alignment[] alignments = alignPanel.getAlignments();
				JScrollPane compPane = gataParams.getGataPlotter().getCompScrollPane();
				JSplitPane innerSplitPane = gataParams.getGataPlotter().getInnerSplitPane();
				AlignPanel alignPanel = gataParams.getAlignPanel();
				boolean vis;
				//show comparative sequence and connecting lines
				if (name.equals("Show")){
					vis = true;
					if (innerSplitPane!=null) innerSplitPane.setDividerLocation(innerDividerLoc);
				}
				//hide em
				else {
					vis = false;
					if (innerSplitPane!=null) innerDividerLoc = innerSplitPane.getDividerLocation();
				}
				//run thru and set boolean 
				for (int i=alignments.length-1; i>=0; i--) alignments[i].setCompVisible(vis);
				alignPanel.setDrawCompBox(vis);
				if (compPane!=null)compPane.setVisible(vis);
				//remake visible alignments
				alignPanel.makeVisableAlignments();
				//repaint the panel
				alignPanel.repaint();
			}
		}
	}
	/**helper method, 0 = linear, 1 = log10, 2 = lnx, 3= reset to solid 4= log2*/				
	public void colorizeGenericGlyphs(int colorizer){
		if (gffRefPres) refAnnotation.colorizeTrackedGenericGlyphs(colorizer);
		if (gffCompPres) compAnnotation.colorizeTrackedGenericGlyphs(colorizer);
	}
	
	public void redrawAlignPanel(){
		//gataParams.clearLines();  might be needed!!!!!!!
		alignPanel.makeAllShapes(gataParams.getWidth(),gataParams.getHeight());	
		alignPanel.repaint();
	}
	/**basically draws annotation twice, once to estimate the size, the second time for real*/
	public void redrawAnnotationPanels(){
		if (gffRefPres) reallyRedrawAnnotationPanel(refParams);
		if (gffCompPres) reallyRedrawAnnotationPanel(compParams);
	}
	public void reallyRedrawAnnotationPanel(AnnoSpecParams asp){
		Annotation anno = asp.getAnnotation();
		anno.drawDryRunAnnotation();		
		int smallestY = asp.getSmallestGlyphPanelY();
		int biggestY= asp.getBiggestGlyphPanelY();
		asp.setAnnoPanelHeight(biggestY-smallestY);
		asp.setReferenceY(asp.getReferenceY()-smallestY+200); //for buffer
		asp.setScaleRulerY(asp.getScaleRulerY()-smallestY+200); //for buffer
		anno.redrawAnnotation();		
	}

	private class MenuAction extends AbstractAction {
		private String path;
		private String name;
		public MenuAction(String path, String name) {
			super(name);
			this.name = name.trim();
			this.path = path.trim();
		}
						
		public void actionPerformed(ActionEvent event) {
			//System.out.println("Menu Path- "+path+ "  Name- "+ text + "  (path split on :)");
			String[] paths = path.split(":");
			if (paths.length==1){
				//File
				if(path.equals("File")){
					if(name.equals("Open New Alignment")) new GATAPlotter();
					else if(name.equals("Close Alignment")) {
						new GATAPlotter();
						Console console = gataParams.getConsole();
						if (console!=null) console.dispose();
						ToolsFrame tf = gataParams.getToolsFrameRef();
						if(tf!=null)tf.dispose();
						dispose();
					}
					else if(name.equals("Save GATAPlotter Settings")) {
						plotterPreferences.fetchAndAssignGATAParams(gataParams);
						GATAUtil.saveObject(new File("GATAPlotterPreferences"), plotterPreferences);
					}
					
					else if(name.equals("Save GATAPlot Image")) {
						//get place to save it, text, pixel width
						new SaveImageForm(gataParams);
					}
					
					//quit gataplotter
					else System.exit(0);
				}
				else if(path.equals("Alignment")){
					if (name.equals("Window Width")) {
						//fetch old number
						int oldNum = gataParams.getWidth();
						//fetch new number
						int newNum =GATAUtil.makeIntInputDialog("Alignment Window Width","Enter a new width.",oldNum);
						//check against old
						if (oldNum!=newNum){
							//clearlines reset width
							gataParams.setWidth(newNum);
							alignPanel.setAlignFrameWidth(newNum); //needed to force zoom update
							gataParams.setZoomedWidth(newNum);
							//redraw plotter
							int h = gataParams.getHeight();
							alignPanel.makeAllShapes(newNum,h);
							alignPanel.setPreferredSize(new Dimension(newNum,h+200));
							//redraw Annotation
							if (gffRefPres) refAnnotation.redrawAnnotation();	
							if (gffCompPres) compAnnotation.redrawAnnotation();
							//resize frame
							pack();
							gataParams.getAlignPanel().repaint();
							//reset location of all windows
						}
					} 
					else if (name.equals("Window Height")){
						//fetch old number
						int oldNum = gataParams.getHeight();
						//fetch new number
						int newNum =GATAUtil.makeIntInputDialog("Alignment Window Height","Enter a new height.",oldNum);
						//check against old
						if (oldNum!=newNum){
							//clearlines reset height
							gataParams.setHeight(newNum);
							//redraw plotter
							int w = gataParams.getWidth();
							gataParams.getAlignPanel().makeAllShapes(w,newNum);
							gataParams.getAlignPanel().setPreferredSize(new Dimension(w,newNum+400)); //400 pix buffer
							//redraw Annotation
							if (gffRefPres) refAnnotation.redrawAnnotation();	
							if (gffCompPres) compAnnotation.redrawAnnotation();
							//resize frame
							pack();
							gataParams.getAlignPanel().repaint();
							//reset location of all windows
						}
					}
					
					else if (name.equals("Dist Btw Top & Ref Seq")) {
						double oldNum = gataParams.getA();
						double newNum =GATAUtil.makeDoublePercentInputDialog("DNA Glyph Spacing", "Enter a new percent of panel height ratio for the distance between\n"
							+"    the top of the panel and the top of the reference sequence",oldNum);
						if (oldNum!=newNum){
							gataParams.setA(newNum);
							redrawAlignPanel();
						}
					}
						 
					else if (name.equals("Dist Btw Top & Comp Seq"))  {
						double oldNum = gataParams.getE();
						double newNum =GATAUtil.makeDoublePercentInputDialog("DNA Glyph Spacing", "Enter a new percent of panel height ratio for the distance between\n"
							+"    the top of the panel and the top of the comparative sequence.",oldNum);
						if (oldNum!=newNum){
							gataParams.setE(newNum);
							redrawAlignPanel();
						}
					}
					
					else if (name.equals("Thickness of Seq Bars")) {
						double oldNum = gataParams.getD();
						double newNum =GATAUtil.makeDoublePercentInputDialog("DNA Glyph Thickness", "Enter a new percent of panel height ratio for the thickness of\n"
							+"    the DNA Sequence Glyphs.",oldNum);
						if (oldNum!=newNum){
							gataParams.setD(newNum);
							redrawAlignPanel();
						}
					}
				
					else if (name.equals("Dist Btw Left & Seq Bar")) {
						double oldNum = gataParams.getB();
						double newNum =GATAUtil.makeDoublePercentInputDialog("DNA Glyph Spacing", "Enter a new percent of panel width ratio for the distance between\n"
							+"    the left side of the panel and the left side of the largest DNA glyph.",oldNum);
						if (oldNum!=newNum){
							gataParams.setB(newNum);
							redrawAlignPanel();
						}
					}
					
					else if (name.equals("Width of Largest Seq Bar")) {
						double oldNum = gataParams.getC();
						double newNum =GATAUtil.makeDoublePercentInputDialog("Width DNA Glyph", "Enter a new percent of panel width ratio for \n"
							+"    the width of the largest DNA glyph",oldNum);
						if (oldNum!=newNum){
							gataParams.setC(newNum);
							redrawAlignPanel();
						}
					}
					else if (name.equals("Set Nucleotides Per Pixel")) {
						double oldNum = gataParams.getNtPerPixel();
						double newNum =GATAUtil.makeDoubleInputDialog("Nucleotides Per Pixel", "Set a new scale: nucleotides per pixel.",oldNum);
						if (oldNum!=newNum){
							int newWidth = (int)Math.round((oldNum*(double)gataParams.getWidth())/newNum);
							//clearlines reset width
							gataParams.setWidth(newWidth);
							alignPanel.setAlignFrameWidth(newWidth); //needed to force zoom update
							gataParams.setZoomedWidth(newWidth);
							//redraw plotter
							int h = gataParams.getHeight();
							alignPanel.makeAllShapes(newWidth,h);
							alignPanel.setPreferredSize(new Dimension(newWidth,h+200));
							//redraw Annotation
							if (gffRefPres) refAnnotation.redrawAnnotation();	
							if (gffCompPres) compAnnotation.redrawAnnotation();
							//resize frame
							pack();
							gataParams.getAlignPanel().repaint();
						}
					}
					
					else if (name.equals("Fetch Conserved Sequences")) {
					  if (conservedSeqs ==null) conservedSeqs = new ConservedSeqs(gataParams);
					  conservedSeqs.generateConservedSeqs();
					  conservedSeqs.printSeqsToConsole();
					}					
					else if (name.equals("GATAligner Parameters")) {
						SimpleTextWindow t = new SimpleTextWindow(
											"GATAligner Parameters",450,300,gataParams.getAlignParams().toHTMLString());
					}
				}
				else if (path.equals("Annotation")){
					if (name.equals("Pix Btw Labels & Trans Groups")) {
						int oldNum = (int)(gataParams.getPixBtwTransGrpLabel());
						int newNum =GATAUtil.makeIntInputDialog("Annotation Glyph Spacing", "Enter the number of pixels to place between\n"
							+"    the labels and the first Transcriptional Group.",oldNum);
						if (oldNum!=newNum){
							gataParams.setPixBtwTransGrpLabel(newNum);
							redrawAnnotationPanels();
						}
					}
					
					else if (name.equals("Pix Btw Trans Groups")){
							int oldNum = (int)(gataParams.getPixBtwTransGrps());
							int newNum =GATAUtil.makeIntInputDialog("Annotation Glyph Spacing", "Enter the number of pixels to place between\n"
								+"    Transcriptional Groups.",oldNum);
							if (oldNum!=newNum){
								gataParams.setPixBtwTransGrps(newNum);
								redrawAnnotationPanels();
						}
					}					
					else if (name.equals("Pix Btw Trans Groups & DNA")){
						int oldNum = (int)(gataParams.getPixBtwTransGrpsGene());
						int newNum =GATAUtil.makeIntInputDialog("Annotation Glyph Spacing", "Enter the number of pixels to place between\n"
							+"    the last Transcriptional Group and DNA glyph.",oldNum);
						if (oldNum!=newNum){
							gataParams.setPixBtwTransGrpsGene(newNum);
							redrawAnnotationPanels();
						}
					}			
					
					else if (name.equals("Pix Btw Protein & RNA")) {
						int oldNum = (int)(gataParams.getPixBtwProtRNA());
						int newNum =GATAUtil.makeIntInputDialog("Annotation Glyph Spacing", "Enter the number of pixels to place between\n"
							+"    Protein glyphs and RNA glyphs.",oldNum);
						if (oldNum!=newNum){
							gataParams.setPixBtwProtRNA(newNum);
							redrawAnnotationPanels();
						}
					}			
	

					else if (name.equals("Pix Btw DNA & Tracks")) {
						int oldNum = (int)(gataParams.getPixBtwDNATracks());
						int newNum =GATAUtil.makeIntInputDialog("Annotation Glyph Spacing", "Enter the number of pixels to place between\n"
							+"    DNA glyphs and Tracks.",oldNum);
						if (oldNum!=newNum){
							gataParams.setPixBtwDNATracks(newNum);
							redrawAnnotationPanels();
						}
					}	
					else if (name.equals("Pix Btw Labels & Tracks")) {
						int oldNum = (int)(gataParams.getPixBtwLabelTrack());
						int newNum =GATAUtil.makeIntInputDialog("Annotation Glyph Spacing", "Enter the number of pixels to place between\n"
							+"    Track Labels and Tracks.",oldNum);
						if (oldNum!=newNum){
							gataParams.setPixBtwLabelTrack(newNum);
							redrawAnnotationPanels();
						}
					}	
				
					else if (name.equals("Pix Btw Tracks")) {
						int oldNum = (int)(gataParams.getPixBtwTracks());
						int newNum =GATAUtil.makeIntInputDialog("Annotation Glyph Spacing", "Enter the number of pixels to place between\n"
							+"Tracks.",oldNum);
						if (oldNum!=newNum){
							gataParams.setPixBtwTracks(newNum);
							redrawAnnotationPanels();
						}
					}	
					
					else if (name.equals("Thin Line Thickness")) {
						int oldNum = (int)(gataParams.getPixThinLineThickness());
						int newNum =GATAUtil.makeIntInputDialog("Glyph Line Thickness", "Enter the number of pixels to make glyph\n"
							+"    thin lines (DNA, RNA, Protein)",oldNum);
						if (oldNum!=newNum){
							gataParams.setPixThinLineThickness(newNum);
							redrawAnnotationPanels();
						}
					}
						
					else if (name.equals("Normal Line Thickness")) {
						int oldNum = (int)(gataParams.getPixNormalLineThickness());
						int newNum =GATAUtil.makeIntInputDialog("Glyph Line Thickness", "Enter the number of pixels to make glyph\n"
							+"     normal lines (DNA coding segments).",oldNum);
						if (oldNum!=newNum){
							gataParams.setPixNormalLineThickness(newNum);
							redrawAnnotationPanels();
						}
					}
						
					else if (name.equals("Fat Line Thickness")) {
						int oldNum = (int)(gataParams.getPixFatLineThickness());
						int newNum =GATAUtil.makeIntInputDialog("Glyph Line Thickness", "Enter the number of pixels to make glyph\n"
							+"     fat lines (RNA exons).",oldNum);
						if (oldNum!=newNum){
							gataParams.setPixFatLineThickness(newNum);
							redrawAnnotationPanels();
						}
					}	
	
					else if (name.equals("Arrow Thickness")) {
						int oldNum = (int)(gataParams.getPixHeightArrows());
						int newNum =GATAUtil.makeIntInputDialog("Glyph Line Thickness", "Enter the number of pixels to reset\n"
							+"    the height of Arrows and Butts.",oldNum);
						if (oldNum!=newNum){
							gataParams.setPixHeightArrows(newNum);
							redrawAnnotationPanels();
						}
					}	
					else if (name.equals("Background Color")) {
						Color oldColor = gataParams.getBackgroundColor();
						Color newColor = JColorChooser.showDialog(gataFrame,"Annotation Window Background Color", oldColor);
						if (oldColor.equals(newColor)==false){
							gataParams.setBackgroundColor(newColor);
							if (gffRefPres) {
								refGlyphPanel.setBackground(newColor);
								refGlyphPanel.repaint();
							} 
							if (gffCompPres) {
								compGlyphPanel.setBackground(newColor);
								compGlyphPanel.repaint();
							} 
						}
					}		
					else if (name.equals("Label Color")) {
						Color oldColor = gataParams.getLabelColor();
						Color newColor = JColorChooser.showDialog(gataFrame,"Label Color", oldColor);
						if (oldColor.equals(newColor)==false){
							gataParams.setLabelColor(newColor);
							if (gffRefPres) {
								refGlyphPanel.setDrawParams();
								refGlyphPanel.repaint();
							} 
							if (gffCompPres) {
								compGlyphPanel.setDrawParams();
								compGlyphPanel.repaint();
							} 
						}
					}	
				}
				else if (path.equals("Windows")){
					if (name.equals("Show All")) {
						setWindowsVisibility(true);
					}
					else if (name.equals("Hide All")){
						setWindowsVisibility(false);
					}
					else if (name.equals("Redraw Using Defaults")) {
						gataParams.getPlotterPreferences().setPlotterDefaults(gataParams);
						redrawAlignPanel();
						if(gffRefPres){
							if (refParams.getAnyGenericGlyphs()) {
								refAnnotation.makeGenericGlyphs();
								refAnnotation.setTrackDefaults();
								} 
						}
						if(gffCompPres){
							if (compParams.getAnyGenericGlyphs()) {
								compAnnotation.makeGenericGlyphs();
								compAnnotation.setTrackDefaults();
								} 
						}
						redrawAnnotationPanels();

					}
				}
				else { // in documentation pull html files
					if (documentWindow==null) documentWindow =new TextWindow(name,600,600,"file:./documents/"+name.replaceAll("\\s","")+".html");
					else documentWindow.loadNewURL("file:./documents/"+name.replaceAll("\\s","")+".html");	
				}
			}
			else{
				if (paths[1].equals("DNA")){
					Color oldColor = gataParams.getDNAColor();
					Color newColor = JColorChooser.showDialog(gataFrame,"Annotation DNA Color", oldColor);
					if (newColor!=null){
						if (oldColor.equals(newColor)==false){
							gataParams.setDNAColor(newColor);
							if (gffRefPres) {
								System.out.println("resetting DNA color ref"); 
								refGlyphPanel.setDrawParams();
								refGlyphPanel.repaint();
							} 
							if (gffCompPres) {
								System.out.println("resetting DNA color comp");
								compGlyphPanel.setDrawParams();
								compGlyphPanel.repaint();
							} 
						}
					}
				}	
				else if (paths[1].equals("RNA")){
					Color oldColor = gataParams.getRNAColor();
					Color newColor = JColorChooser.showDialog(gataFrame,"Annotation RNA Color", oldColor);
					if (newColor!=null){
						if (oldColor.equals(newColor)==false){
							gataParams.setRNAColor(newColor);
							if (gffRefPres) {
								refGlyphPanel.setDrawParams();
								refGlyphPanel.repaint();
							} 
							if (gffCompPres) {
								compGlyphPanel.setDrawParams();
								compGlyphPanel.repaint();
							} 
						}
					}
				}	
				else if (paths[1].equals("Protein")){
					Color oldColor = gataParams.getProteinColor();
					Color newColor = JColorChooser.showDialog(gataFrame,"Annotation Protein Color", oldColor);
					if (newColor!=null){
						if (oldColor.equals(newColor)==false){
							gataParams.setProteinColor(newColor);
							if (gffRefPres) {
								refGlyphPanel.setDrawParams();
								refGlyphPanel.repaint();
							} 
							if (gffCompPres) {
								compGlyphPanel.setDrawParams();
								compGlyphPanel.repaint();
							} 
						}
					}
				}	
				else if (paths[1].startsWith("Scale Ruler")){
					//set which params should be used
					if (paths[1].endsWith("Ref")) genParams = refParams;
					else genParams = compParams;
					if (name.equals("Color")){
						Color oldColor = genParams.getScaleLineColor();
						Color newColor = JColorChooser.showDialog(gataFrame,"Annotation Scale Bar Ruler Color", oldColor);
						if (newColor!=null){
							if (oldColor.equals(newColor)==false){
								genParams.setScaleLineColor(newColor);
								GlyphPanel pan = genParams.getGlyphPanel();
									pan.setDrawParams();
									pan.repaint();
							}
						}
					}
					else { //they want to change the location
						int oldNum = (int)genParams.getScaleRulerY();
						int newNum =GATAUtil.makeIntInputDialog("Scale Ruler Location", "Enter a new Y coordinate to draw the Scale Ruler.",oldNum);
						if (oldNum!=newNum){
							genParams.setScaleRulerY((double)newNum);
							reallyRedrawAnnotationPanel(genParams);
						}
					}
				}
				//modify all tracks
				else if (paths[1].startsWith("Tracks")){
					System.out.println("in Tracks, want to change "+name);
					float defaultThickness = gataParams.getDefaultTrackThickness();
					Object ob =JOptionPane.showInputDialog(null,"Enter a thickness in pixels for all generic tracks",
						"Set All Track Thickness",JOptionPane.PLAIN_MESSAGE,null,null,Integer.toString((int)defaultThickness));
					if (ob != null) { // if user hits cancel or closes the window
						float newNum;
						try{
							newNum = Float.parseFloat((String) ob);
						}
						catch (NumberFormatException e){
							newNum = defaultThickness;
						}
		
						if (gffRefPres){
							float[] tts = refParams.getTrackThickness();
							GATAUtil.setFloatArray(tts, newNum);
							refParams.setTrackThickness(tts);
							}
						if (gffCompPres) {
							float[] tts = compParams.getTrackThickness();
							GATAUtil.setFloatArray(tts, newNum);
							compParams.setTrackThickness(tts);
						}
						gataParams.setDefaultTrackThickness(newNum);
						redrawAnnotationPanels();
					}	
				}
				//into individual tracks	
				else {
					//get track number and ref/comp
					String[] items = paths[1].split(" "); //RefTrk or CompTrk is in items[0]
					int trackNum = Integer.parseInt(items[1]) - 1;
					
					//set which params should be used
					if (items[0].equals("RefTrk")) genParams = refParams;
					else genParams = compParams;

					//get whether they want to change the color or thickness
					if(name.equals("Color")){
						GenericGlyph[][] tgg = genParams.getTrackedGenericGlyphs();
						Color baseTrackColor = tgg[trackNum][0].getColor();
						Color baseColor = new Color(baseTrackColor.getRed(), baseTrackColor.getGreen(), baseTrackColor.getBlue());
						Color newColor = JColorChooser.showDialog(gataFrame,paths[1]+" "+paths[2]+" Color", baseColor);
						if (newColor!=null){
							if (baseColor.equals(newColor)==false){ //selected a new color
								//run through tgg's and reset to new color but using old alpha
								int len = tgg[trackNum].length;
								for (int i=0; i<len; i++){
									Color oldColor = tgg[trackNum][i].getColor();
									int alpha = oldColor.getAlpha();
									Color replacementColor = new Color (newColor.getRed(),newColor.getGreen(),newColor.getBlue(),alpha);
									tgg[trackNum][i].setColor(replacementColor);
								}
								genParams.getGlyphPanel().repaint();
							}
						}
					}
					//by default they want to modify track thickness
					else {
						float[] tts = genParams.getTrackThickness();
						int oldNum = (int)(tts[trackNum]);
						int newNum =GATAUtil.makeIntInputDialog(paths[1]+" "+paths[2]+" Thickness", "Enter a track thickness in pixels for "+paths[1]+" "+paths[2]+".",oldNum);
						if (oldNum!=newNum){
							tts[trackNum]= (float)newNum;
							genParams.setTrackThickness(tts);
							reallyRedrawAnnotationPanel(genParams);
						}
						
					}
				} 
			}	
		}
	}
	public void setAlignPanel(AlignPanel panel) {
		alignPanel = panel;
	}
	public void setCompGlyphPanel(GlyphPanel panel) {
		compGlyphPanel = panel;
	}
	public void setRefGlyphPanel(GlyphPanel panel) {
		refGlyphPanel = panel;
	}
	
	/**Add new windows/panels here to enable them to be shown or hidden in sync with the main frame*/
	public void setWindowsVisibility(boolean setting){
		if (gataParams.getConsole()!=null)gataParams.getConsole().setVisible(setting);
		if (gataParams.getToolsFrameRef()!=null)gataParams.getToolsFrameRef().setVisible(setting);
		if (conservedSeqs!= null && conservedSeqs.getConsole()!=null) conservedSeqs.getConsole().setVisible(setting);
		if (documentWindow!= null) documentWindow.setVisible(setting);				
	}	
	
	/**Used to watch for window closing or hiding events, syncs up other windows*/
	private class WindowStateMonitor extends WindowAdapter{
		public void windowClosing(WindowEvent e){
			if (gataParams.getConsole()!=null)gataParams.getConsole().dispose();
			if (gataParams.getToolsFrameRef()!=null)gataParams.getToolsFrameRef().dispose();
			if (conservedSeqs!= null && conservedSeqs.getConsole()!=null) conservedSeqs.getConsole().dispose();
			if (documentWindow!= null) documentWindow.dispose();			
		}
		public void windowIconified(WindowEvent e){
			setWindowsVisibility(false);
		}
		public void windowDeiconified(WindowEvent e){
			setWindowsVisibility(true);
		}
		
	}

}