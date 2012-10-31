package gata.geneGlyphs;

import gata.main.*;

import java.awt.*;

import java.util.*;

import util.bio.annotation.*;
import util.bio.parsers.gff.GadFlyGffExtractor;

/**
 * @author nix
 * Contains methods and fields to construct annotation glyphs.
 */
public class Annotation {
	private GadFlyGffExtractor gff;
	private GeneGroup[] geneGrps;
	private GeneRepGlyph[] geneRepGlyphs;
	private GenericFeature[] genericFeatures;
	private GenericGlyph[][] trackedGenericGlyphs;
	private int numberOfTracks = 0;
	private GlyphPanel glyphPanel;
	private AnnoSpecParams annoParams;
	private GATAParams gataParams;
	private ScaleBar scaleBar;
	private int trackConverter = 3; //number used to scale track scores

	public Annotation(AnnoSpecParams annoParams, GATAParams gataParams) {
		this.annoParams = annoParams;
		this.gataParams = gataParams;
		//make GadFlyGFFExtractor- Parses, extracts and builds Gene feature representations from a GFF file
		gff = new GadFlyGffExtractor(annoParams.getGffFile(), 0, 1000000000);		
		ArrayList alGG = gff.getGeneGrps();
		geneGrps = new GeneGroup[alGG.size()];
		alGG.toArray(geneGrps);
		Arrays.sort(geneGrps);
		//Make glyphs from the gff object,call once
		makeGlyphs();
		//DrawGlyphs, call everytime you want to redraw after a user change
		drawGlyphs();
		//make generic feature glyphs (might not be any) these are the non standard user added items
		if (gff.isGenericFeaturesFound()) {
			annoParams.setAnyGenericGlyphs(true);
			assignGenericFeatureTracks(); //call once
			makeGenericGlyphs(); //call once
			setTrackDefaults(); //call once
			initializeTrackLabels(); //call once
			drawGenericGlyphs(); //call each time things change, need remake if you hide the generic glyphs
		}
		else annoParams.setBiggestGlyphPanelY((int)annoParams.getScaleRulerY()+10);

		//make scalebar and scaleRuler
		scaleBar = new ScaleBar(annoParams, gataParams);
		
		//scaleBar.makeScaleBar(); //call each time the ntPerPix changes
		scaleBar.makeScaleRuler(); //ditto

		//make panel
		glyphPanel = new GlyphPanel(annoParams, gataParams);
		glyphPanel.setBackground(gataParams.getBackgroundColor());
		//fetch all the draw parameters from GlyphParams, call everytime a user setting changes
		glyphPanel.setDrawParams();
		//converts all of the arraylists and some other stuff to real line arrays, call after makingGlyphs and drawGlyphs
		glyphPanel.setLine2DArrays(); //converts the ArrayList of lines to []
	}
	
	public void redrawAnnotation(){
		annoParams.clearLines();
		glyphPanel.setDrawParams();
		glyphPanel.setBackground(gataParams.getBackgroundColor());
		drawGlyphs();
		if (gff.isGenericFeaturesFound()) drawGenericGlyphs();
		//if (annoParams.isScaleBarVis()) scaleBar.makeScaleBar(); 
		if (annoParams.isScaleRulerVis()) scaleBar.makeScaleRuler();
		glyphPanel.setLine2DArrays();
		glyphPanel.setPreferredSize(new Dimension(gataParams.getWidth(), annoParams.getAnnoPanelHeight()+400)); //400pix buffer
		glyphPanel.repaint();
	}
	public void drawDryRunAnnotation(){
		glyphPanel.setDrawParams();
		glyphPanel.setTrackParams();
		drawGlyphs();
		if (gff.isGenericFeaturesFound()) drawGenericGlyphs();
		//if (annoParams.isScaleBarVis()) scaleBar.makeScaleBar(); 
		if (annoParams.isScaleRulerVis()) scaleBar.makeScaleRuler();
		glyphPanel.setLine2DArrays();
	}
	
	public void setAnnotationBackground(){
		glyphPanel.setBackground(gataParams.getBackgroundColor());
	}
	
	/**Call each time parameters effecting GenericGlyphs change*/
	public void drawGenericGlyphs() {
		//wipe clean trackLabels, additive ArrayLists
		initializeTrackLabels();
		//set yCoor
		float yCoor =
			(float) (annoParams.getReferenceY() + gataParams.getPixBtwDNATracks());
		//get thickness of DNA rep and add to yCoor
		double arrow = gataParams.getPixHeightArrows();
		double normalLine = gataParams.getPixNormalLineThickness() / 2;
		if (arrow > normalLine)
			yCoor += arrow;
		else
			yCoor += normalLine;
		//call draw using yCoor
			boolean[] trackVis = annoParams.getTrackVis();
			boolean[] trackLabelVis = annoParams.getTrackLabelVis();
			float[] trackThickness = annoParams.getTrackThickness();
			double pixBtwTracks = gataParams.getPixBtwTracks();
			double letterHeight = gataParams.getHeightLabels();
			double pixBtwLabelTrack = gataParams.getPixBtwLabelTrack();
			int len = trackedGenericGlyphs.length;
			float pixPerNt = (float) gataParams.getPixPerNt();
			float ntAtXZero = (float) annoParams.getNtAtXZero();
			ArrayList[] trackLabels = annoParams.getTrackLabels();
		//run through all the tracks
		for (int i = 0; i < len; i++) {
			//check if track is visible
			if (trackVis[i]) {			
			int len2 = trackedGenericGlyphs[i].length;
				//draw each generic glyph in the track
				for (int j = 0; j < len2; j++) {
					trackedGenericGlyphs[i][j].draw(
						yCoor,
						ntAtXZero,
						trackThickness[i],
						(float) arrow,
						pixPerNt,
						trackLabels[i],
						(float) letterHeight,
						(float) pixBtwLabelTrack,
						trackLabelVis[i]);
				}
				//calculate thickness based on whether trackLabel is vis
				if (trackLabelVis[i])
					yCoor += letterHeight;
				if (trackLabelVis[i] && trackVis[i])
					yCoor += pixBtwLabelTrack;
				if (trackVis[i])
					yCoor += arrow + trackThickness[i];
				yCoor += pixBtwTracks;
			}
		}
		annoParams.setBiggestGlyphPanelY((int)yCoor);
	}

	/**This makes a sorted array of GenericGlyphs broken down by track*/
	public void makeGenericGlyphs() {
		//sort by track
		Arrays.sort(genericFeatures);
		ArrayList ggs = new ArrayList();
		trackedGenericGlyphs = new GenericGlyph[numberOfTracks][];
		int len = genericFeatures.length;
		int trackCounter = 0;
		for (int i = 0; i < len; i++) {
			boolean flag = true;
			while (flag) {
				//add new GenericGlyphs
				if (trackCounter == genericFeatures[i].getTrackNumber()) {
					ggs.add(new GenericGlyph(genericFeatures[i]));
					flag = false;
					if (i == len - 1) {
						GenericGlyph[] gg = new GenericGlyph[ggs.size()];
						ggs.toArray(gg);
						trackedGenericGlyphs[trackCounter] = gg;
					}

				} else {
					GenericGlyph[] gg = new GenericGlyph[ggs.size()];
					ggs.toArray(gg);
					trackedGenericGlyphs[trackCounter] = gg;
					ggs.clear();
					trackCounter++;
				}
			}
		}
		annoParams.setTrackedGenericGlyphs(trackedGenericGlyphs);
	}
	
	/**This method converts scores in generic glyphs to changes in alpha, the opacity 
	 * number for the glyph color.  The lower the score the lower the alpha.  The trickie business is
	 * using different scales since sometimes the scores are in log or lnx units.  70 is
	 * the lowest alpha allowed, everything else normalized from 70 to 255*/
	public void colorizeTrackedGenericGlyphs(int converter){
	//		0 = linear, 1 = log10, 2 = lnx, 3= reset to solid, 4= log2
		trackConverter = converter;
		int len = trackedGenericGlyphs.length;
		//calculate min and max score for each track
		for (int i=0; i<len; i++){
			double minScore = trackedGenericGlyphs[i][0].getScore();
			double maxScore = minScore;
			int size = trackedGenericGlyphs[i].length;
			for (int j=1; j<size; j++){
				double score = trackedGenericGlyphs[i][j].getScore();
				if (score<minScore) minScore = score;
				if (score>maxScore) maxScore = score;
			}
			//skip if no scores found
			if (minScore ==0 && maxScore==0) continue;
			//possibly recalculate min max based on log scaling
			if (converter ==1){ //log10
				minScore = Math.pow(10,minScore);
				maxScore = Math.pow(10,maxScore);
			}
			else if (converter ==4){ //log2
				minScore = Math.pow(2,minScore);
				maxScore = Math.pow(2,maxScore);
			}
			else if (converter==2){ //lnx
				minScore = Math.exp(minScore);
				maxScore = Math.exp(maxScore);
			}
			double rangeNum = 185/(maxScore-minScore);
			//if max and min are the same, ie only one annotation glyph then don't scale it
			if (maxScore==minScore)converter=3; 
			
			//reset color alpha based on score; 255 opaque,  0 transparent
			//alpha = (185/range)(score-minScore) + 70; to keep even lowest hits visible
			for (int k=0; k<size; k++){
				//calculate alpha based on score and converter
				int alpha;
				if (converter==3) alpha = 255; //reset no scaling
				else {
					double score = trackedGenericGlyphs[i][k].getScore();
					if (converter==1) score = Math.pow(10,score);
					else if (converter==4) score = Math.pow(2,score);
					else if (converter==2) score = Math.exp(score);
					alpha = (int)Math.round(rangeNum*(score-minScore) + 70);
				}
				//reset color	
				Color oldColor = trackedGenericGlyphs[i][k].getColor();
				Color newColor = new Color(oldColor.getRed(), oldColor.getGreen(), oldColor.getBlue(), alpha);
				trackedGenericGlyphs[i][k].setColor(newColor);
			}			
		}
	}
	
	public void setTrackDefaults(){
		boolean[] trackVis = new boolean[numberOfTracks];
		boolean[] labelVis = new boolean[numberOfTracks];
		float num = gataParams.getDefaultTrackThickness();
		float[] nums = new float[numberOfTracks];
		for (int i = 0; i < numberOfTracks; i++){
			trackVis[i] = true;
			labelVis[i] = true;
			nums[i] = num;
		}
		annoParams.setTrackVis(trackVis);
		annoParams.setTrackLabelVis(labelVis);
		annoParams.setTrackThickness(nums);	
		//initialize colors
		setInitialTrackColors();
	}
	
	/**Use to initialize all generic glyphs color field*/
	public void setInitialTrackColors(){
		int len = trackedGenericGlyphs.length;
		Color[] trackColors = new Color[len];
		for (int i=0; i<len; i++){
			int size = trackedGenericGlyphs[i].length;
			trackColors[i]=GATAUtil.COLORS[i];
			for (int j=0; j<size; j++){
				trackedGenericGlyphs[i][j].setColor(GATAUtil.COLORS[i]);	
			}
		}
	}

	public void setDefaultTrackThickness(){
		float num = gataParams.getDefaultTrackThickness();
		float[] nums = new float[numberOfTracks];
		for (int i = 0; i < numberOfTracks; i++){
			nums[i] = num;
		}
		annoParams.setTrackThickness(nums);	
	}
	
	public void initializeTrackLabels() {
		ArrayList[] ip = new ArrayList[numberOfTracks];
		for (int i = 0; i < numberOfTracks; i++)
			ip[i] = new ArrayList();
		annoParams.setTrackLabels(ip);
	}
	public void assignGenericFeatureTracks() {
		LinkedHashSet hash = gff.getGenericFeatureHash();
		numberOfTracks = hash.size();
		//convert ArrayList of GenericFeatures to an Array
		ArrayList gfsAL = gff.getGenericFeatures();
		int len = gfsAL.size();
		genericFeatures = new GenericFeature[len];
		gfsAL.toArray(genericFeatures);
		//covert hash to array
		String[] featureTypes = new String[numberOfTracks];
		hash.toArray(featureTypes);
		int min = 0;
		//assign a track number to each GenericFeature, 0,1,2...
		for (int j = 0; j < len; j++) {
			String genericType = genericFeatures[j].getFeatureType();
			boolean flag = true;
			while (flag) {
				String hashALType = featureTypes[min];
				if (genericType.equalsIgnoreCase(hashALType)) {
					genericFeatures[j].setTrackNumber(min);
					flag = false;
				} else {
					min++;
					if (min == numberOfTracks)
						min = 0;
					flag = true;
				}
			}
		}
	}

	/**Call each time the scale changes or user alters a glyph viewablility*/
	public void drawGlyphs() {
		double yCoordinate = annoParams.getReferenceY();
		double yAdj = yCoordinate;
		double smallestY =yCoordinate; //used to determine max size of glyphPanel
		int[] oldGlyphSS = new int[] { 0, 0 };	
		for (int i = geneRepGlyphs.length - 1; i >= 0; i--) {
			//check if overlap, if so alter yCoordinate
			int[] newGlyphSS = {geneRepGlyphs[i].getStart(), geneRepGlyphs[i].getEnd()};
			if ((newGlyphSS[1] >= oldGlyphSS[0]
				&& newGlyphSS[1] <= oldGlyphSS[1])
				|| //overlap
			 (
					newGlyphSS[0] >= oldGlyphSS[0]
						&& newGlyphSS[0] <= oldGlyphSS[1])
				|| //overlap
			 (
					newGlyphSS[0] <= oldGlyphSS[0]
						&& newGlyphSS[1] >= oldGlyphSS[1])
				|| //contained within
			 (
					newGlyphSS[0] >= oldGlyphSS[0]
						&& newGlyphSS[1] <= oldGlyphSS[1])) {//contained within
							

				geneRepGlyphs[i].draw(yAdj);
				//check if new extends past old, if so set old == new
				if (newGlyphSS[0]<oldGlyphSS[0]) {oldGlyphSS = newGlyphSS;
					yAdj = geneRepGlyphs[i].getNonOverlapYCoor();
				}
				 else {
				 	//possibly set smallest
				 	double checker =  geneRepGlyphs[i].getNonOverlapYCoor();
					if (checker<smallestY) smallestY=checker;
				 }
				
			} else {				
				
				geneRepGlyphs[i].draw(yCoordinate);
				oldGlyphSS = newGlyphSS;
				yAdj = geneRepGlyphs[i].getNonOverlapYCoor();
			}
			
			if (yAdj<smallestY) smallestY=yAdj;
		}
		annoParams.setSmallestGlyphPanelY((int)smallestY);
	}
	/**
	 * Call once to make the glyphs used in rendering the gene annotation.*/
	public void makeGlyphs() {
		//for each gene group make glyphs, saving Shape lines in the GlyphParams
		int len = geneGrps.length;

		//make arrays of GeneRepGlyps and TransGrpGlyphs
		geneRepGlyphs = new GeneRepGlyph[len];
		ArrayList transGrpGlyphs = new ArrayList();
		for (int i = len - 1; i >= 0; i--) {
			//build GeneGrpGlyphs
			geneRepGlyphs[i] =
				new GeneRepGlyph(geneGrps[i].getGeneRep(), gataParams, annoParams);
			//build TransGrpGlyphs
			TransGroup[] transGrps = geneGrps[i].getTransGrps();
			for (int j = transGrps.length - 1; j >= 0; j--) {
				TransGrpGlyph tgg = new TransGrpGlyph(transGrps[j], gataParams, annoParams);
				transGrpGlyphs.add(tgg);
			}
			geneRepGlyphs[i].setTransGrpGlyphs(
				(ArrayList) transGrpGlyphs.clone());
			transGrpGlyphs.clear();
		}
		//set GeneRepGlyph ref in params object
		annoParams.setGeneGlyphs(geneRepGlyphs);
	}
	public GlyphPanel getGlyphPanel() {
		return glyphPanel;
	}
	public int getTrackConverter() {
		return trackConverter;
	}

}
