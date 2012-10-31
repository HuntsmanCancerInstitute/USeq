package gata.geneGlyphs;

import gata.main.*;

import java.awt.geom.*;
import java.util.*;


import java.awt.*;

import util.bio.annotation.*;

/**
 * @author Nix
 * Contains info related to drawing a GenericFeature, non standard gff annotation
 */
public class GenericGlyph {
	//fields
	private int start;
	private int end;
	private String name;
	private int orientation;
	private GeneralPath path;
	private int track;
	private ArrayList[] trackLabels;
	private Rectangle2D.Double boundingBox;
	private double score;
	private Color color;
	private GenericFeature genericFeature;
	
	//constructor
	public GenericGlyph(GenericFeature gf) {
		genericFeature = gf;
		start = gf.getStart();
		end = gf.getEnd();
		track = gf.getTrackNumber();
		orientation = gf.getOrientation();
		name = gf.getName();
		score = gf.getScore();
	}
	public String toString(){
		return (
		"SeqName: "+name+" "+start+"-"+end+
		"\nFeature: "+genericFeature.getFeatureType()+
		"\nScore: "+score+
		"\nStrand: "+genericFeature.getStrand()+
		"\nAttributes: "+genericFeature.getAttributes()+
		"\n\n"
		);
	}
	
	public void draw(
		float yCoord,
		float ntAtXZero,
		float trackThickness,
		float arrowHeight,
		float pixPerNt,
		ArrayList trackLabels,
		float letterHeight,
		float pixBtwLabelTrack,
		boolean trackLabelVis) {
		
		float halfTrackThickness = trackThickness/2;
		double startY = yCoord;
		
		//make label, draw first since it is on top
		if (trackLabelVis) {
			trackLabels.addAll(
				GATAUtil.makeLabel(name, start, end, pixPerNt, ntAtXZero, yCoord+5));
			yCoord += letterHeight + pixBtwLabelTrack;
		}
		//type line(s) based on orientation
			yCoord += halfTrackThickness;
			float startX = (start - ntAtXZero) * pixPerNt;
			float stopX = (end+1 - ntAtXZero) * pixPerNt;
		if (orientation == 0) {
			path = new GeneralPath();
			path.moveTo(startX, yCoord);
			path.lineTo(stopX, yCoord);
		} 
		else {
			path =
			GATAUtil.makeArrowLine(
					start,end,
					1,
					yCoord,
					orientation,
					pixPerNt,
					ntAtXZero,
					arrowHeight,
					trackThickness);
		}
		//make bounding box
		boundingBox = new Rectangle2D.Double(startX-2, startY, stopX-startX+4, (yCoord+halfTrackThickness)-startY+2);
	}
	public GeneralPath getPath() {
		return path;
	}
	public Rectangle2D.Double getBoundingBox() {
		return boundingBox;
	}
	public Color getColor() {
		return color;
	}
	public double getScore() {
		return score;
	}
	public void setColor(Color color) {
		this.color = color;
	}

	public GenericFeature getGenericFeature() {
		return genericFeature;
	}
}
