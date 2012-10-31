package gata.geneGlyphs;


import gata.main.*;

import java.util.*;
import java.awt.geom.*;

import util.bio.annotation.*;

/**
 * @author Nix
 * Shape representation for a GeneRep- the DNA component of genomic annotation.
 */
public class GeneRepGlyph {
	//fields
	private GeneRep geneRep;
	private GATAParams gataParams;
	private AnnoSpecParams annoParams;
	private double yCoordinate;
	private double modYCoor;
	private int orientation;
	private int start;
	private int end;
	private String label;
	private double ntAtXZero;
	private double pixPerNt;
	private double transGrpThickness;
	private ArrayList transGrpGlyphs = new ArrayList();
	private double letterHeight;
	private Rectangle2D.Double boundingBox;
	private ArrayList thinLines = new ArrayList();
	private ArrayList thickLines = new ArrayList();
	private double arrowHeight;
	private double halfNormalLineThickness;

	public GeneRepGlyph(GeneRep geneRep, GATAParams gataParams, AnnoSpecParams annoParams) {
		this.geneRep = geneRep;
		this.gataParams = gataParams;
		this.annoParams = annoParams;
		GeneGroup geneGrp = geneRep.getGeneGrp();
		orientation = geneGrp.getOrientation();
		start = geneGrp.getStart();
		end = geneGrp.getEnd();
		label = geneGrp.getLabel();
	}

	//primary methods
	public String toString(){
		return geneRep.getGeneGrp().toString();
	}
	
	public void draw(double yCoordinate) {
		//set scale
		ntAtXZero = annoParams.getNtAtXZero();	
		pixPerNt = gataParams.getPixPerNt();
		//set y coordinate adjustors
		this.yCoordinate = yCoordinate;
		arrowHeight = gataParams.getPixHeightArrows();
		halfNormalLineThickness =
			(double) gataParams.getPixNormalLineThickness() / 2;
		//reset adjusters
		modYCoor = yCoordinate;
		transGrpThickness =0;
		//		calculate transGrpThickness
		if (gataParams.isProteinVis() && gataParams.isRNAVis()) {
			transGrpThickness = (double) gataParams.getPixFatLineThickness()
				+ arrowHeight
				+ gataParams.getPixBtwProtRNA()
				+ gataParams.getPixBtwTransGrps();
		} else {
			if (gataParams.isProteinVis())
				transGrpThickness = arrowHeight
					+ gataParams.getPixBtwTransGrps();
			if (gataParams.isRNAVis())
				transGrpThickness
					= (double) gataParams.getPixFatLineThickness()
					+ gataParams.getPixBtwTransGrps();
		}
		//draw DNA
		if (gataParams.isDNAVis()) {
			//gene arrow line with coding Segs
			thinLines =
				GATAUtil.makeArrowLine(
					start,end,
					1,
					yCoordinate,
					orientation,
					pixPerNt,
					ntAtXZero,
					arrowHeight);
			//coding segments		
			thickLines =
				GATAUtil.makeSegments(
					geneRep.getCodingSegments(),
					yCoordinate,
					pixPerNt,
					ntAtXZero);		
			annoParams.addThinDNALines(thinLines);
			annoParams.addNormalDNALines(thickLines);

			modYCoor -= arrowHeight
				+ halfNormalLineThickness
				+ gataParams.getPixBtwTransGrpsGene();
			//adjust yCoord to reflect DNA line thickness
			if (arrowHeight>halfNormalLineThickness) yCoordinate +=arrowHeight;
			else yCoordinate += halfNormalLineThickness;		
		}
		//draw transGrpGlyphs
		drawTransGrpGlyphs();
		//make label last since it must be ontop
		if (gataParams.areLabelsVis()) {
			modYCoor -= gataParams.getPixBtwTransGrpLabel();
			annoParams.addLabel(GATAUtil.makeLabel(
				label,
				start, end,
				pixPerNt,
				ntAtXZero,
				modYCoor));
			modYCoor -= gataParams.getHeightLabels();
		}
		//add bounding box, don't add to glyphParams since one will need to associate it with this object
		boundingBox =
		GATAUtil.makeRectangle(start, end, yCoordinate+6, modYCoor, pixPerNt, ntAtXZero);
	}
	//draw glyphs
	public void drawTransGrpGlyphs() {
		for (int i = transGrpGlyphs.size() - 1; i >= 0; i--) {
			((TransGrpGlyph) transGrpGlyphs.get(i)).draw(modYCoor);
			modYCoor -= transGrpThickness;
		}
	}
	public void setTransGrpGlyphs(ArrayList transGrpGlyphs) {
		this.transGrpGlyphs = transGrpGlyphs;
	}
	public double getNonOverlapYCoor() {
		return modYCoor;
	}
	public Rectangle2D.Double getBoundingBox() {
		return boundingBox;
	}

	public int getEnd() {
		return end;
	}
	public int getStart() {
		return start;
	}
}
