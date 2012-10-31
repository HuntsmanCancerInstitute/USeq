package gata.geneGlyphs;


import gata.main.*;

import java.util.*;

import util.bio.annotation.*;
/**
 * @author Nix
 * Container class holding components that make up a TransGroup (shapes for DNA, exons, prot....)
 */
public class TransGrpGlyph {
	//fields
	private GATAParams gataParams;
	AnnoSpecParams annoParams;
	private TransGroup transGrp;
	private ExonIntron[] exons;
	private Transcript transcript;
	private int RNAStart;
	private int RNAEnd;
	private int orientation;
	private Translation prot; //might be null
	private double ntAtXZero;
	private double pixPerNt;
	private double yTop = 0; //set to 0 not null

	//ArrayLists
	private ArrayList thinLinesProt = new ArrayList();
	private ArrayList dottedLinesProt = new ArrayList();
	private ArrayList thickLinesRNA = new ArrayList();

	public TransGrpGlyph(TransGroup transGrp, GATAParams gataParams, AnnoSpecParams annoParams) {
		this.transGrp = transGrp;
		this.gataParams = gataParams;
		this.annoParams = annoParams;
		exons = transGrp.getExons();
		transcript = transGrp.getTranscript();
		RNAStart = transcript.getStart();
		RNAEnd = transcript.getEnd();
		orientation = transcript.getOrientation();
		prot = transGrp.getTranslation();	
	}
	//primary methods
	public void draw(double yCoord) {
		
		//set constants for this particular draw
		pixPerNt = gataParams.getPixPerNt();
		ntAtXZero = annoParams.getNtAtXZero();
		double interval = 0;
		//RNA
		if (gataParams.isRNAVis()) {
			interval += gataParams.getPixBtwProtRNA()+(gataParams.getFatLineThicknes() / 2);
			thickLinesRNA =
			GATAUtil.makeSegmentsInt(
					transGrp.getExtractedExons(),
					yCoord,
					pixPerNt,
					ntAtXZero);
			//add thin base line for RNA (Line2D.double) and a List of thicklines
			annoParams.addRNALines(
			GATAUtil.makeLine(RNAStart, RNAEnd, yCoord, pixPerNt, ntAtXZero),
				thickLinesRNA);
		}
		yTop = yCoord - interval;

		//Protein
		if (prot != null && gataParams.isProteinVis()) {
			dottedLinesProt =
			GATAUtil.makeSegments(
					prot.getIntrons(),
					yTop,
					pixPerNt,
			ntAtXZero);
			thinLinesProt =
			GATAUtil.makeSegments(
					prot.getCodingSegments(),
					yTop,
					pixPerNt,
			ntAtXZero);
			//add arrow and butt lines
			double arrowHeight = gataParams.getPixHeightArrows();
			thinLinesProt.addAll(
			GATAUtil.makeArrowAndButt(
					prot.getStart(),
					prot.getEnd(),
					-1,
					yTop,
					orientation,
					pixPerNt,
			ntAtXZero,
			arrowHeight));
			//set yTop
			yTop -= arrowHeight;
			//add lines
			annoParams.addProteinLines(thinLinesProt, dottedLinesProt);
		}
	}
	/**Top most y coordinate for this transGroup, including transGrp spacer, start
		 * next Shape here for no overlap*/
	public double getYTop() {
		return yTop;
	}
}
