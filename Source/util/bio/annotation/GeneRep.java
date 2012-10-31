package util.bio.annotation;

import java.util.*;
import java.io.*;

/**
 * A GeneRep resentation contains a genomic sequence designating a gene broken into various segments.
 This takes the most conservative approach.  Any annotation that says something is coding is shown as coding.
 Exons are maximized, introns and UTRs minimalized (actually UTR refers to RNA, not genomic).
 Different translations/transcriptions will give various versions of these segments.
 There is no orientation for these features! Use the GeneGroup to get the real start/stop and thus the orientation
 and strand.*/

public class GeneRep implements Serializable{
	//Fields
	private GeneGroup geneGrpRef; 
	private ArrayList codingSegments = new ArrayList();
	private ArrayList nonCodingSegments = new ArrayList();
	private int startATGPosition = -1;
	private int[] fivePrimeNonCodingRegion;
	private int[] threePrimeNonCodingRegion;
	private boolean translationFlag = true; //true if it exists, false otherwise; some features ie RNAs/psuedo genes don't translate
	
	//Constructor
	public GeneRep(TransGroup[] transGrps, GeneGroup geneGrpRef) {
		this.geneGrpRef = geneGrpRef;
		//get all coding and non coding segs
		ArrayList allCodingSegs = new ArrayList();
		int smallestStart = 1000000000;
		int biggestEnd = 0;
		int len = transGrps.length;
		for (int i = 0; i < len; i++) {
			//coding
			Translation t = transGrps[i].getTranslation();
			if (t == null) {
				translationFlag = false;
				break;
			}
			allCodingSegs.addAll(t.getCodingSegments());
			//lengths
			Transcript transcript = transGrps[i].getTranscript();
			if (transcript.getStart()<smallestStart) smallestStart = transcript.getStart();
			if (transcript.getEnd()>biggestEnd) biggestEnd = transcript.getEnd();
		}
		//join segments to get coding Segs
		if (translationFlag){
			codingSegments = GeneGroup.joinSegments(allCodingSegs);
			codingSegments = GeneGroup.sortSegments(codingSegments);
			//make 5'UTR, 3'UTR, and noncoding introns
			nonCodingSegments = GeneGroup.fetchIntrons(codingSegments);
			//fetch start of 1st exon
			int start1stCodingSeg = ((int[])codingSegments.get(0))[0];
			int endLastCodingSeg = ((int[])codingSegments.get(codingSegments.size()-1))[1];	
			int orientation = geneGrpRef.getOrientation();
			if (orientation == 1) {
				if (smallestStart == start1stCodingSeg) fivePrimeNonCodingRegion = new int[]{smallestStart, start1stCodingSeg};
				else fivePrimeNonCodingRegion = new int[]{smallestStart, start1stCodingSeg-1};
				if (endLastCodingSeg == biggestEnd) threePrimeNonCodingRegion = new int[]{endLastCodingSeg, biggestEnd};
				else threePrimeNonCodingRegion = new int[]{endLastCodingSeg+1, biggestEnd};
			}
			else if (orientation == -1) {
				if (smallestStart == start1stCodingSeg)threePrimeNonCodingRegion = new int[]{smallestStart, start1stCodingSeg};
				else threePrimeNonCodingRegion = new int[]{smallestStart, start1stCodingSeg-1};
				if (endLastCodingSeg == biggestEnd)fivePrimeNonCodingRegion = new int[]{endLastCodingSeg, biggestEnd};
				else fivePrimeNonCodingRegion = new int[]{endLastCodingSeg+1, biggestEnd};
			}
		}
	}
	
	//primary methods
	public String toString() {
		StringBuffer sb = new StringBuffer("");
		if (translationFlag == true) {
			sb.append("\nCoding Segment(s):");
			sb.append(GeneGroup.ALToString(codingSegments));
			sb.append("\nNon Coding Intronic Segment(s):");
			sb.append(GeneGroup.ALToString(nonCodingSegments));
			sb.append("\n5' Non Coding GenomicRegion: ");
			sb.append(fivePrimeNonCodingRegion[0]);
			sb.append("-");
			sb.append(fivePrimeNonCodingRegion[1]);
			sb.append("\n3' Non Coding GenomicRegion: ");
			sb.append(threePrimeNonCodingRegion[0]);
			sb.append("-");
			sb.append(threePrimeNonCodingRegion[1]);
			sb.append("\n");
		}
		return sb.toString();
	}
	//getters
	public GeneGroup getGeneGrp() {
		return geneGrpRef;
	}
	public boolean isTranslationFlag() {
		return translationFlag;
	}
	/**Conservative estimate of coding segments. Returns biggest from all translations.*/
	public ArrayList getCodingSegments() {
		return codingSegments;
	}
	/**Conservative estimate of non coding segments between coding exons. Not introns!*/
	public ArrayList getNonCodingSegments() {
		return nonCodingSegments;
	}
	/**Conservative estimate returning smallest non coding region. Not the UTR!*/
	public int[] getFivePrimeNonCodingRegion() {
		return fivePrimeNonCodingRegion;
	}
	/**Conservative estimate returning smallest non coding region. Not the UTR!*/
	public int[] getThreePrimeNonCodingRegion() {
		return threePrimeNonCodingRegion;
	}
	/**Conservative estimate of the start ATG position, from the longest transcript*/
	public int getStartATGPosition() {
		if (startATGPosition !=-1) return startATGPosition;
		if (codingSegments.size() == 0){
			if (geneGrpRef.getOrientation() == 1) startATGPosition = geneGrpRef.getStart();
			else startATGPosition = geneGrpRef.getEnd();
		}
		else{
			if (geneGrpRef.getOrientation() == 1) startATGPosition = ((int[])codingSegments.get(0))[0];
			else startATGPosition = ((int[])codingSegments.get(codingSegments.size()-1))[1];
		}
		return startATGPosition;
	}
}
