package util.bio.annotation;

import java.util.*;

import util.bio.parsers.gff.*;
import util.gen.*;

/**
 * Class to hold user defined gff items like CRMs, enhancers, unpublished single track notation.
 */
public class GenericFeature implements Comparable {
	//fields
	private int orientation;
	private int trackNumber;
	private String featureType;
	private int start;
	private int end;
	private double score;
	private String strand;
	private String attributes;

	
	/**using to sort the GenericFeatures onto different tracks
	attempts to assign from the key value of the attributes section ie text=CG666;
	failing that it assigns from the SeqName gff descriptor (often this is used to
	hold the chromosome text, not the feature text*/
	private String name;

	public GenericFeature(GffFeature gffFeature) {
		//this.gffFeature = gffFeature;
		//orientation
		strand = gffFeature.getStrand();
		orientation = GeneGroup.convertPlusToNumOrientation(strand);	
		featureType = gffFeature.getFeature();
		score = gffFeature.getScore();
		start = gffFeature.getStart();
		end = gffFeature.getEnd();
		attributes = gffFeature.getAttsString();
		
		//text
		HashMap hash = gffFeature.getAttsHash();
		name = (String) hash.get("text");
		if (Misc.isEmpty(name))
			name = gffFeature.getSeqName();
	}
	
	public GenericFeature(Gff3Feature f){
		orientation = GeneGroup.convertPlusToNumOrientation(f.getStrand());
		featureType = f.getType();
		score = f.getScore();
		start = f.getStart();
		end = f.getEnd();
		strand = f.getStrand();
		attributes = f.getAttributes();
		name = f.getId();
		if (Misc.isEmpty(name)) name = f.getSeqId();
	}
	
	
	public int compareTo(Object other) {
		GenericFeature otherGF = (GenericFeature) other;
		if (trackNumber>otherGF.trackNumber)
			return 1;
		if (trackNumber<otherGF.trackNumber)
			return -1;
		return 0;
	}
	public String getFeatureType() {
		return featureType;
	}
	public int getOrientation() {
		return orientation;
	}
	public int getTrackNumber() {
		return trackNumber;
	}
	public void setTrackNumber(int i) {
		trackNumber = i;
	}
	public String getName(){
		return name;
	}
	public double getScore() {
		return score;
	}
	public String getAttributes() {
		return attributes;
	}
	public String getStrand() {
		return strand;
	}
	public int getEnd() {
		return end;
	}
	public int getStart() {
		return start;
	}
}
