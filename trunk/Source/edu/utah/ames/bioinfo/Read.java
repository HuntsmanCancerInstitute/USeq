package edu.utah.ames.bioinfo;

public class Read {

	//fields from input file
	private String displayName;
	private String geneName;
	private String chr;
	private String strand;
	private String start;
	private String stop;
	private String totalBPs;
	private String pValue;
	private String FDR;
	private String varCorLg2R;
	private String spliceChiPVal;
	private String spliceMaxLg2R;
	private String spliceMaxExon;
	
	//index fields
	public static final int displayNameIndex = 0;
	public static final int geneNameIndex = 1;
	public static final int chrIndex = 2;
	public static final int strandIndex = 3;
	public static final int startIndex = 4;
	public static final int stopIndex = 5;
	public static final int totalBPsIndex = 6;
	public static final int pValueIndex = 7;
	public static final int FDRIndex = 8;
	public static final int varCorLg2RIndex = 9;
	public static final int spliceChiPValIndex = 10;
	public static final int spliceMaxLg2RIndex = 11;
	public static final int spliceMaxExonIndex = 12;

	//constructor
	public Read(String[] dataValue) {
		
		//assign variables to their indexed positions
		displayName = dataValue[displayNameIndex];
		geneName = dataValue[geneNameIndex];
		chr = dataValue[chrIndex];
		strand = dataValue[strandIndex];
		start = dataValue[startIndex];
		stop = dataValue[stopIndex];
		totalBPs = dataValue[totalBPsIndex];
		pValue = dataValue[pValueIndex];
		FDR = dataValue[FDRIndex];
		varCorLg2R = dataValue[varCorLg2RIndex];
		spliceChiPVal = dataValue[spliceChiPValIndex];
		spliceMaxLg2R = dataValue[spliceMaxLg2RIndex];
		spliceMaxExon = dataValue[spliceMaxExonIndex];
	}
	
	public String getGeneName() {
		return geneName;
	}

	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}

	public String getDisplayName() {
		return displayName;
	}

	public void setDisplayName(String displayName) {
		this.displayName = displayName;
	}

	public String getChr() {
		return chr;
	}

	public String getStrand() {
		return strand;
	}

	public String getStart() {
		return start;
	}

	public String getStop() {
		return stop;
	}

	public String getTotalBPs() {
		return totalBPs;
	}

	public String getpValue() {
		return pValue;
	}

	public String getFDR() {
		return FDR;
	}

	public String getVarCorLg2R() {
		return varCorLg2R;
	}

	public String getSpliceChiPVal() {
		return spliceChiPVal;
	}

	public String getSpliceMaxLg2R() {
		return spliceMaxLg2R;
	}

	public String getSpliceMaxExon() {
		return spliceMaxExon;
	}

	public void setChr(String chr) {
		this.chr = chr;
	}

	public void setStrand(String strand) {
		this.strand = strand;
	}

	public void setStart(String start) {
		this.start = start;
	}

	public void setStop(String stop) {
		this.stop = stop;
	}

	public void setTotalBPs(String totalBPs) {
		this.totalBPs = totalBPs;
	}

	public void setpValue(String pValue) {
		this.pValue = pValue;
	}

	public void setFDR(String fDR) {
		FDR = fDR;
	}

	public void setVarCorLg2R(String varCorLg2R) {
		this.varCorLg2R = varCorLg2R;
	}

	public void setSpliceChiPVal(String spliceChiPVal) {
		this.spliceChiPVal = spliceChiPVal;
	}

	public void setSpliceMaxLg2R(String spliceMaxLg2R) {
		this.spliceMaxLg2R = spliceMaxLg2R;
	}

	public void setSpliceMaxExon(String spliceMaxExon) {
		this.spliceMaxExon = spliceMaxExon;
	}
	
	
}
