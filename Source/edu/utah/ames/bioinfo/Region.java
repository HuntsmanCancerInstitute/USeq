package edu.utah.ames.bioinfo;

/**
 * 
 * @author darren.ames@hci.utah.edu
 *
 */

public class Region {
	
	//fields from input bed file report
	private String chromNum;
	private String start;
	private String stop;
	
	//index fields
	public static final int chromNumIndex = 0;
	public static final int startIndex = 1;
	public static final int stopIndex = 2; 
	
	//constructor
	public Region(String[] dataValue) {
		chromNum = dataValue[chromNumIndex];
		start = dataValue[startIndex];
		stop = dataValue[stopIndex];
	}

	public String getChromNum() {
		return chromNum;
	}

	public void setChromNum(String chromNum) {
		this.chromNum = chromNum;
	}

	public String getStart() {
		return String.valueOf(start);
	}

	public void setStart(String start) {
		this.start = start;
	}

	public String getStop() {
		return String.valueOf(stop);
	}

	public void setStop(String stop) {
		this.stop = stop;
	}
}
