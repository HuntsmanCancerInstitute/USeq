package util.bio.wrappers;

/**
 * BL2Seq result container.
 *
 */
public class BL2SeqHit {
	
	//fields 
	private int rawScore;
	private int startSeq1;
	private int stopSeq1;
	private String seq1;
	private int startSeq2;
	private int stopSeq2;
	private String seq2;
	private int orientation; //0 for +/+, 1 for +/-
	
	public int getOrientation() {
		return orientation;
	}
	public void setOrientation(int orientation) {
		this.orientation = orientation;
	}
	public int getRawScore() {
		return rawScore;
	}
	public void setRawScore(int rawScore) {
		this.rawScore = rawScore;
	}
	public String getSeq1() {
		return seq1;
	}
	public void setSeq1(String seq1) {
		this.seq1 = seq1;
	}
	public String getSeq2() {
		return seq2;
	}
	public void setSeq2(String seq2) {
		this.seq2 = seq2;
	}
	public int getStartSeq1() {
		return startSeq1;
	}
	public void setStartSeq1(int startSeq1) {
		this.startSeq1 = startSeq1;
	}
	public int getStartSeq2() {
		return startSeq2;
	}
	public void setStartSeq2(int startSeq2) {
		this.startSeq2 = startSeq2;
	}
	public int getStopSeq1() {
		return stopSeq1;
	}
	public void setStopSeq1(int stopSeq1) {
		this.stopSeq1 = stopSeq1;
	}
	public int getStopSeq2() {
		return stopSeq2;
	}
	public void setStopSeq2(int stopSeq2) {
		this.stopSeq2 = stopSeq2;
	}
}
