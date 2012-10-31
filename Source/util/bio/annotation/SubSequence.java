package util.bio.annotation;

public class SubSequence {
	private String sequence;
	private int start;
	private int end;
	
	public SubSequence (String sequence, int start, int end){
		this.sequence = sequence;
		this.start = start;
		this.end = end;
	}
	
	public String getName(){
		return start+"-"+end;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}
}
