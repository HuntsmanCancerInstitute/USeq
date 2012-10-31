package edu.utah.seq.methylation454;

/**Sequence aligned to a reference.*/
public class Read {
	
	//fields
	private String sequence;
	private String strand;
	private int alignmentStart;
	private int alignmentEnd;
	
	//constructor
	/** example of a header and a seqLine
	 * >ETK72AO01CGP6K /amp="I14" /cId="CON_1" /id="ETK72AO01CGP6K" /qCd="24-45" /lCp="5" /rCd="24-45" /rRg="24-56" /mFc="60" /str="1"
	 * 23	AA-T-AA-T-AAA-AA-TTTT-A-G-A-TA-TG	1452
	 * */
	public Read (String header, String seqLine){
		//parse from header the strand
		if (header.indexOf("str=\"1") != -1) strand = "+";
		else strand = "-";
		//pull start of aligment and sequence
		String[] tokens = seqLine.split("\\s+");
		alignmentStart = Integer.parseInt(tokens[0]);
		sequence = tokens[1];
		alignmentEnd = alignmentStart + sequence.length();
	}

	public int getAlignmentEnd() {
		return alignmentEnd;
	}

	public void setAlignmentEnd(int alignmentEnd) {
		this.alignmentEnd = alignmentEnd;
	}

	public int getAlignmentStart() {
		return alignmentStart;
	}

	public void setAlignmentStart(int alignmentStart) {
		this.alignmentStart = alignmentStart;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public String getStrand() {
		return strand;
	}

	public void setStrand(String strand) {
		this.strand = strand;
	}
	
}
