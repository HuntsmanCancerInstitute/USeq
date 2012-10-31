package util.bio.wrappers;

/**
 * Information related to a blast results.
 *
 */
public class BlastResult {
	//fields
	private String matchId;
	private String matchName;
	private String matchLength;
	private String score;
	private String expect;
	private String identities;
	private String gaps;
	private String strand;
	private String[] alignment;
	
	
	//methods
	public String toString(){
		StringBuffer sb = new StringBuffer();
		sb.append("matchId "+matchId);
		sb.append("\n");
		sb.append("matchName "+matchName);
		sb.append("\n");
		sb.append("matchLength "+matchLength);
		sb.append("\n");
		sb.append("score "+score);
		sb.append("\n");
		sb.append("expect "+expect);
		sb.append("\n");
		sb.append("identities "+identities);
		sb.append("\n");
		sb.append("gaps "+gaps);
		sb.append("\n");
		sb.append("strand "+strand);
		sb.append("\n");
		int num = alignment.length;
		for (int i=0; i<num; i++){
			sb.append(alignment[i]);
			sb.append("\n");
		}
		return sb.toString();
	}
	//getters setters
	public String[] getAlignment() {
		return alignment;
	}
	public void setAlignment(String[] alignment) {
		this.alignment = alignment;
	}
	public String getExpect() {
		return expect;
	}
	public void setExpect(String expect) {
		this.expect = expect;
	}
	public String getGaps() {
		return gaps;
	}
	public void setGaps(String gaps) {
		this.gaps = gaps;
	}
	public String getIdentities() {
		return identities;
	}
	public void setIdentities(String identities) {
		this.identities = identities;
	}
	public String getMatchLength() {
		return matchLength;
	}
	public void setMatchLength(String matchLength) {
		this.matchLength = matchLength;
	}
	public String getMatchId() {
		return matchId;
	}
	public void setMatchId(String matchId) {
		this.matchId = matchId;
	}
	public String getScore() {
		return score;
	}
	public void setScore(String score) {
		this.score = score;
	}
	public String getStrand() {
		return strand;
	}
	public void setStrand(String strand) {
		this.strand = strand;
	}
	public String getMatchName() {
		return matchName;
	}
	public void setMatchName(String matchName) {
		this.matchName = matchName;
	}
}
