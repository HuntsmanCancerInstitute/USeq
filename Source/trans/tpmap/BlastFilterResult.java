package trans.tpmap;

/**
 * Container for a blast result: chromosome, numMatches, exactMatches, oneBPMismatches, startExactMatch.
 */
public class BlastFilterResult {
	//fields
	private String chromosome;
	private int matches;
	private int exactMatches;
	private int oneBPMisMatches;
	private int startExactMatch;
	
	public BlastFilterResult(String chromosome, int matches, int exactMatches, int oneBPMismatches, int startExactMatch){
		this.chromosome = chromosome;
		this.matches= matches;
		this.exactMatches= exactMatches;
		this.oneBPMisMatches = oneBPMismatches;
		this.startExactMatch=startExactMatch;
	}

	public String getChromosome() {
		return chromosome;
	}
	public int getExactMatches() {
		return exactMatches;
	}
	public int getMatches() {
		return matches;
	}
	public int getOneBPMisMatches() {
		return oneBPMisMatches;
	}
	public int getStartExactMatch() {
		return startExactMatch;
	}
}
