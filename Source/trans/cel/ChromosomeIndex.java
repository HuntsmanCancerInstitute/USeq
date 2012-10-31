package trans.cel;

/**Use to keep track of intensities in a array, their chromosome text, and start stop indexes.*/
public class ChromosomeIndex {

	private String chromosome;
	private int startIndex;
	private int stopIndex;
	
	public ChromosomeIndex(String chromosome, int startIndex, int stopIndex){
		this.chromosome = chromosome;
		this.startIndex = startIndex;
		this.stopIndex = stopIndex;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int getStartIndex() {
		return startIndex;
	}

	public int getStopIndex() {
		return stopIndex;
	}
	
}
