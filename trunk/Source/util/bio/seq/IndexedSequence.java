package util.bio.seq;
import java.io.*;
import java.util.regex.*;

public class IndexedSequence implements Comparable{
	
	//fields
	private String chromosome;
	private int start;
	private int stop;
	private String sequence;
	private File binarySeqFile;
	public static final Pattern nameParser = Pattern.compile("(.+)_(.+)-(.+)");
	

	public IndexedSequence (File binarySeqFile){
		this.binarySeqFile = binarySeqFile;
		Matcher mat = nameParser.matcher(binarySeqFile.getName());
		if (mat.matches() == false) System.out.println("\tWARNING: cannot parse binary sequence file info for "+binarySeqFile+", skipping!");
		chromosome = mat.group(1);
		start = Integer.parseInt(mat.group(2));
		stop = Integer.parseInt(mat.group(3));
	}

	/**Assumes starts are unique.*/
	public int compareTo(Object obj){
		IndexedSequence other = (IndexedSequence) obj;
		if (other.start < start) return 1;
		return -1;
	}
	
	public boolean contained (int basePosition){
		if (start <= basePosition && basePosition <= stop) return true;
		return false;
	}

	public String getChromosome() {
		return chromosome;
	}

	public String getSequence() {
		if (sequence == null) sequence = Seq.readBinarySequence(binarySeqFile);
		return sequence;
	}

	public int getStart() {
		return start;
	}

	public int getStop() {
		return stop;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

}
