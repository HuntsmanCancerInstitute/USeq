package util.bio.seq.splice;
import java.util.*;

public class Transcript {
	ArrayList<Exon> exons = new ArrayList<Exon>();
	
	public Transcript(Exon firstExon){
		exons.add(firstExon);
	}

	public ArrayList<Exon> getExons() {
		return exons;
	}
}
