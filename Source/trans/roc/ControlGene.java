package trans.roc;
import java.io.File;
import util.gen.*;
import java.util.*;

public class ControlGene{
	
	//fields
	private String name;
	private String chromosome;
	private String strand;
	private int[] exonStarts;
	private int[] exonEnds;
	private int[][] exons;
	private String notes;
	private float[] intensities;
	
	//constructor 
	public ControlGene(String line){
		//#text	chrom	strand	exonCount	exonStarts	exonEnds   notes
		String[] tokens = line.split("\\t");
		if (tokens.length != 6) System.out.println("Problem parsing -> "+line);
		else{
			name = tokens[0];
			chromosome = tokens[1];
			strand = tokens[2];
			exonStarts = Num.parseInts(tokens[3].split(","));
			exonEnds = Num.parseInts(tokens[4].split(","));
			notes = tokens[5];
		}
	}
	
	/**Looks to see if point is contained in any of the exons, ends are inclusive.*/
	public boolean containsPoint (int bp){
		for (int i=0; i<exons.length; i++){
			if (bp >= exons[i][0] && bp <= exons[i][1]) return true;
		}
		return false;
	}
	
	/**Subtracts the endOffSet from the stop and only creates an exon if 1+ stop-start > 0*/
	public void makeExons(int endOffSet){
		if (exonStarts.length != exonEnds.length) System.out.print("Number of starts and ends differ!");
		else {
			ArrayList al = new ArrayList(exonStarts.length);
			for (int i=0; i< exonStarts.length; i++){
				int[] test = {exonStarts[i], exonEnds[i]-endOffSet};
				if ((1+test[1]- test[0]) > 0) al.add(test);
			}
			exons = Num.arrayList2IntArrayArray(al);
			//Misc.printArray(exons);
		}
	}
	
	/**Fires makeExons on each CG within the ControlGene[].*/
	public static void makeExons(ControlGene[] cgs, int endOffSet){
		for (int i=0; i< cgs.length; i++){
			cgs[i].makeExons(endOffSet);
		}
	}
	
	/**parses a control gene file
	 * ex: #text	chrom	strand	exonStarts	exonEnds	notes
	 * ex: NM_001035267	chr12	+	54796640,54797239,54797532,	54796850,54797262,54797882,	Z12962_Homologue to yeast ribosomal protein L41_H
	 * */
	public static ControlGene[] parseControlGeneFile( File file){
		String[] lines = IO.loadFileIntoStringArray(file);
		ArrayList al = new ArrayList(lines.length);
		for (int i=0; i<lines.length; i++){
			if (lines[i].startsWith("#") == false) al.add(new ControlGene(lines[i]));
		}
		ControlGene[] cg = new ControlGene[al.size()];
		al.toArray(cg);
		return cg;
	}

	public float[] getIntensities() {
		return intensities;
	}

	public void setIntensities(float[] intensities) {
		this.intensities = intensities;
	}

	public String getNotes() {
		return notes;
	}

	public void setNotes(String notes) {
		this.notes = notes;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int[] getExonEnds() {
		return exonEnds;
	}

	public int[] getExonStarts() {
		return exonStarts;
	}

	public String getName() {
		return name;
	}

	public String getStrand() {
		return strand;
	}
}