package edu.utah.seq.parsers.mpileup.concordance;

import edu.utah.seq.parsers.mpileup.concordance.ConcordanceChunk.ParsedSample;
import util.gen.Histogram;
import util.gen.IO;
import util.gen.Num;

public 	class Similarity implements Comparable<Similarity>{
	private double numMatchesAB = 0;
	private double numMisMatchesAB = 0;
	private double numMatchesBA = 0;
	private double numMisMatchesBA = 0;
	private double minAFForHom = 0;
	private double minAFForMatch = 0;
	private int sampleA;
	private int sampleB;
	private double maxMatch = 0;
	private Histogram sampleAMisMatchAFs = new Histogram(0,1,101);
	private Histogram sampleBMisMatchAFs = new Histogram(0,1,101);
	
	public void add(Similarity other) throws Exception {
		this.numMatchesAB += other.numMatchesAB;
		this.numMisMatchesAB += other.numMisMatchesAB;
		this.numMatchesBA += other.numMatchesBA;
		this.numMisMatchesBA += other.numMisMatchesBA;
		this.sampleAMisMatchAFs.addCounts(other.sampleAMisMatchAFs);
		this.sampleBMisMatchAFs.addCounts(other.sampleBMisMatchAFs);
	}
	
	Similarity(int a, int b, double minAFForHom, double minAFForMatch){
		sampleA = a;
		sampleB = b;
		this.minAFForHom = minAFForHom;
		this.minAFForMatch = minAFForMatch;
	}

	public String toString(String[] sampleNames){
		double fracMatchAB = (double)numMatchesAB/(double)(numMatchesAB+ numMisMatchesAB);
		String ab = sampleNames[sampleA]+"->"+ sampleNames[sampleB]+"\t"+ Num.formatNumberMinOne(fracMatchAB, 3)+" ("+(int)numMatchesAB+"/"+(int)(numMatchesAB+numMisMatchesAB)+")\n";
		double fracMatchBA = (double)numMatchesBA/(double)(numMatchesBA+ numMisMatchesBA);
		String ba = sampleNames[sampleB]+"->"+ sampleNames[sampleA]+"\t"+ Num.formatNumberMinOne(fracMatchBA, 3)+" ("+(int)numMatchesBA+"/"+(int)(numMatchesBA+numMisMatchesBA)+")\n";
		return ab+ba;
	}
	
	public void toStringMismatch(String[] sampleNames) {
		IO.pl(sampleNames[sampleA] +" <-> "+ sampleNames[sampleB]);
		IO.pl(sampleNames[sampleA] +" mismatches ");
		sampleAMisMatchAFs.printScaledHistogram();
		IO.pl(sampleNames[sampleB] +" mismatches ");
		sampleBMisMatchAFs.printScaledHistogram();
	}
	
	public void calculateMaxMatch(){
		double fracMatchAB = (double)numMatchesAB/(double)(numMatchesAB+ numMisMatchesAB);
		double fracMatchBA = (double)numMatchesBA/(double)(numMatchesBA+ numMisMatchesBA);
		if (fracMatchAB >= fracMatchBA) maxMatch = fracMatchAB;
		else maxMatch = fracMatchBA;
	}

	/**Returns whether a proper mismatch was observed.
	 * @param commonVarChecker */
	public boolean contrast(ParsedSample[] ps){
		//check if both pass DP and INDEL
		if (ps[sampleA].passDpIndel == false || ps[sampleB].passDpIndel == false) return false;
		
		//is A homozygous?
		boolean misMatch = false;
		if (ps[sampleA].maxAFIndex[0] >= minAFForHom){
			int baseIndex = (int)ps[sampleA].maxAFIndex[1];
			double afB = ps[sampleB].ms.getAlleleFreq(baseIndex);
			//increment similarity
			if (afB >= minAFForMatch) numMatchesAB++;
			else {
				numMisMatchesAB++;
				sampleBMisMatchAFs.count(afB);
				misMatch = true;
			}
		}

		//is B homozygous?
		if (ps[sampleB].maxAFIndex[0] >= minAFForHom){
			int baseIndex = (int)ps[sampleB].maxAFIndex[1];
			double afA = ps[sampleA].ms.getAlleleFreq(baseIndex);
			//increment similarity
			if (afA >= minAFForMatch) numMatchesBA++;
			else {
				numMisMatchesBA++;
				sampleAMisMatchAFs.count(afA);
				misMatch = true;
			}
		}	
		return misMatch;
	}

	public int compareTo(Similarity o) {
		if (this.maxMatch > o.maxMatch) return -1;
		if (this.maxMatch < o.maxMatch) return 1;
		return 0;
	}







}
