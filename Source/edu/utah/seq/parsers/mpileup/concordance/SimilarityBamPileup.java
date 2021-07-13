package edu.utah.seq.parsers.mpileup.concordance;

import java.io.IOException;

import edu.utah.seq.parsers.mpileup.concordance.ConcordanceChunkBamPileup.ParsedSample;
import util.gen.Histogram;
import util.gen.IO;
import util.gen.Num;

public 	class SimilarityBamPileup implements Comparable<SimilarityBamPileup>{
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
	private double fracMatchAB = -1;
	private double fracMatchBA = -1;
	private double minSim = 0;
	private String toIgnoreForSim;
	
	public static final char[] indexBases = new char[] {'G', 'A', 'T', 'C'};
	
	public void add(SimilarityBamPileup other) throws Exception {
		this.numMatchesAB += other.numMatchesAB;
		this.numMisMatchesAB += other.numMisMatchesAB;
		this.numMatchesBA += other.numMatchesBA;
		this.numMisMatchesBA += other.numMisMatchesBA;
		this.sampleAMisMatchAFs.addCounts(other.sampleAMisMatchAFs);
		this.sampleBMisMatchAFs.addCounts(other.sampleBMisMatchAFs);
	}
	
	public SimilarityBamPileup(int a, int b, SampleConcordance bc){ //double minAFForHom, double minAFForMatch, double minSim, String toIgnoreForSim){
		sampleA = a;
		sampleB = b;
		minAFForHom = bc.getMinAFForHom();
		minAFForMatch = bc.getMinAFForMatch();
		minSim = bc.getMinFracSim();
		toIgnoreForSim = bc.getToIgnoreForCall();
	}
	
	public String passSimilarity(String[] sampleNames){
		if (fracMatchBA == -1) fracMatchBA = (double)numMatchesBA/(double)(numMatchesBA+ numMisMatchesBA);
		if (fracMatchAB == -1) fracMatchAB = (double)numMatchesAB/(double)(numMatchesAB+ numMisMatchesAB);
		boolean pass = true;
		if (fracMatchBA < minSim) pass = false;
		if (fracMatchAB < minSim) pass = false;
		if (pass == false){
			if (sampleNames[sampleA].contains(toIgnoreForSim) || sampleNames[sampleB].contains(toIgnoreForSim)) return "WARNING";
			return "FAIL";
		}
		return "PASS";
	}

	public String toString(String[] sampleNames){
		if (fracMatchAB == -1) fracMatchAB = (double)numMatchesAB/(double)(numMatchesAB+ numMisMatchesAB);
		String ab = sampleNames[sampleA]+"->"+ sampleNames[sampleB]+"\t"+ Num.formatNumberMinOne(fracMatchAB, 3)+" ("+(int)numMatchesAB+"/"+(int)(numMatchesAB+numMisMatchesAB)+")\n";
		if (fracMatchBA == -1) fracMatchBA = (double)numMatchesBA/(double)(numMatchesBA+ numMisMatchesBA);
		String ba = sampleNames[sampleB]+"->"+ sampleNames[sampleA]+"\t"+ Num.formatNumberMinOne(fracMatchBA, 3)+" ("+(int)numMatchesBA+"/"+(int)(numMatchesBA+numMisMatchesBA)+")\n";
		return ab+ba+passSimilarity(sampleNames)+"\n";
	}
	
	public String toStringShort(String[] sampleNames){
		double fracMatchAB = (double)numMatchesAB/(double)(numMatchesAB+ numMisMatchesAB);
		String ab = sampleNames[sampleA]+"->"+ sampleNames[sampleB]+" "+ Num.formatNumberMinOne(fracMatchAB, 3)+"("+(int)numMatchesAB+"/"+(int)(numMatchesAB+numMisMatchesAB)+")";
		double fracMatchBA = (double)numMatchesBA/(double)(numMatchesBA+ numMisMatchesBA);
		String ba = Num.formatNumberMinOne(fracMatchBA, 3)+"("+(int)numMatchesBA+"/"+(int)(numMatchesBA+numMisMatchesBA)+")";
		return ab+" "+ba+" "+ passSimilarity(sampleNames);
	}
	
	public void toStringMismatch(String[] sampleNames) {
		IO.pl(sampleNames[sampleA] +" <-> "+ sampleNames[sampleB]);
		IO.pl(sampleNames[sampleA] +" mismatches ");
		sampleAMisMatchAFs.printScaledHistogram();
		IO.pl(sampleNames[sampleB] +" mismatches ");
		sampleBMisMatchAFs.printScaledHistogram();
	}
	
	public void calculateMaxMatch(){
		if (fracMatchBA == -1) fracMatchBA = (double)numMatchesBA/(double)(numMatchesBA+ numMisMatchesBA);
		if (fracMatchAB == -1) fracMatchAB = (double)numMatchesAB/(double)(numMatchesAB+ numMisMatchesAB);
		if (fracMatchAB >= fracMatchBA) maxMatch = fracMatchAB;
		else maxMatch = fracMatchBA;
	}

	/**Returns whether a proper mismatch was observed.
	 * @param commonVarChecker 
	 * @throws IOException */
	public boolean contrast(ParsedSample[] ps) throws IOException{
		//check if both pass DP and INDEL
		if (ps[sampleA].passDpIndel == false || ps[sampleB].passDpIndel == false) return false;
		
		//is A homozygous?
		boolean misMatch = false;
		if (ps[sampleA].maxAFIndex[0] >= minAFForHom){
			int baseIndex = (int)ps[sampleA].maxAFIndex[1];
			double afB = ps[sampleB].ms.getSnvAlleleFreq(indexBases[baseIndex]);
			
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
			double afA = ps[sampleA].ms.getSnvAlleleFreq(indexBases[baseIndex]);
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

	public int compareTo(SimilarityBamPileup o) {
		if (this.maxMatch > o.maxMatch) return -1;
		if (this.maxMatch < o.maxMatch) return 1;
		return 0;
	}







}
