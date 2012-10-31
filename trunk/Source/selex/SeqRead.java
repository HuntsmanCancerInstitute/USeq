package selex;
import java.util.*;
import java.util.regex.*;

import util.gen.*;

/**Container for information related to a particular sequencing read.*/
public class SeqRead {
	
	//fields
	private String seq;
	private String seqCap;
	private String qualCap;
	private SelexParams sp;
	private SeqFiles seqFile;
	private int[] qual;
	private int cutOff;
	private int lenMin;
	private int lenMax;
	private int expectedSize;
	private String cutSite;
	private SubSeq[] subSeqs;
	private int numSubs = 0;
	private int numSubsPassScore = 0;
	private int numSubsPassScoreLength = 0;
	private int numExpectedSubSeqs;
	
	//constructor
	public SeqRead(
			String sequence,
			String seqComment,
			String quality,
			String qualComment,
			SelexParams aSP,
			SeqFiles aSF) {
		
		//assign vars
		seq = sequence.toUpperCase();
		seqCap = seqComment;
		qualCap = qualComment;
		sp = aSP;
		seqFile = aSF;
		
		//get and assign sp values
		cutOff = sp.getScoreCutOff();
		lenMin = sp.getLenMin();
		lenMax = sp.getLenMax();
		expectedSize = sp.getExpectedSize();
		cutSite = sp.getRestSite().toUpperCase();
		
		//process quality text and make qual array
		String[] qs = quality.toUpperCase().split("\\s");
		int len = qs.length;
		//check that seq and qual are the same size
		if (len != seq.length()) {
			System.out.println("");
			System.out.println(
			"\n\nFatal Problem! The number of quality scores and sequence letters, in this sequence read, are different!\nDid you reverse the inputs? -s sequence file -q quality file...\n\n");
			System.out.println("File: " + seqFile.getSeqFile());
			System.out.println("Seq Comment: >" + seqCap);
			System.out.println("Qual Comment: >" + qualCap);
			System.out.println("");
			System.exit(0);
		}
		qual = new int[len];
		for (int i = 0; i < len; i++) {
			qual[i] = Integer.parseInt(qs[i]);
		}
		//make an array of SubSeqs and print report
		sp.printSave("***************************************************************************\n");
		sp.printSave("***************************Processing New SeqRead**************************\n");
		boolean subsMade = makeSubSeqs();
		if (subsMade) checkEstimateNumSubSeqs();
		printSeqRead();
		sp.incNumSubs(numSubs, numSubsPassScore, numSubsPassScoreLength);
	}
	public boolean checkEstimateNumSubSeqs(){
		double expected = (double)numExpectedSubSeqs;
		if (numExpectedSubSeqs<1) expected = 0.000000001;
		double ratio = (double)numSubsPassScore/expected;
		System.out.println("\n\n******* short ratio "+ratio+"\n\n");			
		if (ratio < 0.70) {
			sp.incrementNumShortInsertSeqReads();
			return false;
		}
		return true;	 
	}
	
	public void printSeqRead() {
		sp.printSave(
				"\n**********SeqRead Report**********\n"
				+ "Seq Comment: >"
				+ seqCap
				+ "\nSequence: "
				+ seq
				+ "\nQual Comment: >"
				+ qualCap
				+ "\nQuality Scores: "
				+ Misc.intArrayToString(qual," ")
				+ "\n# of sub seqs found in this sequence read: "
				+ numSubs
				+ "\n# of sub seqs with sufficient quality scores: "
				+ numSubsPassScore
				+ "\n# of sub seqs also meeting length boundaries: "
				+ numSubsPassScoreLength
				+"\n# of expected sub seqs: "+
				+ numExpectedSubSeqs
				+ "\n\n\n");
	}
	
	public boolean makeSubSeqs() {
		//process SeqRead
		//check quality
		boolean qualCheck;
		sp.printSave("Checking overall quality of sequence... \n");
		if (qualCheck = checkSeqQuality(qual, 100, cutOff)) {
			sp.printSave("   Looks OK\n");
		} else {
			sp.printSave("   Looks poor, mostlikely a failed reaction.\n");
			sp.incrementNumBadSeqReads();
		}
		//estimate number of expected sub sequences
		sp.printSave("Estimating the number of sub sequences...\n");
		int size = expectedSize+cutSite.length();
		numExpectedSubSeqs = getNumStretchesGoodQuality(qual, cutOff, size);
		if (numExpectedSubSeqs>1) numExpectedSubSeqs -= 2; //minus two for ragged ends
		//save expected for average expects statistic
		if (numExpectedSubSeqs>0)sp.addExpectedSubSeq(numExpectedSubSeqs);
		//look for cut sites to parse off of
		sp.printSave("Looking for the cut sites....\n");
		ArrayList ends = findEnds();
		if (ends.size() < 2) {
			sp.printSave(
			"Too few cut sites! Skipping this sequence read.\n");
			if (qualCheck == true)
				sp.incrementNumNoInsertSeqReads();
			return false;
		}
		//extract sub sequences
		sp.printSave("Extracting raw sub sequences....\n");
		int[][] extractedSegs = extractSegs(ends);
		numSubs = extractedSegs.length;
		//clear low quality subs
		sp.printSave("Removing low quality scoring sub sequences....\n");
		ArrayList segs = removeLowQualitySegs(extractedSegs);
		numSubsPassScore = segs.size() / 2; //use to create histagram
		if (numSubsPassScore == 0) {
			sp.printSave(
			"No sub sequences were found with sufficient quality! Skipping this sequence read.\n");
			return false;
		}	
		//clear odd sized subs
		sp.printSave("Removing odd sized sub sequences....\n");
		ArrayList goodSegs = removeOddSizedSegs(segs);
		numSubsPassScoreLength = goodSegs.size() / 2;
		//save for average statistic
		if (numExpectedSubSeqs>0){
			sp.addObsExpRatio((double)numSubsPassScoreLength/numExpectedSubSeqs);
			sp.addObservedSubSeq(numSubsPassScoreLength);
		}
		if (numSubsPassScoreLength == 0) {
			sp.printSave(
			"No sub sequences were found with correct lengths! Skipping this sequence read.\n");
			return false;
		}
		//make SubSeq objects
		sp.printSave("Making sub sequence objects....\n\n");
		int lenSegs = goodSegs.size();
		subSeqs = new SubSeq[lenSegs / 2];
		int j = 0;
		int lenCutSite = sp.getRestSite().length();
		for (int i = 0; i < lenSegs; i += 2) {
			int start = ((Integer) goodSegs.get(i)).intValue() + lenCutSite;
			int stop = ((Integer) goodSegs.get(i + 1)).intValue() - lenCutSite;
			String subSeq = seq.substring(start, stop + 1);
			sp.appendSub(subSeq);
			//get array of qScores
			int lenFrag = 1 + stop - start;
			int[] qScores = new int[lenFrag];
			for (int k = 0; k < lenFrag; k++) {
				qScores[k] = qual[k + start];
			}
			subSeqs[j] = new SubSeq(subSeq, qScores, sp, this);
			j++;
		}
		return true;
	}
	
	//primary methods
	public ArrayList findEnds() {
		//returns an array of Integers where the int refers to the index of the last base in
		//  the pattern ie "ccGATCcccc" looking for GATC would return 5
		ArrayList ends = new ArrayList();
		Pattern pat = Pattern.compile(cutSite);
		Matcher x = pat.matcher(seq);
		while (x.find()) {
			ends.add(new Integer(x.end() - 1));
		}
		return ends;
	}
	
	public int[][] extractSegs(ArrayList ends) {
		//returns an array of start/stop indexes demarking the REsite-oligo-REsite segment
		int len = ends.size();
		int subSeqs = len - 1;
		int[][] segs = new int[subSeqs][2];
		int lenCutSite = cutSite.length() - 1;
		for (int i = 0; i < subSeqs; i++) {
			int start = ((Integer) ends.get(i)).intValue() - lenCutSite;
			int stop = ((Integer) ends.get(i + 1)).intValue();
			segs[i][0] = start;
			segs[i][1] = stop;
			sp.printSave(
					"  raw sub seq: " + seq.substring(start, stop + 1) + "\n");
		}
		return segs;
	}
	
	public ArrayList removeLowQualitySegs(int[][] segs) {
		//scans all the quality scores for each base between start and stop
		//returns those >= the cut off score
		int len = segs.length;
		ArrayList goodSegs = new ArrayList(len * 2);
		for (int i = 0; i < len; i++) {
			//System.out.println("*************");
			int start = segs[i][0];
			int stop = segs[i][1];
			//scan qual array for low scoring scores
			boolean good = true;
			for (int j = start; j < stop + 1; j++) {
				//System.out.println("S: "+seq.charAt(j)+ " Q: "+qual[j]);
				if (qual[j] < cutOff) {
					good = false;
					//System.out.println("Bad, breaking");
					break;
				}
			}
			if (good) {
				goodSegs.add(new Integer(start));
				goodSegs.add(new Integer(stop));
			}
		}
		goodSegs.trimToSize();
		return goodSegs;
	}
	
	public ArrayList removeOddSizedSegs(ArrayList segs) {
		//removes segs less than or greater than size limits
		int len = segs.size();
		ArrayList goodSegs = new ArrayList(len);
		int sizer = cutSite.length() * 2 - 1;
		for (int i = 0; i < len; i += 2) {
			int start = ((Integer) segs.get(i)).intValue();
			int stop = ((Integer) segs.get(i + 1)).intValue();
			int size = stop - start - sizer;
			sp.appendLength(size);
			//record the size of all good quality subSequences
			
			if (size >= lenMin && size <= lenMax) {
				goodSegs.add(segs.get(i));
				goodSegs.add(segs.get(i + 1));
			}
		}
		goodSegs.trimToSize();
		return goodSegs;
	}
	/**Examines an int[] for a continuous stretch of scores*/
	public static boolean checkSeqQuality(
			int[] qualityScores,
			int minNumHighQualBases,
			int minQualScore) {
		int len = qualityScores.length;
		int counter = 0;
		boolean flag = false;
		for (int i = 0; i < len; i++) {
			if (counter == minNumHighQualBases) {
				flag = true;
				break;
			}
			if (qualityScores[i] >= minQualScore)
				counter++;
			else
				counter = 0;
		}
		return flag;
	}
	/**Counts the number of stretches of scores above a cut off in an int[].  Useful for estimating the 
	 * number of concatinated selex oligos that should be found in a sequence read*/
	public static int getNumStretchesGoodQuality(
			int[] qualityScores, int qualCutOff, int seqLength) {
		int len = qualityScores.length;
		int counter = 1;
		int stretchCounter = 0;
		for (int i = 0; i < len; i++) {
			if (counter == seqLength) {
				stretchCounter++;
				counter = 1;
			}
			if (qualityScores[i] >= qualCutOff)
				counter++;
			else
				counter = 1;
		}
		return stretchCounter;
	}
}
