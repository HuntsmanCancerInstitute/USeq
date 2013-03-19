package edu.utah.seq.base;

import java.util.Arrays;
import java.util.List;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;
import util.gen.Num;

public class VariantAlignment{
	private int centerIndex;
	private int startIndex;
	private int endIndex;
	private SAMRecord alignment;
	private String leftSeq;
	private String centerSeq;
	private String rightSeq;
	private String indexStrand = null;
	private float[] weightedCorrelationsACGT;
	private double[] likelihoodACGTRatios;
	private String classifiedCenterSeq = null;
	private BaseClassifier baseClassifier;
	private boolean changedToN = false;
	private boolean changedToDiff = false;

	//constructor
	public VariantAlignment(BaseClassifier baseClassifier, SAMRecord alignment, int centerReferenceBasePosition, int nMerAdder){
		this.baseClassifier = baseClassifier;
		this.alignment = alignment;

		//find read index of variant, between 0 and 100 for 101bp
		centerIndex = findReadIndex(centerReferenceBasePosition,alignment);

		if (centerIndex == -1) {
			//System.err.println("\nCenterBase is not in this alignment? Deletion? Skipping\n"+alignment.getSAMString()+centerReferenceBasePosition);
			//findReadIndexDebug(centerReferenceBasePosition, alignment);
		}
		else {
			//define start stop
			startIndex = centerIndex - nMerAdder;
			if (startIndex < 0) startIndex = 0;
			endIndex = centerIndex + nMerAdder +1;
			if (endIndex > alignment.getReadLength()) endIndex = alignment.getReadLength();

			String seq = alignment.getReadString();
			//get real read sequence at center, must make this a new string since the seq can be changed
			centerSeq = new String(seq.substring(centerIndex, centerIndex+1));
			//watch out for N's
			if (centerSeq.equals("N")) return;
			leftSeq = alignment.getReadString().substring(startIndex, centerIndex);
			rightSeq = alignment.getReadString().substring(centerIndex+1, endIndex); 

			//get strand
			String strand = "+";
			if (alignment.getReadNegativeStrandFlag()) strand = "-";
			indexStrand = centerIndex + strand;
		}

	}
	public String getCenterIndexBase(){
		return alignment.getReadString().substring(centerIndex, centerIndex+1);
	}
	public byte getCenterIndexBaseQuality(){
		return alignment.getBaseQualities()[centerIndex];
	}

	/**Returns the actual Illumina called nmer sequence*/
	public String getVariantSeq(){
		return leftSeq+ centerSeq +rightSeq;
	}
	/**Returns a concat of the left and right sequence absent the middle base.*/
	public String getVariantSeqFlanks(){
		return leftSeq+rightSeq;
	}
	/**Returns variant sequences with ACGT at the center position.*/
	public String[] getVariantSeqs(){
		String[] seqs = new String[4];
		//A
		seqs[0] = leftSeq+"A"+rightSeq;
		//C
		seqs[1] = leftSeq+"C"+rightSeq;
		//G
		seqs[2] = leftSeq+"G"+rightSeq;
		//T
		seqs[3] = leftSeq+"T"+rightSeq;
		return seqs;
	}

	/**Returns a concatenate of the intensities for the ACGT channels.*/
	public double[] fetchACGTIntensities(){
		//fetch intensities
		short[] aIntensities = alignment.getSignedShortArrayAttribute("IA");
		short[] cIntensities = alignment.getSignedShortArrayAttribute("IC");
		short[] gIntensities = alignment.getSignedShortArrayAttribute("IG");
		short[] tIntensities = alignment.getSignedShortArrayAttribute("IT"); 
		//if negative strand must reverse the intensity arrays, these then represent the complement sequence 
		if (alignment.getReadNegativeStrandFlag()){
			aIntensities = Num.reverse(aIntensities);
			cIntensities = Num.reverse(cIntensities);
			gIntensities = Num.reverse(gIntensities);
			tIntensities = Num.reverse(tIntensities); 
		}
		//build vector
		double[] vals = new double[4 * (endIndex-startIndex)];
		//for each position in the nMer
		int counter = 0;
		//A
		for (int i=startIndex; i<endIndex; i++) vals[counter++] = aIntensities[i];
		//C
		for (int i=startIndex; i<endIndex; i++) vals[counter++] = cIntensities[i];
		//G
		for (int i=startIndex; i<endIndex; i++) vals[counter++] = gIntensities[i];
		//T
		for (int i=startIndex; i<endIndex; i++) vals[counter++] = tIntensities[i];

		return vals;
	}


	public static int findReadIndex(int centerBase, SAMRecord a) {
		List<AlignmentBlock> blocks = a.getAlignmentBlocks();
		for (AlignmentBlock ab: blocks){
			int start = ab.getReferenceStart()-1;
			int stop = ab.getLength() + start;
			//contained?
			if (centerBase >= start && centerBase< stop){
				int diff = centerBase - start;
				return (ab.getReadStart() -1) + diff;
			}
		}
		return -1;
	}
	
	public static int findReadIndexDebug(int centerBase, SAMRecord a) {
		List<AlignmentBlock> blocks = a.getAlignmentBlocks();
		for (AlignmentBlock ab: blocks){
			int start = ab.getReferenceStart()-1;
			int stop = ab.getLength() + start;
			System.out.println("\tBlock "+start+" "+stop);
			//contained?
			if (centerBase >= start && centerBase< stop){
				int diff = centerBase - start;
				return (ab.getReadStart() -1) + diff;
			}
		}
		return -1;
	}
	
	
	
	public int getCenterIndex() {
		return centerIndex;
	}
	public int getStartIndex() {
		return startIndex;
	}
	public int getEndIndex() {
		return endIndex;
	}
	public SAMRecord getAlignment() {
		return alignment;
	}
	public String getIndexStrand() {
		return indexStrand;
	}
	public float[] getWeightedCorrelationsACGT() {
		return weightedCorrelationsACGT;
	}
	public void setWeightedCorrelationsACGT(float[] weightedCorrelationsACGT) {
		this.weightedCorrelationsACGT = weightedCorrelationsACGT;
	}
	
	public void callCenterBaseWithCorrelation() {
		//sort the correlations
		float[] sortedCorr = new float[weightedCorrelationsACGT.length];
		System.arraycopy(weightedCorrelationsACGT, 0, sortedCorr, 0, sortedCorr.length);
		Arrays.sort(sortedCorr);

		//does it pass minimum correlation
		if (sortedCorr[3] < baseClassifier.getMinimumCorrelation()) {
			baseClassifier.incrementNumVarAlignFailingMinCor();
			baseClassifier.incrementNumVarAlignChanged2N();
			changedToN = true;
			classifiedCenterSeq = "N";
		}

		//does it pass minimum diff
		else if ( (sortedCorr[3] - sortedCorr[2]) < baseClassifier.getMinimumCorrelationDifference()) {
			baseClassifier.incrementNumVarAlignFailingMinCorDiff();
			baseClassifier.incrementNumVarAlignChanged2N();
			changedToN = true;
			classifiedCenterSeq = "N";
		}

		//assign base, simple method.
		else {
			int index = 0;
			float biggest = weightedCorrelationsACGT[0];
			for (int i=1; i< 4; i++){
				if (weightedCorrelationsACGT[i] > biggest) {
					biggest = weightedCorrelationsACGT[i];
					index = i;
				}
			}
			if (index == 0) classifiedCenterSeq = "A";
			else if (index == 1) classifiedCenterSeq = "C";
			else if (index == 2) classifiedCenterSeq = "G";
			else classifiedCenterSeq = "T";

			//did it change?
			if (classifiedCenterSeq.equals(centerSeq) == false) {
				baseClassifier.incrementNumVarAlignChanged2Novel();
				changedToDiff = true;
			}
		}
		//add reference
		if (changedToN || changedToDiff) {
			String seqOrder = "1";
			if (alignment.getSecondOfPairFlag()) seqOrder = "2";
			baseClassifier.getChangedAlignments().put(alignment.getReadName()+seqOrder, alignment);
			changeSAMRecord();
		}
	}
	
	public void callCenterBaseWithLikelihood() {
		//sort the correlations
		double[] sorted = new double[likelihoodACGTRatios.length];
		System.arraycopy(likelihoodACGTRatios, 0, sorted, 0, sorted.length);
		Arrays.sort(sorted);

		//does it pass minimum likelihood?
		if (sorted[3] < baseClassifier.getMinimumLikelihood()) {
			baseClassifier.incrementNumVarAlignFailingMinCor();
			baseClassifier.incrementNumVarAlignChanged2N();
			changedToN = true;
			classifiedCenterSeq = "N";
		}

		//does it pass minimum diff
		else if ( (sorted[3] - sorted[2]) < baseClassifier.getMinimumCorrelationDifference()) {
			baseClassifier.incrementNumVarAlignFailingMinCorDiff();
			baseClassifier.incrementNumVarAlignChanged2N();
			changedToN = true;
			classifiedCenterSeq = "N";
		}

		//assign base, simple method.
		else {
			int index = 0;
			double biggest = likelihoodACGTRatios[0];
			for (int i=1; i< 4; i++){
				if (likelihoodACGTRatios[i] > biggest) {
					biggest = likelihoodACGTRatios[i];
					index = i;
				}
			}
			if (index == 0) classifiedCenterSeq = "A";
			else if (index == 1) classifiedCenterSeq = "C";
			else if (index == 2) classifiedCenterSeq = "G";
			else classifiedCenterSeq = "T";

			//did it change?
			if (classifiedCenterSeq.equals(this.centerSeq) == false) {
				baseClassifier.incrementNumVarAlignChanged2Novel();
				changedToDiff = true;
			}
		}
		//add reference
		if (changedToN || changedToDiff) {
			String seqOrder = "1";
			if (alignment.getSecondOfPairFlag()) seqOrder = "2";
			baseClassifier.getChangedAlignments().put(alignment.getReadName()+seqOrder, alignment);
			changeSAMRecord();
		}
	}
	
	public void changeSAMRecord(){
		//add new base or just N?
		String toChange2 = "N";
		if (baseClassifier.isSetRecalledBasesToN() == false) toChange2 = classifiedCenterSeq;
		
		//modify sequence 
		String seq = alignment.getReadString();
		seq = seq.substring(0, centerIndex) + toChange2 + seq.substring(centerIndex+1, seq.length());
		alignment.setReadString(seq);
	}
	
	public String getClassifiedCenterSeq() {
		return classifiedCenterSeq;
	}
	public boolean isChanged() {
		if (changedToN || changedToDiff) return true;
		return false;
	}
	public String getCenterSeq() {
		return centerSeq;
	}
	public void setLikelihoodACGTRatios(double[] likelihoodACGTRatios) {
		this.likelihoodACGTRatios = likelihoodACGTRatios;
	}
	public double[] getLikelihoodACGTRatios() {
		return likelihoodACGTRatios;
	}
	public boolean isChangedToN() {
		return changedToN;
	}
	public boolean isChangedToDiff() {
		return changedToDiff;
	}


}