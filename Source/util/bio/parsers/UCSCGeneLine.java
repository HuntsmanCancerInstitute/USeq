package util.bio.parsers;
import util.bio.annotation.*;
import util.gen.*;
import java.util.*;
import edu.utah.seq.analysis.*;

public class UCSCGeneLine {
	//fields
	private String displayName = null;
	private String name;
	private String chrom;
	private String strand;
	private String notes;
	private int txStart;
	private int txEnd;
	private int cdsStart;
	private int cdsEnd;
	private ExonIntron[] exons = null;
	private ExonIntron[] introns = null;
	private int tss;
	private float[] scores;
	private SpliceJunction[] spliceJunctions = null;
	private float[][] exonCounts;
	private float[][] treatmentExonCounts;
	private float[][] controlExonCounts;
	//diff exp
	private float log2Ratio = 0;
	private float pValue = 0;
	private float fdr = 0;
	private float maxAbsLog2Ratio = 0;
	private float maxPValue = 0;
	//diff splicing
	private float splicingLog2Ratio = 0;
	private float splicingPValue = 0;
	private int splicingExon = 0;
	private StringBuilder text;
	private boolean flagged = false;
	//all condition edgeR ANOVA
	private float fdrEdgeR = 0;
	


	//constructors
	public UCSCGeneLine (){};
	public UCSCGeneLine (String name, Coordinate corr){
		this.displayName = name;
		chrom = corr.getChromosome();
		txStart = cdsStart = corr.getStart();
		txEnd = cdsEnd = corr.getStop();
		exons = new ExonIntron[1];
		exons[0] = new ExonIntron(txStart, txEnd);
		strand = ".";
	}
	
	public UCSCGeneLine (Bed bed){
		displayName = bed.getName();
		chrom = bed.getChromosome();
		txStart = cdsStart = bed.getStart();
		txEnd = cdsEnd = bed.getStop();
		exons = new ExonIntron[1];
		exons[0] = new ExonIntron(txStart, txEnd);
		strand = ""+bed.getStrand();
		scores = new float[]{(float)bed.getScore()};
	}
	
	public UCSCGeneLine (String[] tokens, int numToSubtractFromEnd, boolean displayNamePresent){
		try {
			int i = 0;
			if (displayNamePresent) {
				i = 1;
				displayName = tokens[0];
			}
			name = tokens[0+i];
			chrom = tokens[1+i];
			strand = tokens[2+i];
			txStart = Integer.parseInt(tokens[3+i]);
			txEnd = Integer.parseInt(tokens[4+i]);
			cdsStart = Integer.parseInt(tokens[5+i]);
			cdsEnd = Integer.parseInt(tokens[6+i]);
			txEnd-= numToSubtractFromEnd;
			cdsEnd-= numToSubtractFromEnd;

			//calc tss
			if (strand.equals("-")) tss = txEnd;
			else tss = txStart;

			//make exons
			int exonCount = Integer.parseInt(tokens[7+i]);
			if (exonCount > 0){
				exons = new ExonIntron[exonCount];
				int[] starts = Num.stringArrayToInts(tokens[8+i],",");
				int[] stops = Num.stringArrayToInts(tokens[9+i],",");		
				if (starts.length != exonCount || stops.length != exonCount) {
					Misc.printExit("\nError: problem parsing exons from "+Misc.stringArrayToString(tokens, "\t"));
				}
				for (int j=0; j< exonCount; j++){
					int stop = stops[j]-numToSubtractFromEnd;
					if (stop <= starts[j] ) Misc.printExit("\nError: One of your exon stops is <= its start! Correct and restart.  See ->\n"+Misc.stringArrayToString(tokens, "\t"));
					exons[j] = new ExonIntron(starts[j], stops[j]-numToSubtractFromEnd);
				}
				if (exonCount > 1) Arrays.sort(exons);
				//check that last exon stop is == txEnd, this is a hack to catch bad UCSCLines!  
				int lastExonEnd = exons[exonCount-1].getEnd();
				if (lastExonEnd != txEnd) txEnd = lastExonEnd;
				if (cdsEnd > lastExonEnd) cdsEnd = lastExonEnd;
			}

			//make introns
			if (exonCount > 1){
				//find smallest start and size
				int smallest = exons[0].getStart();
				int size = exons[exons.length-1].getEnd() - smallest;
				//make boolean array to mask
				boolean[] masked = new boolean[size+1];
				//for each exon mask boolean
				for (int k=0; k< exons.length; k++){
					int end = exons[k].getEnd() - smallest;
					int start = exons[k].getStart() - smallest;
					for (int j= start; j< end; j++){
						masked[j] = true;
					}
				}
				//fetch blocks of false
				int[][] startStop = ExportIntergenicRegions.fetchFalseBlocks(masked,0,0);
				//add smallest to starts and ends to return to real coordinates and make introns
				if (startStop.length>0){
					introns = new ExonIntron[startStop.length];
					for (int k=0; k< startStop.length; k++){
						introns[k] = new ExonIntron(startStop[k][0]+smallest, startStop[k][1]+smallest+1);
					}
				}
			}

			//any notes?
			if (tokens.length > 10+i){
				StringBuffer sb = new StringBuffer(tokens[10]);
				for (int j=11+i; j< tokens.length; j++){
					sb.append("\t");
					sb.append(tokens[j]);
				}
				notes = sb.toString();
			}
		} catch (Exception e){
			System.out.println("\nError parsing UCSCGeneLine "+Misc.stringArrayToString(tokens, " "));
			e.printStackTrace();
		}
	}
	
	public void zeroNullScores(){
		scores = null;
		exonCounts = null;
		treatmentExonCounts = null;
		controlExonCounts = null;
		log2Ratio = 0;
		pValue = 0;
		fdr = 0;
		splicingLog2Ratio =0;
		splicingPValue = 0;
		splicingExon = 0;
	}
	
	public UCSCGeneLine partialCloneForSpliceJunctionGeneration(){
		UCSCGeneLine clone = new UCSCGeneLine();
		clone.displayName = new String (this.displayName);
		clone.name = "merge";
		clone.chrom = this.chrom;
		clone.strand = this.strand;
		clone.exons = ExonIntron.clone(exons);
		clone.txStart = this.txStart;
		clone.txEnd = this.txEnd;
		clone.cdsStart = this.cdsStart;
		clone.cdsEnd = this.cdsEnd;
		return clone;
	}
	
	public String toStringBedFormat(){
		StringBuilder sb = new StringBuilder();
		sb.append(chrom);
		sb.append("\t");
		sb.append(txStart);
		sb.append("\t");
		sb.append(txEnd);
		sb.append("\t");
		sb.append(name);
		sb.append("\t");
		sb.append(scores[0]);
		sb.append("\t");
		sb.append(strand);
		return sb.toString();
	}
	
	public String getNames(String divider){
		if (displayName == null) return name;
		return name+ divider+ displayName;
	}

	public String scoresToString(){
		StringBuilder sb = new StringBuilder();
		if (displayName != null) {
			sb.append(displayName); sb.append("\t");
		}
		sb.append(name); sb.append("\t");
		sb.append(chrom); sb.append("\t");
		sb.append(txStart); sb.append("\t");
		sb.append(txEnd); 
		//scores
		if (scores != null){
			for (int i=0; i< scores.length; i++){
				sb.append("\t");
				sb.append(scores[i]);
			}
		}
		return sb.toString();
	}
	
	public String  coordinates(){
		StringBuilder sb = new StringBuilder();
		sb.append(chrom); sb.append("\t");
		sb.append(strand); sb.append("\t");
		sb.append(txStart); sb.append("\t");
		sb.append(txEnd); 
		return sb.toString();
	}
	
	/*Assumes interbase coordinates**/
	public boolean intersects(int start, int stop){
		if (stop < this.txStart || start >= this.txEnd) return false;
		return true;
	}
	/**Returns null if no overlap, otherwise the start and stop coordinates of the overlap.
	 * Assumes interbase coordinates*/
	public int[] fetchOverlap (int start, int stop) {
		if (intersects(start,stop) == false) return null;
		int beginning = 0;
		int ending = 0;
		//define the start
		if (this.txStart < start) beginning = start;
		else beginning = this.txStart;
		//define the end
		if (this.txEnd> stop) ending = stop;
		else ending = this.txEnd;
		return new int[]{beginning, ending};
	}
	
	public String exonsToString(){
		StringBuilder sb = new StringBuilder();
		String chrTab = chrom +"\t";
		for (int i=0; i<exons.length; i++){
			sb.append(chrTab);
			sb.append(exons[i].getStart());
			sb.append("\t");
			sb.append(exons[i].getEnd());
			sb.append("\n");
		}
		return sb.toString();
	}

	public String simpleToString(){
		StringBuffer sb = new StringBuffer();
		if (displayName != null) {
			sb.append(displayName); sb.append("\t");
		}
		sb.append(name); sb.append("\t");
		sb.append(chrom); sb.append("\t");
		sb.append(txStart); sb.append("\t");
		sb.append(txEnd); sb.append("\t");
		sb.append(strand); sb.append("\t");
		sb.append(tss); sb.append("\t");
		return sb.toString();
	}

	public String toString(){
		StringBuffer sb = new StringBuffer();
		if (displayName != null) {
			sb.append(displayName); sb.append("\t");
		}
		sb.append(name); sb.append("\t");
		sb.append(chrom); sb.append("\t");
		sb.append(strand); sb.append("\t");
		sb.append(txStart); sb.append("\t");
		sb.append(txEnd); sb.append("\t");
		sb.append(cdsStart); sb.append("\t");
		sb.append(cdsEnd); sb.append("\t");
		if (exons != null){
			sb.append(exons.length);sb.append("\t");
			sb.append(exons[0].getStart());
			for (int i=1; i<exons.length; i++){
				sb.append(",");
				sb.append(exons[i].getStart());
			}
			sb.append("\t");
			sb.append(exons[0].getEnd());
			for (int i=1; i<exons.length; i++){
				sb.append(",");
				sb.append(exons[i].getEnd());
			}
		}
		return sb.toString();
	}
	
	public String toStringBed12(int score0to1000){
		StringBuilder sb = new StringBuilder();
		//chrom
		sb.append(chrom);
		sb.append("\t");
		//start
		sb.append(this.txStart);
		sb.append("\t");
		//stop
		sb.append(this.txEnd);
		sb.append("\t");
		//name
		if (displayName != null && displayName.equals(name) == false){
			sb.append(displayName);
			sb.append(" ");
			sb.append(name);
		}
		else sb.append(name);
		sb.append("\t");
		//score 0 - 1000
		if (score0to1000 < 0 || score0to1000 > 1000) return null;
		sb.append(score0to1000);
		sb.append("\t");
		//strand
		sb.append(strand);
		sb.append("\t");
		//thick start (coding seq)
		sb.append(this.cdsStart);
		sb.append("\t");
		//thick end (coding seq)
		sb.append(this.cdsEnd);
		sb.append("\t");
		//rgb string 255,0,0
		sb.append(0);
		sb.append("\t");
		//number exons
		sb.append(exons.length);
		sb.append("\t");
		//exon sizes 321,9,1222
		////add first
		sb.append((exons[0].getEnd()-exons[0].getStart()));
		////add remainder
		for (int i=1; i< exons.length; i++){
			sb.append(",");
			sb.append((exons[i].getEnd()-exons[i].getStart()));
		}
		sb.append("\t");
		//exon starts relative to chromStart 0,234,9949
		////add first
		sb.append((exons[0].getStart()-txStart));
		////add remainder
		for (int i=1; i< exons.length; i++){
			sb.append(",");
			sb.append((exons[i].getStart()-txStart));
		}
		sb.append("\t");
		return sb.toString();
	}
	
	public String toStringBed12Float(float scoreFloat, String customName){
		StringBuilder sb = new StringBuilder();
		//chrom
		sb.append(chrom);
		sb.append("\t");
		//start
		sb.append(this.txStart);
		sb.append("\t");
		//stop
		sb.append(this.txEnd);
		sb.append("\t");
		//name
		sb.append(customName);
		sb.append("\t");
		sb.append(scoreFloat);
		sb.append("\t");
		//strand
		sb.append(strand);
		sb.append("\t");
		//thick start (coding seq)
		sb.append(this.cdsStart);
		sb.append("\t");
		//thick end (coding seq)
		sb.append(this.cdsEnd);
		sb.append("\t");
		//rgb string 255,0,0
		sb.append(0);
		sb.append("\t");
		//number exons
		sb.append(exons.length);
		sb.append("\t");
		//exon sizes 321,9,1222
		////add first
		sb.append((exons[0].getEnd()-exons[0].getStart()));
		////add remainder
		for (int i=1; i< exons.length; i++){
			sb.append(",");
			sb.append((exons[i].getEnd()-exons[i].getStart()));
		}
		sb.append("\t");
		//exon starts relative to chromStart 0,234,9949
		////add first
		sb.append((exons[0].getStart()-txStart));
		////add remainder
		for (int i=1; i< exons.length; i++){
			sb.append(",");
			sb.append((exons[i].getStart()-txStart));
		}
		sb.append("\t");
		return sb.toString();
	}
	
	public String toStringAll(){
		StringBuffer sb = new StringBuffer();
		if (displayName != null) {
			sb.append(displayName); sb.append("\t");
		}
		sb.append(name); sb.append("\t");
		sb.append(chrom); sb.append("\t");
		sb.append(strand); sb.append("\t");
		sb.append(txStart); sb.append("\t");
		sb.append(txEnd); sb.append("\t");
		sb.append(cdsStart); sb.append("\t");
		sb.append(cdsEnd); sb.append("\t");
		if (exons != null){
			sb.append(exons.length);sb.append("\t");
			sb.append(exons[0].getStart());
			for (int i=1; i<exons.length; i++){
				sb.append(",");
				sb.append(exons[i].getStart());
			}
			sb.append("\t");
			sb.append(exons[0].getEnd());
			for (int i=1; i<exons.length; i++){
				sb.append(",");
				sb.append(exons[i].getEnd());
			}
			sb.append("\t");
		}		
		if (introns != null){
			sb.append(introns.length);sb.append("\t");
			sb.append(introns[0].getStart());
			for (int i=1; i<introns.length; i++){
				sb.append(",");
				sb.append(introns[i].getStart());
			}
			sb.append("\t");
			sb.append(introns[0].getEnd());
			for (int i=1; i<introns.length; i++){
				sb.append(",");
				sb.append(introns[i].getEnd());
			}
		}
		else sb.append("\t");
		
		if (notes != null) {
			sb.append("\t");
			sb.append(notes);
		}
		return sb.toString();
	}
	
	public int getTotalExonicBasePairs(){
		if (exons == null) return 0;
		int totalExonicBP = 0;
		for (int k=0; k< exons.length; k++){
			int start = exons[k].getStart();
			int stop = exons[k].getEnd();
			totalExonicBP+= (stop - start);
		}
		return totalExonicBP;
	}
	
	public int getLength(){
		return txEnd-txStart;
	}

	public int getCdsEnd() {
		return cdsEnd;
	}

	public int getCdsStart() {
		return cdsStart;
	}

	public String getChrom() {
		return chrom;
	}

	public ExonIntron[] getExons() {
		return exons;
	}

	public String getName() {
		return name;
	}

	public String getStrand() {
		return strand;
	}

	public int getTxEnd() {
		return txEnd;
	}

	public int getTxStart() {
		return txStart;
	}

	public String getNotes() {
		return notes;
	}

	public void setNotes(String notes) {
		this.notes = notes;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getTss() {
		return tss;
	}

	public float[] getScores() {
		return scores;
	}

	public void setScores(float[] scores) {
		this.scores = scores;
	}

	public String getDisplayName() {
		return displayName;
	}
	
	public String getDisplayNameThenName() {
		if (displayName != null) return displayName;
		return name;
	}

	public void setDisplayName(String displayName) {
		this.displayName = displayName;
	}

	public void setCdsEnd(int cdsEnd) {
		this.cdsEnd = cdsEnd;
	}

	public void setCdsStart(int cdsStart) {
		this.cdsStart = cdsStart;
	}

	public void setExons(ExonIntron[] exons) {
		this.exons = exons;
	}

	public void setTss(int tss) {
		this.tss = tss;
	}

	public void setTxEnd(int txEnd) {
		this.txEnd = txEnd;
	}

	public void setTxStart(int txStart) {
		this.txStart = txStart;
	}

	public ExonIntron[] getIntrons() {
		return introns;
	}

	public void setIntrons(ExonIntron[] introns) {
		this.introns = introns;
	}

	public SpliceJunction[] getSpliceJunctions() {
		return spliceJunctions;
	}

	public void setSpliceJunctions(SpliceJunction[] spliceJunctions) {
		this.spliceJunctions = spliceJunctions;
	}

	public float[][] getExonCounts() {
		return exonCounts;
	}

	public void setExonCounts(float[][] exonCounts) {
		this.exonCounts = exonCounts;
	}

	public void setChrom(String chrom) {
		this.chrom = chrom;
	}
	/**float[replicas][exonCounts]*/
	public float[][] getTreatmentExonCounts() {
		return treatmentExonCounts;
	}
	/**float[replicas][exonCounts]*/
	public void setTreatmentExonCounts(float[][] treatmentExonCounts) {
		this.treatmentExonCounts = treatmentExonCounts;
	}
	/**float[replicas][exonCounts]*/
	public float[][] getControlExonCounts() {
		return controlExonCounts;
	}
	/**float[replicas][exonCounts]*/
	public void setControlExonCounts(float[][] controlExonCounts) {
		this.controlExonCounts = controlExonCounts;
	}
	public float getLog2Ratio() {
		return log2Ratio;
	}
	public void setLog2Ratio(float log2Ratio) {
		this.log2Ratio = log2Ratio;
	}
	public StringBuilder getText() {
		return text;
	}
	public void setText(StringBuilder text) {
		this.text = text;
	}
	public void setpValue(float pValue) {
		this.pValue = pValue;
	}
	public float getpValue() {
		return pValue;
	}
	public float getSplicingLog2Ratio() {
		return splicingLog2Ratio;
	}
	public void setSplicingLog2Ratio(float splicingLog2Ratio) {
		this.splicingLog2Ratio = splicingLog2Ratio;
	}
	public float getSplicingPValue() {
		return splicingPValue;
	}
	public void setSplicingPValue(float splicingPValue) {
		this.splicingPValue = splicingPValue;
	}
	public float getFdr() {
		return fdr;
	}
	public void setFdr(float fdr) {
		this.fdr = fdr;
	}
	public float getMaxAbsLog2Ratio() {
		return maxAbsLog2Ratio;
	}
	public void setMaxAbsLog2Ratio(float maxAbsLog2Ratio) {
		this.maxAbsLog2Ratio = maxAbsLog2Ratio;
	}
	public float getMaxPValue() {
		return maxPValue;
	}
	public void setMaxPValue(float maxPValue) {
		this.maxPValue = maxPValue;
	}
	public int getSplicingExon() {
		return splicingExon;
	}
	public void setSplicingExon(int splicingExon) {
		this.splicingExon = splicingExon;
	}
	public boolean isFlagged() {
		return flagged;
	}
	public void setFlagged(boolean flaggedGene) {
		this.flagged = flaggedGene;
	}
	public float getFdrEdgeR() {
		return fdrEdgeR;
	}
	public void setFdrEdgeR(float fdrEdgeR) {
		this.fdrEdgeR = fdrEdgeR;
	}
}
