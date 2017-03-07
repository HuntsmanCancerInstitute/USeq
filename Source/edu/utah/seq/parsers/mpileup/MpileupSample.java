package edu.utah.seq.parsers.mpileup;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import util.bio.seq.Seq;
import util.gen.Misc;
import util.gen.Num;

public class MpileupSample {

	//fields
	private int mpileupReadCount = 0;
	private String baseSymbols = "";
	private String qualitySymbols = "";
	private char[] maskedBaseCalls = new char[0];
	private MpileupLine record;
	private int minBaseQuality;
	
	//mpileup only marks the upstream base as having an insertion, not the downstream, IGV marks both
	private int insertions = 0;
	private int deletions = 0;
	private int[] forwardGATC = new int[4];
	private int[] reverseGATC = new int[4];
	private int poorQualBases = 0;
	private int readCoverageAll = 0;
	private int readCoverageForwardBases = 0;
	private int readCoverageReverseBases = 0;
	private boolean pass = true;
	
	public static final int G_INDEX = 0;
	public static final int A_INDEX = 1;
	public static final int T_INDEX = 2;
	public static final int C_INDEX = 3;
	
	//Pattern fields
	public static final Pattern comma = Pattern.compile(",");
	public static final Pattern period = Pattern.compile("\\.");
	public static final Pattern indel = Pattern.compile("(\\+|-)([0-9]+)");
	public static final Pattern start = Pattern.compile("\\^");
	public static final Pattern doubleHats = Pattern.compile("\\^\\^");
	public static final Pattern odd = Pattern.compile("[\\$<>]");
	
	public MpileupSample(String depth, String baseSymbols, String qualitySymbols, MpileupLine record) {
		mpileupReadCount = Integer.parseInt(depth);
		this.baseSymbols = baseSymbols;
		this.qualitySymbols = qualitySymbols;
		this.record = record;
		minBaseQuality = record.getMinBaseQuality();
		
		if (mpileupReadCount !=0 ){
			char[] maskedBases = maskBaseSymbols();
			scanBases(maskedBases);
		}
	}
	
	public MpileupSample() {}

	/**Appends a comma delimited list counts of observed base observations that pass the minBaseQual: GATC forward, GATC reverse, ins, del, as well as the num that failed minBaseQual so 11 values*/
	public void appendCounts(StringBuilder sb, boolean combineStrands){
		if (combineStrands){
			for (int i=0; i< forwardGATC.length; i++){
				sb.append((forwardGATC[i]+reverseGATC[i]));
				sb.append(",");
			}
		}
		else {
			for (int c: forwardGATC){
				sb.append(c);
				sb.append(",");
			}
			for (int c: reverseGATC){
				sb.append(c);
				sb.append(",");
			}
		}
		sb.append(insertions);
		sb.append(",");
		sb.append(deletions);
		sb.append(",");
		sb.append(poorQualBases);
		sb.append(",");
		sb.append((getNonRefBaseCounts()+ insertions + deletions));
		sb.append(",");
		sb.append(getAlleleFreqNonRefPlusIndels());
	}
	
	public void debug(){
		System.out.println("LN: "+record.getLine());
		System.out.println("BS: "+baseSymbols);
		System.out.println("MC: "+new String(maskedBaseCalls));
		System.out.println("QS: "+qualitySymbols);
		System.out.print("QV: "); Misc.printArray(Seq.convertSangerQualityScores(qualitySymbols));
		System.out.println("INS: "+insertions);
		System.out.println("DEL: "+deletions);
		System.out.println("PQ: "+poorQualBases);
		System.out.print("F: " ); Misc.printArray(forwardGATC);
		System.out.print("R: " ); Misc.printArray(reverseGATC);
		int total = deletions +poorQualBases+ readCoverageForwardBases + readCoverageReverseBases;
		System.out.println("TOT: "+total +" from rec "+mpileupReadCount);
	}
	
	public void printMinInfo(){
		System.out.println(record.getChr()+" "+(record.getZeroPos()+1)+" "+readCoverageAll+" "+ getAlleleFreqNonRefPlusIndels());
	}
	
	
	
	public void scanBases(char[] bases){
		//fetch qualities
		int[] qualScores = Seq.convertSangerQualityScores(qualitySymbols);
		
		int qsIndex = 0;
		for (int i=0; i< bases.length; i++){
			switch (bases[i]){
			case 'M': break; //skip it's masked
			case 'g': count(reverseGATC, G_INDEX, qualScores[qsIndex++]); break;
			case 'G': count(forwardGATC, G_INDEX, qualScores[qsIndex++]); break;
			case 'a': count(reverseGATC, A_INDEX, qualScores[qsIndex++]); break;
			case 'A': count(forwardGATC, A_INDEX, qualScores[qsIndex++]); break;
			case 't': count(reverseGATC, T_INDEX, qualScores[qsIndex++]); break;
			case 'T': count(forwardGATC, T_INDEX, qualScores[qsIndex++]); break;
			case 'c': count(reverseGATC, C_INDEX, qualScores[qsIndex++]); break;
			case 'C': count(forwardGATC, C_INDEX, qualScores[qsIndex++]); break;
			case 'n': poorQualBases++; qsIndex++; break;
			case 'N': poorQualBases++; qsIndex++; break;
			case '*': countDeletion(qualScores[qsIndex++]); break;
			default: Misc.printExit("\nProblem, unknown base symbol! " + bases[i] +" in "+ record.getLine());
			}
		}
		//check lengths
		if ((qualScores.length) != qsIndex) {
			debug();
			Misc.printErrAndExit("Quality score to masked base mismatch! "+record.getLine());
		}
		//compare with record
		readCoverageForwardBases = Num.sumIntArray(forwardGATC);
		readCoverageReverseBases = Num.sumIntArray(reverseGATC);
		
		if (mpileupReadCount != (poorQualBases+ deletions+ readCoverageForwardBases+ readCoverageReverseBases)) {
			debug();
			Misc.printErrAndExit("Read count mismatch! "+record.getLine());
		}
		
		readCoverageAll = deletions + insertions + readCoverageForwardBases + readCoverageReverseBases; 
	} 
	
	/**Modifies the first sample by adding the counts from the subsequent, note! the Strings and char[] aren't joined, just the counts.
	 * Don't use the first sample after calling this method.*/
	public static MpileupSample mergeSampleCounts(MpileupSample[] s){
		MpileupSample m = s[0];
		for (int i=1; i< s.length; i++) {
			if (s[i].mpileupReadCount == 0) continue;
			m.mpileupReadCount+= s[i].mpileupReadCount;
			m.insertions+= s[i].insertions;
			m.deletions+= s[i].deletions;
			for (int x=0; x< m.forwardGATC.length; x++){
				m.forwardGATC[x] += s[i].forwardGATC[x];
				m.reverseGATC[x] += s[i].reverseGATC[x];
			}
			m.poorQualBases+= s[i].poorQualBases;
			m.readCoverageAll+= s[i].readCoverageAll;
			m.readCoverageForwardBases+= s[i].readCoverageForwardBases;
			m.readCoverageReverseBases+= s[i].readCoverageReverseBases;
		}
		return m;
	}
	
	public double getAlleleFreqPoorQual(){
		return ((double)(poorQualBases))/((double)readCoverageAll+ poorQualBases);
	}
	
	public double getAlleleFreqINDEL(){
		return ((double)(insertions+deletions))/((double)readCoverageAll);
	}
	
	public double getAlleleFreq(int index){
		return ((double)(forwardGATC[index]+reverseGATC[index])) / ((double)(readCoverageForwardBases+ readCoverageReverseBases));
	}
	
	public double getAlleleFreqForward(int index){
		return ((double)(forwardGATC[index])) / ((double)(readCoverageForwardBases));
	}
	
	public double getAlleleFreqReverse(int index){
		return ((double)(reverseGATC[index])) / ((double)(readCoverageReverseBases));
	}
	
	public double getAlleleFreqNonRef(){
		double nonRefCounts = getNonRefBaseCounts();
		return nonRefCounts/((double)(readCoverageForwardBases+ readCoverageReverseBases));
	}
	
	public double getAlleleFreqNonRefPlusIndels(){
		double nonRefCounts = getNonRefBaseCounts()+ insertions + deletions;
		return nonRefCounts/((double)(readCoverageAll));
	}
	
	public int getNonRefBaseCounts() {
		char ref = record.getRef().charAt(0);
		int index = -1;
		if (ref == 'G') index = G_INDEX;
		else if (ref == 'A') index = A_INDEX;
		else if (ref == 'T') index = T_INDEX;
		else if (ref == 'C') index = C_INDEX;
		return readCoverageForwardBases - forwardGATC[index] + readCoverageReverseBases - reverseGATC[index];
	}

	private void countDeletion(int baseQuality) {
		if (baseQuality < minBaseQuality) poorQualBases++;
		else deletions++;
	}

	private void count(int[] baseCounts, int index, int baseQuality) {
		if (baseQuality < minBaseQuality) poorQualBases++;
		else baseCounts[index]++;
	}

	private char[] maskBaseSymbols(){
		String ref = record.getRef();
		
		//replace . and , with Ref or ref
		String y = period.matcher(baseSymbols).replaceAll(ref);
		y = comma.matcher(y).replaceAll(ref.toLowerCase());
		
		//replace all ^^ with MM
		y= doubleHats.matcher(y).replaceAll("MM");

		maskedBaseCalls = y.toCharArray();

		//mask indels
		Matcher mat = indel.matcher(y);
		while (mat.find()){
			int startMask=mat.start();
			int endMask = mat.end()+ Integer.parseInt(mat.group(2));
			if (mat.group(1).equals("+")) insertions++;
			//don't count deletions, just mask
			for (int i=startMask; i< endMask; i++) maskedBaseCalls[i] = 'M';
		}

		//mask ^ and base afterward
		mat = start.matcher(y);
		while (mat.find()){
			maskedBaseCalls[mat.end()] = 'M';
			maskedBaseCalls[mat.end()-1] = 'M';
		}

		//mask $ > <
		mat = odd.matcher(y);
		while (mat.find())maskedBaseCalls[mat.start()] = 'M';

		return maskedBaseCalls;
	}
	
	public double findMaxSnvAF(){
		double maxAF = Double.MIN_NORMAL;
		char ref = record.getRef().charAt(0);
		int[] indexesToScan = null;
		if (ref == 'G') indexesToScan = new int[]{A_INDEX, T_INDEX, C_INDEX};
		else if (ref == 'A') indexesToScan = new int[]{G_INDEX, T_INDEX, C_INDEX};
		else if (ref == 'T') indexesToScan = new int[]{G_INDEX,A_INDEX,C_INDEX};
		else if (ref == 'C') indexesToScan = new int[]{G_INDEX,A_INDEX, T_INDEX};
		for (int index: indexesToScan){
			double af = getAlleleFreq(index);
			if (af> maxAF) maxAF = af;
		}
		return maxAF;
	}

	public int getReadCoverageAll() {
		return readCoverageAll;
	}

	public int getReadCoverageForwardBases() {
		return readCoverageForwardBases;
	}

	public int getReadCoverageReverseBases() {
		return readCoverageReverseBases;
	}

	public boolean isPass() {
		return pass;
	}

	public void setPass(boolean pass) {
		this.pass = pass;
	}

	public MpileupLine getRecord() {
		return record;
	}

	public int[] getForwardGATC() {
		return forwardGATC;
	}

	public int[] getReverseGATC() {
		return reverseGATC;
	}

	public int getMpileupReadCount() {
		return mpileupReadCount;
	}

	public int getInsertions() {
		return insertions;
	}

	public int getDeletions() {
		return deletions;
	}


}
