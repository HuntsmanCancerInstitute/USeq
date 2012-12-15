package edu.utah.seq.data.sam;

/**Class to deal with the horrendous bitwise flags used by SAM.
 * @author davidnix*/
public class SamAlignmentFlags {

	private boolean partOfAPairedAlignment = false;
	private boolean aProperPairedAlignment = false;
	private boolean unmapped = false;
	private boolean mateUnMapped = false;
	private boolean reverseStrand = false;
	private boolean mateReverseStrand = false;
	private boolean firstPair = false;
	private boolean secondPair = false;
	private boolean notAPrimaryAlignment = false;
	private boolean failedQC = false;
	private boolean aDuplicate = false;
	
	private short flags;
	private boolean flagResetNeeded;
	
	//constructor
	public SamAlignmentFlags(short flags){
		setBooleanFlags(flags);
	}
	public SamAlignmentFlags() {
	}
	
	//methods

	/*Called to set all booleans upon construction*/
	private void setBooleanFlags (short flags){
		this.flags = flags;
		flagResetNeeded = false;
		partOfAPairedAlignment = isPartOfAPairedAlignment();
		aProperPairedAlignment = isAProperPairedAlignment();
		unmapped = isUnmapped();
		mateUnMapped = isMateUnMapped();
		reverseStrand = isReverseStrand();
		mateReverseStrand = isMateReverseStrand();
		firstPair = isFirstPair();
		secondPair = isSecondPair();
		notAPrimaryAlignment = isNotAPrimaryAlignment();
		failedQC = failedQC();
		aDuplicate = isADuplicate();
		
	}

	/*Will be called to reset the flags value if any of the booleans have changed.*/
	private void resetFlags(){
		int count = 0;
		if (partOfAPairedAlignment) count += Math.pow(2, 0);
		if (aProperPairedAlignment) count += Math.pow(2, 1);
		if (unmapped) count += Math.pow(2, 2);
		if (mateUnMapped) count += Math.pow(2, 3);
		if (reverseStrand) count += Math.pow(2, 4);
		if (mateReverseStrand) count += Math.pow(2, 5);
		if (firstPair) count += Math.pow(2, 6);
		if (secondPair) count += Math.pow(2, 7);
		if (notAPrimaryAlignment) count += Math.pow(2, 8);
		if (failedQC) count += Math.pow(2, 9);
		if (aDuplicate) count += Math.pow(2, 10);
		flags = (short) count;
		flagResetNeeded = false;
		
	}
	
	/**
	 * Index	Index^2	MethodName	Flag	DescriptionFromSamSpec
	 * 0	1	isReadPartOfAPairedAlignment	0x0001	the read is paired in sequencing, no matter whether it is mapped in a pair 
	 * 1	2	isReadAProperPairedAlignment	0x0002	the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1 
	 * 2	4	isQueryUnmapped	0x0004	the query sequence itself is unmapped 
	 * 3	8	isMateUnMapped	0x0008	the mate is unmapped  
	 * 4	16	isQueryReverseStrand	0x0010	strand of the query (false for forward; true for reverse strand) 
	 * 5	32	isMateReverseStrand	0x0020	strand of the mate 
	 * 6	64	isReadFirstPair	0x0040	the read is the Þrst read in a pair 
	 * 7	128	isReadSecondPair	0x0080	the read is the second read in a pair 
	 * 8	256	isAlignmentNotPrimary	0x0100	the alignment is not primary (a read having split hits may have multiple primary alignment records) 
	 * 9	512	doesReadFailVendorQC	0x0200	the read fails platform/vendor quality checks 
	 * 10	1024	isReadADuplicate	0x0400	the read is either a PCR duplicate or an optical duplicate 
	 */
	boolean testBitwiseFlag(int testValue){
		if (flagResetNeeded) resetFlags();
		return ((flags & testValue) == testValue);
	}
	/**The read is paired in sequencing, no matter whether it is mapped in a pair*/
	public boolean isPartOfAPairedAlignment(){
		return testBitwiseFlag(1);
	}
	/**The read is mapped in a proper pair (depends on the protocol, normally inferred during alignment)*/
	public boolean isAProperPairedAlignment(){
		return testBitwiseFlag(2);
	}
	/**The query sequence itself is unmapped*/
	public boolean isUnmapped(){
		return testBitwiseFlag(4);
	}
	/**The mate is unmapped*/
	public boolean isMateUnMapped(){
		return testBitwiseFlag(8);
	}
	/**Strand of the query (false for forward; true for reverse strand)*/
	public boolean isReverseStrand(){
		return testBitwiseFlag(16);
	}
	/**Strand of the mate*/
	public boolean isMateReverseStrand(){
		return testBitwiseFlag(32);
	}
	/**The read is the Þrst read in a pair*/
	public boolean isFirstPair(){
		return testBitwiseFlag(64);
	}
	/**The read is the second read in a pair*/
	public boolean isSecondPair(){
		return testBitwiseFlag(128);
	}
	/**The alignment is not primary (a read having split hits may have multiple primary alignment records)*/
	public boolean isNotAPrimaryAlignment(){
		return testBitwiseFlag(256);
	}
	/**The read fails platform/vendor quality checks*/
	public boolean failedQC(){
		return testBitwiseFlag(512);
	}
	/**The read is either a PCR duplicate or an optical duplicate*/
	public boolean isADuplicate(){
		return testBitwiseFlag(1024);
	}

	public void printFlags(){
		if (flagResetNeeded) resetFlags();
		System.out.println(flags+"\tflags");
		System.out.println(isPartOfAPairedAlignment()+"\tisPartOfAPairedAlignment()");
		System.out.println(isAProperPairedAlignment()+"\tisAProperPairedAlignment()");
		System.out.println(isUnmapped()+"\tisUnmapped()");
		System.out.println(isMateUnMapped()+"\tisMateUnMapped()");
		System.out.println(isReverseStrand()+"\tisReverseStrand()");
		System.out.println(isMateReverseStrand()+"\tisMateReverseStrand()");
		System.out.println(isFirstPair()+"\tisFirstPair()");
		System.out.println(isSecondPair()+"\tisSecondPair()");
		System.out.println(isNotAPrimaryAlignment()+"\tisNotAPrimaryAlignment()");
		System.out.println(failedQC()+"\tfailedQC()");
		System.out.println(isADuplicate()+"\tisADuplicate()");

	}

	//setters
	public void setPartOfAPairedAlignment(boolean partOfAPairedAlignment) {
		this.partOfAPairedAlignment = partOfAPairedAlignment;
		flagResetNeeded = true;
	}

	public void setaProperPairedAlignment(boolean aProperPairedAlignment) {
		this.aProperPairedAlignment = aProperPairedAlignment;
		flagResetNeeded = true;
	}

	public void setUnmapped(boolean unmapped) {
		this.unmapped = unmapped;
		flagResetNeeded = true;
	}

	public void setMateUnMapped(boolean mateUnMapped) {
		this.mateUnMapped = mateUnMapped;
		flagResetNeeded = true;
	}

	public void setReverseStrand(boolean reverseStrand) {
		this.reverseStrand = reverseStrand;
		flagResetNeeded = true;
	}

	public void setMateReverseStrand(boolean mateReverseStrand) {
		this.mateReverseStrand = mateReverseStrand;
		flagResetNeeded = true;
	}

	public void setFirstPair(boolean firstPair) {
		this.firstPair = firstPair;
		flagResetNeeded = true;
	}

	public void setSecondPair(boolean secondPair) {
		this.secondPair = secondPair;
		flagResetNeeded = true;
	}

	public void setNotAPrimaryAlignment(boolean notAPrimaryAlignment) {
		this.notAPrimaryAlignment = notAPrimaryAlignment;
		flagResetNeeded = true;
	}

	public void setFailedQC(boolean failedQC) {
		this.failedQC = failedQC;
		flagResetNeeded = true;
	}

	public void setaDuplicate(boolean aDuplicate) {
		this.aDuplicate = aDuplicate;
		flagResetNeeded = true;
	}

	public void setFlags(short flags) {
		setBooleanFlags(flags);
	}


	public short getFlags() {
		if (flagResetNeeded) resetFlags();
		return flags;
	}


	
}
