package edu.utah.seq.data.sam;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.bio.seq.Seq;
import util.gen.Misc;

public class SamLayout{
	private char[] seq;
	private int[] qual;
	private char[] call;
	private Pattern CIGAR_SUB = Pattern.compile("(\\d+)([MDIN])");

	public SamLayout(int size){
		seq = new char[size];
		qual = new int[size];
		call = new char[size];
	}

	public SamLayout(char[] seq, int[] qual, char[] call) {
		this.seq = seq;
		this.qual = qual;
		this.call = call;
	}

	public void print(){
		printArray(seq);
		printArray(qual);
		printArray(call);
	}

	public void layoutCigar(int start, SamAlignment sam){
		//for each cigar block in first, looking for MDIN s
		Matcher mat = CIGAR_SUB.matcher(sam.getCigar());
		int position = sam.getPosition() - start;
		int index = 0;
		int layoutIndex = position;
		char[] samSeq = sam.getSequence().toCharArray();
		int[] samQual = Seq.convertScores(sam.getQualities());
		while (mat.find()){
			String cCall = mat.group(2);
			int numberBases = Integer.parseInt(mat.group(1));
			//System.out.println("\t\t"+mat.group());
			//a match
			if (cCall.equals("M")) {
				//layout Ms
				for (int i = 0; i< numberBases; i++){
					seq[layoutIndex] = samSeq[index];
					qual[layoutIndex] = samQual[index];
					call[layoutIndex] = 'M';
					index++;
					layoutIndex++;
				}
			}
			else if (cCall.equals("N")) {
				//layout Ns
				for (int i = 0; i< numberBases; i++){
					seq[layoutIndex] = 'N';
					qual[layoutIndex] = -1;
					call[layoutIndex] = 'N';
					layoutIndex++;
				}
			}
			else if (cCall.equals("D")) {
				//layout Ds, deletion in the read
				for (int i = 0; i< numberBases; i++){
					seq[layoutIndex] = 'D';
					qual[layoutIndex] = -2;
					call[layoutIndex] = 'D';
					layoutIndex++;
				}
			}
			else if (cCall.equals("I")) {
				//layout Is
				for (int i = 0; i< numberBases; i++){
					seq[layoutIndex] = samSeq[index];
					qual[layoutIndex] = samQual[index];
					call[layoutIndex] = 'I';
					index++;
					layoutIndex++;
				}
			}
			else Misc.printErrAndExit("\nError: unsupported CIGAR! Should never hit this, contact admin!\n");

		}
	}
	public static int[] deleteInt(int [] bad, int index){
		int[] good = new int[bad.length];
		System.arraycopy(bad, 0, good, 0, index);
		System.arraycopy(bad, index+1, good, index, bad.length-index-1);
		return good;
	}

	public static char[] deleteChar(char [] bad, int index){
		char[] good = new char[bad.length];
		System.arraycopy(bad, 0, good, 0, index);
		System.arraycopy(bad, index+1, good, index, bad.length-index-1);
		return good;
	}

	/**Joins the pair.  Returns null if it couldn't.*/
	public static SamLayout mergeLayouts(SamLayout f, SamLayout s, int minDiffQualScore, double minimumInFrameMismatch) {
		
		char[] seq = new char[f.seq.length];
		int[] qual = new int[seq.length];
		char[] call = new char[seq.length];
		boolean baseRemoved = false;
		for (int i=0; i< f.seq.length; i++){
			//do calls differ? and are not empty
			if (f.call[i] != s.call[i] && f.call[i] != '\u0000' && s.call[i] != '\u0000'){

				//are one of the calls an I and the other M
				if ((f.call[i]=='M' && s.call[i]=='I') || (f.call[i]=='I' && s.call[i]=='M')){
					//are downstream bases in frame?
					if (inframeDownstreamSeqI(f,s,i,minimumInFrameMismatch)){
						call[i] = 'I';
						if (f.call[i]=='I'){
							seq[i] = f.seq[i];
							qual[i] = f.qual[i];
						}
						else{
							seq[i] = s.seq[i];
							qual[i] = s.qual[i];
						}
						continue;
					}
					//not in frame so delete the insertion
					else{
						//has this been modified before if so give up
						if (baseRemoved) return null;
						//nope delete it
						if (f.call[i]=='I'){
							f.seq = deleteChar(f.seq, i);
							f.qual = deleteInt(f.qual, i);
							f.call = deleteChar(f.call, i);
						}
						else{
							s.seq = deleteChar(s.seq, i);
							s.qual = deleteInt(s.qual, i);
							s.call = deleteChar(s.call, i);
						}
						i--;
						baseRemoved = true;
						continue;
					}
				}
				
				//are one of the calls an N and the other a M?
				else if ((f.call[i]=='M' && s.call[i]=='N') || (f.call[i]=='N' && s.call[i]=='M')){
					call[i] = 'M';
					if (f.call[i]=='M'){
						seq[i] = f.seq[i];
						qual[i] = f.qual[i];
					}
					else{
						seq[i] = s.seq[i];
						qual[i] = s.qual[i];
					}
					continue;
				}
				
				// deletion?
				else if ((f.call[i]=='M' && s.call[i]=='D') || (f.call[i]=='D' && s.call[i]=='M')){
					//are downstream bases in frame?
					if (inframeDownstreamSeqD(f,s,i,minimumInFrameMismatch)){
						call[i] = 'D';
						if (f.call[i]=='D'){
							seq[i] = f.seq[i];
							qual[i] = f.qual[i];
						}
						else{
							seq[i] = s.seq[i];
							qual[i] = s.qual[i];
						}
						continue;
					}
					//nope not in frame so delete it
					else{
						//has this been modified before? if so give up!
						if (baseRemoved) return null;
						if (f.call[i]=='D'){
							f.seq = deleteChar(f.seq, i);
							f.qual = deleteInt(f.qual, i);
							f.call = deleteChar(f.call, i);
						}
						else{
							s.seq = deleteChar(s.seq, i);
							s.qual = deleteInt(s.qual, i);
							s.call = deleteChar(s.call, i);
						}
						i--;
						baseRemoved = true;
						continue;
					}
				}
				
				//something weird 
				else return null;

			}
			else call[i] = f.call[i];
			
			
			//check bases for differences
			
			//same base?
			if (f.seq[i] == s.seq[i]){
				seq[i] = f.seq[i];
				qual[i] = larger(f.qual[i], s.qual[i]);
			}
			//missing first
			else if (f.seq[i] == '\u0000'){
				seq[i] = s.seq[i];
				qual[i] = s.qual[i];
				call[i] = s.call[i];
			}
			//missing second
			else if (s.seq[i] == '\u0000'){
				seq[i] = f.seq[i];
				qual[i] = f.qual[i];
				call[i] = f.call[i];
			}
			//nope different base
			else {
				//assign base and initial quality
				if (f.qual[i] >= s.qual[i]) {
					seq[i] = f.seq[i];
					qual[i] = f.qual[i];
				}
				else {
					seq[i] = s.seq[i];
					qual[i] = s.qual[i];
				}
				//check score difference
				int diffScores = Math.abs(f.qual[i] - s.qual[i]);
				//too little difference to make clear call?
				if (diffScores < minDiffQualScore) qual[i] = 0;
			}
		}

		SamLayout sl = new SamLayout(seq, qual, call);
		return sl;
	}
	
	/**Returns the number of overlapping bases, non-overlapping bases, alignmentSize.*/
	public static int[] countOverlappingBases(SamLayout f, SamLayout s){
		int numNonOverlappingBases = 0;
		int numOverlappingBases = 0;
		int numBases = 0;
		int length = f.call.length;
		for (int i=0; i< length; i++){
			boolean fBase = (f.call[i] == 'M' || f.call[i] == 'I');
			boolean sBase = (s.call[i] == 'M' || s.call[i] == 'I');
			boolean gBase = (f.call[i] == '\u0000' && s.call[i] == '\u0000');
			if (fBase==true || sBase==true || gBase==true){
				numBases++;
				if (fBase == true && sBase == true) numOverlappingBases++;
				else if (fBase == true || sBase == true) numNonOverlappingBases++;
			}
			
		}
		return new int[]{numOverlappingBases, numNonOverlappingBases, numBases};
	}
	
	public String fetchCigar(){
		StringBuilder cigar = new StringBuilder();
		char currentCall = call[0];
		int numCurrentCall = 1;
		for (int i=1; i< call.length; i++){
			//has call changed?
			if (call[i] != currentCall){
				cigar.append(numCurrentCall);
				if (currentCall == '\u0000') cigar.append('N');
				else cigar.append(currentCall);
				currentCall = call[i];
				numCurrentCall = 1;
			}
			else numCurrentCall++;
		}
		//add last so long as it isn't blank
		if (currentCall != '\u0000'){
			cigar.append(numCurrentCall);
			cigar.append(currentCall);
		}
		
		return cigar.toString();
	}
	
	public void setSequenceAndQualities(SamAlignment mergedSam) {
		StringBuilder seqSB = new StringBuilder();
		StringBuilder qualSB = new StringBuilder();
		for (int i=0; i<call.length; i++){
			if (call[i] != 'N' && call[i] != 'D' &&  call[i] != '\u0000' ){
				seqSB.append(seq[i]);
				qualSB.append(Seq.ORDERED_ASCII[qual[i]]);
			}
		}
		mergedSam.setSequence(seqSB.toString());
		mergedSam.setQualities(qualSB.toString());
	}

	public static boolean inframeDownstreamSeqI (SamLayout f, SamLayout s, int startIndex, double minimumInFrameMismatch){
		double numMatch = 0;
		double numMisMatch = 0;
		for (int i=startIndex; i<f.seq.length; i++){
			//are both seqs present
			if (f.seq[i] != '\u0000' && s.seq[i] != '\u0000'){
				if (f.seq[i] == s.seq[i]) numMatch++;
				else numMisMatch++;
			}
		}
		double total = numMatch+numMisMatch;
		double fractionMismatch = numMisMatch/total;
		//System.err.println("\nInframe countsI "+fractionMismatch+" "+numMisMatch+" "+numMatch);
		if (total == 0 || fractionMismatch <= minimumInFrameMismatch) return true;
		return false;
	}

	public static boolean inframeDownstreamSeqD (SamLayout f, SamLayout s, int startIndex, double minimumInFrameMismatch){
		double numMatch = 0;
		double numMisMatch = 0;
		for (int i=startIndex; i<f.seq.length; i++){
			//are both seqs present and not deletions
			if (f.seq[i] != '\u0000' && s.seq[i] != '\u0000' && f.seq[i] != 'D' && s.seq[i] != 'D'){
				if (f.seq[i] == s.seq[i]) numMatch++;
				else numMisMatch++;
			}
		}
		double total = numMatch+numMisMatch;
		double fractionMismatch = numMisMatch/total;
		//System.err.println("\nInframe countsD "+fractionMismatch+" "+numMisMatch+" "+numMatch);
		if (total == 0 || fractionMismatch <= minimumInFrameMismatch) return true;
		return false;
	}

	public static int larger(int i, int j) {
		if (i>j) return i;
		return j;
	}

	/**Prints a char[] to System.out*/
	public static void printArray(char[] array){
		int len = array.length;
		for (int i=0; i<len; i++) {
			if (array[i] == '\u0000') System.out.print("- ");
			else System.out.print(array[i]+" ");
		}
		System.out.println();
	}

	/**Prints a char[] to System.out*/
	public static void printArray(int[] array){
		int len = array.length;
		for (int i=0; i<len; i++) {
			System.out.print(array[i]+" ");
		}
		System.out.println();
	}


}
