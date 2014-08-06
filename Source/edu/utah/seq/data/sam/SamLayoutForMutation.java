package edu.utah.seq.data.sam;

import java.util.regex.Matcher;
import edu.utah.seq.vcf.VCFRecord;
import util.bio.seq.Seq;
import util.gen.Misc;
import util.gen.Num;

public class SamLayoutForMutation{
	private char[] seq;
	private int[] qual;
	private char[] call;
	private int[] pos;
	private SamAlignment sam;

	public SamLayoutForMutation(SamAlignment sam){
		this.sam = sam;
		int size = SamAlignment.countLengthOfCIGAR(sam.getCigar());
		seq = new char[size];
		qual = new int[size];
		call = new char[size];
		pos = new int[size];
		layoutCigar(sam);
	}

	public void print(){
		printArray(seq);
		printArray(qual);
		printArray(call);
		printArray(pos);
	}
	
	public static void main (String[] args){
		SamAlignment s;
		try {
			//s = new SamAlignment("HWI-ST1124:106:C15APACXX:1:1103:12835:181960\t147\t10\t42928657\t70\t12M9D88M\t=\t42928620\t-146\tTTTTGAATGAAGCCTTGGCATTATAGTTGTACTCACGATTTTAGTGTCATTGTAAGATGAGGCTGGCTGGCTATGCTGTCATCAAGGAATATTGTCGAAC\t;;<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<9<<<;<<;;;\tMD:Z:12^G0^A0^G0^T0^C0^T0^A0^T0^C79G8\tPG:Z:novoalign\tAM:i:70\tNM:i:10\tSM:i:70\tPQ:i:225\tUQ:i:124\tAS:i:124", true);
			//s = new SamAlignment("HWI-ST1124:106:C15APACXX:1:1103:18430:24840\t147\t10\t132250\t70\t29H71M\t=\t132029\t-292\tCACCCCAGCAGATTGCTCAGCGGGTCTCCATCCGCTCTCACACTGTGTCCTGTGAAAGACGCTGGGCACCA\t<:299(;;:<<82:<;487+/<;6)<<8<835(6/):1(4(,/(9:6(:;7,1<:995,<;<<<::<5++;\tMD:Z:32T2A14C20\tPG:Z:novoalign\tAM:i:70\tNM:i:3\tSM:i:70\tPQ:i:69\tUQ:i:60\tAS:i:60",true);
			//s = new SamAlignment("HWI-ST1124:106:C15APACXX:1:1101:14173:112386\t147\t10\t1041944\t70\t80M3I17M\t=\t1041826\t-215\tATAAATATTGCCAAAAATTTAGTGGACAAGTAAGGTACTTTTTTCTGTAGTGTCTTTTAAGTTATCGAAAGTGGAAAAAAAACTAATATTTTTTATTTTT\t<;<<<;<<<<<<<<<==<<<<<;<<<=<;<<<<<<;:<;;<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<;<<<<<<<<<<<<<<<<<<<<;;;;;\tMD:Z:74T22\tPG:Z:novoalign\tAM:i:1\tNM:i:4\tSM:i:1\tPQ:i:96\tUQ:i:88\tAS:i:88",true);
			//s = new SamAlignment("HWI-ST1124:106:C15APACXX:1:1103:5502:56189\t83\t10\t357440\t70\t96M4S\t=\t357435\t-101\tGGGAGGCNCGGGTGTGATGGACAGGAAAGTCCCCTACCTGCTCGGCCAGGGTCCACAGCCTGTGAGGCTACAGCAGCAGAAGGGTGAGGACGAGAAGAGA\t<<;<<<:#;<;<<<<<<<<<<<<<<<<;;8<<<<;;<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<;<<;:;<<<<<<<<<;<;<;4<<<<<<<;:;\tMD:Z:7A60C27\tPG:Z:novoalign\tAM:i:70\tNM:i:2\tSM:i:70\tPQ:i:105\tUQ:i:66\tAS:i:66", true);
			//s = new SamAlignment("HWI-ST1124:106:C15APACXX:1:1104:6073:105070\t83\t22\t16416318\t26\t8H5S87M\t=\t16416197\t-208\tTTATCTTTCCTTAAAAAAGAAATGTTTTAATCCATCACATTTTTCTTCCCTTCCCTTTAGTTTTTGATAAATGATAAAAATGAGCCAGTTAT\t<;<<;<<<;78;4--;;;;;;<<<<;;<<<<<<;9;:;:;<<<<<<<<<<8*<;;89;<<<:<<:<::<;9;9;<<;<<<;;7<<<:902:8\tMD:Z:46A40\tPG:Z:novoalign\tAM:i:70\tNM:i:1\tSM:i:70\tPQ:i:64\tUQ:i:60\tAS:i:60", true);
			s= new SamAlignment("HWI-ST1124:106:C15APACXX:1:1104:11418:109939\t83\t22\t25601197\t70\t3S97M\t=\t25601123\t-171\tCAGATCANAAGCTGCATCTGTTTGAGAACCCAGCTTTCAGTGGCCGCAAGATGGAGATAGTGGATGATGACGTGCCCAGCCTGTGGGCTCATGGCTTCCA\t<<<<;<2#<<<<<<<<<<<<<<<<;<;:<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<;\tMD:Z:4C92\tPG:Z:novoalign\tAM:i:70\tNM:i:1\tSM:i:70\tPQ:i:66UQ:i:36\tAS:i:36", true);
			System.out.println("Ori\t"+s);
			SamLayoutForMutation sl = new SamLayoutForMutation(s);
			sl.print();
			//sl.changeBase(1042017, '*');
			//sl.print();
			sl.setCigar();
			sl.setSequenceAndQualities();
			sl.setAlignmentPosition();
			s.removeMD();
			System.out.println("Mod\t"+s);
			//sl.print();
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
	}
	
	/**Sets new CIGAR, Seq, Qual, Position, and removed MD field if present.*/
	public void setChangesInSam(){
		setCigar();
		setSequenceAndQualities();
		setAlignmentPosition();
		sam.removeMD();
	}
	
	public boolean changeBase(int bpPosition, char reference, char mutation){
		int index = Num.findSmallestIndexToValue(pos, bpPosition);
		if (index < 0 || index >= seq.length) return false;
		//check for reference, skip if not, represents a bad base call or mis alignment artifact
		if (seq[index] != reference) return false;
		seq[index] = mutation;
		return true;
	}
	
	public int findBaseIndex(int bpPosition){
		int index = Num.findSmallestIndexToValue(pos, bpPosition);
		if (index < 0 || index >= seq.length) return -1;
		return index;
	}
	
	public void changeInsertion(VCFRecord vcf) {
		System.out.println("Swapping Insertion");
		print();
		System.exit(0);
	}

	public void layoutCigar(SamAlignment sam){
		//for each cigar block in first, looking for MDIN s
		Matcher mat = SamAlignment.CIGAR_SUB.matcher(sam.getCigar());
		int position = sam.getUnclippedStart();
		int index = 0;
		int layoutIndex = 0;
		char[] samSeq = sam.getSequence().toUpperCase().toCharArray();
		int[] samQual = Seq.convertScores(sam.getQualities());
		while (mat.find()){
			String cCall = mat.group(2);
			int numberBases = Integer.parseInt(mat.group(1));
			//a match
			if (cCall.equals("M")) {
				//layout Ms
				for (int i = 0; i< numberBases; i++){
					seq[layoutIndex] = samSeq[index];
					qual[layoutIndex] = samQual[index];
					call[layoutIndex] = 'M';
					pos[layoutIndex] = position;
					index++;
					layoutIndex++;
					position++;
				}
			}
			else if (cCall.equals("N")) {
				//layout Ns
				for (int i = 0; i< numberBases; i++){
					seq[layoutIndex] = 'N';
					qual[layoutIndex] = -1;
					call[layoutIndex] = 'N';
					pos[layoutIndex] = position;
					position++;
					layoutIndex++;
				}
			}
			else if (cCall.equals("D")) {
				//layout Ds, deletion in the read
				for (int i = 0; i< numberBases; i++){
					seq[layoutIndex] = 'D';
					qual[layoutIndex] = -1;
					call[layoutIndex] = 'D';
					pos[layoutIndex] = position;
					position++;
					layoutIndex++;
					
				}
			}
			else if (cCall.equals("I")) {
				//layout Is
				for (int i = 0; i< numberBases; i++){
					seq[layoutIndex] = samSeq[index];
					qual[layoutIndex] = samQual[index];
					call[layoutIndex] = 'I';
					pos[layoutIndex] = position;
					index++;
					layoutIndex++;
				}
			}
			else if (cCall.equals("S")){
				for (int i = 0; i< numberBases; i++){
					seq[layoutIndex] = samSeq[index];
					qual[layoutIndex] = samQual[index];
					call[layoutIndex] = 'S';
					pos[layoutIndex] = position;
					index++;
					layoutIndex++;
					position++;
				}
			}
			else if (cCall.equals("H")){
				for (int i = 0; i< numberBases; i++){
					seq[layoutIndex] = 'H';
					qual[layoutIndex] = -1;
					call[layoutIndex] = 'H';
					pos[layoutIndex] = position;
					position++;
					layoutIndex++;
				}
			}
			else Misc.printErrAndExit("\nError: unsupported CIGAR! Should never hit this, contact admin!\n");

		}
	}
	
	/**Skip first H or S position in sam. Not both though! For example Novoalign puts the position on the first S when H precedes it, eg 9H5S87M.  Seems like a bug! */
	public void setAlignmentPosition(){
		int ap = -1;
		//first call an H?
		char toSkip = ' ';
		if (call[0] == 'H') toSkip = 'H';
		else if (call[0] == 'S') toSkip = 'S';
		for (int i=0; i< pos.length; i++){
			if (call[i] != toSkip) {
				ap = pos[i];
				break;
			}
		}
		sam.setPosition(ap);
	}
	
	public void setCigar(){
		StringBuilder cigar = new StringBuilder();
		char currentCall = call[0];
		int numCurrentCall = 1;
		for (int i=1; i< call.length; i++){
			//has call changed?
			if (call[i] != currentCall){
				cigar.append(numCurrentCall);
				cigar.append(currentCall);
				currentCall = call[i];
				numCurrentCall = 1;
			}
			else numCurrentCall++;
		}
		//add last
		cigar.append(numCurrentCall);
		cigar.append(currentCall);
		sam.setCigar(cigar.toString());;
	}
	
	public void setSequenceAndQualities() {
		StringBuilder seqSB = new StringBuilder();
		StringBuilder qualSB = new StringBuilder();
		for (int i=0; i<call.length; i++){
			//add S, M, I
			if (call[i] == 'M' || call[i] == 'S' || call[i] == 'I'){
				seqSB.append(seq[i]);
				qualSB.append(Seq.ORDERED_ASCII[qual[i]]);
			}
		}
		sam.setSequence(seqSB.toString());
		sam.setQualities(qualSB.toString());
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

	public SamAlignment getSam() {
		return sam;
	}

	public char[] getSeq() {
		return seq;
	}

	public int[] getQual() {
		return qual;
	}

	public char[] getCall() {
		return call;
	}

	public int[] getPos() {
		return pos;
	}

}
