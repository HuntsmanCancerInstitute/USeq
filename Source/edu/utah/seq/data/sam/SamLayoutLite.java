package edu.utah.seq.data.sam;

import htsjdk.samtools.SAMRecord;
import java.util.regex.Matcher;
import util.gen.Misc;
import util.gen.Num;

public class SamLayoutLite{
	private char[] seq;
	private char[] qual;
	private char[] call;
	private int[] pos;
	private SamAlignment samAligment;
	private SAMRecord samRecord;
	
	public SamLayoutLite(SamAlignment samAlignment){
		this.samAligment = samAlignment;
		int size = SamAlignment.countLengthOfCIGAR(samAligment.getCigar());
		seq = new char[size];
		qual = new char[size];
		call = new char[size];
		pos = new int[size];
		layoutCigar(samAligment.getCigar(), samAligment.getUnclippedStart(), samAligment.getSequence(), samAligment.getQualities());
	}
	
	public SamLayoutLite(SAMRecord samRecord){
		this.samRecord = samRecord;
		int size = SamAlignment.countLengthOfCIGAR(samRecord.getCigarString());
		seq = new char[size];
		qual = new char[size];
		call = new char[size];
		pos = new int[size];
		layoutCigar(samRecord.getCigarString(), samRecord.getUnclippedStart()-1, samRecord.getReadString(), samRecord.getBaseQualityString());
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
			SamLayoutLite sl = new SamLayoutLite(s);
			sl.print();
		} catch (Exception e) {
			e.printStackTrace();
		} 
	}
	

	
	public boolean changeBase(int bpPosition, char reference, char mutation, boolean warn){
//System.out.println("SNV: "+bpPosition+" "+reference+" "+mutation);		
		int index = Num.findSmallestIndexToValue(pos, bpPosition);
		if (index < 0 || index >= seq.length) {
			if (warn) System.err.println("WARNING index outside bounds of seq for snv "+index);
			return false;
		}
		//check for reference, skip if not, represents a bad base call or mis alignment artifact
		if (seq[index] != reference) {
			if (warn) System.err.println("WARNING ref does not match aligned base "+seq[index]);
			return false;
		}
		//check that the call is an M 
		if (call[index] != 'M') {
			if (warn) System.err.println("WARNING not an M base "+call[index]);
			return false;
		}
		seq[index] = mutation;
		return true;
	}
	
	/**This replaces the refBases with the altBases and changes the difference to D but does not actually modify the layout length.*/
	public boolean markDeletion(int startPos, String refBases, String altBases, boolean warn) {
		//System.out.println("DELETION: "+startPos+" "+refBases+" "+altBases);
		int index = Num.findSmallestIndexToValue(pos, startPos);
		
		if (index < 0 || index >= seq.length) {
			//sometimes the starting position of the deletion is before the start of the alignment so delete one until ref is empty
			if (index == -1 && refBases.length()!=0){
				if (warn) System.err.println("WARNING index outside bounds of seq for deletion, attempting 1bp advance for "+startPos);
				//add one to pos
				int pos = startPos+1;
				String ref = refBases.substring(1);
				String alt = "";
				if (altBases.length()!=0) alt = altBases.substring(1);
				return markDeletion(pos, ref, alt, warn);
			}
			return false;
		}
		
		int endIndex = refBases.length() + index;
		if (endIndex > seq.length) endIndex = seq.length;
		
		//change seq
		int stop = index + altBases.length();
		int counter = 0;
		boolean changed=false;
		for (int i= index; i< stop; i++){
			char a = altBases.charAt(counter++);
			if (seq[i] != a) {
				seq[i] = a;
				changed = true;
			}
		}
		
		//change for difference, note the pos[] is not updated!
		for (int i=stop; i< endIndex; i++){
			call[i] = 'D';
			seq[i] = 'D';
			qual[i] = '!';
			changed = true;
		}
		return changed;
	}
	
	/**This swaps the alt bases for the refBases and expands the layout arrays, except for the pos[] which isn't updated.*/
	public boolean insertBases(int startPos, String refBases, String altBases, boolean warn) {
//System.out.println("INSERTION: "+startPos+" "+refBases+" "+altBases);
		int index = Num.findSmallestIndexToValue(pos, startPos);
		if (index < 0 || index >= seq.length) {
			if (warn) System.err.println("WARNING index outside bounds of seq for insertion "+index);
			//attempt to trunk ref?
			if (refBases.length() >1){
				if (warn) System.err.println("Clipping refBases");
				refBases = refBases.substring(1);
				altBases = altBases.substring(1);
				startPos++;
				return insertBases(startPos, refBases, altBases, warn);
			}
			return false;
		}

		//for each refBase, swap it with the alt base, alt length always bigger, could be a compound snv insertion  e.g.  CC -> TACGG
		int endIndex = refBases.length() + index;
		if (endIndex > seq.length) endIndex = seq.length;
		int counter=0;
		boolean changed = false;
		for (int i= index; i< endIndex; i++) {
			char a = altBases.charAt(counter++);
			if (seq[i] != a){
				changed = true;
				seq[i] = a;
			} 
		}
		
		int indexToInsert=index+refBases.length();
		String altToInsert= altBases.substring(refBases.length());
//System.out.println("Inserting "+altToInsert+" at "+indexToInsert);
		
		//watch out for making an insertion past the last base
		if (indexToInsert >= seq.length) return changed;
		
		//insert it, note the pos[] is not updated!
		seq = Misc.insertString(seq, altToInsert, indexToInsert);
		String q30s = Misc.concatinateRepeats("?", altToInsert.length());
		qual =  Misc.insertString(qual, q30s, indexToInsert);
		String iis = Misc.concatinateRepeats("I", altToInsert.length());
		call = Misc.insertString(call, iis, indexToInsert);
		
		return true;
	}
	
	
	
	

	public void layoutCigar(String cigar, int unclippedStart, String sequ, String quali){
		//for each cigar block in first, looking for MDIN s
		Matcher mat = SamAlignment.CIGAR_SUB.matcher(cigar);
		int position = unclippedStart;
		int index = 0;
		int layoutIndex = 0;
		char[] samSeq = sequ.toUpperCase().toCharArray();
		char[] samQual = quali.toCharArray();
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
					qual[layoutIndex] = '!';
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
					qual[layoutIndex] = '!';
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
					qual[layoutIndex] = '!';
					call[layoutIndex] = 'H';
					pos[layoutIndex] = position;
					position++;
					layoutIndex++;
				}
			}
			else Misc.printErrAndExit("\nError: unsupported CIGAR! Should never hit this, contact admin!\n" +cigar);
		}
	}
	
	public String[] getSequenceAndQualtities(){
		StringBuilder s = new StringBuilder();
		StringBuilder q = new StringBuilder();
		//walk and add in S M I but not D H N
		for (int i=0; i< call.length; i++){
			if (call[i]=='M' || call[i]=='S' || call[i]=='I'){
				s.append(seq[i]);
				q.append(qual[i]);
			}
		}
		return new String[]{s.toString(), q.toString()};
	}
	
	/**Returns qualities where the CIGAR is M or I*/
	public String getMIQualities(){
		StringBuilder q = new StringBuilder();
		//walk and add in M I but not D H N S
		for (int i=0; i< call.length; i++){
			if (call[i]=='M' || call[i]=='I') q.append(qual[i]);
		}
		return q.toString();
	}
	
	public String getSequence(){
		StringBuilder s = new StringBuilder();
		//walk and add in S M I but not D H N
		for (int i=0; i< call.length; i++){
			if (call[i]=='M' || call[i]=='S' || call[i]=='I'){
				s.append(seq[i]);
			}
		}
		return s.toString();
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
