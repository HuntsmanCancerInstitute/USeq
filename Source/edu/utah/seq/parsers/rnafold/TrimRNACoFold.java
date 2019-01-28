package edu.utah.seq.parsers.rnafold;

import java.io.IOException;
import java.util.regex.Pattern;

import util.gen.IO;
import util.gen.Misc;

public class TrimRNACoFold {

	//fields
	private String seq = "CUGAGAUGUCGUGUGUGU&ACACUACAUAUUUAUUGUUUAUUUAUACCUGAUAAUGCUGCUUUAAUCACACACGAGGUAUCGCUUCUACUA";
	private String fold = ".(((.((.(((((((((.&............(((((((...........)))))))..........))))))))).)).))).........";
	private int gapDotCount = 3;
	private Pattern and = Pattern.compile("&");
	
	public static void main (String[] args) throws IOException{
		new TrimRNACoFold(args);
	}
	
	public TrimRNACoFold(String[] args) throws IOException{
		//check lengths
		if (seq.length() != fold.length()) throw new IOException("Mismatched length between "+seq+" and "+fold);
		
		//split
		String[] seqs = and.split(seq);
		String[] folds = and.split(fold);
		
		//trim
		for (int i=0; i< seqs.length; i++){
			IO.pl();
			IO.pl(seqs[i] +"\n"+ folds[i]);
			String[] trimmed = trimEdges(seqs[i], folds[i]);
			IO.pl(trimmed[0] +"\n"+ trimmed[1]);
			IO.pl();
			String[] trimInternal = trimInternal(trimmed[0], trimmed[1]);
			
		}
		
	}
	
	private String[] trimInternal(String seq, String fold) {
		int left = 1;
		int right = findNextRightBracket(fold);
		
		if (right == 3) {
			IO.pl("internal no trim "+ left+" "+right);
			return new String[]{seq, fold};
		}
		IO.pl("trim "+ left+" "+right);
		String trimmedSeq = seq.substring(left, right);
		String trimmedFold = fold.substring(left, right);
		IO.pl(trimmedSeq+"\n"+trimmedFold);
		return new String[]{trimmedSeq, trimmedFold};
	}

	private int findNextRightBracket(String fold) {
		int index = 2;
		boolean rightFound = false;
		for (int i=2; i< fold.length(); i++){
			if (fold.charAt(i) == ')') {
				rightFound = true;
				index = i;
			}
			else if (rightFound) break;
			
		}
		return ++index;
	}

	private String[] trimEdges(String seq, String fold){
		int left = findLeftIndex(fold);
		int right = findRightIndex(fold);
		//IO.pl(left+" "+right);
		String trimmedSeq = seq.substring(left, right);
		String trimmedFold = fold.substring(left, right);
		return new String[]{trimmedSeq, trimmedFold};
	}

	private int findLeftIndex(String fold) {
		int leftIndex = 0;
		for (int i=0; i<fold.length(); i++){
			if (fold.charAt(i) == '.') leftIndex = i;
			else break;
		}
		return leftIndex;
	}
	
	private int findRightIndex(String fold) {
		int rightIndex = 0;
		for (int i=fold.length()-1; i>=0; i--){
			if (fold.charAt(i)=='.') rightIndex = i;
			else break;
		}
		return ++rightIndex;
	}
	
	
}
