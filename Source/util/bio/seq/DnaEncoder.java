package util.bio.seq;

import java.util.HashMap;

import util.gen.IO;

public class DnaEncoder {

	public static void main(String[] args) {
		// 4 codons produce 64 possible sequences
		String letters = "qwertyuioplkjhgfdsazxcvbnmMNBVCXZASDFGHJKLPOIUYTREWQ1234567890\t-";
		IO.pl(letters.length());

		HashMap<Character,String> letterToMer = new HashMap<Character, String>();
		HashMap<String, Character> merToLetter = new HashMap<String, Character>();
		
		String[] bases = new String[] {"g","a","t","c"};
		int letterCounter = 0;
		toBreak:
		for (int i=0; i< 4; i++) {
			for (int j=0; j< 4; j++) {
				for (int k=0; k< 4; k++) {
					for (int l=0; l< 4; l++) {
						StringBuilder mer = new StringBuilder(bases[i]);
						mer.append(bases[j]);
						mer.append(bases[k]);
						mer.append(bases[l]);
						Character c = letters.charAt(letterCounter++);
						String dna = mer.toString();
						letterToMer.put(c,dna);
						merToLetter.put(dna, c);
						if (letters.length() == letterCounter) break toBreak;
					}
				}
			}
		}
		IO.pl(letterToMer);
		IO.pl(merToLetter);

	}

}
