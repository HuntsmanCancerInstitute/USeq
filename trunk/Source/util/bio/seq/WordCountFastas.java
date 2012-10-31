package util.bio.seq;
import java.io.*;
import java.util.*;
import util.gen.*;
import util.apps.MergeRegions;
import util.bio.parsers.*;

public class WordCountFastas {

	public static void main (String[] args){
		File fastaFile = new File (args[0]);
		MultiFastaParser fastas = new MultiFastaParser(fastaFile);
		String[] seqs = fastas.getSeqs();
		String[] names = fastas.getNames();

		int sizeOfWord = Integer.parseInt(args[1]);
		int maxNumber = Integer.parseInt(args[2]);

		HashMap<String, Integer> seqCounts = new HashMap<String, Integer>();


		for (int i=0; i< seqs.length; i++){
			int numWords = seqs[i].length()- sizeOfWord;
			for (int j=0; j< numWords; j++){
				int end = j+ sizeOfWord;
				String word = seqs[i].substring(j, end);
				Integer count = seqCounts.get(word);
				int num = 1;
				if (count != null) num = count.intValue() + 1;
				seqCounts.put(word, new Integer(num));
			}
		}

		//find bad words
		Iterator<String> it = seqCounts.keySet().iterator();
		HashSet<String> badWords = new HashSet<String>();
		while (it.hasNext()){
			String word = it.next();
			int count = seqCounts.get(word).intValue();
			if (count >= maxNumber) {
				badWords.add(word);
				System.out.println(word + "\t" + count);
			}
		}

		//make mask and print
		/*
		try{
			File results = new File (args[0]+".bed");
			PrintWriter out = new PrintWriter (new FileWriter(results));
			for (int i=0; i< seqs.length; i++){
				int numWords = seqs[i].length()- sizeOfWord;
				boolean[] mask = new boolean[seqs[i].length()];
				for (int j=0; j< numWords; j++){
					int stop = j+ sizeOfWord;
					String word = seqs[i].substring(j, stop);
					if (badWords.contains(word)){
						for (int k=j; k< stop; k++) mask[k] = true;
					}
				}
				MergeRegions.print(names[i], mask, out);
			}
		}catch (Exception e){
			e.printStackTrace();
		}*/


		//print good sequences
		if (maxNumber ==0) System.exit(0);
		int numGoodSeqs = 0;
		for (int i=0; i< seqs.length; i++){
			int numWords = seqs[i].length()- sizeOfWord;
			boolean printMe = true;
			for (int j=0; j< numWords; j++){
				int end = j+ sizeOfWord;
				String word = seqs[i].substring(j, end);
				if (badWords.contains(word)) {
					printMe = false;
					break;
				}
			}
			if (printMe) {
				System.out.println(names[i]);
				System.out.println(seqs[i]);
				numGoodSeqs++;
			}
		}
		
		System.out.println(seqs.length+"\tNumber starting seqs");
		System.out.println(numGoodSeqs+"\tNumber filtered seqs");
		 


	}

}
