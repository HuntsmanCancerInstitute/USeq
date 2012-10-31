package util.bio.converters;
import java.io.*;
import java.util.*;

import util.bio.seq.*;
import util.gen.*;

/**
 * Counts and parses the bases in a multiple alignment where each text is broken into three sections
 * like gatcc GGCCAATT cgatcca, ccat GGCCAATT cctttt, like the output from meme
 * the central portion is aligned, the left and right flanks are not.
 */
public class AlignedBaseCounter {
	public static void main(String[] args) {
		if (args.length!=4){
			System.out.println("\nTo count the bases, and parse cut meme results oligos, enter 4 params:\n"+
					"\t1) reverse complement motif? enter yes or no\n" +
					"\t2) left cut number? (ie 5)\n" +
					"\t3) right right number? (ie 19)\n" +
					"\t4) full path file text, just seq lines from meme report.\n");
			System.exit(0);
		}
	//parse params
	File file = new File(args[3]);
	boolean reverseComplement = args[0].equalsIgnoreCase("yes");
	int leftTweak=0;
	int rightTweak=0;
	try {
		leftTweak = Integer.parseInt(args[1]);
		rightTweak = Integer.parseInt(args[2]);
	}catch(NumberFormatException e){
		System.out.println("\nEnter integers for left and right teaking!\n");
		System.exit(0);
	}
		
	//read in all seqs
	String line;
	String[] oligos = null;
	try {
		if (file.exists()==false) {
			System.out.println("\nCannot find your file!\n");
			System.exit(0);
		}
		ArrayList oligosAL = new ArrayList();
		BufferedReader in = new BufferedReader(new FileReader(file));
		String[] tokens;
		while ((line = in.readLine()) !=null) {
			line = line.trim();                     //kill whitespace and test if exists
			if (line.length() == 0) continue;       //skip blank lines
			line = line.replaceAll("\\.","n");	  //replace . with n
			tokens = line.split("\\s+");
			if (tokens.length==7) line = tokens[4]+"\t"+tokens[5]+"\t"+tokens[6];
			else line = tokens[4]+"\t"+tokens[5]+"\tn";
			oligosAL.add(line);
		}
		in.close();
		int len = oligosAL.size();
		oligos = new String[len];
		oligosAL.toArray(oligos);
	}
	catch (IOException e) {
		e.printStackTrace();
	}
	
	//get max sizes of each section
	int len = oligos.length;
	int left =0;
	int right=0;
	for (int i=0; i<len; i++){
		String[] sections = oligos[i].split("\\t");
		if (left<sections[0].length()) left = sections[0].length();
		if (right<sections[2].length()) right = sections[2].length();
	}
	

	
	//fill empty spaces with n's
	System.out.println("\nUntrimmed oligos for a seq logo:");
	String[] ns = {"", "n", "nn", "nnn", "nnnn", "nnnnn", "nnnnnn", "nnnnnnn", "nnnnnnnn", "nnnnnnnnn", "nnnnnnnnnn","nnnnnnnnnnn","nnnnnnnnnnnn","nnnnnnnnnnnnn","nnnnnnnnnnnnnn","nnnnnnnnnnnnnnn","nnnnnnnnnnnnnnnn"};
	for (int i=0; i<len; i++){
		String[] sections = oligos[i].split("\\t");
			//left
			if (sections[0].length()<left) sections[0] = ns[left-sections[0].length()]+sections[0];
			if (sections[2].length()<right) sections[2] = sections[2] + ns[right-sections[2].length()];
			oligos[i]= sections[0]+sections[1]+sections[2];
			
			if (reverseComplement) oligos[i]= Seq.reverseComplementDNA(oligos[i]);
			System.out.println(oligos[i]);
			
			oligos[i]= (String)oligos[i].subSequence(leftTweak,rightTweak); //tweak to limit
		}
	
	//count
	int maxSize = oligos[0].length();
	int[] g = new int[maxSize];
	int[] a = new int[maxSize];
	int[] t = new int[maxSize];
	int[] c = new int[maxSize];
	char Aye = 'a';
	char Gee = 'g';
	char Cee = 'c';
	char Tea = 't';
	
	for (int i=0; i<len; i++){
		//walk across each oligo
		String oligo = oligos[i].toLowerCase();
		for (int j=0; j<maxSize; j++){
			char test = oligo.charAt(j);
			if (test==Aye) a[j]++;
			else if (test==Gee) g[j]++;
			else if (test==Cee) c[j]++;
			else if (test==Tea) t[j]++;
		}
	}
	
	//print the oligos
	System.out.println ("*** Parsed Oligos ***\nFasta Format...\n");
	for (int i=0; i<len; i++){
		System.out.println(">oligo "+(i+1)+"\n"+oligos[i]);
		//System.out.println(">oligo "+(i+1)+"\n"+oligos[i].subSequence(3,17));
	}
	
	System.out.println("\n*** Count Matrix ***\n");
	
	System.out.println("TransFac Style...");
	System.out.println("NA\tMotif \nXX \nDE\tMotif \nXX \nBF\tMotif \nXX \nP0\tA\tC\tG\tT");
	
	for (int i=0; i<maxSize; i++){
		System.out.println("0"+i+"\t"+a[i]+"\t"+c[i]+"\t"+g[i]+"\t"+t[i]+"\t X");
	}
	
	
	}
	
}
