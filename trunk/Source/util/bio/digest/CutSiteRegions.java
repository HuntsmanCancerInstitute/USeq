package util.bio.digest;

import util.bio.seq.*;
import util.bio.parsers.*;
import util.gen.*;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class CutSiteRegions {

	//fields
	Pattern cutSite;
	HashMap<String, File> fastas;
	int radius;
	int halfLength;
	String chromosomeName;
	String chromosomeSequence;
	PrintWriter out;
	File bedFile;
	
	//constructors
	public CutSiteRegions (String[] args){
		
		//set fields
		fastas = Seq.fetchChromosomeFastaFileHashMap(new File (args[0]));
		cutSite = Pattern.compile(args[1], Pattern.CASE_INSENSITIVE);
		halfLength = Math.round(((float)args[1].length())/2.0f);
		radius = Integer.parseInt(args[2]);
		bedFile = new File (args[3]);
		
		//for each chromosome
		try {
			out = new PrintWriter (new FileWriter (bedFile));
			for (String seq: fastas.keySet()){
				chromosomeName = seq;
				MultiFastaParser mfp = new MultiFastaParser(fastas.get(seq));
				chromosomeSequence = mfp.getSeqs()[0];
				System.out.println(chromosomeName+"\t0\t"+chromosomeSequence.length());
				makeRegions();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	//methods
	public void makeRegions(){
		Matcher m = cutSite.matcher(chromosomeSequence);
		int counter =0;
		while (m.find()) {
				int center = m.start() + halfLength;
				int start = center - radius;
				if (start < 0) start = 0;
				int end = center + radius; 
				out.println(chromosomeName +"\t"+start+"\t"+end);
				if (counter++ < 5) System.out.println("\tFound: "+m.group());
		}
	}
	
	
	public static void main(String[] args) {
		new CutSiteRegions(args);

	}

}
