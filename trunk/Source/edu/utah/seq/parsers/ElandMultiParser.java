
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;
import util.bio.parsers.*;
import edu.utah.seq.data.*;
import util.gen.*;

import java.util.*;

/**Parses an Eland xxx.eland_multi.txt file tabulating hits to each fasta entry.  Good for scoring hits to a transcriptome where
 * every fasta entry represents a different gene.
 * @author david.nix@hci.utah.edu 
 **/
public class ElandMultiParser {
	//fields
	private File[] dataFiles;
	private File resultsFile;
	private HashMap<String, Multi> map = new HashMap<String, Multi>();
	private int numNM = 0;
	private int numQC = 0;
	private int numRM = 0;
	private int numTM = 0;
	private int numAligningReads;
	private Pattern tab = Pattern.compile("\\t+");
	private Pattern comma = Pattern.compile(",");
	private Pattern divider = Pattern.compile("/");

	//constructors
	public ElandMultiParser(String[] args){
		processArgs(args);
		
		//for each dataFile parse and load HashMap
		System.out.println("Parsing and loading ...");
		for (int i=0; i< dataFiles.length; i++) {
			System.out.println("\t"+dataFiles[i].getName());
			parseAndLoad (dataFiles[i]);
		}
		
		//print results
		printOutMap();
		
		//print out summary info
		System.out.println("\nNumber of sequences...");
		System.out.println("\tSucessfully aligned\t"+numAligningReads);
		System.out.println("\tToo many alignments\t"+numTM);
		System.out.println("\tFailed QC\t"+numQC);
		System.out.println("\tNo match\t"+numNM);
		System.out.println("\tRepeat masked\t"+numRM);
		System.out.println("\nDone!\n");
	}
	
	public void printOutMap(){
		try {
			Iterator<String> it = map.keySet().iterator();
			PrintWriter out = new PrintWriter ( new FileWriter (resultsFile));
			out.println("Name\tFO\tF1\tF2\tR0\tR1\tR2\t#_Aligned="+numAligningReads);
			while (it.hasNext()){
				String name = it.next();
				out.println(name+"\t"+map.get(name));
			}
			out.close();
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	/*>HWI-EAS240_SR18_204L8AAXX:8:1:766:360	TAGGTTGTCTAAAAATA	2:6:40	t2.fa/uc004eal.1_chrX:71280758-71298324:1241F1,t4.fa/uc002rro.1_chr2:39329925-39517723:9704F1,10699F1,t4.fa/uc002rrp.1_chr2:39329925-39517723:9458F1,10453F1,t4.fa/uc003dsl.1_chr3:99592199-99593165:218R1,t1.fa/uc001uzy.1_chr13:44809303-44813297:215F0,t1.fa/uc001uzz.1_chr13:44809303-44813297:215F0*/
	public void parseAndLoad(File multiFile){
		String line = null;
		try{
			BufferedReader in;
			if (multiFile.getName().endsWith(".zip"))in = IO.fetchReaderOnZippedFile(multiFile);
			else in = new BufferedReader ( new FileReader (multiFile));
			while ((line = in.readLine()) !=null){
				line = line.trim();
				if (line.length() == 0 || line.startsWith("#")) continue;
				String[] items = tab.split(line);
				if (items.length < 3) continue;
				//what kind of match
				if (items[2].contains(":")){
					//alignments found, too many?
					if (items.length == 3) numTM++;
					else{
						//good alignments found, parse
						numAligningReads++;
						loadHits(items[3]);
					}
				}
				else if (items[2].equals("NM")) numNM++;
				else if (items[2].equals("QC")) numQC++;
				else if (items[2].equals("RM")) numRM++;
				else Misc.printErrAndExit("\nProblem parsing match type from "+line);
			}
			in.close();
		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nProblem parsing and loading \n\t-> "+multiFile+"\n\t-> "+line);
		}
	}

	/*
	 t2.fa/uc004eal.1_chrX:71280758-71298324:1241F1
		t4.fa/uc002rro.1_chr2:39329925-39517723:9704F1
		10699F1
		t4.fa/uc002rrp.1_chr2:39329925-39517723:9458F1
		10453F1
		t4.fa/uc003dsl.1_chr3:99592199-99593165:218R1
	t1.fa/uc001uzy.1_chr13:44809303-44813297:215F0
	t1.fa/uc001uzz.1_chr13:44809303-44813297:215F0
	 */
	public void loadHits(String alignments){
		String[] hits = comma.split(alignments);
		String name = null;
		for (int i=0; i< hits.length; i++){
			//any dividers?
			int dividerIndex = hits[i].indexOf('/');
			if (dividerIndex != -1){
				//split on last ':'
				int index = hits[i].lastIndexOf(':');
				name = hits[i].substring(dividerIndex+1, index);
			}
			//parse match type
			String matchType = hits[i].substring(hits[i].length()-2);
			//find and modify Multi object
			Multi multi;
			if (map.containsKey(name)) multi = map.get(name);
			else {
				multi = new Multi();
				map.put(name, multi);
			}
			//modify
			if (matchType.equals("F0")) multi.f0++;
			else if (matchType.equals("F1")) multi.f1++;
			else if (matchType.equals("F2")) multi.f2++;
			else if (matchType.equals("R0")) multi.r0++;
			else if (matchType.equals("R1")) multi.r1++;
			else if (matchType.equals("R2")) multi.r2++;
			else Misc.printErrAndExit("\nError parsing match type! for "+Misc.stringArrayToString(hits, ","));
		}
	}

	private class Multi {
		int f0 = 0;
		int f1 = 0;
		int f2 = 0;
		int r0 = 0;
		int r1 = 0;
		int r2 = 0;
		StringBuilder seqs = new StringBuilder();

		public String toString(){
			return f0 +"\t"+ f1 +"\t"+ f2 +"\t"+ r0 +"\t"+ r1 +"\t"+ r2 +"\t"+ seqs;
		}

	}

	/**Returns number space (percent%)*/
	public static String formatFraction(int number, int total, int numberDecimals){
		double fraction = 100.0 *(double)number / (double)total ;
		String fractionString = Num.formatNumber(fraction, numberDecimals);
		return fractionString+"%\t"+number;
	}


	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new ElandMultiParser(args);
	}		

	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		File forExtraction = null;
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': forExtraction = new File(args[i+1]); i++; break;
					case 'r': resultsFile = new File (args[i+1]); i++; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}

		//pull files
		File[][] tot = new File[2][];
		tot[0] = IO.extractFiles(forExtraction,"eland_multi.txt");
		tot[1] = IO.extractFiles(forExtraction,"eland_multi.txt.zip");
		dataFiles = IO.collapseFileArray(tot);
		if (dataFiles == null || dataFiles.length ==0 || dataFiles[0].canRead() == false) Misc.printExit("\nError: cannot find your xxx_eland_multi.txt(.zip) file(s)!\n");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Eland Multi Parser: October 2008                        **\n" +
				"**************************************************************************************\n" +
				"Parses an Eland xxx.eland_multi.txt alignment file tabulating hits to each fasta entry.\n" +
				"Good for scoring hits to a transcriptome where every fasta entry represents a\n" +
				"different gene.\n\n" +

				"-f The full path directory/file text of your xxx.eland_multi.txt(.zip) file(s). Files\n" +
				"      will be merged.\n" +
				"-r Full path file text for saving the results.\n"+

				"\n\nExample: java -Xmx1500M -jar pathToUSeq/Apps/ElandMultipParser -f \n" +
				"      /data/MultiFiles/ -r /data/transcriptomeResults.xls \n\n" +

		"**************************************************************************************\n");

	}	

}
