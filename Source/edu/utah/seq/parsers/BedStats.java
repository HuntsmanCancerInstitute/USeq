
package edu.utah.seq.parsers;

import java.io.*;
import java.util.regex.*;

import edu.utah.seq.data.*;
import edu.utah.seq.useq.data.RegionScoreText;
import util.gen.*;

import java.util.*;


import util.bio.annotation.*;
import util.bio.seq.Seq;

/**Calculates various statistics on the length and sequence composition of the reads in bed file format.
 * @author david.nix@hci.utah.edu 
 **/
public class BedStats {
	
	//fields
	private File[] bedFiles;
	private int[][] countData;
	private Histogram histogram = new Histogram(10,100, 100-10);
	private int trim = 0;
	private int secondBase = -1;
	private boolean reverseComplementSequence = false;
	
	
	//constructors
	public BedStats(String[] args){
		long startTime = System.currentTimeMillis();
		processArgs(args);
		histogram.setTrimLabelsToSingleInteger(true);
		
		//make array to hold data, GATCN
		countData = new int[100][10];
		for (int i=0; i< countData.length; i++){
			countData[i] = new int[10];
		}
		
		loadCountData();
		
		//make histogram of data
		System.out.println("A histogram of read lengths:");
		histogram.printScaledHistogram();
		System.out.println((int)histogram.getTotalBinCounts()+" total reads\n");
		
		//print base stats
		printBaseUsage();
		
		//finish and calc run time
		double diffTime = ((double)(System.currentTimeMillis() -startTime))/1000;
		System.out.println("\nDone! "+Math.round(diffTime)+" seconds\n");
	}
	
	public void printBaseUsage(){
		System.out.println("f= first 5' base; l= last 3' base");
		System.out.println("ReadLength\t#fG\t#fA\t#fT\t#fC\t#fN\t#lG\t#lA\t#lT\t#lC\t#lN\t%fG\t%fA\t%fT\t%fC\t%fN\t%lG\t%lA\t%lT\t%lC\t%lN");
		for (int i=0; i< countData.length; i++){
			int[] gatcn = countData[i];
			double totalCounts = Num.sumIntArray(gatcn);
			if (totalCounts !=0){
				System.out.print(i);
				for (int x=0; x< gatcn.length; x++){
					System.out.print("\t");
					System.out.print(gatcn[x]);
				}
				double totalFirst = 0;
				for (int x=0; x< 5; x++) totalFirst+= gatcn[x];
				double totalSecond = 0;
				for (int x=5; x< 10; x++) totalSecond+= gatcn[x];
				for (int x=0; x< 5; x++){
					System.out.print("\t");
					System.out.print(Num.formatPercentOneFraction((double)gatcn[x]/totalFirst));
				}
				for (int x=5; x< 10; x++){
					System.out.print("\t");
					System.out.print(Num.formatPercentOneFraction((double)gatcn[x]/totalSecond));
				}
				System.out.println();
			}
		}
	}
	
	public void loadCountData(){
		for (int i=0; i< bedFiles.length; i++){
			System.out.println("\t"+bedFiles[i]);
			String line = null;
			PrintWriter out = null;
			try {
				if (trim !=0){
					String name = Misc.removeExtension(bedFiles[i].getName());
					File trimmed = new File(bedFiles[i].getParentFile(), name+"_trimmed"+trim+".bed");
					out = new PrintWriter (new FileWriter (trimmed));
				}
				BufferedReader in = IO.fetchBufferedReader(bedFiles[i]);
				String[] tokens;
				Pattern tab = Pattern.compile("\\t");
				String tabString = "\t";
				while ((line=in.readLine()) != null){
					//break by tab
					tokens = tab.split(line);
					//split on first semicolon
					String seq;
					int index = tokens[3].indexOf(";");
					if (index == -1){
						//System.err.println("Error: bad line, skipping -> "+line);
						//continue;
						seq = tokens[3].toLowerCase();
					}
					//read sequence, not reference sequence
					else seq = tokens[3].substring(0, index).toLowerCase();
					if (reverseComplementSequence) seq = Seq.reverseComplementDNA(seq);
					
					//stats
					int seqLength = seq.length();
					histogram.count(seqLength);
					//set 5' base, gatcn
					int[] bases = countData[seqLength];
					char firstBase = seq.charAt(0);
					if (firstBase == 'g') bases[0]++;
					else if (firstBase == 'a') bases[1]++;
					else if (firstBase == 't') bases[2]++;
					else if (firstBase == 'c') bases[3]++;
					else bases[4]++; //n
					//last base
					char lastBase;
					if (secondBase == -1) lastBase = seq.charAt(seqLength-1);
					else lastBase = seq.charAt(secondBase);
					if (lastBase == 'g') bases[5]++;
					else if (lastBase == 'a') bases[6]++;
					else if (lastBase == 't') bases[7]++;
					else if (lastBase == 'c') bases[8]++;
					else bases[9]++; //n
					//write out?
					if (out !=null) {
						int start = Integer.parseInt(tokens[1]);
						int stop = Integer.parseInt(tokens[2]);
						if (tokens.length !=6 || tokens[5].equals("+")) stop = start + trim;
						else start = stop - trim;
						seq = seq.substring(0, trim);
						StringBuilder sb = new StringBuilder();
						//chrom
						sb.append(tokens[0]); sb.append(tabString);
						//start
						sb.append(start);sb.append(tabString);
						//stop
						sb.append(stop); sb.append(tabString);
						//name
						sb.append(seq); sb.append(tabString);
						//score
						if (tokens.length>=5) sb.append(tokens[4]); sb.append(tabString);
						//strand
						if (tokens.length ==6) sb.append(tokens[5]);
						sb.append("\n");
						out.print(sb);
					}
				}
				if (out != null) out.close();
				in.close();
			} catch (Exception e) {
				System.err.println("\nBad line -> "+line);
				e.printStackTrace();
			}	
		}
	}
	

	
	public static void main(String[] args) {
		if (args.length ==0){
			printDocs();
			System.exit(0);
		}
		new BedStats(args);
	}		
	
	/**This method will process each argument and assign new variables*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		System.out.println("\n"+IO.fetchUSeqVersion()+" Arguments: "+Misc.stringArrayToString(args, " ")+"\n");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'b': bedFiles = IO.extractFiles(new File(args[++i])); break;
					case 't': trim = Integer.parseInt(args[++i]); break;
					case 's': secondBase = Integer.parseInt(args[++i]); break;
					case 'r': reverseComplementSequence = true; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					Misc.printExit("\nSorry, something doesn't look right with this parameter: -"+test+"\n");
				}
			}
		}
		if (bedFiles == null || bedFiles.length ==0) Misc.printErrAndExit("Cannot find your xxx.bed(.gz/.zip OK) files?");
	}	

	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                                BedStats: June 2010                               **\n" +
				"**************************************************************************************\n" +
				"Calculates several statistics on bed files where the name column contains a short read\n" +
				"sequence. This includes a read length distribution and frequencies of the 1st and last\n" +
				"bps. Can also trim your read to a particular length. \n" +

				"\nOptions:\n"+
				"-b Full path file name for your alignment bed file or directory containing such. The\n" +
				"       name column should contain your just you sequence or seq;qual .\n" +
				"-t Trim the 3' ends of your reads to the indicated length, defaults to not trimming.\n"+
				"-s Calculate base frequencies for the given 0 indexed base instead of the last base.\n"+
				"-r Reverse complement sequences before calculating stats and trimming.\n"+

				"\nExample: java -Xmx1500M -jar pathToUSeq/Apps/BedStats -b /Res/ex1.bed.gz -s 9 -t 10\n\n" +

		"**************************************************************************************\n");

	}	

}
