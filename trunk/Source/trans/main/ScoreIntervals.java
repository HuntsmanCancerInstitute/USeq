package trans.main;
import meme.*;
import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.*;

import util.bio.parsers.*;
import util.gen.*;

/**
 * Scores Intervals for the presence of transcription factor binding sites.
 */
public class ScoreIntervals {
	//fields
	private Interval[] intervals;
	private File flySeqDir;
	private File fastaSeqsFile;
	private File[] intervalFiles;
	private String[] bindingSites;	//aligned and trimmed
	private String chromName;
	private String chromSeq;
	private MultiFastaParser fastaParser;
	private MotifScanner motifScanner;
	private double cutOff = -1;
	private int[] numberMotifHits;	//use to make a histogram, index is the number of hits, int[] is the number of these found, int[100] is actually 100 or more.
	private int[] sizeOfIntervals;	//ditto
	private int totalNumberHits;
	private float totalNumHitsPerLength;
	private int clusterWindow =350;	//size in bp of window used to count motif hits
	private boolean scoreBestWindow = false;
	private boolean printAveHitsPerKb = false;
	
	public static void main(String[] args) {
		if (args.length < 3){
			printDocs();
			System.exit(0);
		}
		new ScoreIntervals(args);
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"**************************************************************************************\n" +
				"**                          Score Intervals: Jan 2005                               **\n" +
				"**************************************************************************************\n" +
				"SI scores serialized Interval[] arrays for the presence of transcription factor\n" +
				"binding sites. Use the following options:\n\n" +
				
				"-f Full path file text for the Interval[], if a directory is specified, all files\n" +
				"      within will be processed. (Required)\n" +
				"-g The full path directory text to the split genomic sequences (i.e. chr2L.fasta, \n"+
				"      chr3R.fasta...).  The file prefix names should match the chromosome names in the\n"+
				"      Intervals (ie chr2L, chr3R...). (Required)\n" +
				"-t Full path file text for the fasta file containing aligned trimmed examples of\n"+
				"      transcription factor binding sites.  A log likelihood position specific\n"+
				"      probability matrix will be generated from these sequences and used to scan all\n"+
				"      the Intervals for hits to the matrix.\n"+
				"-s Score cut off for the matrix. Defaults to the score of the lowest scoring sequence\n"+
				"      used in making the LLPSPM.\n"+
				"-w A window size in bp for calculating the maximum binding cluster, defaults to 350bp\n" +
				"-i Score the best average intensity window instead of the full interval.\n"+
				"-a Just print average hits per kb for best windows in interval array.\n"+
				"\n" +
				"Example: java -jar pathTo/T2/Apps/ScoreIntervals -f /my/affy/intervals/ -g\n" +
				"      /my/affy/DmelSeqs/ -t /my/affy/zeste.fasta -s 4.9 -w 375 -i\n"+
				"\n" +
			
		"**************************************************************************************\n");
	}
	/**This method will process each argument and assign any new varibles*/
	public void processArgs(String[] args){
		File directory = null;
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'f': directory = new File(args[i+1]); i++; break;
					case 'g': flySeqDir = new File(args[i+1]); i++; break;
					case 't': fastaSeqsFile = new File(args[i+1]); i++; break;
					case 's': cutOff = Double.parseDouble(args[i+1]); i++; break;
					case 'w': clusterWindow =Integer.parseInt(args[i+1]); i++; break;
					case 'i': scoreBestWindow=true; break;
					case 'a': printAveHitsPerKb=true; scoreBestWindow=true; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter request: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		//check to see if they entered required params
		if (directory==null || directory.exists() == false){
			System.out.println("\nEnter a serialized Interval[] array file or directory!\n");
			System.exit(0);
		}
		if (flySeqDir ==null || flySeqDir.exists()==false){
			System.out.println("\nEnter the full path text for the directory  containing split genomic sequences " +
					"(i.e. chr2L.fasta, chr3R.fasta...).  The file prefex names should match the chromosome names " +
			"in the Intervals (ie chr2L, chr3R...)\n");
			System.exit(0);
		}
		if (fastaSeqsFile == null || fastaSeqsFile.exists()==false){
			System.out.println("\nEnter a full path file text for the fasta file containing aligned trimmed examples of transcription factor binding sites.\n");
			System.exit(0);
		}
		
		//get files to process
		intervalFiles = IO.extractFiles(directory);	
		Arrays.sort(intervalFiles);
	}
	
	
	
	public ScoreIntervals(String[] args){
		processArgs(args);
		
		//get binding sites
		fastaParser = new MultiFastaParser(fastaSeqsFile);
		bindingSites = fastaParser.getSeqs();
		
		//make a MotifScanner
		double[][] PSPM = MemeMotif.makePSPM(bindingSites);
		double[][] llPSPM = MemeMotif.makeLLPSPM(PSPM, 0.25, 0.25, 0.25, 0.25);
		motifScanner = new MotifScanner(llPSPM);
		
		if (cutOff == -1) cutOff = motifScanner.findLowestScoringSeq(bindingSites);
		
		//for each interval file
		if (printAveHitsPerKb) System.out.println("Name\t#Intervals\tAveBestWindowHitsPerKb");
		
		for (int i=0; i<intervalFiles.length; i++){
			//get intervals
			intervals = (Interval[])IO.fetchObject(intervalFiles[i]);
			//sort intervals by chromosome to speed lookup
			HashMap chroms = new HashMap();
			int chromNumber = 0;
			for (int j = intervals.length-1; j>=0; j--){
				if (chroms.containsKey(intervals[j].getChromosome()) == false) chroms.put( intervals[j].getChromosome(), new Integer(chromNumber++) );
			}
			for (int j = intervals.length-1; j>=0; j--){
				intervals[j].setSortBy(((Integer)chroms.get(intervals[j].getChromosome())).intValue());
			}
			Arrays.sort(intervals);
			if (printAveHitsPerKb == false) System.out.println("Scanning "+intervalFiles[i].getName());
			
			//scan intervals, this updates each Interval obj with lots o stuff
			scanIntervals();
			
			if (printAveHitsPerKb){
				double[] hitsPerKb = new double[intervals.length];
				for (int j = intervals.length-1; j>=0; j--){
					double size = intervals[j].getBestWindow().getStartLastOligo()+25 - intervals[j].getBestWindow().getStart1stOligo();
					hitsPerKb[j] = 1000 * ((double)intervals[j].getNumberMotifHits()) / size;
				}
				
				System.out.println(intervalFiles[i].getName()+"\t"+intervals.length+"\t"+Num.mean(hitsPerKb));
			}
			else {
				//save modified intervals
				IO.saveObject(new File(intervalFiles[i].toString()+"Scr"),intervals);
				
				//print stats
				String word = "Interval";
				if (scoreBestWindow) word = "Best Window";
				int ave = Math.round(100.0f*(float)totalNumberHits/(float)intervals.length);
				System.out.println("\n"+totalNumberHits+ " Total number of hits above a cut off of "+cutOff);
				System.out.println(ave+"% of "+word+"s with one or more hits.  See histogram below.");
				
				System.out.println(totalNumHitsPerLength + " Total number of hits divided by "+word+" lengths.");
				float aveF = 100.0f*totalNumHitsPerLength/(float)intervals.length;
				System.out.println(aveF+"% of "+word+"s with one or more hits, using length penalty.");
				
			}
		}
	}
	
	
	public void scanIntervals(){
		totalNumberHits = 0;
		totalNumHitsPerLength = 0;
		MotifHit[] hits = new MotifHit[0];
		numberMotifHits = new int[101];
		sizeOfIntervals = new int[11];
		int numHits = 0;
		int seqLength = 0;
		int maxCluster;
		//initialize params
		String sequence;
		setGenomicSequenceParams(intervals[0]);
		int numIntervals = intervals.length;
		for (int i=0; i<numIntervals; i++){
			//check chromosome
			if (chromName.equals(intervals[i].getChromosome()) == false) setGenomicSequenceParams(intervals[i]);
			//get sequence, by some lucky fluke the Start and End return the exact sequence from 3.1, be sure to check
			//only holds when using the affy bpmap coordinates, when using to blast based coordinates subtract one 
			//new option, score window instead
			if (scoreBestWindow){
				sequence = new String(chromSeq.substring(intervals[i].getBestWindow().getStart1stOligo()-1, 
						intervals[i].getSizeOfOligoMinusOne()+intervals[i].getBestWindow().getStartLastOligo()-1));
			}
			else {				
				sequence = new String(chromSeq.substring(intervals[i].getStart1stOligo()-1, 
						intervals[i].getSizeOfOligoMinusOne()+intervals[i].getStartLastOligo()-1));
			}
			seqLength = sequence.length();
			hits = motifScanner.scoreSequence(cutOff,sequence);
			maxCluster = MotifScanner.calculateMaxCluster(hits, clusterWindow);
			
			//increment hit counter array
			numHits = hits.length;
			if (numHits>=100)numberMotifHits[100] +=1;
			else numberMotifHits[numHits] +=1;
			if (numHits !=0) {
				totalNumHitsPerLength += (float)numHits/(float)seqLength;
				totalNumberHits++;
			}
			
			//rescann full interval? need base scores
			if (scoreBestWindow){
				//re assign full sequence
				sequence = new String(chromSeq.substring(intervals[i].getStart1stOligo()-1, 
						intervals[i].getSizeOfOligoMinusOne()+intervals[i].getStartLastOligo()-1));
				//rescan to get base scores
				motifScanner.scoreSequence(cutOff,sequence);
			}
			
			//set fields in interval, always set full interval sequence
			intervals[i].setSequence(sequence);
			intervals[i].setNumberMotifHits(numHits);
			intervals[i].setMaxCluster(maxCluster);
			intervals[i].setBaseScores(motifScanner.getBaseScores());
			intervals[i].setBestWindowScored(scoreBestWindow);
			
			//save interval lengths
			seqLength = Math.round((float)intervals[i].getSequence().length()/1000.0f);
			if (seqLength>=10) sizeOfIntervals[10]+=1;
			else sizeOfIntervals[seqLength]+=1;
		}
	}
	
	/**Loads up the fasta file for a particular chromosome and sets some params.*/
	public void setGenomicSequenceParams(Interval interval){
		chromName = interval.getChromosome();
		if (printAveHitsPerKb == false) System.out.println("\tProcessing intervals from chromosome "+chromName);
		File chromFastaFile = new File (flySeqDir,chromName+".fasta");
		fastaParser.parseIt(chromFastaFile);
		chromSeq = fastaParser.getSeqs()[0];
	}
	
}
