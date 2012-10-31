package trans.roc;

import meme.MemeMotif;
import meme.MotifHit;
import meme.MotifScanner;

import java.io.*;
import java.util.*;
import java.util.regex.*;

import trans.anno.BindingRegion;
import trans.main.Window;
import util.bio.annotation.GeneGroup;
import util.bio.parsers.MultiFastaParser;
import util.gen.*;



/**
 * Scores windows for hits to a transcription factor PWM, for {@link RocWindow} filtering.
 */
public class PWMScoreWindows {
	//fields
	private String bpmapFile; 
	private String infoFile;
	private MultiFastaParser fastaParser;
	private MotifScanner motifScanner;
	private String[] bindingSites;	//aligned and trimmed
	private File flySeqDir;
	private File fastaSeqsFile;
	private String chromName;
	private String chromSeq;
	private double cutOff = 0;
	private int start = 1;
	private boolean startSelected = false;
	
	public PWMScoreWindows(String[] args){
		processArgs(args);
		
		System.out.println("Score Cut Off: "+cutOff);
		
		//get binding sites
		fastaParser = new MultiFastaParser(fastaSeqsFile);
		bindingSites = fastaParser.getSeqs();
		
		//make a MotifScanner
		double[][] PSPM = MemeMotif.makePSPM(bindingSites);
		double[][] llPSPM = MemeMotif.makeLLPSPM(PSPM, 0.25, 0.25, 0.25, 0.25);
		motifScanner = new MotifScanner(llPSPM);
		
		//get .bpmap file info
		ArrayList info = (ArrayList)IO.fetchObject(new File(infoFile));
		
		//for each chromosome, safe to parallel process on cluster
		for (int i=start; i<info.size(); i+=4){
			ArrayList allWindows = new ArrayList();
			
			//load chromosome sequence
			chromName = (String)info.get(i);
			setGenomicSequenceParams();
			
			//fetch testWindows to scan, the int[][] made in the WindowMaker
			System.out.println("\nTesting chromosome: "+chromName);
			int[][] testWindows = (int[][])IO.fetchObject(new File(bpmapFile+chromName+"Win"));
			
			//fetch bp positions to use in converting window indexes to real base pairs
			int[] positions = (int[])IO.fetchObject(new File(bpmapFile+chromName));
			
			int numWindows = testWindows.length;
			double numWinWithHits = 0;
			double totalHitDensity = 0;
			for (int j=0; j<numWindows; j++){
				//fetch sequence	
				int[] startStop = testWindows[j];
				String windowSequence = chromSeq.substring(positions[startStop[0]], positions[startStop[1]]+ 24);
				
				//score against matrix				
				double seqLength = windowSequence.length();
				MotifHit[] hits  = motifScanner.scoreSequence(cutOff,windowSequence);
				double hitsPerBP = ((double)hits.length)/seqLength;
				
				//record totals
				if (hitsPerBP !=0 ){
					numWinWithHits++;
					totalHitDensity += hitsPerBP;
				}
				//make windows
				RocWindow win = new RocWindow(chromName, positions[startStop[0]], positions[startStop[1]]+25, hitsPerBP);
				allWindows.add(win);
				
			}
			//print sum for chrom
			System.out.println("Total overlaping hit density   = "+totalHitDensity);
			System.out.println("Total number Windows with hits = "+numWinWithHits);
			System.out.println("Divided = "+totalHitDensity/numWinWithHits);
			IO.saveObject(new File(IO.getFullPathName(fastaSeqsFile)+chromName+"Win"), allWindows);
			if (startSelected) System.exit(0);
		}	
	}
	

	
	
	public static void main(String[] args) {
		if (args.length<4){
			printDocs();
			System.exit(0);
		}	
		new PWMScoreWindows(args);
	}		
	
	/**This method will process each argument and assign new varibles*/
	public void processArgs(String[] args){
		Pattern pat = Pattern.compile("-[a-z]");
		for (int i = 0; i<args.length; i++){
			String lcArg = args[i].toLowerCase();
			Matcher mat = pat.matcher(lcArg);
			if (mat.matches()){
				char test = args[i].charAt(1);
				try{
					switch (test){
					case 'c': flySeqDir = new File(args[i+1]); i++; break;
					case 'f': fastaSeqsFile = new File(args[i+1]); i++; break;
					case 's': cutOff = Double.parseDouble(args[i+1]); i++; break;
					case 'b': bpmapFile = new File(args[i+1],"bpmap.DupFree").toString(); i++; break;
					case 't': start = Integer.parseInt(args[i+1]); i++; startSelected = true; break;
					case 'h': printDocs(); System.exit(0);
					default: System.out.println("\nProblem, unknown option! " + mat.group());
					}
				}
				catch (Exception e){
					System.out.print("\nSorry, something doesn't look right with this parameter: -"+test);
					System.out.println();
					System.exit(0);
				}
			}
		}
		//check if directories are directories
		infoFile = bpmapFile+"Info";
		if (new File(infoFile).exists()==false){
			System.out.println("\nCheck your BPMapFiles directory, it appears to be incomplete!\n");
			System.exit(1);
		}
	}	
	
	/**Loads up the fasta file for a particular chromosome and sets some params.*/
	public void setGenomicSequenceParams(){
		File chromFastaFile = new File (flySeqDir,chromName+".fasta");
		fastaParser.parseIt(chromFastaFile);
		chromSeq = fastaParser.getSeqs()[0];
	}
	
	public static void printDocs(){
		System.out.println("\n" +
				"-t Start\n"+
				"-b BPMap\n"+
				"-c fly sequence directory\n" +
				"-s score cut off for pwm\n" +
				"-f fasta sequences for PWM\n" +
				"-t start at diff chroms 1, 5, 9, 13, 17, 21\n"
				
				);		
	}	
	
}
