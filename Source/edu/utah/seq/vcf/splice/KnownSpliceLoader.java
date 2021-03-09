package edu.utah.seq.vcf.splice;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.regex.Matcher;
import edu.utah.seq.its.Interval1D;
import edu.utah.seq.its.IntervalST;
import edu.utah.seq.mes.MaxEntScanScore3;
import edu.utah.seq.mes.MaxEntScanScore5;
import edu.utah.seq.vcf.xml.foundation.SimpleVcf;
import htsjdk.samtools.reference.ReferenceSequence;
import util.bio.annotation.ExonIntron;
import util.bio.parsers.UCSCGeneLine;
import util.bio.seq.Seq;
import util.gen.Gzipper;
import util.gen.Histogram;
import util.gen.Misc;
import util.gen.Num;


/**This is a loader for splice junction info for use with the VCFSpliceAnnotator app.*/
public class KnownSpliceLoader implements Runnable {

	//fields
	private KnownSpliceJunctionScanner vsa;
	private MaxEntScanScore3 score3;
	private MaxEntScanScore5 score5;
	private boolean failed = false;
	private String workingSequence;
	private Gzipper bedOut;
	private String chromName;
	private UCSCGeneLine[] transcripts;
	private Histogram spliceHistogram5;
	private Histogram spliceHistogram3;

	//counters to upload to main thread
	private int numUniqueJunctions = 0;

	//constructor
	public KnownSpliceLoader(KnownSpliceJunctionScanner vsa, String chromName, File tempBed) {
		try {
			this.vsa = vsa;
			this.chromName = chromName;

			//start up mes
			score5 = new MaxEntScanScore5(vsa.getSpliceModelDirectory());
			score3 = new MaxEntScanScore3(vsa.getSpliceModelDirectory());

			//make histograms
			spliceHistogram5 = new Histogram(-14, 14, 100);
			spliceHistogram3 = new Histogram(-14, 14, 100);

			//results writer, MUST close at end of loader life!!!!!!!!!
			bedOut = new Gzipper (tempBed);

		} catch (IOException e) {
			try { if (bedOut!= null) bedOut.close(); } catch (IOException e1) {}
			e.printStackTrace();
			failed = true;
			tempBed.delete();
		}
	}


	private void loadChromosomeData() {
		try {

			//find and load sequence
			ReferenceSequence p = vsa.getFasta().getSequence(chromName);
			if (p == null ) throw new IOException ("\n\nFailed to find or load a fasta sequence for '"+chromName+"', aborting.\n");
			workingSequence = new String(p.getBases());
			workingSequence = workingSequence.toUpperCase();

			//fetch transcripts
			transcripts = vsa.getChromTranscripts().get(chromName);
			System.out.println("\tAnnotating: "+chromName+"\tLen: "+workingSequence.length()+"\t Trans: "+transcripts.length);


		} catch (Exception e){
			e.printStackTrace();
			Misc.printErrAndExit("\nDo your chromosome names in the gene table and indexed fasta match? Watch for 'chr'.\n");
		}
	}

	public void run() {	
		try {

			loadChromosomeData();
			scoreKnownSplices();
			bedOut.close();
			workingSequence = null;
			
		} catch (Exception e) {
			failed = true;
			System.err.println("\nError: problem annotating splices for "+chromName );
			bedOut.getGzipFile().delete();
			bedOut.closeNoException();
			bedOut = null;
			e.printStackTrace();
		}
	}


	private  void scoreKnownSplices() throws IOException {

		int lenMinOne = workingSequence.length() - 1;
		HashSet<String> scored5Plus = new HashSet<String>();
		HashSet<String> scored3Plus = new HashSet<String>();
		HashSet<String> scored5Minus = new HashSet<String>();
		HashSet<String> scored3Minus = new HashSet<String>();


		//for each transcript
		for (UCSCGeneLine g: transcripts){
			ExonIntron[] introns = g.getIntrons();
			if (introns == null) continue;
			boolean plusStrand = g.getStrand().equals("+");
			//mask both 5' and 3' junctions
			int startJunc = 0;
			int endJunc = 0;
			for (int i=0; i< introns.length; i++){
				//plus strand 
				if (plusStrand){
					//5'
					startJunc = introns[i].getStart()-3;     
					endJunc = introns[i].getStart()+6;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					String coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored5Plus.contains(coor) == false){
						scored5Plus.add(coor);
						String seq = workingSequence.substring(startJunc, endJunc);
						double score = score5.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE) {
							spliceHistogram5.count(score);
							bedOut.println(chromName+"\t"+startJunc+"\t"+endJunc+"\t5_"+seq+"\t"+score+"\t+");
							numUniqueJunctions++;
						}
					}
					//3'
					startJunc = introns[i].getEnd()-20;
					endJunc = introns[i].getEnd()+3;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					coor = startJunc+"-"+endJunc;
					if (scored3Plus.contains(coor) == false){
						scored3Plus.add(coor);
						String seq = workingSequence.substring(startJunc, endJunc);
						double score = score3.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE){
							spliceHistogram3.count(score);
							//System.out.println("+ 3\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);
							bedOut.println(chromName+"\t"+startJunc+"\t"+endJunc+"\t3_"+seq+"\t"+score+"\t+");
							numUniqueJunctions++;
						}
						
					}
				}
				//minus strand
				else {
					//3'
					startJunc = introns[i].getStart()-3;
					endJunc = introns[i].getStart()+20;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					String coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored3Minus.contains(coor) == false){
						scored3Minus.add(coor);
						String seq = workingSequence.substring(startJunc, endJunc);
						seq = Seq.reverseComplementDNA(seq);
						double score = score3.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE) {
							spliceHistogram3.count(score);
							//System.out.println("- 3\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);
							bedOut.println(chromName+"\t"+startJunc+"\t"+endJunc+"\t3_"+seq+"\t"+score+"\t-");
							numUniqueJunctions++;
						}
						
					}

					//5'
					startJunc = introns[i].getEnd()-6;
					endJunc = introns[i].getEnd()+3;
					if (startJunc < 0) startJunc = 0;
					if (endJunc > lenMinOne) endJunc = lenMinOne;
					coor = startJunc+"-"+endJunc;
					//add to histogram of known splice scores?
					if (scored5Minus.contains(coor) == false){
						scored5Minus.add(coor);
						String seq = workingSequence.substring(startJunc, endJunc);
						seq = Seq.reverseComplementDNA(seq);
						double score = score5.scoreSequenceWithChecks(seq);
						if (score != Double.MIN_VALUE) {
							spliceHistogram5.count(score);
							//System.out.println("- 5\t"+startJunc+"\t"+endJunc+"\t"+seq+"\t"+score);
							bedOut.println(chromName+"\t"+startJunc+"\t"+endJunc+"\t5_"+seq+"\t"+score+"\t-");
							numUniqueJunctions++;
						}
					}
				}
			}
		}
	}






	public boolean isFailed() {
		return failed;
	}

	public Gzipper getBedOut() {
		return bedOut;
	}


	public Histogram getSpliceHistogram5() {
		return spliceHistogram5;
	}


	public Histogram getSpliceHistogram3() {
		return spliceHistogram3;
	}


	public int getNumUniqueJunctions() {
		return numUniqueJunctions;
	}


}
