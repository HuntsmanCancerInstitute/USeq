package util.bio.wrappers;

import java.util.*;
import java.io.*;
import java.util.regex.*;

import util.gen.*;

/**
 * This is a bit of a hack from the SimpleBlast, expand as needed.
 */
public class BL2Seq {
	//fields
	private String fullPathToBL2seq; //full path filename to bl2seqprogram
	private String[] rawUnparsedResults;
	private String[] commandArray;
	private BL2SeqHit[] hits;
	
	//defaults, add more as needed and modify commandArray
	private double rawScoreCutOff = 10;
	private double e = 0.1;	//expectation cut off
	private int W = 7;		//word size
	private int E = -1;		//cost to extend a gap, use positive numbers -1 is for default behavior
	
	public BL2Seq(String fullPathToBL2seq){
		this.fullPathToBL2seq = fullPathToBL2seq;
	}
	
	//public BL2SeqHit[] blastSequences(){}
	
	/**Uses BLAST bl2seq program to align two sequences found in the seq files, won't work on windows?*/
	public BL2SeqHit[] blastFASTAFiles(File seq1File, File seq2File){
		//convert to String[], cannot use a text since exec uses as stringTokenizer to bust up command, if spaces exist in the folder or  file names then this creates havoc
		commandArray = new String[]{fullPathToBL2seq,"-p","blastn","-i", IO.getFullPathName(seq1File),"-j", 
				IO.getFullPathName(seq2File),"-e", e+"","-W",W+"", "-E", E+""};
		
		ArrayList dataArrayList = new ArrayList(500);
		try {
			//System.out.println("launching blast...");
			Runtime rt = Runtime.getRuntime();
			//rt.traceInstructions(true); //for debugging
			//rt.traceMethodCalls(true); //for debugging
			Process p = rt.exec(commandArray);
			BufferedReader data = new BufferedReader(new InputStreamReader(p.getInputStream()));
			//BufferedReader data = new BufferedReader(new InputStreamReader(p.getErrorStream())); //for debugging
			String line;
			while ((line = data.readLine()) != null){
				dataArrayList.add(line);
				//System.out.println("X: "+line);
			}
			data.close();   //this is close/null stuff is needed to invoke the garbage collector
			p.waitFor();
			p=null;
			data = null;
			rt = null;
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		dataArrayList.trimToSize();
		rawUnparsedResults = new String[dataArrayList.size()];
		dataArrayList.toArray(rawUnparsedResults);
		
		//parseBlast results
		return parseBl2SeqResults(rawUnparsedResults);
	}
	
	/**This takes an array of BLAST bl2seq results and parses the info into an array*/
	public BL2SeqHit[] parseBl2SeqResults(String[] data) {
		
		//make some precompiled matchers
		Pattern patScore = Pattern.compile("Score = .+\\((\\d+)\\)");
		Pattern patSSS = Pattern.compile("(\\d+)\\s+([a-zA-Z-]+)\\s+(\\d+)");
		Pattern patLambda = Pattern.compile("^Lambda.+");
		
		int size = data.length;
		ArrayList parsedData = new ArrayList(size/8);
		try {
			for (int i=6; i< size; i++){
				Matcher w = patScore.matcher(data[i]);
				Matcher l = patLambda.matcher(data[i]);
				
				if (w.find()){
					//get score
					int score = new Integer((w.group(1))).intValue();
					//get orientation
					i+=2;
					String o = data[i].substring(17,18);
					int ori;
					if (o.equals("P")) ori =0;
					else ori=1;
					
					//get everything else
					i+=3;
					String seq1Ln = data[i];
					i+=2;
					String seq2Ln = data[i];
					
					//extract start stop seq
					ArrayList datSeq1 = parseLine(seq1Ln, patSSS);
					ArrayList datSeq2 = parseLine(seq2Ln, patSSS);
					
					//look ahead to see if there are any more lines
					boolean inSeq = true;
					while (inSeq){
						i+=3;
						if (data[i].startsWith("Query")){
							//get additional lines
							ArrayList new1 = parseLine(data[i], patSSS);
							i+=2;
							ArrayList new2 = parseLine(data[i], patSSS);
							
							//join with previous seqs and change stop numbers
							datSeq1.set(2,((String)datSeq1.get(2)).concat((String)new1.get(2)));
							datSeq2.set(2,((String)datSeq2.get(2)).concat((String)new2.get(2)));
							datSeq1.set(1, new1.get(1));
							datSeq2.set(1, new2.get(1));
						}
						else inSeq = false;
					}
					//add data
					if (score >=rawScoreCutOff){
						BL2SeqHit hit = new BL2SeqHit();
						hit.setRawScore(score);
						hit.setStartSeq1( ( (Integer)datSeq1.get(0) ).intValue());
						hit.setStopSeq1( ( (Integer)datSeq1.get(1) ).intValue());
						hit.setSeq1(  (String)datSeq1.get(2) );
						hit.setStartSeq2( ( (Integer)datSeq2.get(0) ).intValue());
						hit.setStopSeq2( ( (Integer)datSeq2.get(1) ).intValue());
						hit.setSeq2(  (String)datSeq2.get(2) );
						hit.setOrientation(ori);
						parsedData.add(hit);
					}
				}
			}
		} catch (NumberFormatException n){
			n.printStackTrace();
		}
		hits = new BL2SeqHit[parsedData.size()];
		parsedData.toArray(hits);
		return hits;
	}
	
	public static ArrayList parseLine(String line, Pattern pat){
		//returns start, stop, sequence
		ArrayList dat = new ArrayList(3);
		Matcher mat = pat.matcher(line);
		mat.find();
		dat.add(new Integer(mat.group(1)));
		dat.add(new Integer (mat.group(3)));
		dat.add(mat.group(2));      
		return dat;        
	}
	public double getExpectation() {
		return e;
	}
	public void setExpectation(double e) {
		this.e = e;
	}
	public String getFullPathToBL2seq() {
		return fullPathToBL2seq;
	}
	public void setFullPathToBL2seq(String fullPathToBL2seq) {
		this.fullPathToBL2seq = fullPathToBL2seq;
	}
	public double getRawScoreCutOff() {
		return rawScoreCutOff;
	}
	public void setRawScoreCutOff(double rawScoreCutOff) {
		this.rawScoreCutOff = rawScoreCutOff;
	}
	public String[] getRawUnparsedResults() {
		return rawUnparsedResults;
	}
	public void setRawUnparsedResults(String[] rawUnparsedResults) {
		this.rawUnparsedResults = rawUnparsedResults;
	}
	public int getW() {
		return W;
	}
	public void setW(int w) {
		W = w;
	}
	public String[] getCommandArray() {
		return commandArray;
	}
	public BL2SeqHit[] getHits() {
		return hits;
	}
	public void setGapExtensionCost(int e) {
		E = e;
	}
}
