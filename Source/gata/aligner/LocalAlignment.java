package gata.aligner;

import gata.main.*;

import java.io.*;
import java.util.*;


/**
 * @author nix
 * Similar to but not exactly like Alignment obj, represents the Local Alignment returned by
 * Blast.  It's saved and made availible to GATAPlotter.*/
public class LocalAlignment implements Serializable{
    //fields
    private String seq1;
    private int startSeq1;
    private int stopSeq1;
    private String seq2;
    private int startSeq2;
    private int stopSeq2;
    private int ori;        //0 is forward, 1 reverse-complement
    private int score;
	private double bitScore;
	private String expect;
    private String oriWord = "+/+";
    private static AlignParams ap;
    private static int match=0;
    private static int misMatch=0;
    private static int gapCreate=0;
    private static int gapExtend=0;
    private static int window=0;
    private static int minScore=0;
	private final double LAMBDA;
	private final double K;
	private final double EFF_M;
	private final double EFF_N;
    
    //constructor
    public LocalAlignment(String seq_1, int startSeq_1, int stopSeq_1, String seq_2,
    int startSeq_2,int stopSeq_2, int ori_, AlignParams ap_, int score_) {        
        seq1=seq_1;
        this.startSeq1=startSeq_1;
		this.stopSeq1=stopSeq_1;
		this.seq2=seq_2;
		this.startSeq2=startSeq_2;
		this.stopSeq2=stopSeq_2;
		this.ori=ori_;
		ap=ap_;
		this.score=score_;
		match = ap.getMATCH();
        misMatch = ap.getMISMATCH();
        gapCreate = ap.getGAP_CREATE();
        gapExtend = ap.getGAP_EXT();
        window = ap.getWIN_SIZE();
        minScore = ap.getMIN_SCORE();
        if (ori==1) oriWord = "+/-";
		//calculate bitScore and E
		LAMBDA = ap.getLAMBDA();
		K = ap.getK();
		EFF_M = ap.getEFF_M();
		EFF_N = ap.getEFF_N();
		expect = GATAUtil.calculateEValue(LAMBDA,K,score,EFF_N,EFF_M);
		bitScore = GATAUtil.convertRawScoreToBit(LAMBDA,K,score);
		        
    }
    //getter methods
    public AlignParams getAP(){return ap;}
    public int getOri(){return ori;}
    public int getLAScore(){return score;}
    
    //primary methods
    public void printLocalAlignment(){
        System.out.println("Local Alignment");
        System.out.println("        " + startSeq1 + GATAUtil.spaces(seq1,startSeq1) + stopSeq1);
        System.out.println("refseq: " + seq1);
        System.out.println("        " + GATAUtil.genDashes(seq1, seq2));
        System.out.println("cmpseq: " + seq2);
        System.out.println("        " + startSeq2 + GATAUtil.spaces(seq2,startSeq2) + stopSeq2);
		System.out.println(" Score: " + bitScore+"("+score+") bits,  Expect: "+expect+",  Ori: " + oriWord);
        System.out.println();
    }
    
    public String getLocalAlignString(int subStart, int subLength){
        //find relative position in local alignment
        int firstBase = subStart - startSeq1+1;  //number of real bases to move in
        int len = seq1.length();
        int counter = 0;
        int capStart;
        for (capStart=0; capStart<len; capStart++){
            if (counter == firstBase) break;
            if ((seq1.substring(capStart,capStart+1)).equals("-")) continue;
            counter++;
        }
        
        capStart--;
        int capLen = capStart+subLength;
        String dashes = GATAUtil.genDashes(seq1, seq2);
        //modify substrings, capitalizing or | -> *
        String modSeq1 = seq1.substring(0,capStart)+
                        (seq1.substring(capStart,capLen)).toUpperCase()+
                        seq1.substring(capLen);
        String modSeq2 = seq2.substring(0,capStart)+
                        (seq2.substring(capStart,capLen)).toUpperCase()+
                        seq2.substring(capLen);
        String modDashes = dashes.substring(0,capStart)+
                        (dashes.substring(capStart,capLen)).replace('|','*')+
                        dashes.substring(capLen);
        
        return  
                "Local Alignment"+
                "\n        " + startSeq1 + GATAUtil.spaces(seq1,startSeq1) + stopSeq1+
                "\nrefseq: " + modSeq1+
                "\n        " + modDashes+
                "\ncmpseq: " + modSeq2+
                "\n        " + startSeq2 + GATAUtil.spaces(seq2,startSeq2) + stopSeq2+
				"\nScore: " + bitScore+"("+score+") bits,  Expect: "+expect+",  Ori: " + oriWord+"\n";
    }
    
    
    public ArrayList processLocalAlignment(){
        ArrayList subAligns;
        
        //run window across alignment
        int[][] windows = makeWindowsFromAlign();
        
        //optimize score by trimming ends
        ArrayList optiAligns = optimizeScores(windows);
        //check to see if there are any optiAligns (survived the score cut off)
        if (optiAligns.size()==0){
            return subAligns=null;
        }
        
        //compare each window, if same score across entire overlap then fuse
        //only run if more than one subalignment
        if(optiAligns.size()>1){
            optiAligns = mergeSubAligns(optiAligns);
        }
        
        //remove those alignments with the same or lesser score entirely contained within another
        //only run if more than one subalignment
        if (optiAligns.size()>1){
            optiAligns = clearDups(optiAligns);
        }
        //make alignment objects
        subAligns = makeAlignmentObjects(optiAligns);
        return subAligns;
    }
    
    public int[][] makeWindowsFromAlign() {
        //runs a sliding window across an alignment returning the start,stop indexes
        int len = seq1.length()-(window-1);
        //check to see if local alignment is bigger than window size, just return loc align
        if (len<1){
            int[][] subAligns = {{0, seq1.length()-1}};
            return subAligns;
        }
        int[][] subAligns = new int[(len)][];
        for (int i=0; i<len; i++){
            //save window (start, stop)
            subAligns[i]= new int[2];
            subAligns[i][0]=i;
            subAligns[i][1]=i+window-1;
        }
        return subAligns;
    }
    
    public ArrayList optimizeScores(int[][] aligns){
        //trims back either stop attempting to optimize scores
        int len = aligns.length;
        ArrayList optiAligns = new ArrayList(len);
        //System.out.println("optimizingScores\n");
        for (int i=0; i<len; i++){
            String subSeq1 = seq1.substring(aligns[i][0], aligns[i][1]+1);
            String subSeq2 = seq2.substring(aligns[i][0], aligns[i][1]+1);
            //System.out.println("\nold nums: "+aligns[i][0]+" "+aligns[i][1]);
            int[] x = GATAUtil.optimizeAlignmentScore(subSeq1, subSeq2, match, misMatch, gapCreate, gapExtend, minScore);
            //adjust start and stop to match index numbers,check score first
            if (x[2]>=minScore) {
                x[0] = aligns[i][0]+x[0];
                x[1] = aligns[i][0]+x[1];
                optiAligns.add(x);
                //System.out.println("adding new nums: "+x[0]+" "+x[1]+" "+x[2]);
            }
        }
        optiAligns.trimToSize();
        return optiAligns;
    }
    
    public ArrayList mergeSubAligns(ArrayList winds){
        //merges alignments that are adjacent and with the same score
        //returns an Arraylist containing int[] 0 =start, 1=stop, 2=score
        //note, assumes that the last base is aligned
        int len = winds.size();
        ArrayList mergedWinds = new ArrayList(len);
        boolean end = false;
        for (int i=0; i<len; i++){
            if (end) break;
            int[] s1 = (int[])winds.get(i);
            for (int j=i+1; j<len; j++){
                int[] s2 = (int[])winds.get(j);
                //if scores are the same and 1bp apart then merge
                if (s1[2]==s2[2] && s1[1]==s2[1]-1){
                    s1[1]++;
                    if (len == j+1){
                        mergedWinds.add(s1); //add if last one
                        end = true;
                    }
                }
                //otherwise save and go forward
                else {
                    mergedWinds.add(s1);
                    i=j-1;
                    break;
                }
            }
        }
        mergedWinds.trimToSize();
        return mergedWinds;
    }
    
    public ArrayList clearDups(ArrayList aligns){
        //removes subAligns that are the same score or less and entirely contained within a preceeding alignment
        int len = aligns.size();
        ArrayList nonDups = new ArrayList(len);
        if (len==0){
            nonDups.trimToSize();
            return nonDups;
        }
        //get first tester
        int[] align1 = (int[])aligns.get(0);
        for (int i=1; i<len; i++){
            //get new alignment to compare (0=start, 1=stop, 2=score)
            int[] align2 = (int[])aligns.get(i);

            //check whether contained within
            //align1 is bigger
            if (align1[2]>=align2[2] && align1[0]<=align2[0] && align1[1]>=align2[1]){
            }
            //align2 is bigger
            else if (align1[2]<=align2[2] && align1[0]>=align2[0] && align1[1]<=align2[1]){
                //skip align1 set align2 as align 1
                align1=new int[]{align2[0],align2[1],align2[2]};
            }
            // not correct score or not contained within
            else{
                // save align 1, set align2 as align1
                nonDups.add(align1);
                align1=new int[]{align2[0],align2[1],align2[2]};
            }
        }
        
        //add last one
        nonDups.add(align1);
        nonDups.trimToSize();
        return nonDups;
    }
    
    public ArrayList makeAlignmentObjects(ArrayList alignNums){
        
        //takes int array and converts it into actual objects
        int len = alignNums.size();
        ArrayList aligns = new ArrayList(len);
        
        //get ref number arrays
        int[] seq1Nums = GATAUtil.getBaseNumArray(seq1,startSeq1);
        int[] seq2Nums;
        if (ori==0) seq2Nums = GATAUtil.getBaseNumArray(seq2,startSeq2); //forward
        //reverse
        else {
            StringBuffer sb = new StringBuffer(seq2);
            sb.reverse();
            String revSeq2 = new String(sb);
            seq2Nums = GATAUtil.getBaseNumArray(revSeq2,stopSeq2);
            //reverse the array
            int len2 = seq2Nums.length;
            int[] revSeq2Nums = new int[len2];
            for (int i=0; i<len2; i++){
                revSeq2Nums[len2-i-1] = seq2Nums[i];
            }
            seq2Nums = revSeq2Nums;
        }
        
        //make Aligns
        for (int i=0; i<len; i++){
            int[] subAlign = (int[])alignNums.get(i);
            int start1= seq1Nums[subAlign[0]];
            int stop1= seq1Nums[subAlign[1]];
            String subSeq1 = seq1.substring(subAlign[0], subAlign[1]+1);
            
            int start2 = seq2Nums[subAlign[0]];
            int stop2 = seq2Nums[subAlign[1]];
            
            String subSeq2 = seq2.substring(subAlign[0], subAlign[1]+1);
            
            aligns.add(new Alignment(subSeq1, start1, stop1, subSeq2, start2, stop2, subAlign[2], ori, ap, this));
        }
        return aligns;
    }
}