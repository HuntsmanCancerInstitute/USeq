package selex;
import util.gen.*;

/**Holds information related to a particular parsed SELEX oligo.*/
public class SubSeq {
     
    //parameters
    String seq;         //the actual sequence without flanking cutsites
    String subSeq;      //seq flanked by cutSites
    int length;         //length of seq
    int[] qualScores;   //array of quality scores for each base in subSeq
    SeqRead seqRead;    //reference to parent SeqRead
    SelexParams sp;
    
    //constructor
    public SubSeq(String sequence, int[] qualityScores, SelexParams SP, SeqRead srObjectRef) {
        seq = sequence;
        qualScores = qualityScores;
        seqRead = srObjectRef;
        sp=SP;
        length = seq.length();
        printSubSeq();
    }
    
    //getter methods
    public String getSeq(){return seq;}
    public String getSubSeq(){
        if (subSeq==null){
            String cutSite = sp.getRestSite();
            subSeq = cutSite+seq+cutSite;
        }
        return subSeq;
    }
    public int getLength(){return length;}
    public int[] getQualScores(){return qualScores;}
    public SeqRead getSeqRead(){return seqRead;}
    
    
    //make getter method to pull full path info and use as id
    
    
    
    //primary methods
    public void printSubSeq() {
       sp.printSave(
       "*****SubSeq*****"+
       "\nSeq: "+seq+
       "\nQuality Scores: "+ Misc.intArrayToString(qualScores, " ")+
       "\nLength: "+length+"\n");
    }
    
}