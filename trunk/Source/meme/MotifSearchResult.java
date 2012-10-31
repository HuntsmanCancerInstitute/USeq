package meme;

import java.util.*;

/**Holds the results of searching a multi-fasta seq file with a particular motif. */
public class MotifSearchResult implements Comparable{
    
    //fields
    private MemeMotif memeMotif;  //reference to MemeMotif used to generate the search result
    //use this to retrieve all the items needed, MemeMotif has a ref to the MemeParser
    //which contains many more fields
    private double memeLLPSPMCutOffScore;  //used to calculate whether a seq contains the motif
    private double addOneLLPSPMCutOffScore;     //ditto
    private String searchedSeqFile; //file containing sequences searched for the presence of this motif
    private ArrayList hitsMemeLLPSPM; //arraylist containing an array of MotifHits for each seq in the file found to have hits
    private ArrayList hitsAddOneLLPSPM; //ditto
    private int numHitsMemeLLPSPM; //number of hits to this motif
    private int numHitsAddOneLLPSPM; //ditto
    private int numSeqsSearched; //total number of sequences in the fasta file searched
    
    //fields for compareTo method
    private double motifSig;  //signature representing the motif, file used to make the motif. the motif number, ie 2.1  (file 2 . motif 1)
    private int fileNumber;  //from hash in MemeResults
    
    //constructor
    public MotifSearchResult(MemeMotif mm, String file, int numberOfSeqsSearched){
        memeMotif = mm;
        searchedSeqFile = file;
        numSeqsSearched = numberOfSeqsSearched;
    }
    
    public int compareTo(Object otherObject){
        //sort first by motif then by file searched
        MotifSearchResult other = (MotifSearchResult)otherObject;
        if (motifSig<other.motifSig) return -1;
        if (motifSig>other.motifSig) return 1;
        //if equal order by file searched
        if (fileNumber<other.fileNumber) return -1;
        if (fileNumber>other.fileNumber) return 1;
        return 0;
    }
    
    //setter methods
    public void setMotifSig(double data){motifSig = data;}
    public void setFileNumber(int data){fileNumber = data;}
    public void setMemeLLPSPMCutOffScore(double data){memeLLPSPMCutOffScore = data;}
    public void setAddOneLLPSPMCutOffScore(double data){addOneLLPSPMCutOffScore = data;}
    public void setHitsMemeLLPSPM(ArrayList data){hitsMemeLLPSPM = data;}
    public void setHitsAddOneLLPSPM(ArrayList data){hitsAddOneLLPSPM = data;}
    public void setNumHitsMemeLLPSPM(int data) {numHitsMemeLLPSPM = data;}
    public void setNumHitsAddOneLLPSPM(int data) {numHitsAddOneLLPSPM = data;}
    
    //getter methods
    public int getFileNumber(){return fileNumber;}
    public double getMotifSig(){return motifSig;}
    public MemeMotif getMemeMotif(){return memeMotif;}
    public double getMemeLLPSPMCutOffScore(){return memeLLPSPMCutOffScore;}
    public double getAddOneLLPSPMCutOffScore(){return addOneLLPSPMCutOffScore;}
    public String getSearchedSeqFile(){return searchedSeqFile;}
    public ArrayList getHitsMemeLLPSPM(){return hitsMemeLLPSPM;}
    public ArrayList getHitsAddOneLLPSPM(){return hitsAddOneLLPSPM;}
    public int getNumHitsMemeLLPSPM() {return numHitsMemeLLPSPM;}
    public int getNumHitsAddOneLLPSPM() {return numHitsAddOneLLPSPM;}
    public int getnumSeqsSearched(){return numSeqsSearched;}
}